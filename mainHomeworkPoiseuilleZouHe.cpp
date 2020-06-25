//////////////////////////////////////////////////////////////
//
//  Contributors: 
//	1) Sebastien Leclaire (sebastien.leclaire@polymtl.ca)
//  2) ...
//
//////////////////////////////////////////////////////////////

#include "../AF_VTK/AF_VTK.h"
#include "../AF_VTK/AF_VTK.cpp"

#include "../CPP_UTILITY_FUNCTIONS/CPP_UTILITY_FUNCTIONS.h"
#include "../CPP_UTILITY_FUNCTIONS/CPP_UTILITY_FUNCTIONS.cpp"

#include "../AF_FILE_IO/AF_FILE_IO.h"
#include "../AF_FILE_IO/AF_FILE_IO.cpp"

#include "../AF_POISEUILLE2D/AF_POISEUILLE2D.h"
#include "../AF_POISEUILLE2D/AF_POISEUILLE2D.cpp"

#include "../tinyxml2_edit/tinyxml2_edit.h"

#include <arrayfire.h>

#include <stdio.h>

#include <cstdio>
#include <cstdlib>
#include <iostream> // std::cout, std::fixed
#include <fstream>
#include <iomanip>  // std::setprecision
#include <algorithm> // min, max
#include <limits>

using namespace af;

typedef double T;
#define T_BACKEND f64
#define AF_BACKEND_COMPUTE AF_BACKEND_CPU
//#define AF_BACKEND_COMPUTE AF_BACKEND_OPENCL
//#define AF_BACKEND_COMPUTE AF_BACKEND_CUDA

#undef min
#undef max 

const T     EPS = std::numeric_limits<T>::epsilon();
const long long int INTMIN = std::numeric_limits<long long int>::min();

template<typename T>
T get_latticeWeightD2Q9_NSWC(const T alpha, int &nq, array& weight, array& phiWeight, array& varPhiWeight, array& zeroVelocityWeight, const dtype arrayType)
{
	nq = 9;

	T weightHost[] = { (T)4 / (T)9,
	(T)1 / (T)9, (T)1 / (T)36, (T)1 / (T)9, (T)1 / (T)36,
	(T)1 / (T)9, (T)1 / (T)36, (T)1 / (T)9, (T)1 / (T)36 };
	array weightDevice(1, nq, weightHost);
	weight = weightDevice;

	T phiWeightHost[] = { (T)0,
		(T)1 / (T)5, (T)1 / (T)20, (T)1 / (T)5, (T)1 / (T)20,
		(T)1 / (T)5, (T)1 / (T)20, (T)1 / (T)5, (T)1 / (T)20 };
	array phiWeightDevice(1, nq, phiWeightHost);
	phiWeight = phiWeightDevice.as(arrayType);

	T varPhiWeightHost[] = { (T)1,
		-(T)1 / (T)5, -(T)1 / (T)20, -(T)1 / (T)5, -(T)1 / (T)20,
		-(T)1 / (T)5, -(T)1 / (T)20, -(T)1 / (T)5, -(T)1 / (T)20 };
	array varPhiWeightDevice(1, nq, varPhiWeightHost);
	varPhiWeight = varPhiWeightDevice.as(arrayType);

	zeroVelocityWeight = phiWeight + varPhiWeight * alpha;

	T soundSpeedWeight = (T)3 / (T)5;
	T cs = std::sqrt(soundSpeedWeight*((T)1 - alpha));

	return cs;
}

template<typename T>
void get_latticeVelocityD2Q9(array& cx, array& cy, array& cNorm, const dtype arrayType)
{
	T cxHost[] = { (T)0, (T)0, (T)1, (T)1, (T)1, (T)0, (T)-1, (T)-1, (T)-1 };
	array cxDevice(1, 9, cxHost);
	cx = cxDevice;

	T cyHost[] = { (T)0, (T)1, (T)1, (T)0, (T)-1, (T)-1, (T)-1, (T)0, (T)1 };
	array cyDevice(1, 9, cyHost);
	cy = cyDevice;

	cNorm = sqrt(cx * cx + cy * cy); cNorm(0) = (T)1;
}

template<typename T>
void get_bbIndexD2Q9(array& bbIndex, const dtype arrayType)
{
	int bbIndexHost[] = { 0, 5, 6, 7, 8, 1, 2, 3, 4 };
	array bbIndexDevice(1, 9, bbIndexHost);
	bbIndex = bbIndexDevice;
}

int main(int argc, char *argv[])
{
	af::setBackend(AF_BACKEND_COMPUTE);
	af::info();

	int npre = loadXMLvariableINT16("mainHomeworkPoiseuilleZouHe_inputs.xml", "npre");
	int ns = loadXMLvariableINT16("mainHomeworkPoiseuilleZouHe_inputs.xml", "ns");
	int nc = loadXMLvariableINT16("mainHomeworkPoiseuilleZouHe_inputs.xml", "nc");
	int isForceStartScratch = loadXMLvariableINT16("mainHomeworkPoiseuilleZouHe_inputs.xml", "isForceStartScratch");

	long long int nt = loadXMLvariableINT64("mainHomeworkPoiseuilleZouHe_inputs.xml", "nt");
	int nx = loadXMLvariableINT16("mainHomeworkPoiseuilleZouHe_inputs.xml", "nx");
	int ny = loadXMLvariableINT16("mainHomeworkPoiseuilleZouHe_inputs.xml", "ny");

	T rho0 = loadXMLvariableDBL("mainHomeworkPoiseuilleZouHe_inputs.xml", "rho0");
	T deltaUx = loadXMLvariableDBL("mainHomeworkPoiseuilleZouHe_inputs.xml", "deltaUx");

	T nu = loadXMLvariableDBL("mainHomeworkPoiseuilleZouHe_inputs.xml", "nu");
	
	T omega = ((T)2 / ((T)6 * nu + (T)1));

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set lattice
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int nq;
	T alpha = (T)4 / (T)9;
	T cs;

	array weight, phiWeight, varPhiWeight;
	array zeroVelocityWeight;
	array cx, cy, cNorm;
	array bbIndex;

	cs = get_latticeWeightD2Q9_NSWC<T>(alpha, nq,
		weight, phiWeight, varPhiWeight, zeroVelocityWeight, T_BACKEND);

	get_latticeVelocityD2Q9<T>(cx, cy, cNorm, T_BACKEND);

	get_bbIndexD2Q9<T>(bbIndex, T_BACKEND);
	
	int *cxIntHost = cx.as(s32).host<int>();
	int *cyIntHost = cy.as(s32).host<int>();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set geometry
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	array whereIsSolid = setWhereIsSolid<T>(nx, ny, T_BACKEND);

	const long long int nSize = nx * ny;

	array isFluid = constant(1, nx, ny, T_BACKEND);
	isFluid(whereIsSolid) = 0; isFluid.eval();
	array isSolid = abs(isFluid - 1); isSolid.eval();
	 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialize BB array 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	array Ndist_isSolid = constant(0, whereIsSolid.elements(), bbIndex.elements(), T_BACKEND);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialize working array 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	array rho;
	array invRho;

	rho				 = constant(0, nx, ny, T_BACKEND);
	invRho           = constant(0, nSize, T_BACKEND);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Declare distribution function
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	array Ndist = constant(0, nSize, nq, T_BACKEND);

	long long int iT;
	array iTgpu = fileToArray<long long int>("iTgpu");
	if (!(iTgpu.scalar<long long int>() == INTMIN) && isForceStartScratch == 0)
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Loading data from previous simulation
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		iT = iTgpu.scalar<long long int>();
		std::cout << "Load solution time step : " << iT << "...";

		Ndist = fileToArray<T>("Ndist");

		std::cout << " Done." << std::endl;
	}
	else
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Start solution from scratch
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		std::cout << "Start solution from scratch." << std::endl;
		iTgpu = constant(0, 1, 1, 1, 1, s64);
		iT    = 0;

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialize density distribution
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for (int iPop = 0; iPop < nq; iPop++)
			Ndist(span, iPop) = constant(rho0, nSize, T_BACKEND)* tile(weight(iPop), nSize, 1);
				
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start main loop
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	timer start1;
	long long int iTstart = iT;
	while (iT < iTstart + nt)
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// BB BC
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		Ndist_isSolid = Ndist(whereIsSolid, bbIndex);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// MACRO
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		rho = sum(Ndist, 1);
		invRho = (T)1 / rho ;

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// VERBOSE
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (iT == iTstart + npre)
			start1 = timer::start();

		if (iT % nc == 0) {

			std::cout << "Time step : " << iT << std::endl;

			T totMass = sum<T>(rho);

			std::cout << "totMass: " << std::setprecision(16) << totMass << std::endl;

			std::string filename = "rho_" + dbl2str(iT);
			array tmp = moddims(rho, nx, ny);
			arrayToVTK<T>(filename, "rho", tmp, true);

			array ux = sum(Ndist * tile(cx, nSize, 1), 1) * invRho; ux.eval();
			filename = "ux_" + dbl2str(iT);
			tmp = moddims(ux, nx, ny);
			arrayToVTK<T>(filename, "ux", tmp, true);
		}
		 
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// SINGLE-PHASE COLLISION
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		{
			array workHyperVolume = Ndist;

			array ux = sum(workHyperVolume * tile(cx, nSize, 1), 1) * invRho + deltaUx;
			array uy = sum(workHyperVolume * tile(cy, nSize, 1), 1) * invRho;
					
			workHyperVolume = tile(ux, 1, nq) * tile(cx, nSize, 1) + tile(uy, 1, nq) * tile(cy, nSize, 1);
			workHyperVolume = ((T)3 + (T)9 / (T)2 * workHyperVolume) * workHyperVolume;
			workHyperVolume -= tile((T)3 / (T)2 * (ux * ux + uy * uy), 1, nq);
			workHyperVolume *= tile(weight, nSize, 1);
			workHyperVolume += tile(zeroVelocityWeight, nSize, 1);
			workHyperVolume *= tile(rho, 1, nq);

			Ndist = Ndist + (workHyperVolume - Ndist) * omega;
			Ndist.eval();						
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// BB BC
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		Ndist(whereIsSolid, span) = Ndist_isSolid;
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// STREAMING
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Streaming
		Ndist = moddims(Ndist, nx, ny, nq);
		for (int iPop = 1; iPop < nq; iPop++)
			Ndist(span, span, iPop) = shift(Ndist(span, span, iPop), cxIntHost[iPop], cyIntHost[iPop]);
		Ndist = moddims(Ndist, nSize, nq);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// New time step and save solution
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		iT++;
		iTgpu += 1;
		
		if (iT % ns == 0) {

			std::cout << "Save solution time step : " << iT << "...";

			arrayToFile<T>("Ndist", Ndist);

			arrayToFile<long long int>("iTgpu", iTgpu);

			std::cout << " Done." << std::endl;
		}
		//std::getchar();
	}
	af::sync;

	printf("Elapsed seconds: %g\n", timer::stop(start1));
	printf("MLUPS: %g\n", nSize*(nt-npre)/((T)1000000*timer::stop(start1)) );

	{
		std::cout << "Time step : " << iT << std::endl;

		T totMass = sum<T>(rho);

		std::cout << "totMass: " << std::setprecision(16) << totMass << std::endl;

		std::string filename = "rho_" + dbl2str(iT);
		array tmp = moddims(rho, nx, ny);
		arrayToFile<T>(filename, tmp);	
	}

	std::cout << "Save solution time step : " << iT << "...";
	arrayToFile<T>("Ndist", Ndist);

	arrayToFile<long long int>("iTgpu", iTgpu);

	std::cout << " Done." << std::endl;
	
	return 0;
}
