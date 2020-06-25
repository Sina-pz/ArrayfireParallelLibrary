#!/bin/sh

# The first the script is execute, it may be need to execute first
# "sed -i 's/\r//' make.sh" 
# to convert endofline to unix format

mkdir build
cd build
cmake ..
make

cd ..
mkdir workingDirectory
cd build
mv mainHomeworkPoiseuilleZouHe ../workingDirectory/
  
cd ..
cp mainHomeworkPoiseuilleZouHe_inputs.xml ./workingDirectory/ 
