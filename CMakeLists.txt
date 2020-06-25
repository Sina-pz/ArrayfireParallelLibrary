cmake_minimum_required(VERSION 3.8)

project(mainHomeworkPoiseuille)

set(CMAKE_CXX_STANDARD 11)

find_package(ArrayFire REQUIRED)

add_executable(mainHomeworkPoiseuilleZouHe mainHomeworkPoiseuilleZouHe.cpp
../AF_FILE_IO/AF_FILE_IO.h
../AF_FILE_IO/AF_FILE_IO.cpp
../AF_POISEUILLE2D/AF_POISEUILLE2D.h
../AF_POISEUILLE2D/AF_POISEUILLE2D.cpp
../AF_VTK/AF_VTK.h
../AF_VTK/AF_VTK.cpp
../CPP_UTILITY_FUNCTIONS/CPP_UTILITY_FUNCTIONS.h
../CPP_UTILITY_FUNCTIONS/CPP_UTILITY_FUNCTIONS.cpp
../tinyxml2_edit/tinyxml2_edit.h
../tinyxml2_edit/tinyxml2_edit.cpp
../tinyxml2_edit/tinyxml2_customFunctions.cpp
)

# To use Unified backend, do the following.
# Unified backend lets you choose the backend at runtime
target_link_libraries(mainHomeworkPoiseuilleZouHe ArrayFire::af)

