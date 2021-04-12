cmake_minimum_required(VERSION 3.10)

# set the project name
project(SymPiFE)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Set MFEM_HOME and MFEM_BUILD here
set(MFEM_HOME "/u/04/chonel1/unix/Softs/mfem")
set(MFEM_BUILD "mpi-omp-build4.2")

include_directories(include)
add_subdirectory(src)
