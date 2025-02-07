cmake_minimum_required(VERSION 3.15)
project(FortranProject LANGUAGES Fortran)

# Specify the Fortran standard (for gfortran)
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_FLAGS "-Wall -Wextra")

# Set the output directory for the executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)   
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)   

set(FORTRAN_MOD_DIR ${CMAKE_BINARY_DIR}/mod)
file(MAKE_DIRECTORY ${FORTRAN_MOD_DIR}) 

# Specify the output directory for Fortran modules
set(CMAKE_Fortran_MODULE_DIRECTORY ${FORTRAN_MOD_DIR})
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fcheck=all -fbacktrace -finit-real=snan -Wall")
# Detect the Fortran compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "Using GNU Fortran (gfortran)")
    #set(CMAKE_Fortran_FLAGS "-O2 -g -pg -fbacktrace -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS "-O3 -fbacktrace -Wall -Wextra")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    message(STATUS "Using Intel Fortran (ifx)")
    set(CMAKE_Fortran_FLAGS "-O2 -g -traceback -qopt-report=3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(STATUS "Using Intel Fortran (ifort)")
    set(CMAKE_Fortran_FLAGS "-O2 -g -traceback -prof-gen=srcpos")
else()
endif()

enable_testing()

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)
