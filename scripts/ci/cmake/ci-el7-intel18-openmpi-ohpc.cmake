# Client maintainer: chuck.atkins@kitware.com

find_package(EnvModules REQUIRED)

env_module(purge)
env_module(load intel)
env_module(load py2-numpy)
env_module(load openmpi3)
env_module(load phdf5)
env_module(load py2-mpi4py)

set(ENV{CC}  icc)
set(ENV{CXX} icpc)
set(ENV{FC}  ifort)
set(ENV{CFLAGS} -Werror)
set(ENV{CXXFLAGS} -Werror)
set(ENV{FFLAGS} "-warn errors")

set(ENV{CMAKE_PREFIX_PATH} "/opt/libfabric/1.6.0:$ENV{CMAKE_PREFIX_PATH}")

set(dashboard_cache [[
ADIOS2_USE_BZip2:BOOL=ON
ADIOS2_USE_Blosc:BOOL=ON
ADIOS2_USE_DataMan:BOOL=ON
ADIOS2_USE_Fortran:BOOL=ON
ADIOS2_USE_HDF5:BOOL=ON
ADIOS2_USE_MPI:BOOL=ON
ADIOS2_USE_Python:BOOL=ON
ADIOS2_USE_ZeroMQ:STRING=ON

MPIEXEC_MAX_NUMPROCS:STRING=8
PYTHON_EXECUTABLE:FILEPATH=/usr/bin/python2.7
]])

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
list(APPEND CTEST_UPDATE_NOTES_FILES "${CMAKE_CURRENT_LIST_FILE}")
include(${CMAKE_CURRENT_LIST_DIR}/ci-common.cmake)
