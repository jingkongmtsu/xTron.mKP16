# - Find NETLIB math library (BLAS+LAPACK)
#
# This module defines the following variables:
#
#   NETLIB_FOUND         : True if we are able to find both blas and lapack libraries
#   NETLIB_LIBRARIES     : the blas and lapack library to link against.

# firstly assume that netlib is not found
set(NETLIB_FOUND FALSE)

# the NETLIB_ROOT must be defined
# which is ensured by the perl configure file
set(NETLIB_ROOT $ENV{NETLIB_ROOT})
if(NOT NETLIB_ROOT)
	message(${NETLIB_ROOT})
	message(FATAL_ERROR "NETLIB_ROOT is not defined in findNETLIB.cmake")
endif(NOT NETLIB_ROOT)

# blas library
find_library(NETLIB_BLAS_PATH blas PATHS ${NETLIB_ROOT} NO_DEFAULT_PATH)
if (NOT NETLIB_BLAS_PATH)
	message(${NETLIB_BLAS_PATH})
	message(FATAL_ERROR "NETLIB_BLAS can not be found")
endif(NOT NETLIB_BLAS_PATH)

# lapack library
find_library(NETLIB_LAPACK_PATH lapack PATHS ${NETLIB_ROOT} NO_DEFAULT_PATH)
if (NOT NETLIB_LAPACK_PATH)
	message(${NETLIB_LAPACK_PATH})
	message(FATAL_ERROR "NETLIB_LAPACK can not be found")
endif(NOT NETLIB_LAPACK_PATH)

# now add in all of these
set(NETLIB_LIBRARIES ${NETLIB_BLAS_PATH} ${NETLIB_LAPACK_PATH})
#message(${NETLIB_LIBRARIES})

# finally reset the status of NETLIB
set(NETLIB_FOUND TRUE)

