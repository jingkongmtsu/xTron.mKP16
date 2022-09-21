#
# add in macros to the program
#

# if we use mkl turn on the macro
if (BLAS_VENDOR STREQUAL "mkl")
	add_definitions( -DUSE_MKL )
endif (BLAS_VENDOR STREQUAL "mkl")

if (BLAS_VENDOR STREQUAL "netlib")
	add_definitions( -DUSE_GENERAL_BLAS -DUSE_GENERAL_LAPACK )
endif (BLAS_VENDOR STREQUAL "netlib")

# let's see whether we invork GDM
if (WITH_GDM_LIB STREQUAL "1")
	add_definitions( -DWITH_GDM )
endif (WITH_GDM_LIB STREQUAL "1")

# do we add in debug choice?
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions( -DDEBUG )
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

# does the program use the float rather than double?
if(WITH_SINGLE_PRECISION)
	add_definitions( -DWITH_SINGLE_PRECISION )
endif(WITH_SINGLE_PRECISION)

# if we compile code for MIC, we need to define some macro
if (WITH_INTEL_PHI_COPROC STREQUAL "1")
	add_definitions( -DEMUL_WITH_MIC )
endif()

# add macro for testing
if(TARGET_MODULE STREQUAL "DXCINTS")
	if (WITH_INTEL_PHI_COPROC STREQUAL "1")
		add_definitions( -DESPINTS_TIMING )
	endif (WITH_INTEL_PHI_COPROC STREQUAL "1")
endif(TARGET_MODULE STREQUAL "DXCINTS")

