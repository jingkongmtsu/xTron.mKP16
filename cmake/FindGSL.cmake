# - Find GNU GSL library
#
# This module defines the following variables:
#
#   GSL_FOUND            : True if GSL_INCLUDE_DIR are found
#   GSL_INCLUDE_DIRS     : set when GSL_INCLUDE_DIR found
#   GSL_LIBRARIES        : the library to link against.

# firstly assume that GSL is not found
set(GSL_FOUND FALSE)

# the GSL_ROOT may be defined or not
# if it's defined, we will search the path
# which is ensured by the perl configure file
set(GSLROOT $ENV{GSLROOT})
if(GSLROOT)

	# Find include dir
	find_path(GSL_INCLUDE_DIR "gsl/gsl_vector.h" PATHS ${GSLROOT}/include)
	if (NOT GSL_INCLUDE_DIR)
		message(FATAL_ERROR "gsl/gsl_vector.h can not be found in ${GSLROOT}/include")
	endif (NOT GSL_INCLUDE_DIR)
	set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})

	# find the library
	find_library(GSL_CORE_PATH gsl PATHS ${GSLROOT}/lib)
	if (NOT GSL_CORE_PATH)
		message(FATAL_ERROR "libgsl.so can not be found in ${GSLROOT}/lib")
	else (NOT GSL_CORE_PATH)
		set(GSL_LIBRARIES ${GSL_CORE_PATH}) 
	endif(NOT GSL_CORE_PATH)

else(GSLROOT)

	# Find include dir
	find_path(GSL_INCLUDE_DIR "gsl/gsl_vector.h")
	if (NOT GSL_INCLUDE_DIR)
		message(FATAL_ERROR "gsl/gsl_vector.h can not be found, please check whether you installed gnu gsl")
	else (NOT GSL_INCLUDE_DIR)
		set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
	endif (NOT GSL_INCLUDE_DIR)

	# now the GSLROOT is not defined
	# so we assume that the library is in the LD_LIBRARY_PATH
	# we will check it now
	# the head file is assumed to autometically loaded
	find_library(GSL_CORE_PATH gsl PATHS ENV LD_LIBRARY_PATH)
	if (NOT GSL_CORE_PATH)
		message("we need to search libgsl.so in ENV LD_LIBRARY_PATH since you did not define $GSLROOT")
		message("however, it's not found in LD_LIBRARY_PATH. Please check whether you installed GSL")
		message(FATAL_ERROR "libgsl.so can not be found, please install gsl library")
	else (NOT GSL_CORE_PATH)
		set(GSL_LIBRARIES ${GSL_CORE_PATH}) 
	endif(NOT GSL_CORE_PATH)
endif(GSLROOT)

# debugging print
#message(${GSL_LIBRARIES})

# finally reset the status of GSL
set(GSL_FOUND TRUE)

