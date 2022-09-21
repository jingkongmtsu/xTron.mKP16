# - Find Intel MKL
# Find the MKL libraries
# here for MKL library, we always use the dynamic library 
# as well as the multi-threading version. The only thing
# need to be checked here is the integer size (64 or 32?)
#
# limitations:
# 1  now it's only for unix;
# 2  only support emt64/ia32 two cases;
# 3  only support mkl version 11
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.

# firstly assume that MKL is not found
set(MKL_FOUND FALSE)

# check that whether the MKLROOT is set
if (NOT MKLROOT)
	message(${MKLROOT})
	message(FATAL_ERROR "MKLROOT is not defined in findMKL.cmake")
endif(NOT MKLROOT)

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h PATHS ${MKLROOT}/include)
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})

# check whether we use 64 bit integer or 32 bit integer in MKL?
set(MKL_INTEGER_64 0)
if (CMAKE_SIZEOF_INT STREQUAL "8")
	set(MKL_INTEGER_64 1)
endif(CMAKE_SIZEOF_INT STREQUAL "8")
if (USE_INTEGER_8)
	set(MKL_INTEGER_64 1)
endif(USE_INTEGER_8)

# get mkl lib path check the system type
set(MKL_86_64 0)
if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
	set(MKL_86_64 1)
endif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")

# set the searching path for lib
if(MKL_86_64)
	if(EXISTS ${MKLROOT}/lib/intel64)
		set(MKL_LIBRARY_PATH ${MKLROOT}/lib/intel64)
	else()
		message(FATAL_ERROR "MKL_LIBRARY_PATH can not be set for x86_64")
	endif()
else()
	if(EXISTS ${MKLROOT}/lib/ia32)
		set(MKL_LIBRARY_PATH ${MKLROOT}/lib/ia32)
	else()
		message(FATAL_ERROR "MKL_LIBRARY_PATH can not be set for ia32")
	endif()
endif()

# now find the library
# the linking of library is suggested from "link advisor" by mkl

# it's major version is 11
if(MKL_MAJOR_VERSION STREQUAL "11" OR MKL_MAJOR_VERSION STREQUAL "2017" OR MKL_MAJOR_VERSION STREQUAL "2018")  

	# mkl core
	find_library(MKL_CORE_PATH mkl_core  PATHS ${MKL_LIBRARY_PATH})
	if (NOT MKL_CORE_PATH)
		message(${MKL_CORE_PATH})
		message(FATAL_ERROR "MKL core library can not be found")
	endif()

	# integer interface lib
	if (MKL_86_64)
		if (MKL_INTEGER_64)
			find_library(MKL_INT_INTERFACE_PATH mkl_intel_ilp64 PATHS ${MKL_LIBRARY_PATH})
		else()
			find_library(MKL_INT_INTERFACE_PATH mkl_intel_lp64 PATHS ${MKL_LIBRARY_PATH})
		endif()
	else()
		find_library(MKL_INT_INTERFACE_PATH mkl_intel PATHS ${MKL_LIBRARY_PATH})
	endif()
	if (NOT MKL_INT_INTERFACE_PATH)
		message(${MKL_INT_INTERFACE_PATH})
		message(FATAL_ERROR "MKL integer interface library can not be found")
	endif()

	# threading
	find_library(MKL_THREAD_PATH mkl_intel_thread PATHS ${MKL_LIBRARY_PATH})
	if (NOT MKL_THREAD_PATH)
		message(${MKL_THREAD_PATH})
		message(FATAL_ERROR "MKL threads library can not be found")
	endif()

	# sequential library
	find_library(MKL_SEQ_PATH mkl_sequential PATHS ${MKL_LIBRARY_PATH})
	if (NOT MKL_SEQ_PATH)
		message(${MKL_SEQ_PATH})
		message(FATAL_ERROR "MKL sequential library can not be found")
	endif()

	# openmp
	# this one is intel's library
	# it's not in mkl's directory
	# so we search across the LD_LIBRARY_PATH
	find_library(MKL_OPENMP_PATH NAMES iomp5 PATHS ENV LD_LIBRARY_PATH) 
	if (NOT MKL_OPENMP_PATH)
		message(${MKL_OPENMP_PATH})
		message(FATAL_ERROR "MKL openmp library can not be found")
	endif()

	#additionaly libpthread 
	#search across the LD_LIBRARY_PATH
	find_library(MKL_GEN_PATH1 pthread PATHS ENV LD_LIBRARY_PATH)
	if (NOT MKL_GEN_PATH1)
		message(${MKL_GEN_PATH1})
		message(FATAL_ERROR "pthread library for MKL can not be found")
	endif()

	#additionaly libm
	#search across the LD_LIBRARY_PATH
	find_library(MKL_GEN_PATH2 m PATHS ENV LD_LIBRARY_PATH)
	if (NOT MKL_GEN_PATH2)
		message(${MKL_GEN_PATH2})
		message(FATAL_ERROR "libm library for MKL can not be found")
	endif()

	# if not use intel compiler, we need to add something -dl
	# this is called Dynamically Loaded (DL) Libraries
	#search across the LD_LIBRARY_PATH
	if(NOT (CMAKE_CXX_COMPILER_ID STREQUAL "Intel"))
		find_library(MKL_DL_PATH dl PATHS ENV LD_LIBRARY_PATH)
		if (NOT MKL_DL_PATH)
			message(${MKL_DL_PATH})
			message(FATAL_ERROR "libdl library for MKL can not be found")
		endif()
	endif()

	# now add in all of these
	if (USE_SEQUENTIAL_MKL_LIB STREQUAL "0")
		if (MKL_DL_PATH)
			set(MKL_LIBRARIES ${MKL_INT_INTERFACE_PATH} ${MKL_CORE_PATH} ${MKL_THREAD_PATH} 
				${MKL_OPENMP_PATH} ${MKL_GEN_PATH1} ${MKL_GEN_PATH2} ${MKL_DL_PATH})
		else()
			set(MKL_LIBRARIES ${MKL_INT_INTERFACE_PATH} ${MKL_CORE_PATH} ${MKL_THREAD_PATH} 
				${MKL_OPENMP_PATH} ${MKL_GEN_PATH1} ${MKL_GEN_PATH2})
		endif()
	else()
		if (MKL_DL_PATH)
			set(MKL_LIBRARIES ${MKL_INT_INTERFACE_PATH} ${MKL_SEQ_PATH} ${MKL_CORE_PATH} 
				${MKL_GEN_PATH1} ${MKL_GEN_PATH2} ${MKL_DL_PATH})
		else()
			set(MKL_LIBRARIES ${MKL_INT_INTERFACE_PATH} ${MKL_SEQ_PATH} ${MKL_CORE_PATH} 
				${MKL_GEN_PATH1} ${MKL_GEN_PATH2})
		endif()
	endif()
	#message(${MKL_LIBRARIES})

	# finally reset the status of MKL
	set(MKL_FOUND TRUE)
else()
	message(FATAL_ERROR "MKL version is not supported!")
endif()

