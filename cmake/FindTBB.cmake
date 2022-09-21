# - Find Intel TBB
# Find the TBB libraries
# here for TBB library, we always use the dynamic library 
# we get it from the LD_LIBRARY_PATH
#
# This module defines the following variables:
#
#   TBB_FOUND            : True if TBB_INCLUDE_DIR are found
#   TBB_INCLUDE_DIRS     : set when TBB_INCLUDE_DIR found
#   TBB_LIBRARIES        : the library to link against.

# firstly assume that TBB is not found
set(TBB_FOUND FALSE)

# the TBBROOT must be defined
# which is ensured by the perl configure file
set(TBBROOT $ENV{TBBROOT})
if(NOT TBBROOT)
	message(${TBBROOT})
	message(FATAL_ERROR "TBBROOT is not defined in findTBB.cmake")
endif(NOT TBBROOT)

# Find include dir
find_path(TBB_INCLUDE_DIR "tbb.h" PATHS ${TBBROOT}/include/tbb)
set(TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIR})

# TBB core
find_library(TBB_CORE_PATH tbb PATHS ENV LD_LIBRARY_PATH)
if (NOT TBB_CORE_PATH)
	message(${TBB_CORE_PATH})
	message(FATAL_ERROR "libtbb.so can not be found")
endif(NOT TBB_CORE_PATH)

# TBB malloc
find_library(TBB_MALLOC_PATH tbbmalloc PATHS ENV LD_LIBRARY_PATH)
if (NOT TBB_MALLOC_PATH)
	message(${TBB_MALLOC_PATH})
	message(FATAL_ERROR "libtbbmalloc.so can not be found")
endif(NOT TBB_MALLOC_PATH)

# now add in all of these
set(TBB_LIBRARIES ${TBB_CORE_PATH} ${TBB_MALLOC_PATH})
#message(${TBB_LIBRARIES})

# finally reset the status of TBB
set(TBB_FOUND TRUE)

