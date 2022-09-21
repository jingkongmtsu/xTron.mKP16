#
# set the library for the whole program
#
set(THIS_PROJECT_LIBS "")

#
# boost library
# the require version must start from 1.47
# currently, the following libraries in boost are needed:
# 1  filesystem and system;
#
set(Boost_USE_STATIC_LIBS     OFF)
set(Boost_USE_MULTITHREADED   OFF)
set(Boost_USE_STATIC_RUNTIME  OFF)
find_package(Boost 1.46.0 COMPONENTS thread filesystem system timer)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	set(THIS_PROJECT_LIBS ${Boost_LIBRARIES} ${THIS_PROJECT_LIBS})
	
	# output message
	message("boost libs: " ${Boost_LIBRARIES})
	message("boost include dir: " ${Boost_INCLUDE_DIRS})
	message("          ")
else()
	message(FATAL_ERROR "boost library not found")
endif()

#
# TBB library
#
# reset the module path
set(CMAKE_MODULE_PATH_SAVED ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# find tbb
find_package(TBB)
if (TBB_FOUND)
	include_directories(${TBB_INCLUDE_DIRS})
	set(THIS_PROJECT_LIBS ${TBB_LIBRARIES} ${THIS_PROJECT_LIBS})

	# output message
	message("tbb libs: " ${TBB_LIBRARIES})
	message("tbb include dir: " ${TBB_INCLUDE_DIRS})
	message("          ")
else()
	message(FATAL_ERROR "tbb library was not found")
endif()

# unset the module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_SAVED})
unset(CMAKE_MODULE_PATH_SAVED)

#
# now we need to find out the math library
# blas, lapack etc.
#
if(BLAS_VENDOR STREQUAL "netlib")

	# reset the module path
	set(CMAKE_MODULE_PATH_SAVED ${CMAKE_MODULE_PATH})
	set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

	# find the blas lib
	find_package(NETLIB)
	if (NETLIB_FOUND)
		set(THIS_PROJECT_LIBS ${NETLIB_LIBRARIES}  ${THIS_PROJECT_LIBS})
		message("netlib libs: " ${NETLIB_LIBRARIES})
		message("          ")
	else()
		message(FATAL_ERROR "BLAS and LAPACK library from netlib was not found")
	endif()

	# unset the module path
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_SAVED})
	unset(CMAKE_MODULE_PATH_SAVED)

elseif (BLAS_VENDOR STREQUAL "mkl")

	# reset the module path
	set(CMAKE_MODULE_PATH_SAVED ${CMAKE_MODULE_PATH})
	set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

	# find mkl
	find_package(MKL)
	if (MKL_FOUND)
		include_directories(${MKL_INCLUDE_DIRS})
		set(THIS_PROJECT_LIBS ${MKL_LIBRARIES}  ${THIS_PROJECT_LIBS})

		# output message
		message("mkl libs: " ${MKL_LIBRARIES})
		message("mkl include dir: " ${MKL_INCLUDE_DIRS})
		message("          ")
	else()
		message(FATAL_ERROR "mkl library was not found")
	endif()

	# unset the module path
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_SAVED})
	unset(CMAKE_MODULE_PATH_SAVED)

endif()


# GNU GSL library
# right now GSL library is not used anymore in SCF module

# reset the module path
set(CMAKE_MODULE_PATH_SAVED ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# find GSL
find_package(GSL)
if (GSL_FOUND)
	include_directories(${GSL_INCLUDE_DIRS})
	set(THIS_PROJECT_LIBS ${GSL_LIBRARIES} ${THIS_PROJECT_LIBS})

	# output message
	message("gnu gsl libs: " ${GSL_LIBRARIES})
	message("gnu gsl include dir: " ${GSL_INCLUDE_DIRS})
	message("          ")
else()
	message(FATAL_ERROR "gnu gsl library was not found")
endif()

# unset the module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_SAVED})
unset(CMAKE_MODULE_PATH_SAVED)


#
# if this is Intel compiler, we need to 
# remove the fortran library ifcore from 
# the implicit linking list
# 
# this is because ifcore is only for serial calculation
# on the multi-threading env it cause the SCF with
# DFT totally wrong
#
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
	list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES ifcore)
	list(APPEND      CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES ifcoremt)
endif()

