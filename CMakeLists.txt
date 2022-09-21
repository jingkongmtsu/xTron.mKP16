################################################
#       heading for the top cmake file         #
################################################
cmake_minimum_required(VERSION 2.8)
project(phoenix CXX C Fortran)
set(BIN_NAME "xtron.exe")

################################################
#           check the type                     #
################################################
INCLUDE (CheckTypeSize)
check_type_size("int"    CMAKE_SIZEOF_INT     BUILTIN_TYPES_ONLY)
message("integer variable is in size of: ${CMAKE_SIZEOF_INT}")
check_type_size("double" CMAKE_SIZEOF_DOUBLE  BUILTIN_TYPES_ONLY)
message("double  variable is in size of: ${CMAKE_SIZEOF_DOUBLE}")
check_type_size("float"  CMAKE_SIZEOF_FLOAT   BUILTIN_TYPES_ONLY)
message("float   variable is in size of: ${CMAKE_SIZEOF_FLOAT}")

################################################
#        set the compiler flags                #
################################################
include("${PROJECT_SOURCE_DIR}/cmake/set_cxx_compiler_flags.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/set_fortran_compiler_flags.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/set_macros.cmake")

################################################
#            set the library                   #
################################################
include("${PROJECT_SOURCE_DIR}/cmake/set_library.cmake")
#The following line was suggested by Lori. For possible future use
#with pybind11.
#find_package(pybind11 2.2.3 EXACT REQUIRED CONFIG)

################################################
#           set the include dir                #
################################################

# first level modules
include_directories ("${PROJECT_SOURCE_DIR}/src/general/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/geom/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/math/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/shell/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/xcfunc/include")

# second level modules
include_directories ("${PROJECT_SOURCE_DIR}/src/result/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/gints/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/xcints/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/vdw/include")

# third level modules
include_directories ("${PROJECT_SOURCE_DIR}/src/scf/include")
include_directories ("${PROJECT_SOURCE_DIR}/src/scf_deriv/include")

# these are the integral engines
include_directories ("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy")
include_directories ("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/first_deriv")
include_directories ("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/second_deriv")

################################################
#           add in each subdir                 #
################################################

# set the list of my libraries for future target linking
# because all of these modules will be linked into *.a form,
# and for linker it will has a sequence for library searching
# therefore we need to put the lower level modules at last,
# and higher level modules to the first
#
# right now the targets are divided into four groups:
#
# 1 MY_LIBRARIES_CXX: it contains the c++ code compiled from each sub directory;
#
# 2 MY_LIBRARIES_FOR: it contains the Fortran code compiled from each sub directory;
#
# 3 MY_EXT_LIBRARIES: libraries imported from outside so that we do not compile here.
# one thing to note that the codes in the MY_EXT_LIBRARIES are not ordinary external
# libraries like boost, mkl etc. The codes are actually contained inside along the package
# and need to be pre-compiled.
#
# 4 MY_LIBRARIES_PHI: codes compiled for running on Intel Phi co-processor
# if we do not use the Phi co-processor, these codes are same to the MY_LIBRARIES_CXX
#
# one thing to note, that the integral libraries is not in the MY_EXT_LIBRARIES although
# it may be imported outside. However, to maintain a correct linking we need to put it
# inside the MY_LIBRARIES_CXX.
#

# set the flag
# we may not have the external libraries
set(HAS_EXT_LIBRARIES OFF)

# the first level of modules
# they are most foundamental ones
add_subdirectory("${PROJECT_SOURCE_DIR}/src/general")
set (MY_LIBRARIES_CXX "libgen")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/geom")
set (MY_LIBRARIES_CXX "geom"      ${MY_LIBRARIES_CXX})
add_subdirectory("${PROJECT_SOURCE_DIR}/src/math")
set (MY_LIBRARIES_CXX "math"      ${MY_LIBRARIES_CXX})
set (MY_LIBRARIES_FOR "math_f")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/shell")
set (MY_LIBRARIES_CXX "shell"     ${MY_LIBRARIES_CXX})
add_subdirectory("${PROJECT_SOURCE_DIR}/src/xcfunc")
set (MY_LIBRARIES_FOR "xcfunc"    ${MY_LIBRARIES_FOR})

# add in the integral engines
# the energy one is added by default
# the derivatives will be added depending on user
#
# on the other hand, such integral library may already 
# available for linking, so we just bypass these compilation
#
if (WITH_EXTERNAL_GINTS_ENGINE_LIB STREQUAL "0")

	# now it's the energy part
	add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/esp")
	set (MY_LIBRARIES_CXX "gints_engine_esp_d0"            ${MY_LIBRARIES_CXX})
	add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p")
	set (MY_LIBRARIES_CXX "gints_engine_mom_p_d0"          ${MY_LIBRARIES_CXX})
	add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/twobodyoverlap")
	set (MY_LIBRARIES_CXX "gints_engine_twobodyoverlap_d0" ${MY_LIBRARIES_CXX})
	add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/nai")
	set (MY_LIBRARIES_CXX "gints_engine_nai_d0"            ${MY_LIBRARIES_CXX})
	add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/kinetic")
	set (MY_LIBRARIES_CXX "gints_engine_kinetic_d0"        ${MY_LIBRARIES_CXX})
	add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/eri")
	set (MY_LIBRARIES_CXX "gints_engine_eri_d0"            ${MY_LIBRARIES_CXX})

	# let's see the gints derivatives
	if(WITH_GINTS_DERIV STREQUAL "1")
		add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/first_deriv/eri")
		set (MY_LIBRARIES_CXX "gints_engine_eri_d1" ${MY_LIBRARIES_CXX})
		add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/second_deriv/eri")
		set (MY_LIBRARIES_CXX "gints_engine_eri_d2" ${MY_LIBRARIES_CXX})
	endif(WITH_GINTS_DERIV STREQUAL "1")

else (WITH_EXTERNAL_GINTS_ENGINE_LIB STREQUAL "0")

	# now we are going to use external lib
	# check the input path
	set(GINTS_LIB_PATH $ENV{GINTS_ENGINE_LIB_PATH})
	if(NOT GINTS_LIB_PATH)
		message(${GINTS_LIB_PATH})
		message(FATAL_ERROR "GINTS_LIB_PATH is not defined so we can not use external gints_engine libraries")
	endif(NOT GINTS_LIB_PATH)

	# now we just try to add the library
	# energy part first
	set (GINTS_OPERATOR_NAMES "esp" "mom_p" "twobodyoverlap" "nai" "kinetic" "eri")
	FOREACH(GINTS_OPER ${GINTS_OPERATOR_NAMES})
		set (GINTS_LIB_NAME "libgints_engine_${GINTS_OPER}_d0")
		ADD_LIBRARY (${GINTS_LIB_NAME} STATIC IMPORTED)
		set (GINTS_LIB_NAME_PATH "${GINTS_LIB_PATH}/${GINTS_LIB_NAME}.a")
		SET_PROPERTY(TARGET ${GINTS_LIB_NAME} PROPERTY IMPORTED_LOCATION ${GINTS_LIB_NAME_PATH})
		set (MY_LIBRARIES_CXX "${GINTS_LIB_NAME}" ${MY_LIBRARIES_CXX})
	ENDFOREACH(GINTS_OPER)

	# now it's the derivatives
	if(WITH_GINTS_DERIV STREQUAL "1")
		set (GINTS_DERIV_NAMES "d1" "d2")
		FOREACH(GINTS_DERIV ${GINTS_DERIV_NAMES})
			set (GINTS_LIB_NAME "libgints_engine_eri_${GINTS_DERIV}")
			ADD_LIBRARY (${GINTS_LIB_NAME} STATIC IMPORTED)
			set (GINTS_LIB_NAME_PATH "${GINTS_LIB_PATH}/${GINTS_LIB_NAME}.a")
			SET_PROPERTY(TARGET ${GINTS_LIB_NAME} PROPERTY IMPORTED_LOCATION ${GINTS_LIB_NAME_PATH})
			set (MY_LIBRARIES_CXX "${GINTS_LIB_NAME}" ${MY_LIBRARIES_CXX})
		ENDFOREACH(GINTS_DERIV)
	endif(WITH_GINTS_DERIV STREQUAL "1")

endif (WITH_EXTERNAL_GINTS_ENGINE_LIB STREQUAL "0")

# the second level of modules
# they are depending on the first level modules
add_subdirectory("${PROJECT_SOURCE_DIR}/src/gints")
set (MY_LIBRARIES_CXX "gints"     ${MY_LIBRARIES_CXX})
add_subdirectory("${PROJECT_SOURCE_DIR}/src/xcints")
set (MY_LIBRARIES_CXX "xcints"     ${MY_LIBRARIES_CXX})
set (MY_LIBRARIES_FOR "xcints_f"   ${MY_LIBRARIES_FOR})
add_subdirectory("${PROJECT_SOURCE_DIR}/src/vdw")
set (MY_LIBRARIES_CXX "vdw"    ${MY_LIBRARIES_CXX})

# the third level of modules
# they are depending on the first and second level modules
add_subdirectory("${PROJECT_SOURCE_DIR}/src/scf")
set (MY_LIBRARIES_CXX "scf"  ${MY_LIBRARIES_CXX})
add_subdirectory("${PROJECT_SOURCE_DIR}/src/result")
set (MY_LIBRARIES_CXX "result" ${MY_LIBRARIES_CXX})
if(WITH_GINTS_DERIV STREQUAL "1")
	add_subdirectory("${PROJECT_SOURCE_DIR}/src/scf_deriv")
	set (MY_LIBRARIES_CXX "scf_deriv"  ${MY_LIBRARIES_CXX})
endif(WITH_GINTS_DERIV STREQUAL "1")

# additional external libraries

# GDM
if (WITH_GDM_LIB STREQUAL "1")

	# now we are going to use libgdm.a
	# check the input path
	set(GDM_LIB_PATH $ENV{GDM_LIB_PATH})
	if(NOT GDM_LIB_PATH)
		message(${GDM_LIB_PATH})
		message(FATAL_ERROR "GDM_LIB_PATH is not defined so we can not use GDM library")
	endif(NOT GDM_LIB_PATH)

	# now add library
	set(HAS_EXT_LIBRARIES ON)
	ADD_LIBRARY (libgdm STATIC IMPORTED)
	set (GDMLIB_PATH "${GDM_LIB_PATH}/libgdm.a")
	SET_PROPERTY(TARGET libgdm PROPERTY IMPORTED_LOCATION ${GDMLIB_PATH})
	set (MY_EXT_LIBRARIES ${MY_EXT_LIBRARIES} "libgdm")
endif (WITH_GDM_LIB STREQUAL "1")

# debug
# print out the libraries linked
message("MY_LIBRARIES_CXX: ${MY_LIBRARIES_CXX}")
message("MY_LIBRARIES_FOR: ${MY_LIBRARIES_FOR}")
if (${HAS_EXT_LIBRARIES})
	message("MY_EXT_LIBRARIES: ${MY_EXT_LIBRARIES}")
endif ()

# now add all of stuff together to form the final one
set (MY_WHOLE_LIBRARIES ${MY_LIBRARIES_CXX})
set (MY_WHOLE_LIBRARIES ${MY_WHOLE_LIBRARIES} ${MY_LIBRARIES_FOR})
if (${HAS_EXT_LIBRARIES})
	set (MY_WHOLE_LIBRARIES ${MY_WHOLE_LIBRARIES} ${MY_EXT_LIBRARIES})
endif()

################################################
#     now we need to add the target module     #
#     so that to create the binary             #
################################################

# add exe
if(TARGET_MODULE STREQUAL "DGINTS")
	add_executable(${BIN_NAME} 
		${PROJECT_SOURCE_DIR}/src/gints/debug/gints2d_test.cpp
		${PROJECT_SOURCE_DIR}/src/gints/debug/gints4d_test.cpp
		${PROJECT_SOURCE_DIR}/src/gints/debug/gints4deriv_test.cpp
		${PROJECT_SOURCE_DIR}/src/gints/debug/integralinfor_test.cpp
		${PROJECT_SOURCE_DIR}/src/gints/debug/shellpair_test.cpp
		${PROJECT_SOURCE_DIR}/src/gints/debug/main.cpp)
elseif(TARGET_MODULE STREQUAL "DXCINTS")
	add_executable(${BIN_NAME} 
		${PROJECT_SOURCE_DIR}/src/xcints/debug/espints_test.cpp
		${PROJECT_SOURCE_DIR}/src/xcints/debug/exrho_test.cpp
		${PROJECT_SOURCE_DIR}/src/xcints/debug/xcmtrx_test.cpp
		${PROJECT_SOURCE_DIR}/src/xcints/debug/main.cpp)
elseif(TARGET_MODULE STREQUAL "DRESULT")
	add_executable(${BIN_NAME} 
		${PROJECT_SOURCE_DIR}/src/result/debug/atomdenmtrx_test.cpp
		${PROJECT_SOURCE_DIR}/src/result/debug/coreguess_test.cpp
		${PROJECT_SOURCE_DIR}/src/result/debug/main.cpp)
elseif(TARGET_MODULE STREQUAL "DSCF")
	add_executable(${BIN_NAME} 
		${PROJECT_SOURCE_DIR}/src/scf/debug/scf_test.cpp  
		${PROJECT_SOURCE_DIR}/src/scf/debug/b05_test.cpp  
		${PROJECT_SOURCE_DIR}/src/scf/debug/main.cpp)
endif(TARGET_MODULE STREQUAL "DGINTS")

# add in release module
if(TARGET_MODULE STREQUAL "SCF")
	add_executable(${BIN_NAME} 
		${PROJECT_SOURCE_DIR}/src/main/scf/sp.cpp
		${PROJECT_SOURCE_DIR}/src/main/scf/frac.cpp
	   ${PROJECT_SOURCE_DIR}/src/main/scf/main.cpp)
endif(TARGET_MODULE STREQUAL "SCF")

# form the target linking library
target_link_libraries(${BIN_NAME} ${MY_WHOLE_LIBRARIES} ${THIS_PROJECT_LIBS})

# set the linking property
# we use the C++ as linker langurage
# also set the linker option according to the compiler
set_target_properties(${BIN_NAME} PROPERTIES LINKER_LANGUAGE CXX)

