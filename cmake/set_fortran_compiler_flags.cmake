#
# set the Fortran compiler options
#
# we note that Fortran code is always related to the basic calculation,
# it will never runs with threading inside and we simply drop  the openmp
# keyword for Fortran code compilation
#

# reset it to empty
set(CMAKE_Fortran_FLAGS_DEBUG   " ")
set(CMAKE_Fortran_FLAGS_RELEASE " ")

# set the compilation flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

	#######################
	#  set the debug mode #
	#######################

	# open the warning sign
	if (COMPILER_WARNING_OFF  STREQUAL "1")
		set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -w")
	else()
		set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall")
	endif()

	# adding in the -g choice
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")

	# optimizations
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0")

	# add in the preprocessing
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp")

	########################
	# set the release mode #
	########################

	# warning
	if (COMPILER_WARNING_OFF  STREQUAL "1")
		set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w")
	else()
		set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -Wall")
	endif()

	# optimizations
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")

	# add in the preprocessing
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

	#######################
	#  set the debug mode #
	#######################

	# open the warning sign
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -w0")

	# adding in the -g choice
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")

	# pre-processing
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fpp")

	# disable optimization
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0")

	########################
	# set the release mode #
	########################

	# warning
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w0")

	# optimizations
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")

	# add in pre-processing
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fpp")

else()
	message(SEND_ERROR "Unrecognized compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
endif()

