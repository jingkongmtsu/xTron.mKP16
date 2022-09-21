#
# set the C++ compiler options
#

# reset the string
set(CMAKE_CXX_FLAGS_DEBUG   " ")
set(CMAKE_CXX_FLAGS_RELEASE " ")

# set the value
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

	#######################
	#  set the debug mode #
	#######################

	# open the warning sign
	if (COMPILER_WARNING_OFF  STREQUAL "1")
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -w")
	else()
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
	endif()

	# adding in the -g choice
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

	# optimizations
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

	# adding openmp support
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp")

	########################
	# set the release mode #
	########################

	# turn on warning
	if (COMPILER_WARNING_OFF  STREQUAL "1")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w")
	else()
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
	endif()


	# optimizations
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

	# adding openmp support
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

	#######################
	#  set the debug mode #
	#######################

	# open the warning sign
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -w0")

	# adding in the -g choice
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

	# optimizations
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

	# adding openmp support
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -qopenmp")

	########################
	# set the release mode #
	########################

	# open the warning sign
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w0")

	# optimization
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

	# adding openmp support
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -qopenmp")

else()
	message(FATAL_ERROR "Unrecognized CXX compiler ID: ${CMAKE_CXX_COMPILER_ID}")
endif()

