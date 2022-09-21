#
#  build the files in scf_deriv folder
#
set(SRC_SCFDERIV_CXX 
	${PROJECT_SOURCE_DIR}/src/scf_deriv/src/scfderiv.cpp
	)
add_library(scf_deriv ${SRC_SCFDERIV_CXX})

