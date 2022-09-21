#
#  build the files in scf folder
#
set(SRC_SCF_CXX 
	${PROJECT_SOURCE_DIR}/src/scf/src/fock.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfparam.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfmacro.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfconvsol.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfdiis.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfeadiis.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/diiscontroller.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfintscontroller.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scf.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfconv.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/wfn.cpp
	${PROJECT_SOURCE_DIR}/src/scf/src/scfenergyconv.cpp
	)
add_library(scf ${SRC_SCF_CXX})

