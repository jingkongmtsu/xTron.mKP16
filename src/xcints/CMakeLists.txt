#
#  build the files in xcints folder
#
set(SRC_XCINTS_CXX 
	${PROJECT_SOURCE_DIR}/src/xcints/src/xcintsinfor.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/atomgrids.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/batchbasis.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/batchfunc.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/batchgrid.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/batchvar.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/batchxcmtrx.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/batchjkmtrx.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/dftbasisderiv1.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/dftbasisderiv2.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/dftbasisderiv3.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/dftbasisderiv4.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/dftmatrix.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/molegrids.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/oddelec.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/sigatombasis.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/xcints.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/xcvar.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/xcvarinfor.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/hirshfeld.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/espints.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/halfjkrho.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/xcenergyinfor.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/sphereInt.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/integralF.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/gaussian.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/coefCalculator.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/hfxhole.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/xcintsprops.cpp
	${PROJECT_SOURCE_DIR}/src/xcints/src/xcintsoper.cpp
	)
add_library(xcints ${SRC_XCINTS_CXX})
 
set(SRC_XCINTS_FOR
	${PROJECT_SOURCE_DIR}/src/xcints/src/lebedev.f
	)
add_library(xcints_f ${SRC_XCINTS_FOR})
 
 
