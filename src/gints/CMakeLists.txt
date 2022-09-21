#
#  build the files in gints folder
#
#  the normal code on host
#
set(SRC_GINTS_CXX 
	${PROJECT_SOURCE_DIR}/src/gints/src/digestutil.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/gints2d.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/gints4d.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/gints4deriv.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/gintsinfor.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/integralinfor.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/oneemtrx.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/oneintformula.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/shellpair.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/sigshellpairinfor.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/denmtrxinfor.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/symmjkdigest.cpp
	${PROJECT_SOURCE_DIR}/src/gints/src/atomfockmtrx.cpp
	)
add_library(gints ${SRC_GINTS_CXX})


