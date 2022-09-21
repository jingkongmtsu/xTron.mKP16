#
#  build the files in general folder
#  list all of cpp files
#
set(SRC_GENERAL_CXX 
	${PROJECT_SOURCE_DIR}/src/general/src/filerw.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/globalinfor.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/parameterparsing.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/textread.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/xcfunc.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/localmemscr.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/histdataman.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/excep.cpp 
	${PROJECT_SOURCE_DIR}/src/general/src/integraljobs.cpp
	)
add_library(libgen ${SRC_GENERAL_CXX})

