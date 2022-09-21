#
#  build the files in result folder
#
set(SRC_RESULT_CXX 
	${PROJECT_SOURCE_DIR}/src/result/src/atomdenmtrx.cpp
	${PROJECT_SOURCE_DIR}/src/result/src/denmtrx2.cpp
	${PROJECT_SOURCE_DIR}/src/result/src/denmtrx.cpp
	${PROJECT_SOURCE_DIR}/src/result/src/fracspininfor.cpp
	${PROJECT_SOURCE_DIR}/src/result/src/mo.cpp
	)
add_library(result ${SRC_RESULT_CXX})

