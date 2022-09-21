#
#  build the files in vdw folder
#
set(SRC_VDW_CXX 
	${PROJECT_SOURCE_DIR}/src/vdw/src/vdwinfor.cpp
	${PROJECT_SOURCE_DIR}/src/vdw/src/xdm.cpp 
	)
add_library(vdw ${SRC_VDW_CXX})

