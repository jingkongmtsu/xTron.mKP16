#
#  build the files in geom folder
#
set(SRC_GEOM_CXX 
	${PROJECT_SOURCE_DIR}/src/geom/src/element.cpp 
	${PROJECT_SOURCE_DIR}/src/geom/src/molecule.cpp 
	${PROJECT_SOURCE_DIR}/src/geom/src/zmatrix.cpp 
	${PROJECT_SOURCE_DIR}/src/geom/src/geomutil.cpp 
	)
add_library(geom ${SRC_GEOM_CXX})

