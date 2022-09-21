#
#  build the files in math folder
#
set(SRC_MATH_CXX 
	"${PROJECT_SOURCE_DIR}/src/math/src/blas1.cpp" 
	"${PROJECT_SOURCE_DIR}/src/math/src/blas.cpp" 
	"${PROJECT_SOURCE_DIR}/src/math/src/lapack.cpp" 
	"${PROJECT_SOURCE_DIR}/src/math/src/matrix.cpp" 
	"${PROJECT_SOURCE_DIR}/src/math/src/blockmatrix.cpp" 
	"${PROJECT_SOURCE_DIR}/src/math/src/blockmatrixlist.cpp" 
	)
add_library(math ${SRC_MATH_CXX})

set(SRC_MATH_FOR 
	"${PROJECT_SOURCE_DIR}/src/math/src/in_place_transpose.f" 
	)
add_library(math_f ${SRC_MATH_FOR})

