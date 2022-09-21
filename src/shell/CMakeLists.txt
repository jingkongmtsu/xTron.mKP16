#
#  build the files in shell folder
#
set(SRC_SHELL_CXX 
	${PROJECT_SOURCE_DIR}/src/shell/src/basissetnames.cpp
	${PROJECT_SOURCE_DIR}/src/shell/src/cptrans.cpp
	${PROJECT_SOURCE_DIR}/src/shell/src/parseshell.cpp
	${PROJECT_SOURCE_DIR}/src/shell/src/purecart.cpp
	${PROJECT_SOURCE_DIR}/src/shell/src/shell.cpp
	${PROJECT_SOURCE_DIR}/src/shell/src/shellsize.cpp
	${PROJECT_SOURCE_DIR}/src/shell/src/shellprop.cpp
	)
add_library(shell ${SRC_SHELL_CXX})

