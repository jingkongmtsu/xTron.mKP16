#
#
# makefile for shell part
#
#
#
NAME                 = matrix_compare
CC		   	         = g++
CFLAGS       			= -w -g -O3 -fopenmp 
FC                   = gfortran 
FFLAGS               = -cpp -g  -O3 
SRC                  = ../../src
INCLUDE              = -I${SRC}/general/include   -I${SRC}/math/include  \
							  -I${SRC}/shell/include -I${SRC}/geom/include -I${TBBROOT}/include  \
							  -I${SRC}/result/include
MACRO                = -DUSE_GENERAL_BLAS #-DWITH_SINGLE_PRECISION 
LIBRARY              = -L${TBBROOT} -ltbb -ltbbmalloc  -L/usr/lib/lapack \
							  -L/usr/lib/libblas -llapack -lblas -lm -lgomp \
							  -lboost_thread -lboost_filesystem -lboost_system
OBJC                 = parseshell.o shellprop.o shell.o  \
							  excep.o textread.o element.o molecule.o geomutil.o zmatrix.o \
							  main.o globalinfor.o parameterparsing.o filerw.o \
							  blas.o blas1.o lapack.o matrix.o  in_place_transpose.o 

ALL: $(OBJC) 
	$(CC) -o $(NAME) $(OBJC) $(LIBRARY)

main.o: main.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# shell module
parseshell.o: ${SRC}/shell/src/parseshell.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

shellprop.o: ${SRC}/shell/src/shellprop.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

shell.o: ${SRC}/shell/src/shell.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# general module
excep.o: ${SRC}/general/src/excep.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

textread.o: ${SRC}/general/src/textread.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

filerw.o: ${SRC}/general/src/filerw.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

globalinfor.o: ${SRC}/general/src/globalinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

parameterparsing.o: ${SRC}/general/src/parameterparsing.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# molecule module
element.o: ${SRC}/geom/src/element.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

molecule.o: ${SRC}/geom/src/molecule.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

zmatrix.o:  ${SRC}/geom/src/zmatrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

geomutil.o: ${SRC}/geom/src/geomutil.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# math module
blas.o: ${SRC}/math/src/blas.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

blas1.o: ${SRC}/math/src/blas1.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

lapack.o: ${SRC}/math/src/lapack.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

matrix.o: ${SRC}/math/src/matrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

in_place_transpose.o: ${SRC}/math/src/in_place_transpose.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@


.PHONY: clean
clean:
	@rm -rf *.o

