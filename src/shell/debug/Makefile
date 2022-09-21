#
#
# makefile for shell part
#
#
#
NAME                 = test
CC		   	         = g++
CFLAGS       			= -Wall -g -fopenmp 
FC                   = gfortran 
FFLAGS               = -cpp -g  
INCLUDE              = -I../../general/include   -I../../math/include  \
							  -I../include -I../../geom/include -I${TBBROOT}/include 
MACRO                = -DDEBUG -DUSE_GENERAL_BLAS #-DWITH_SINGLE_PRECISION 
LIBRARY              = -L${TBBROOT} -ltbb -ltbbmalloc  -L/usr/lib/lapack \
							  -L/usr/lib/libblas -llapack -lblas -lm -lgomp \
							  -lboost_thread -lboost_filesystem -lboost_system
OBJC                 = parseshell.o shellprop.o basissetnames.o shell.o purecart.o shellsize.o \
							  excep.o textread.o element.o molecule.o geomutil.o zmatrix.o \
							  cptrans.o main.o globalinfor.o parameterparsing.o filerw.o \
							  blas.o blas1.o lapack.o matrix.o  in_place_transpose.o \
							  blockmatrix.o cptrans_test.o cptrans_para_test.o

ALL: $(OBJC) 
	$(CC) -o $(NAME) $(OBJC) $(LIBRARY)

main.o: main.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# shell module
parseshell.o: ../src/parseshell.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

shellprop.o: ../src/shellprop.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

shell.o: ../src/shell.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

purecart.o: ../src/purecart.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

cptrans.o: ../src/cptrans.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

shellsize.o: ../src/shellsize.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

basissetnames.o: ../src/basissetnames.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

cptrans_test.o: cptrans_test.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

cptrans_para_test.o: cptrans_para_test.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# general module
excep.o: ../../general/src/excep.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

textread.o: ../../general/src/textread.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

filerw.o: ../../general/src/filerw.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

globalinfor.o: ../../general/src/globalinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

parameterparsing.o: ../../general/src/parameterparsing.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# molecule module
element.o: ../../geom/src/element.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

molecule.o: ../../geom/src/molecule.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

zmatrix.o:  ../../geom/src/zmatrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

geomutil.o: ../../geom/src/geomutil.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# math module
blas.o: ../../math/src/blas.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

blas1.o: ../../math/src/blas1.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

lapack.o: ../../math/src/lapack.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

matrix.o: ../../math/src/matrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

blockmatrix.o: ../../math/src/blockmatrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

blockmatrixlist.o: ../../math/src/blockmatrixlist.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

in_place_transpose.o: ../../math/src/in_place_transpose.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@


.PHONY: clean
clean:
	@rm -rf *.o
	@rm -rf test

