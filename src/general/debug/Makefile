#
# makefile for testing general module
# the boost library is assumed to be in default location, like 
# /usr/lib etc.
#
NAME                 = test
CC		   	         = g++
FC                   = gfortran 
CFLAGS       			= -Wall  -g 
FFLAGS       			= -Wall  -g -cpp
INCLUDE              = -I../include -I../../math/include -I${TBBROOT}/include
MACRO                = -DDEBUG -DUSE_GENERAL_BLAS -DUSE_GENERAL_LAPACK 
OBJC                 = excep.o globalinfor.o xcfunc.o filerw.o textread.o \
							  parameterparsing.o localmemscr.o histdataman.o main.o \
                       blas.o blas1.o lapack.o matrix.o in_place_transpose.o 
LIB                  = -lboost_filesystem -lboost_system -lboost_thread -ltbb -ltbbmalloc \
							  -llapack -lblas -lm -lgomp -lgfortran

ALL: $(OBJC)
	$(CC) -o $(NAME) $(OBJC) $(LIB)

globalinfor.o: ../src/globalinfor.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

textread.o: ../src/textread.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

xcfunc.o: ../src/xcfunc.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

excep.o: ../src/excep.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

filerw.o: ../src/filerw.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

localmemscr.o: ../src/localmemscr.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

histdataman.o: ../src/histdataman.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

parameterparsing.o: ../src/parameterparsing.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

# additional math files for testing histdataman
blas.o: ../../math/src/blas.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

blas1.o: ../../math/src/blas1.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

lapack.o: ../../math/src/lapack.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

matrix.o: ../../math/src/matrix.cpp
	$(CC) $(CFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

in_place_transpose.o: ../../math/src/in_place_transpose.f
	$(FC) $(FFLAGS) $(MACRO) $(INCLUDE) -c  $<   -o $@

.PHONY: clean
clean:
	@rm *.o
	@rm test.out
	@rm test

