#
#
# makefile for the math  module
# TBB library is required
# you must specify the TBBROOT env
#
#
#
NAME                 = test
CC		   	         = g++ 
FC                   = gfortran 
CFLAGS       			= -DDEBUG -Wall -g -fopenmp # -DWITH_SINGLE_PRECISION  
FFLAGS               = -cpp -g  #-DWITH_SINGLE_PRECISION
INCLUDE              = -I../../general/include   -I../include -I${TBBROOT}/include
OBJ                  = blas.o blas1.o lapack.o blas_test.o para_blas_test.o \
							  lapack_test.o para_lapack_test.o test_util.o \
							  matrix.o matrix_test.o para_matrix_test.o  in_place_transpose.o \
							  blockmatrix.o blockmatrix_test.o  blockmatrixlist.o blockmatrixlist_test.o \
							  excep.o main.o
TBB_LIB              = -ltbb -ltbbmalloc -lgfortran -lboost_system -lboost_thread        

# this is for blas and lapack libaray  
# you choose what do you want to include it
include lib_choice/intel.mk

ALL: ${OBJ}
	$(CC) -o $(NAME) ${OBJ} $(MATH_LIB) $(TBB_LIB)

excep.o: ../../general/src/excep.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

matrix.o: ../src/matrix.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

blockmatrix.o: ../src/blockmatrix.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

blockmatrixlist.o: ../src/blockmatrixlist.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

matrix_test.o: matrix_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

matrix_phi_test.o: matrix_phi_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

blockmatrix_test.o: blockmatrix_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

blockmatrixlist_test.o: blockmatrixlist_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

para_matrix_test.o: para_matrix_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

blas.o: ../src/blas.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

blas1.o: ../src/blas1.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

lapack.o: ../src/lapack.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

blas_test.o: blas_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

para_blas_test.o: para_blas_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

lapack_test.o: lapack_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

para_lapack_test.o: para_lapack_test.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

test_util.o: test_util.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

in_place_transpose.o: ../src/in_place_transpose.f
	$(FC) $(FFLAGS) $(MATH_MACRO) $(INCLUDE) -c  $<   -o $@

.PHONY: clean
clean:
	@rm -rf *.o
	@find . -name '*~' -exec rm -rf {} \;

