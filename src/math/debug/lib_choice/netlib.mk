# this is for netlib library and macro definitions
# here we do not use openmp, but we still need to 
# add openmp choice to complete linking
MATH_MACRO   = -DUSE_GENERAL_BLAS -DUSE_GENERAL_LAPACK 
MATH_LIB     = -L/usr/lib/lapack -L/usr/lib/libblas -llapack -lblas -lm -lgomp

