#
# this is for intel's library and macro definitions
# there are several things to remind:
# 1  different version's linking library is different
#    this is based on version 11.0. Please check the on line mkl liking advisor
# 2  do please choose the multi-threading library, so that we can enable openmp
# 3  do you use 32 bit integer or 64 bit integer? This is what you need to make
#    sure. Write a simple c++ problem and use sizeof to check the integer size 
# 4  you have to set the MKLROOT env ahead
# 
MATH_MACRO   = -DUSE_MKL -DOMP_DEBUG  
MATH_LIB     = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core \
					-lmkl_intel_thread -liomp5 -ldl -lpthread -lm  

