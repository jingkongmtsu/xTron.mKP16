#ifndef INPLACETRANSPOSE_H
#define INPLACETRANSPOSE_H
#include "libgen.h"

// here we need to add a interface to link in the Fortran code to do
// in place transpose operation
#define  in_place_transpose  in_place_transpose_
extern "C" {
	extern void in_place_transpose(Double* A, Int* M, Int* N, Int* MN,
		  Int* MOVE, Int* IWRK, Int* IOK);	
}

#endif
