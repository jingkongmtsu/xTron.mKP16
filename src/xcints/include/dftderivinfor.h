/**
 * \file   derivinfor.h
 * \brief  general derivatives information 
 * \author Fenglai Liu and Jing Kong
 */

#ifndef DERIVINFOR_H
#define DERIVINFOR_H
#include "libgen.h"

// defining the number of derivatives information
const UInt N_DERIV_0 = 1; 
const UInt N_DERIV_1 = 3; 
const UInt N_DERIV_2 = 6; 
const UInt N_DERIV_3 = 10;
const UInt N_DERIV_4 = 15; 

// definition of the derivatives orders implicitly used
// in the program
const UInt XC_DERIV_X    = 1;
const UInt XC_DERIV_Y    = 2;
const UInt XC_DERIV_Z    = 3;
const UInt XC_DERIV_XX   = 4;
const UInt XC_DERIV_XY   = 5;
const UInt XC_DERIV_YY   = 6;
const UInt XC_DERIV_XZ   = 7;
const UInt XC_DERIV_YZ   = 8;
const UInt XC_DERIV_ZZ   = 9;
const UInt XC_DERIV_XXX  = 10;
const UInt XC_DERIV_XXY  = 11;
const UInt XC_DERIV_XYY  = 12;
const UInt XC_DERIV_YYY  = 13;
const UInt XC_DERIV_XXZ  = 14;
const UInt XC_DERIV_XYZ  = 15;
const UInt XC_DERIV_YYZ  = 16;
const UInt XC_DERIV_XZZ  = 17;
const UInt XC_DERIV_YZZ  = 18;
const UInt XC_DERIV_ZZZ  = 19;
const UInt XC_DERIV_XXXX = 20;
const UInt XC_DERIV_XXXY = 21;
const UInt XC_DERIV_XXYY = 22;
const UInt XC_DERIV_XYYY = 23;
const UInt XC_DERIV_YYYY = 24;
const UInt XC_DERIV_XXXZ = 25;
const UInt XC_DERIV_XXYZ = 26;
const UInt XC_DERIV_XYYZ = 27;
const UInt XC_DERIV_YYYZ = 28;
const UInt XC_DERIV_XXZZ = 29;
const UInt XC_DERIV_XYZZ = 30;
const UInt XC_DERIV_YYZZ = 31;
const UInt XC_DERIV_XZZZ = 32;
const UInt XC_DERIV_YZZZ = 33;
const UInt XC_DERIV_ZZZZ = 34;

// let's organize the above derivative data into const arrays
const UInt XC_DERIV_ORDER_1[ ] = { XC_DERIV_X, XC_DERIV_Y, XC_DERIV_Z};
const UInt XC_DERIV_ORDER_2[ ] = { XC_DERIV_XX, XC_DERIV_XY, XC_DERIV_YY, 
	XC_DERIV_XZ, XC_DERIV_YZ, XC_DERIV_ZZ};
const UInt XC_DERIV_ORDER_3[ ] = { XC_DERIV_XXX, XC_DERIV_XXY, XC_DERIV_XYY, 
	XC_DERIV_YYY, XC_DERIV_XXZ, XC_DERIV_XYZ, XC_DERIV_YYZ, XC_DERIV_XZZ,
	XC_DERIV_YZZ, XC_DERIV_ZZZ};
const UInt XC_DERIV_ORDER_4[ ] = { XC_DERIV_XXXX, XC_DERIV_XXXY, XC_DERIV_XXYY, 
	XC_DERIV_XYYY, XC_DERIV_YYYY, XC_DERIV_XXXZ, XC_DERIV_XXYZ, XC_DERIV_XYYZ,
	XC_DERIV_YYYZ, XC_DERIV_XXZZ, XC_DERIV_XYZZ, XC_DERIV_YYZZ, XC_DERIV_XZZZ,
	XC_DERIV_YZZZ, XC_DERIV_ZZZZ};

#endif

