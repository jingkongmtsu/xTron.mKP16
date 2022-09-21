#ifndef LIBGEN_H
#define LIBGEN_H

////////////////////////////////////////////////////////////////////////
//          General Head Files Needed for The Whole Program           //
////////////////////////////////////////////////////////////////////////

// common C head files used in the program
#include<cstdlib>
#include<cstdio>

// external math functions 
#include<cmath>

// common C++ head files used in the program
#include<complex>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<vector>
using namespace std;

/**
 * re-define the types here
 * we do not use unsigned type for int, since the common int
 * should be good enough
 */
typedef int           Int;
typedef long double   Double;
typedef complex<long double> DComplex;

////////////////////////////////////////////////////////////////////////
//              General Data Needed for The Whole Program             //
////////////////////////////////////////////////////////////////////////

// general constant used in the program
#define  MINUS_HALF   -0.5E0L
#define  MINUS_ONE    -1.0E0L
#define  MINUS_TWO    -2.0E0L
#define  MINUS_THREE  -3.0E0L
#define  MINUS_FOUR   -4.0E0L
#define  MINUS_FIVE   -5.0E0L
#define  MINUS_SIX    -6.0E0L
#define  MINUS_SEVEN  -7.0E0L
#define  MINUS_EIGHT  -8.0E0L
#define  MINUS_NINE   -9.0E0L
#define  MINUS_TEN    -10.0E0L
#define  ZERO          0.0E0L
#define  ONE           1.0E0L
#define  TWO           2.0E0L
#define  THREE         3.0E0L
#define  FOUR          4.0E0L
#define  FIVE          5.0E0L
#define  SIX           6.0E0L
#define  SEVEN         7.0E0L
#define  EIGHT         8.0E0L
#define  NINE          9.0E0L
#define  TEN           10.0E0L
#define  HALF          0.5E0L

//
// physical and math constants
//
#define  PI              3.14159265358979323846E0L
#define  THRESHOLD_MATH  1.0E-20L

/**
 * number of Cartesian type of basis set functions
 */
inline Int getCartBas(const Int& lmin, const Int& lmax) {
	return ((lmax+1)*(lmax+2)*(lmax+3)-lmin*(lmin+1)*(lmin+2))/6;
};

/**
 * number of spherical type of basis set functions
 */
inline Int getPureBas(const Int& lmin, const Int& lmax) {
	return (lmax+1)*(lmax+1)-lmin*lmin;
};

#endif

