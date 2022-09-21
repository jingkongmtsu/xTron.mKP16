//
// see the comments in the README
// in the upper folder
//
// common C head files used in the program
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<complex>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<vector>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
using namespace std;
using namespace boost::math;
typedef int           Int;
typedef long double   Double;
typedef complex<Double> DComplex;

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
#define  PI                3.14159265358979323846E0L

//
// constant used to measure the difference between two Double value
// since the double accruracy is around 15-16, therefore we take
// its lower bound
//
#define  THRESHOLD_MATH  1.0E-20L

#include "angmomlisttest.h"

Double bionomialCoe(const Double& n, const Double& l);
Double getFactorial(const Double& n); 
Double doubleFac(const Double& n); 

Double getFactorial(const Double& n) {
	return factorial<Double>(n);
}

Double bionomialCoe(const Double& n, const Double& l) {
	return binomial_coefficient<Double>(n, l);
}

Double doubleFac(const Double& n) {
	if (fabs(n)<THRESHOLD_MATH || fabs(n+ONE)<THRESHOLD_MATH) {
		return ONE; //case for n = -1 and n = 0
	}else{
		return double_factorial<Double>(n);
	}
} 


Int main() 
{

	for (Int L=0; L<=30; L++) {
		for (Int absm=0; absm<=L; absm++) {
			Int nbas = (L+1)*(L+2)/2;
			for(Int index=0; index<nbas; index++) {
				Int lx = 0;
				Int ly = 0;
				Int lz = 0;
				getlmn(L,L,index,lx,ly,lz);

				// first method for C
				Double f1 = getFactorial(2*lx)*getFactorial(2*ly)*getFactorial(2*lz)*getFactorial(L);
				Double f2 = getFactorial(lx)*getFactorial(ly)*getFactorial(lz)*getFactorial(2*L);
				Double f3 = getFactorial(L-absm)/getFactorial(L+absm);
				Double C  = sqrt((f1/f2)*f3);	
				C        *= ONE/(getFactorial(L)*pow(TWO,L));

				// second method for C

				// 
				// sqrt((2*lx-1)!!*(2*ly-1)!!*(2*lz-1)!!/(2*L-1)!!)
				//
				Double c1 = ONE;
				for(Int count=0; count<L; count++) {
					Int lx0 = 2*(lx-count)-1;
					if (lx0<1) lx0 = 1;
					Int ly0 = 2*(ly-count)-1;
					if (ly0<1) ly0 = 1;
					Int lz0 = 2*(lz-count)-1;
					if (lz0<1) lz0 = 1;
					Int L0  = 2*(L-count)-1;
					if (L0<1)  L0  = 1;
					Double val = lx0*ly0*lz0;
					val = val/L0;
					c1 = c1*sqrt(val);
				}

				//
				// sqrt((L-|M|)!/(L+|M|)!)
				//
				Double c2 = ONE;
				for(Int count=1; count<=2*absm; count++) {
					Int t = L-absm+count;
					c2 = c2*sqrt(ONE/t);
				}

				//
				// 1/2^l*l!
				//
				Double c3 = ONE;
				for(Int count=1; count<=L; count++) {
					c3 = c3/(TWO*count);
				}

				Double result = c1*c2*c3;
				if (fabs(C-result)>THRESHOLD_MATH) {
					cout << C << " " << result << " difference " << fabs(C-result) << endl;
					cout << "for the lx " << lx << " ly " << ly << " lz " << lz << endl;
					cout << "absm " << absm << " L " << L << endl;
					cout << "result is not correct" << endl;
				}
			}

			// the above test is actually for small coefficients
			// now let's test big number
			// this value could be goes to (2L)! as highest 
			// we try to compute (2L)!/L!
			Double f = getFactorial(L+absm)/getFactorial(L);
			Double c = ONE;
			for(Int count=1; count<=absm; count++) {
				Int t = L+count;
				c = c*t;
			}

			// here the absolute error is meaningless
			// we will see the relative error for it
			// also make the error range larger
			Double error = fabs(f-c)/c;
			if (error>1.0E-16L) {
				cout << f << " " << c << " error " << error;
				cout << " absm " << absm << " L " << L;
				cout << " big number testing result is not correct" << endl;
			}
		}
	}
	return 0;
}
