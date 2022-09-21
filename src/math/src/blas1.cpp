/**
 * \file blas1.cpp
 * \author fenglai liu 
 *
 * functions for BLAS1
 * ALL of functions in BLAS1 has increment of 1 
 */
#include <cstdio>
#include <cmath>
#include "excep.h"
#include "blas1.h"
using namespace excep;
using namespace blas;

// begin the offload attribute region
//#pragma offload_attribute (push, target (mic))

Double blas::vdot(const Double* x, const Double* y, UInt len) 
{

	// initilize
	UInt n = len%10;
	Double result = ZERO;

	// head part
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			result += x[i]*y[i];
		}
		if (len < 10) return result;
	}

	// body part
	for(UInt i=n; i<len; i = i+10) {
		Double r1 = x[i]*y[i]     + x[i+1]*y[i+1] + x[i+2]*y[i+2] + x[i+3]*y[i+3] + x[i+4]*y[i+4];
		Double r2 = x[i+5]*y[i+5] + x[i+6]*y[i+6] + x[i+7]*y[i+7] + x[i+8]*y[i+8] + x[i+9]*y[i+9];
		result += r1;
		result += r2;
	}
	return result;
}

void blas::vaxpy(const Double* x, Double* y,  Double a,  UInt len) 
{
	// initilize
	UInt n = len%10;

	// head part
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			y[i] = a*x[i]+y[i];
		}
		if (len < 10) return;
	}

	// real body
	for(UInt i=n; i<len; i = i+10) {
		y[i] = a*x[i]+y[i];
		y[i+1] = a*x[i+1]+y[i+1];
		y[i+2] = a*x[i+2]+y[i+2];
		y[i+3] = a*x[i+3]+y[i+3];
		y[i+4] = a*x[i+4]+y[i+4];
		y[i+5] = a*x[i+5]+y[i+5];
		y[i+6] = a*x[i+6]+y[i+6];
		y[i+7] = a*x[i+7]+y[i+7];
		y[i+8] = a*x[i+8]+y[i+8];
		y[i+9] = a*x[i+9]+y[i+9];
	}
}

Double blas::maxSearch(const Double* x, const Double* y,  UInt len) 
{
	Double maxVal = ZERO;
	for(UInt i=0; i<len; i++) {
		Double v = fabs(x[i]-y[i]);
		if (v>maxVal) maxVal = v;
	}
	return maxVal;
}

Double blas::maxSearch(const Double* x, UInt len) 
{
	Double maxVal = ZERO;
	for(UInt i=0; i<len; i++) {
		Double v = fabs(x[i]);
		if (v>maxVal) maxVal = v;
	}
	return maxVal;
}

Double blas::rmsd(const Double* A, UInt len)
{
	Double r = ZERO;
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			r += A[i]*A[i];
		}
		if (len < 10) {
			r = sqrt(r/len);
			return r;
		}
	}
	for(UInt i=n; i<len; i = i+10) {
			r += A[i]*A[i];
			r += A[i+1]*A[i+1];
			r += A[i+2]*A[i+2];
			r += A[i+3]*A[i+3];
			r += A[i+4]*A[i+4];
			r += A[i+5]*A[i+5];
			r += A[i+6]*A[i+6];
			r += A[i+7]*A[i+7];
			r += A[i+8]*A[i+8];
			r += A[i+9]*A[i+9];
	}
	r = sqrt(r/len);
	return r;
}

Double blas::rmsd(const Double* x, const Double* y, UInt len)
{
	Double r = ZERO;
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			Double r0 = (x[i]-y[i]);
			r += r0*r0;
		}
		if (len < 10) {
			r = sqrt(r/len);
			return r;
		}
	}
	Double r0 = ZERO;
	for(UInt i=n; i<len; i = i+10) {
		r0 = (x[i]-y[i]);
		r += r0*r0;
		r0 = (x[i+1]-y[i+1]);
		r += r0*r0;
		r0 = (x[i+2]-y[i+2]);
		r += r0*r0;
		r0 = (x[i+3]-y[i+3]);
		r += r0*r0;
		r0 = (x[i+4]-y[i+4]);
		r += r0*r0;
		r0 = (x[i+5]-y[i+5]);
		r += r0*r0;
		r0 = (x[i+6]-y[i+6]);
		r += r0*r0;
		r0 = (x[i+7]-y[i+7]);
		r += r0*r0;
		r0 = (x[i+8]-y[i+8]);
		r += r0*r0;
		r0 = (x[i+9]-y[i+9]);
		r += r0*r0;
	}
	r = sqrt(r/len);
	return r;
}

void blas::vcopy(const Double* x, Double* y,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			y[i] = x[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		y[i]   = x[i];
		y[i+1] = x[i+1];
		y[i+2] = x[i+2];
		y[i+3] = x[i+3];
		y[i+4] = x[i+4];
		y[i+5] = x[i+5];
		y[i+6] = x[i+6];
		y[i+7] = x[i+7];
		y[i+8] = x[i+8];
		y[i+9] = x[i+9];
	}
}

void blas::vabs(Double* x, UInt len)
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			x[i] = fabs(x[i]);
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		x[i]   = fabs(x[i]);
		x[i+1] = fabs(x[i+1]);
		x[i+2] = fabs(x[i+2]);
		x[i+3] = fabs(x[i+3]);
		x[i+4] = fabs(x[i+4]);
		x[i+5] = fabs(x[i+5]);
		x[i+6] = fabs(x[i+6]);
		x[i+7] = fabs(x[i+7]);
		x[i+8] = fabs(x[i+8]);
		x[i+9] = fabs(x[i+9]);
	}
}

void blas::vscal(Double* x,  Double a,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			x[i] = a*x[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		x[i]   = a*x[i];
		x[i+1] = a*x[i+1];
		x[i+2] = a*x[i+2];
		x[i+3] = a*x[i+3];
		x[i+4] = a*x[i+4];
		x[i+5] = a*x[i+5];
		x[i+6] = a*x[i+6];
		x[i+7] = a*x[i+7];
		x[i+8] = a*x[i+8];
		x[i+9] = a*x[i+9];
	}
}

void blas::vset(Double* x,  Double a, UInt len)
{
	UInt n = len%20;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			x[i] = a;
		}
		if (len < 20) return;
	}
	for(UInt i=n; i<len; i = i+20) {
		x[i   ] = a;
		x[i+1 ] = a;
		x[i+2 ] = a;
		x[i+3 ] = a;
		x[i+4 ] = a;
		x[i+5 ] = a;
		x[i+6 ] = a;
		x[i+7 ] = a;
		x[i+8 ] = a;
		x[i+9 ] = a;
		x[i+10] = a;
		x[i+11] = a;
		x[i+12] = a;
		x[i+13] = a;
		x[i+14] = a;
		x[i+15] = a;
		x[i+16] = a;
		x[i+17] = a;
		x[i+18] = a;
		x[i+19] = a;
	}
}	

void blas::vmul(const Double* x,  const Double* y, Double* z,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			z[i] = x[i]*y[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		z[i] = x[i]*y[i];
		z[i+1] = x[i+1]*y[i+1];
		z[i+2] = x[i+2]*y[i+2];
		z[i+3] = x[i+3]*y[i+3];
		z[i+4] = x[i+4]*y[i+4];
		z[i+5] = x[i+5]*y[i+5];
		z[i+6] = x[i+6]*y[i+6];
		z[i+7] = x[i+7]*y[i+7];
		z[i+8] = x[i+8]*y[i+8];
		z[i+9] = x[i+9]*y[i+9];
	}
}

void blas::vadd(const Double* x, const Double* y, Double* z,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			z[i] = x[i]+y[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		z[i] = x[i]+y[i];
		z[i+1] = x[i+1]+y[i+1];
		z[i+2] = x[i+2]+y[i+2];
		z[i+3] = x[i+3]+y[i+3];
		z[i+4] = x[i+4]+y[i+4];
		z[i+5] = x[i+5]+y[i+5];
		z[i+6] = x[i+6]+y[i+6];
		z[i+7] = x[i+7]+y[i+7];
		z[i+8] = x[i+8]+y[i+8];
		z[i+9] = x[i+9]+y[i+9];
	}
}

void blas::vsub(const Double* x,  const Double* y, Double* z,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			z[i] = x[i]-y[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		z[i] = x[i]-y[i];
		z[i+1] = x[i+1]-y[i+1];
		z[i+2] = x[i+2]-y[i+2];
		z[i+3] = x[i+3]-y[i+3];
		z[i+4] = x[i+4]-y[i+4];
		z[i+5] = x[i+5]-y[i+5];
		z[i+6] = x[i+6]-y[i+6];
		z[i+7] = x[i+7]-y[i+7];
		z[i+8] = x[i+8]-y[i+8];
		z[i+9] = x[i+9]-y[i+9];
	}
}

void blas::vinv(const Double* x, Double* y,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
#ifdef DEBUG
			if(fabs(x[i])<THRESHOLD_MATH) {
				string infor = "the x[i] value is too small to be divided";
				Excep excep("blas","vinv",EXCEPTION_DIVIDE_BY_ZERO,infor);
				handleExcep(excep);
			}
#endif
			y[i] = ONE/x[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
#ifdef DEBUG
		if(fabs(x[i])<THRESHOLD_MATH || 
				fabs(x[i+1])<THRESHOLD_MATH ||
				fabs(x[i+2])<THRESHOLD_MATH ||
				fabs(x[i+3])<THRESHOLD_MATH ||
				fabs(x[i+4])<THRESHOLD_MATH ||
				fabs(x[i+5])<THRESHOLD_MATH ||
				fabs(x[i+6])<THRESHOLD_MATH ||
				fabs(x[i+7])<THRESHOLD_MATH ||
				fabs(x[i+8])<THRESHOLD_MATH ||
				fabs(x[i+9])<THRESHOLD_MATH 
				) {
			string infor = "the x[i] value is too small to be divided";
			Excep excep("blas","vinv",EXCEPTION_DIVIDE_BY_ZERO,infor);
			handleExcep(excep);
		}
#endif
		y[i] = ONE/x[i];
		y[i+1] = ONE/x[i+1];
		y[i+2] = ONE/x[i+2];
		y[i+3] = ONE/x[i+3];
		y[i+4] = ONE/x[i+4];
		y[i+5] = ONE/x[i+5];
		y[i+6] = ONE/x[i+6];
		y[i+7] = ONE/x[i+7];
		y[i+8] = ONE/x[i+8];
		y[i+9] = ONE/x[i+9];
	}
}

void blas::vdiv( const Double* x,  const Double* y, Double* z,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
#ifdef DEBUG
			if(fabs(y[i])<THRESHOLD_MATH) {
				string infor = "the y[i] value is too small to be divided";
				Excep excep("blas","vdiv",EXCEPTION_DIVIDE_BY_ZERO,infor);
				handleExcep(excep);
			}
#endif
			z[i] = x[i]/y[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
#ifdef DEBUG
		if(fabs(y[i])<THRESHOLD_MATH || 
				fabs(y[i+1])<THRESHOLD_MATH ||
				fabs(y[i+2])<THRESHOLD_MATH ||
				fabs(y[i+3])<THRESHOLD_MATH ||
				fabs(y[i+4])<THRESHOLD_MATH ||
				fabs(y[i+5])<THRESHOLD_MATH ||
				fabs(y[i+6])<THRESHOLD_MATH ||
				fabs(y[i+7])<THRESHOLD_MATH ||
				fabs(y[i+8])<THRESHOLD_MATH ||
				fabs(y[i+9])<THRESHOLD_MATH 
		  ) {
			string infor = "the y[i] value is too small to be divided";
			Excep excep("blas","vdiv",EXCEPTION_DIVIDE_BY_ZERO,infor);
			handleExcep(excep);
		}
#endif
		z[i] = x[i]/y[i];
		z[i+1] = x[i+1]/y[i+1];
		z[i+2] = x[i+2]/y[i+2];
		z[i+3] = x[i+3]/y[i+3];
		z[i+4] = x[i+4]/y[i+4];
		z[i+5] = x[i+5]/y[i+5];
		z[i+6] = x[i+6]/y[i+6];
		z[i+7] = x[i+7]/y[i+7];
		z[i+8] = x[i+8]/y[i+8];
		z[i+9] = x[i+9]/y[i+9];
	}
}

void blas::v3add(const Double* x,  const Double* y,  const Double* z, 
		Double* w,  UInt len) 
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			w[i] = x[i]+y[i]+z[i]+w[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		w[i]   = x[i]+y[i]+z[i]+w[i];
		w[i+1] = x[i+1]+y[i+1]+z[i+1]+w[i+1];
		w[i+2] = x[i+2]+y[i+2]+z[i+2]+w[i+2];
		w[i+3] = x[i+3]+y[i+3]+z[i+3]+w[i+3];
		w[i+4] = x[i+4]+y[i+4]+z[i+4]+w[i+4];
		w[i+5] = x[i+5]+y[i+5]+z[i+5]+w[i+5];
		w[i+6] = x[i+6]+y[i+6]+z[i+6]+w[i+6];
		w[i+7] = x[i+7]+y[i+7]+z[i+7]+w[i+7];
		w[i+8] = x[i+8]+y[i+8]+z[i+8]+w[i+8];
		w[i+9] = x[i+9]+y[i+9]+z[i+9]+w[i+9];
	}
}

Double blas::vdot3( const Double* x,  const Double* y,  const Double* z, 
		UInt len)
{
	UInt n = len%10;
	Double result = ZERO;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			result += x[i]*y[i]*z[i];
		}
		if (len < 10) return result;
	}
	for(UInt i=n; i<len; i = i+10) {
		Double r= ZERO;
		r += x[i]*y[i]*z[i];
		r += x[i+1]*y[i+1]*z[i+1];
		r += x[i+2]*y[i+2]*z[i+2];
		r += x[i+3]*y[i+3]*z[i+3];
		r += x[i+4]*y[i+4]*z[i+4];
		r += x[i+5]*y[i+5]*z[i+5];
		r += x[i+6]*y[i+6]*z[i+6];
		r += x[i+7]*y[i+7]*z[i+7];
		r += x[i+8]*y[i+8]*z[i+8];
		r += x[i+9]*y[i+9]*z[i+9];
		result += r;
	}
	return result;
}

Double blas::vdot4(const Double* w, const Double* x, const Double* y, 
		const Double* z, UInt len)
{
	UInt n = len%10;
	Double result = ZERO;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			result += w[i]*x[i]*y[i]*z[i];
		}
		if (len < 10) return result;
	}
	for(UInt i=n; i<len; i = i+10) {
		Double r= ZERO;
		r += w[i]*x[i]*y[i]*z[i];
		r += w[i+1]*x[i+1]*y[i+1]*z[i+1];
		r += w[i+2]*x[i+2]*y[i+2]*z[i+2];
		r += w[i+3]*x[i+3]*y[i+3]*z[i+3];
		r += w[i+4]*x[i+4]*y[i+4]*z[i+4];
		r += w[i+5]*x[i+5]*y[i+5]*z[i+5];
		r += w[i+6]*x[i+6]*y[i+6]*z[i+6];
		r += w[i+7]*x[i+7]*y[i+7]*z[i+7];
		r += w[i+8]*x[i+8]*y[i+8]*z[i+8];
		r += w[i+9]*x[i+9]*y[i+9]*z[i+9];
		result += r;
	}
	return result;
}

void blas::vmuladd(const Double* x, const Double* y, Double* z,  UInt len)
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			z[i] += x[i]*y[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		z[i]   += x[i]*y[i];
		z[i+1] += x[i+1]*y[i+1];
		z[i+2] += x[i+2]*y[i+2];
		z[i+3] += x[i+3]*y[i+3];
		z[i+4] += x[i+4]*y[i+4];
		z[i+5] += x[i+5]*y[i+5];
		z[i+6] += x[i+6]*y[i+6];
		z[i+7] += x[i+7]*y[i+7];
		z[i+8] += x[i+8]*y[i+8];
		z[i+9] += x[i+9]*y[i+9];
	}
}

void blas::vaxyaddz( const Double* x,  const Double* y,  Double a, 
		Double* z,  UInt len)
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			z[i] += a*x[i]*y[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		z[i]   += a*x[i]*y[i];
		z[i+1] += a*x[i+1]*y[i+1];
		z[i+2] += a*x[i+2]*y[i+2];
		z[i+3] += a*x[i+3]*y[i+3];
		z[i+4] += a*x[i+4]*y[i+4];
		z[i+5] += a*x[i+5]*y[i+5];
		z[i+6] += a*x[i+6]*y[i+6];
		z[i+7] += a*x[i+7]*y[i+7];
		z[i+8] += a*x[i+8]*y[i+8];
		z[i+9] += a*x[i+9]*y[i+9];
	}
}

void blas::vmulsub(const Double* x, const Double* y, Double* z,  UInt len)
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			z[i] = x[i]*y[i]-z[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		z[i] = x[i]*y[i]-z[i];
		z[i+1] = x[i+1]*y[i+1]-z[i+1];
		z[i+2] = x[i+2]*y[i+2]-z[i+2];
		z[i+3] = x[i+3]*y[i+3]-z[i+3];
		z[i+4] = x[i+4]*y[i+4]-z[i+4];
		z[i+5] = x[i+5]*y[i+5]-z[i+5];
		z[i+6] = x[i+6]*y[i+6]-z[i+6];
		z[i+7] = x[i+7]*y[i+7]-z[i+7];
		z[i+8] = x[i+8]*y[i+8]-z[i+8];
		z[i+9] = x[i+9]*y[i+9]-z[i+9];
	}
}

void blas::vmul3add(const Double* x, const Double* y, const Double* z, Double* w, 
		UInt len)
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			w[i] = x[i]*y[i]*z[i]+w[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		w[i] = x[i]*y[i]*z[i]+w[i];
		w[i+1] = x[i+1]*y[i+1]*z[i+1]+w[i+1];
		w[i+2] = x[i+2]*y[i+2]*z[i+2]+w[i+2];
		w[i+3] = x[i+3]*y[i+3]*z[i+3]+w[i+3];
		w[i+4] = x[i+4]*y[i+4]*z[i+4]+w[i+4];
		w[i+5] = x[i+5]*y[i+5]*z[i+5]+w[i+5];
		w[i+6] = x[i+6]*y[i+6]*z[i+6]+w[i+6];
		w[i+7] = x[i+7]*y[i+7]*z[i+7]+w[i+7];
		w[i+8] = x[i+8]*y[i+8]*z[i+8]+w[i+8];
		w[i+9] = x[i+9]*y[i+9]*z[i+9]+w[i+9];
	}
}

void blas::v3axpy(const Double& c1, const Double& c2, const Double& c3,
		const Double* x, const Double* y, const Double* z, Double* w, UInt len)
{
	UInt n = len%10;
	if (n != 0) {
		for(UInt i=0; i<n; i++) {
			w[i] = c1*x[i]+c2*y[i]+c3*z[i]+w[i];
		}
		if (len < 10) return;
	}
	for(UInt i=n; i<len; i = i+10) {
		w[i] = c1*x[i]+c2*y[i]+c3*z[i]+w[i];
		w[i+1] = c1*x[i+1]+c2*y[i+1]+c3*z[i+1]+w[i+1];
		w[i+2] = c1*x[i+2]+c2*y[i+2]+c3*z[i+2]+w[i+2];
		w[i+3] = c1*x[i+3]+c2*y[i+3]+c3*z[i+3]+w[i+3];
		w[i+4] = c1*x[i+4]+c2*y[i+4]+c3*z[i+4]+w[i+4];
		w[i+5] = c1*x[i+5]+c2*y[i+5]+c3*z[i+5]+w[i+5];
		w[i+6] = c1*x[i+6]+c2*y[i+6]+c3*z[i+6]+w[i+6];
		w[i+7] = c1*x[i+7]+c2*y[i+7]+c3*z[i+7]+w[i+7];
		w[i+8] = c1*x[i+8]+c2*y[i+8]+c3*z[i+8]+w[i+8];
		w[i+9] = c1*x[i+9]+c2*y[i+9]+c3*z[i+9]+w[i+9];
	}
}

// end the offload attribute region
//#pragma offload_attribute (pop)

