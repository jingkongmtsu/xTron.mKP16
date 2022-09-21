/**
 * NOTE:
 * *******************
 * now this file is kept as archive
 * we do not use the GSL functions inside anymore
 * *******************
 *
 * this is the multi-dimensional minimization for SCF in terms of 
 * GSL library
 * \author fenglai liu 
 *
 * for the error functional expression, please see the EDIIS/ADIIS
 * paper, which has been given in the head of scfeadiis.h
 *
 * for the use of GSL multi-dimensional minimizer, we use the sample
 * code provided by GSL documents. Please see the GSL website:
 *
 * http://www.gnu.org/software/gsl/manual/html_node/
 *
 * the section of Multidimensional Minimization
 *
 * additional note:
 *
 * originally we use the gradient based algorithms to solve the 
 * minimization problem for ADIIS/EDIIS. However, this method
 * seems to have a deficiency that the optimized coefficients 
 * may have some zero elements. This may bring the interpolated 
 * Fock matrix back to previous iteration result, then accordingly;
 * the SCF iteration may get stuck in some iterations and 
 * unsuccessfully jump out.
 *
 * Therefore we switch our way to the direct minimization of error 
 * function without gradient algorithm. You can see the GSL documents
 * for more details about the method we choose. We comment out all 
 * of our current gradient methods.
 *
 * in the original codes, we tested our codes with idempotent condition
 * that the error functional is appended with another term:
 * \f$ f = f + \lambda*\lambda*(NE-\sum_{i}\sum_{j}c_{i}c_{j}(D_{i}SD_{j}|S)) \f$
 * 
 * where the new density matrix is expressed as:
 *
 * \f$ D = \sum_{i}c_{i} D_{i} \f$
 *
 * All of Di are density matrices in DIIS history. However, with forcing
 * idempotent condition to new D, the SCF convergence is very bad therefore
 * we finally turned off the idempotent condition requirement.
 */
#include <cmath>
#include <string>
#include <iostream>
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"
#include "libgen.h"
#include "excep.h"
#include "scfeadiis.h"
#include "gslscfminimizer.h"
using namespace excep;
using namespace scfeadiis;
using namespace gslscfminimizer;
using namespace std;

void gslscfminimizer::getCoefs(const gsl_vector *t, gsl_vector *c)
{
	//
	// to remove the condition that \sum_{i}c_{i} = 1
	// and c_{i} >= 0, we can re-define the c_{i} via
	// t_{i}: 
	// c_{i} = t_{i}^{2}/\sum_{k}t_{k}^{2}
	//
	// here we will calculate the ci according to the 
	// definition above
	//
	
	// get the sum of ti
	double sumTi = ZERO;
	UInt nTerms = t->size;
	for(UInt i=0; i<nTerms; i++) {
		double ti  = gsl_vector_get(t,i);
		sumTi += ti*ti;
	}

	// now do the calculation
	for(UInt i=0; i<nTerms; i++) {
		double ti  = gsl_vector_get(t,i);
		double v   = ti*ti/sumTi;
		gsl_vector_set(c,i,v);
	}
}

void gslscfminimizer::getGradCoefs(const gsl_vector *t, gsl_vector *gc)
{
	//
	// to remove the condition that \sum_{i}c_{i} = 1
	// and c_{i} >= 0, we can re-define the c_{i} via
	// t_{i}: 
	// c_{i} = t_{i}*t_{i}/\sum_{k}t_{k}^{2}
	//
	// here we will calculate the gradient of ci 
	// according to the definition above
	//
	// the result of gc is a n*n matrix.
	// for the derivatives, the index of ci is made 
	// in first dimension, the index of tj is made 
	// in second dimension. That is  to say,
	// \parital ci/ \partical tj is in row of i,
	// and in column of j
	//
	
	// get the sum of ti
	double sumTi = ZERO;
	UInt nTerms = t->size;
	for(UInt i=0; i<nTerms; i++) {
		double ti  = gsl_vector_get(t,i);
		sumTi += ti*ti;
	}

	// now do the calculation
	// i is refer to t index, and j is refer to result c index
	for(UInt i=0; i<nTerms; i++) {
		double ti  = gsl_vector_get(t,i);
		for(UInt j=0; j<nTerms; j++) {
			double tj  = gsl_vector_get(t,j);
			double gv  = ZERO;
			if (i == j) {
				// this is for \partial ci/\partial ti 
				gv  = TWO*ti/sumTi - TWO*ti*ti*ti/(sumTi*sumTi);
			}else{
				// this is for \partial cj/\partial ti when i != j
				gv  = -TWO*tj*tj*ti/(sumTi*sumTi);
			}

			// write the value
			UInt index = j+i*nTerms;
			gsl_vector_set(gc,index,gv);
		}
	}
}

double gslscfminimizer::ADIIS_f(const gsl_vector *t, void *params)
{
	// we compute the coefs according to the definition above
	UInt n = t->size;
	gsl_vector * c = gsl_vector_alloc(n);
	getCoefs(t,c);

	// the zeroth term
	double *p = (double *)params;
	UInt offset = 0;
	double f  = p[offset];

	// first order term
	offset += 1;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		double pi = p[offset+i];
		f = f + TWO*ci*pi;
	}

	// now it's the second order term
	offset += n;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		for(UInt j=0; j<n; j++) {
			double cj = gsl_vector_get(c,j);
			UInt index= j+i*n;
			double pji= p[offset+index];
			f = f + ci*cj*pji;
		}
	}

	// now it's the idempoetnt requirement term
	// this is abandoned
	/*
	offset += n*n;
	double v = ZERO;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		for(UInt j=0; j<n; j++) {
			double cj = gsl_vector_get(c,j);
			UInt index= j+i*n;
			double pji= p[offset+index];
			v = v + ci*cj*pji;
		}
	}
	double lambda = gsl_vector_get(t,n);
	f = f + lambda*lambda*(NE-v)*(NE-v);
	*/

	// free the space
	gsl_vector_free (c);

	// now return f
	return f;
}

/*
void gslscfminimizer::ADIIS_df(const gsl_vector *t, void *params, 
       gsl_vector *df)
{
	// we compute the coefs according to the definition above
	UInt n = t->size;
	gsl_vector * c  = gsl_vector_alloc(n);
	gsl_vector * gc = gsl_vector_alloc(n*n);
	getCoefs(t,c);
	getGradCoefs(t,gc);

	// set the parameters
	double *p = (double *)params;

	// the number of total electrons
	//double NE = p[0];

	// now it's first order term
	// j is loop over t index, and i loop over number of terms
	// so i is also the coefs index
	UInt offset = 2;
	for(UInt j=0; j<n; j++) {
		double v = ZERO;
		for(UInt i=0; i<n; i++) {
			UInt index  = i + j*n;
			double gcij = gsl_vector_get(gc,index);
			double pi   = p[offset+i];
			v += TWO*gcij*pi;
		}
		gsl_vector_set(df,j,v);
	}

	//
	// because the idempotent term also contains the second order
	// term. Therefore we need to do some pre-calculation
	// the offset contains:
	// 1   ne
	// 1   total energy
	// n   first order term
	// n*n second order term 
	double pref1 = ZERO;
	UInt offsetIdemMtrx = 2+n+n*n;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		for(UInt j=0; j<n; j++) {
			double cj = gsl_vector_get(c,j);
			UInt index= j+i*n;
			double pji= p[offsetIdemMtrx+index];
			pref1 = pref1 + ci*cj*pji;
		}
	}
	double lam  = gsl_vector_get(t,n);
	double pref = TWO*lam*lam*(NE-pref1);

	// second order term
	// i is loop over ci index, and j loop over cj
	// and k is loop over the t index
	offset = 2+n;
	for(UInt k=0; k<n; k++) {
		double v = ZERO;
		for(UInt j=0; j<n; j++) {
			double cj   = gsl_vector_get(c,j);
			UInt index1 = j+k*n;
			double gcjk = gsl_vector_get(gc,index1);
			for(UInt i=0; i<n; i++) {
				double ci   = gsl_vector_get(c,i);
				UInt index2 = i+k*n;
				double gcik = gsl_vector_get(gc,index2);

				// now it's the derivative term
				// for the second order term
				// in error functional
				UInt index3 = i+j*n;
				double pij = p[offset+index3];
				v += (gcik*cj+ci*gcjk)*pij;

				// this is the derivative term
				// for the idempotent part
				// simply this is (partial (NE-pref1)/ partial tk)
				//pij = p[offsetIdemMtrx+index3];
				//v -= pref*(gcik*cj+ci*gcjk)*pij;
			}
		}

		// write the value back
		double v0 = gsl_vector_get(df,k);
		v = v0 + v;
		gsl_vector_set(df,k,v);
	}

	// now it's derivatives for the lambda
	//double fx = TWO*lam*(NE-pref1)*(NE-pref1);
	//double fx = ZERO;
	//gsl_vector_set(df,n,fx);

	// free the space
	gsl_vector_free(c);
	gsl_vector_free(gc);
}

void gslscfminimizer::ADIIS_fdf(const gsl_vector * t, void * params, 
		double *f, gsl_vector *df)
{
	*f = ADIIS_f(t,params);
	ADIIS_df(t,params,df);
}

*/
double gslscfminimizer::EDIIS_f(const gsl_vector *t, void *params)
{
	// we compute the coefs according to the definition above
	UInt n = t->size;
	gsl_vector * c = gsl_vector_alloc(n);
	getCoefs(t,c);

	// the zeroth term
	// EDIIS does not have zeroth term
	// so we directly go to first order term
	double *p = (double *)params;
	UInt offset = 0;
	double f = ZERO;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		double pi = p[offset+i];
		f = f + ci*pi;
	}

	// now it's the second order term
	// we note that if i == j
	// the second term is 0(pij = 0)
	// thereofre we omit this term
	// 
	// on the other hand, in error 
	// function EDIIS has coefficient
	// of -1/2, and because we only
	// have half loop, we need to multiply 
	// 2 (because j=i+1 to n); therefore
	// the total is 1
	offset += n;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		for(UInt j=i+1; j<n; j++) {
			double cj = gsl_vector_get(c,j);
			UInt index= j+i*n;
			double pij= p[offset+index];
			f = f - ci*cj*pij;
		}
	}

	// now it's the idempoetnt requirement term
	// this is abandoned
	/*
	offset += n*n;
	double v = ZERO;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		for(UInt j=0; j<n; j++) {
			double cj = gsl_vector_get(c,j);
			UInt index= j+i*n;
			double pji= p[offset+index];
			v = v + ci*cj*pji;
		}
	}
	double lambda = gsl_vector_get(t,n);
	f = f + lambda*lambda*(NE-v)*(NE-v);
	*/

	// free the space
	gsl_vector_free (c);

	// now return f
	return f;
}

/*
void gslscfminimizer::EDIIS_df(const gsl_vector *t, void *params, 
       gsl_vector *df)
{
	// we compute the coefs according to the definition above
	UInt n = t->size;
	gsl_vector * c  = gsl_vector_alloc(n);
	gsl_vector * gc = gsl_vector_alloc(n*n);
	getCoefs(t,c);
	getGradCoefs(t,gc);

	// set the parameters
	double *p = (double *)params;

	// the number of total electrons
	//double NE = p[0];

	// now it's first order term
	// j is loop over t index, and i loop over number of terms
	// so i is also the coefs index
	UInt offset = 1;
	for(UInt j=0; j<n; j++) {
		double v = ZERO;
		for(UInt i=0; i<n; i++) {
			UInt index  = i + j*n;
			double gcij = gsl_vector_get(gc,index);
			double pi   = p[offset+i];
			v += gcij*pi;
		}
		gsl_vector_set(df,j,v);
	}

	//
	// because the idempotent term also contains the second order
	// term. Therefore we need to do some pre-calculation
	// the offset contains:
	// 1   ne
	// n   first order term
	// n*n second order term 
	double pref1 = ZERO;
	UInt offsetIdemMtrx = 1+n+n*n;
	for(UInt i=0; i<n; i++) {
		double ci = gsl_vector_get(c,i);
		for(UInt j=0; j<n; j++) {
			double cj = gsl_vector_get(c,j);
			UInt index= j+i*n;
			double pji= p[offsetIdemMtrx+index];
			pref1 = pref1 + ci*cj*pji;
		}
	}
	double lam  = gsl_vector_get(t,n);
	double pref = TWO*lam*lam*(NE-pref1);

	// second order term
	// i is loop over ci index, and j loop over cj
	// and k is loop over the t index
	//
	// here we note that even though the second order
	// term in EDIIS is symmetrical, but the idempotent
	// matrix may not be symmetrical. Therefore here
	// we loop over all of elements
	offset = 1+n;
	for(UInt k=0; k<n; k++) {
		double v = ZERO;
		for(UInt j=0; j<n; j++) {
			double cj   = gsl_vector_get(c,j);
			UInt index1 = j+k*n;
			double gcjk = gsl_vector_get(gc,index1);
			for(UInt i=0; i<n; i++) {
				double ci   = gsl_vector_get(c,i);
				UInt index2 = i+k*n;
				double gcik = gsl_vector_get(gc,index2);

				// now it's the derivative term
				// for the second order term
				// in error functional
				UInt index3 = i+j*n;
				double pij = p[offset+index3];
				v -= HALF*(gcik*cj+ci*gcjk)*pij;

				// this is the derivative term
				// for the idempotent part
				// simply this is (partial (NE-pref1)/ partial tk)
				//pij = p[offsetIdemMtrx+index3];
				//v -= pref*(gcik*cj+ci*gcjk)*pij;
			}
		}

		// write the value back
		double v0 = gsl_vector_get(df,k);
		v = v0 + v;
		gsl_vector_set(df,k,v);
	}

	// now it's derivatives for the lambda
	//double fx = TWO*lam*(NE-pref1)*(NE-pref1);
	//double fx = ZERO;
	//gsl_vector_set(df,n,fx);

	// free the space
	gsl_vector_free(c);
	gsl_vector_free(gc);
}

void gslscfminimizer::EDIIS_fdf(const gsl_vector * t, void * params, 
		double *f, gsl_vector *df)
{
	*f = EDIIS_f(t,params);
	EDIIS_df(t,params,df);
}
*/

void gslscfminimizer::scfMinimize(SCFEADIIS& diis)
{
	// now set the function information
	// for the parameters appearing the EDIIS/ADIIS
	// it should already been solved in diis object
	// also we derive the error here
	// the variable length will be plus the lambda
	//gsl_multimin_function_fdf my_func;
	gsl_multimin_function my_func;
	UInt n = diis.getNTerms();
	my_func.n = n;
	my_func.params = (void*)diis.getParams();
	if (diis.useADIIS()) {
		my_func.f = ADIIS_f;
		//my_func.df = ADIIS_df;
		//my_func.fdf = ADIIS_fdf;
	}else if (diis.useEDIIS()) {
		my_func.f = EDIIS_f;
		//my_func.df = EDIIS_df;
		//my_func.fdf = EDIIS_fdf;
	}else{
		string infor = "job invalid. Only EDIIS/ADIIS is allowed to pass in this function.";
		Excep excep("gslscfminimizer","scfMinimize",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}

	// initilize the variables t
	gsl_vector* t = gsl_vector_alloc(n);
	double initVal = ONE/sqrt(n);
	for(UInt i=0; i<n; i++) {
		gsl_vector_set(t,i,initVal);
	}

	// set the step size
	double stepSize = static_cast<double>(diis.getStepSize());
	gsl_vector* ss = gsl_vector_alloc(n);
	gsl_vector_set_all(ss,stepSize);

	// set the minimization method
	// basically we use the bfgs
	// in the future we may introduce flexibility 
	// to change it
	//double tol     = static_cast<double>(diis.getTol());
	double sizeTol = static_cast<double>(diis.getSizeTol());
	//const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2;
	//gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc (T, n);
	//gsl_multimin_fdfminimizer_set(s,&my_func,t,stepSize,tol);
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T,n);
	gsl_multimin_fminimizer_set(s,&my_func,t,ss);


	// now this is the real step of minimization
	// the default return status is false
	bool success = false;
	UInt iter = 0;
	Int status;
	//double previousEnergy = 0.0E0;
	do {

		// do the job
		iter++;
		//status = gsl_multimin_fdfminimizer_iterate(s);
		status = gsl_multimin_fminimizer_iterate(s);

		// testing the gradient
		//status = gsl_multimin_test_gradient (s->gradient, tol);

		// let's see the gradient
		//double dnrm = ZERO;
		//for(UInt i=0; i<n; i++) {
		//	double v = gsl_vector_get(s->gradient,i);
		//	dnrm += v*v;
		//}
		//dnrm = sqrt(dnrm/n);

		// debugging information
		//double err = gsl_multimin_fdfminimizer_minimum(s);
		double size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size (size,sizeTol);
		double err = gsl_multimin_fminimizer_minimum(s);
		//printf("current error function value %-15.10f, size is %-15.10f for iter %d\n", 
		//		err, size, (Int)iter);

		//
		// testing the energy difference
		// if the energy difference is less than 
		// the tol and this energy is smaller than
		// the previous one, we think it's converged,
		// too
		// this is because the gradient may not really
		// converged to the given tol
		/*
		if (status == GSL_CONTINUE && iter>1) {
			double diff = fabs(err-previousEnergy);
			if (diff<tol) {
				status = GSL_SUCCESS;
			}
		}
		previousEnergy = err;
		*/

		//
		// if the last element becomes very small,
		// we need to adjust the whole ti 
		/*
		for(UInt i=0; i<n; i++) {
			double ti = gsl_vector_get(s->x,i);
			if (ti<0.001E0) {
				gsl_vector_set(s->x,i,0.05E0);
				status = GSL_CONTINUE;
			}
		}
		*/

		// let's see the result
		// first case, if we reach the maximum iterations
		if (status == GSL_SUCCESS) {
			success = true;
			for(UInt i=0; i<n; i++) {
				double v = gsl_vector_get(s->x,i);
				gsl_vector_set(t,i,v);
			}
			break;
		}else if (status == GSL_CONTINUE && iter == diis.getMaxIters()) {

			// we returns
			success = true;
			for(UInt i=0; i<n; i++) {
				double v = gsl_vector_get(s->x,i);
				gsl_vector_set(t,i,v);
			}

			// print out information for users
			printf("current error function value %-15.10f, size is %-15.10f for iter %d\n", 
					err, size, (Int)iter);

			// print out some warning message and break
			string infor = "minimization reaches it's maximum interations, therefore we use approx. results";
			Excep excep("gslscfminimizer","scfMinimize",EXCEPTION_SCFEADIIS_WARNING,infor);
			handleExcep(excep);
			break;
		}else if (status == GSL_ENOPROG ) {
			break;
		}

	} while (status == GSL_CONTINUE && iter<=diis.getMaxIters());

	// write the result back, and free space
	//gsl_multimin_fdfminimizer_free(s);
	gsl_multimin_fminimizer_free(s);
	if (success) {

		// we copy the coefficients code above
		// so that we do not have extra vector constructed
		Double* coefs = diis.coef();

		// get the sum of ti
		double sumTi = ZERO;
		UInt nTerms = n;
		for(UInt i=0; i<nTerms; i++) {
			double ti  = gsl_vector_get(t,i);
			sumTi += ti*ti;
		}

		// now do the calculation
		for(UInt i=0; i<nTerms; i++) {
			double ti  = gsl_vector_get(t,i);
			double v   = ti*ti/sumTi;
			coefs[i]   = v;
		}

		// check the last coefs
		/*
		double lastTi = coefs[nTerms-1];
		if (lastTi<0.001E0) {
			coefs[nTerms-1] = 0.05E0;
			Double sum = ZERO;
			for(UInt i=0; i<nTerms; i++) {
				sum += coefs[i];
			}
			for(UInt i=0; i<nTerms; i++) {
				coefs[i] = coefs[i]/sum;
			}
		}
		*/

	}
	gsl_vector_free(t);
	gsl_vector_free(ss);

	// double check the coefficients
	// each coefs should be >=0 and sum should be 1
	diis.checkCoefs();

	// double check the idempotency requirement
	// if false, let's check the lambda value
	//bool pass = diis.checkIdempotent();
	//if (! pass) {
	//	printf("the lambda value is %-12.10f\n", lambda);
	//}

	// now deal with possible error
	if (! success) {
		string infor;
		if (status == GSL_ENOPROG) {
			infor = "minimization fails and gsl reports error";
		}else {
			infor = "unknown error arises in minimization, please check the code in scfMinimize";
		}
		Excep excep("gslscfminimizer","scfMinimize",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}
}

