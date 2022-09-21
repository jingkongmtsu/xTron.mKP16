/**
 * \file    gslscfminimizer.h
 * \brief   do SCF minimization process through GSL library
 * \author  Fenglai Liu 
 *
 * see the comments in gslscfminimizer.cpp for more information.
 */
#ifndef GSLSCFMINIMIZER_H
#define GSLSCFMINIMIZER_H
#include "libgen.h"
#include "gsl/gsl_vector.h"

namespace scfeadiis {
	class SCFEADIIS;
}

namespace gslscfminimizer {

	using namespace scfeadiis;

	///
	/// to remove the condition that \f$\sum_{i}c_{i} = 1 \f$
	/// and \f$c_{i} >= 0\f$, we can re-define the \f$c_{i}\f$ via
	/// \f$t_{i}\f$:
	///
	/// \f$ c_{i} = t_{i}^{2}/\sum_{k}t_{k}^{2} \f$
	///
	/// t_{i} is the auxiliary variable, c_{i} is the 
	/// real variable 
	///
	void getCoefs(const gsl_vector *t, gsl_vector *c);

	///
	/// this is to derive the derivatives for \f$ c_{i} \f$ 
	/// in terms of \f$ t_{i} \f$:
	/// \f$ gc_{ij} = \frac{\partial c_{i}}{\partial t_{j}} \f$
	/// where i and j goes from 1 to n
	///
	/// n is the dimension of t
	///
	/// the result gc is a matrix, first dimension is on
	/// c, and the second dimension is on t
	///
	void getGradCoefs(const gsl_vector *t, gsl_vector *gc);

	///
	/// defines the ADIIS functional in terms of t and 
	/// parameters 
	///
	double ADIIS_f(const gsl_vector *t, void *params);

	///
	/// define the ADIIS functional derivatives in terms of 
	/// t
	///
	//void ADIIS_df(const gsl_vector *t, void *params, 
	//		gsl_vector *df);

	///
	/// this is a top function that simply put ADIIS_f
	/// and ADIIS_df together
	///
	//void ADIIS_fdf(const gsl_vector * t, void * params, 
	//		double *f, gsl_vector *df);

	///
	/// defines the EDIIS functional in terms of t and 
	/// parameters 
	///
	double EDIIS_f(const gsl_vector *t, void *params);

	///
	/// define the EDIIS functional derivatives in terms of 
	/// t
	///
	//void EDIIS_df(const gsl_vector *t, void *params, 
	//		gsl_vector *df);

	///
	/// this is a top function that simply put EDIIS_f
	/// and EDIIS_df together
	///
	//void EDIIS_fdf(const gsl_vector * t, void * params, 
	//		double *f, gsl_vector *df);

	///
	/// driver function to perform the SCF minimization work
	/// in terms of ADIIS/EDIIS
	///
	void scfMinimize(SCFEADIIS& diis);
}

#endif

