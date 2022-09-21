/**
 * \file    scfeadiis.h
 * \brief   produce the SCF DIIS result in terms of EDIIS/ADIIS procedure
 * \author  Fenglai Liu 
 *
 * please see the scfeadiis.cpp for more notes about EDIIS/ADIIS.
 *
 * EDIIS paper:
 * "A black-box self-consistent field convergence algorithm: One step closer"
 * Kudin, Konstantin N. and Scuseria, Gustavo E. and  Cances Eric
 * The Journal of Chemical Physics, 116, 8255-8261, 2002
 *
 * ADIIS paper:
 * "Accelerating self-consistent field convergence with the augmented Roothaan
 * Hall energy function"
 * Hu, Xiangqian and Yang, Weitao
 * The Journal of Chemical Physics, 132, 054109, 2010
 *
 */
#ifndef SCFEADIIS_H
#define SCFEADIIS_H
#include "libgen.h"
#include "scalablevec.h"

namespace denmtrx {
	class DenMtrx;
}

namespace spinmatrix{
	class SpinMatrix;
}

namespace scfconv {
	class SCFConv;
}

namespace diiscontroller {
	class DIISController;
}

namespace oneemtrx {
	class OneEMtrx;
}

namespace matrix {
	class Mtrx;
}

namespace scfeadiis {

	using namespace denmtrx;
	using namespace spinmatrix;
	using namespace scfconv;
	using namespace diiscontroller; 
	using namespace oneemtrx;
	using namespace matrix;

	/**
	 * this class is used to obtain the SCF extrapolated MO result
	 * in terms of EDIIS/ADIIS procedure
	 *
	 * both EDIIS and ADIIS construct the error functional based on the 
	 * historical density matrix and historical fock matrix just like 
	 * Pulay's DIIS.
	 *
	 * However, both EDIIS and ADIIS optimized the coefficients in the 
	 * convex set - that means the coefficients should satisfy the condition
	 * that c_{i} >=0 and \sum_{i}c_{i} = 1
	 *
	 * therefore, we do not use the Lagrange multiplier method as we did in DIIS.
	 * Here we will use algorithms to directly search the answer. Because the 
	 * exhausted search grows very fast in terms of the coefficient space,
	 * therefore we restrict the space less than 10 so that to make it doable.
	 *
	 * the zero, 1st, 2ed term defined below is related to the terms appearing 
	 * in the error functional. For example, in the paper of ADIIS the error 
	 * functional is defined as:
	 *
	 * \f$ f = E(D_{n}) + 2\sum_{i}^{n}c_{i}(D_{i}-D_{n}|F(D_{n})) + 
	 * \sum_{i=1}^{n}\sum_{j=1}^{n} c_{i}c_{j}(D_{i}-D_{n}|F(D_{j})-F(D_{n}))\f$
	 *
	 * the zero order term is \f$ E(D_{n}) \f$
	 *
	 * the first order term is \f$ (D_{i}-D_{n}|F(D_{n})) \f$
	 *
	 * the second order term is \f$ (D_{i}-D_{n}|F(D_{j})-F(D_{n})) \f$
	 *
	 * user can use EDIIS or ADIIS upon their choice.
	 *
	 */
	class SCFEADIIS {

		private:

			UInt jobType;            ///< whether it's EDIIS or ADIIS?
			///
			/// maximum precision that the result coefficient has
			/// 3 meaning the coefficients is as accurate as up to 0.001
			///
			UInt maxPrec;       
			UInt nTerms;             ///< number of terms in current calculation

			///
			/// data to perform the EDIIS/ADIIS calculation
			/// we do not count in the 0th term in ADIIS
			///
			DoubleVec params;        ///< all of 1th and 2th terms 
			DoubleVec coefs;         ///< the result coefficients 

			///
			/// get the 1st order term
			///
			Double* getT1() {
				return &params[0];
			};

			///
			/// get the 2ed order term
			///
			Double* getT2() {
				UInt offset = nTerms;
				return &params[offset];
			};

			///
			/// get idemPotent requirement matrix
			///
			//double* getIdemMatrix() {
			//	UInt offset = nTerms*nTerms+nTerms+1;
			//	if (useADIIS()) offset += 1;
			//	return &params[offset];
			//};

			/**
			 * get the length of params array
			 * n   first order term
			 * n*n second order term
			 */
			UInt getLenParams() const {
				UInt len = nTerms+nTerms*nTerms;
				return len;
			};

			///
			/// update the error functional terms, see the comments inside the function
			///
			/// \param cont  diis like algorithm control information
			/// \param conv  scfconv, which is the topper class
			/// \param den   current spin polarized density matrix
			/// \param Fock  current spin polarized Fock matrix
			///
			void updateTerms(const DIISController& cont, const SCFConv& conv, 
					const DenMtrx& den, const SpinMatrix& fock);

			///
			/// this is used to compute the idempotent matrix requirement for the result 
			/// density matrix
			//
			/// Note for the idempotent requirement for density matrix. Here for
			/// both EDIIS/ADIIS the new density matrix (as a sum of the density
			/// matrices over the coefficients) is required to be idempotent.
			/// For achieving that purpose, we append the following term to the 
			/// error functional:
			///
			/// \f$ f = \sum_{i}\sum_{j}C_{i}C_{j}{-(D_{i}SD_{j}|S) + N} \f$
			///
			/// When the new density matrix is:
			///
			/// \f$ D = \sum_{i}C_{i} D_{i}\f$
			///
			/// It's easy to see this is just DSD = D combined with tr(DS) = N
			/// If the error function is minimized, we expected that the f 
			/// above should be approaching to zero.
			///
			//void updateIdemMtrx(const SCFConv& conv, const Mtrx& S);

			///
			/// this is the driver function to derive the best coefficients 
			/// who minimize the error functional of ADIIS/EDIIS through
			/// the contraints:
			/// \f$ \sum_{i}c_{i} = 1 \f$ and \f$ c_{i} >=0 \f$
			///
			/// more details please see the comments inside the function
			///
			void findCoefs(bool withExpensiveSolver);

		public:

			/**
			 * constructor
			 * 
			 * initilize the information from the top level class of scfconv
			 */
			SCFEADIIS(const DIISController& cont);

			/**
			 * destructor
			 */
			~SCFEADIIS() { };

			/**
			 * get the length of terms
			 */
			UInt getNTerms() const { return nTerms; };

			/**
			 * do the SCF minimization job so that to derive the result 
			 * coefficients
			 */
			void minimization(const DIISController& cont, const SCFConv& conv, 
					const DenMtrx& den, const SpinMatrix& Fock);

			/**
			 * get the coefficients - in non-constant way
			 */
			Double* coef() {
				return &coefs[0];
			};

			/**
			 * update the outside coefficients
			 */
			void updateCoefs(DoubleVec& coe) const {
				coe = coefs;
			};

			/**
			 * checking the coefs status
			 */
			void checkCoefs() const;

			/**
			 *  checking the idempotent requirement status
			 *
			 *  the below comment goes with the original implementation, which is
			 *  together with use of GSL library. We may need to find a time to
			 *  look into this, too.
			 *
			 *	 Note for the idempotent requirement for density matrix. Here for
			 *	 both EDIIS/ADIIS the new density matrix (as a sum of the density
			 *	 matrices over the coefficients) is required to be idempotent.
			 *	 For achieving that purpose, we append the following term to the 
			 *	 error functional:
			 *	
			 *	 \f$ err = \sum_{i}\sum_{j}C_{i}C_{j}{-(D_{i}SD_{j}|S) + N} \f$
			 *	
			 *	 When the new density matrix is:
			 *	
			 *	 \f$ D = \sum_{i}C_{i} D_{i}\f$
			 *	
			 *	 It's easy to see this is just DSD = D combined with tr(DS) = N
			 *	 If the error function is minimized, we expected that the err
			 *	 above approaching to zero in the optimized coefficients.
			 *
			 *	 This is made by appending err to the original error functional f:
			 *
			 *	 \f$ f^{'} = f + \lambda*\lambda*err^{2}\f$
			 *
			 *	 so that if we have derivatives for f' on lambda, then this term
			 *	 could be zero.
			 *
			 *	 However, in the practical calculation we usually found that the 
			 *	 idempotent condition is not achieved in optimizing the coefficients
			 *	 of c (you can see gslscfminimizer.cpp). As a result, the err is not
			 *	 0 but lambda is going to zero. Therefore the idempotent requirement
			 *	 is failed in practical calculation, and even worse; the ADIIS/EDIIS
			 *	 with idempotent condition behaves much worse than the one without 
			 *	 idempotent condition. Therefore we finally abandoned this requirement.
			 *	 and all of related codes are commented out.
			 */
			//bool checkIdempotent() const;

			/**
			 * print out the parameters as well as coefficients
			 * for debug purpose
			 */
			void print() const;
	};
}

#endif

