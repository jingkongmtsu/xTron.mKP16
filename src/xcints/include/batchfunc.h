/**
 * \file    batchfunc.h
 * \brief   monitoring functional derivatives etc. calculation in batches
 * \author  Fenglai Liu and Jing Kong
 *
 * note
 *
 * for functional derivatives implementation, the user and programmer should be 
 * aware of the small tau and small exchange energy density (see the comments of 
 * function scalFuncDeriv). After you add the functional implementation to the 
 * doFuncDeriv, please be aware that you may need to scale the derivatives by
 * call the function of scalFuncDeriv.
 *
 */
#ifndef BATCHFUNC_H
#define BATCHFUNC_H
#include "libgen.h"
#include "scalablevec.h"
#include "matrix.h"
using namespace matrix;

namespace xcvar {
	class XCVar;
}

namespace xcfunc {
	class XCFunc;
}

namespace batchvar {
	class BatchVar;
}

namespace batchgrid {
	class BatchGrid;
}

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace hirshfeld {
	class Hirshfeld; 
}

namespace denmtrx {
	class DenMtrx;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace xcenergyinfor {
	class XCEnergyInfor;
}

namespace batchfunc {

	using namespace xcvar;
	using namespace xcfunc;
	using namespace batchvar;
	using namespace batchgrid;
	using namespace xcintsinfor;
	using namespace xcenergyinfor; 
	using namespace sigatombasis;
	using namespace hirshfeld; 
	using namespace denmtrx;

	/**
	 * \class   BatchFunc
	 * \brief   functional value and derivatives etc.
	 *
	 * the setting of batchfunc is consistent with the xcfunc class. If only one functional
	 * is given, then no matter whether it's exchange or correlation we will push it to 
	 * the funvals. It could be exchange functional, or only correlation functional. Therefore,
	 * we only provide the member function of getFunVals() and did not distingurish the 
	 * exchange and correlation functional values.
	 */
	class BatchFunc {

		private:

			// data
			UInt ng;                  ///< number of the grids
			UInt nDen;                ///< number of densities actually used
			UInt nAlpha;              ///< number of alpha electrons
		  UInt nBeta;               ///< number of beta electrons	
			UInt nEXFunc;             ///< number of exchange functionals
			UInt nECFunc;             ///< number of correlation functionals
			UInt nFunc;               ///< number of functionals
			UInt nD1FWithGamma;       ///< number of variables for D1F(with gamma version)
			UInt nD1F;                ///< number of variables for D1F(with gradient rho)
			Double tol;               ///< tolerance
			DoubleVec hirWeights;     ///< the possible modified hirshfeld weights
			DoubleVec funvals;        ///< functional values 
			Mtrx d1F;                 ///< 1st order of functional derivatives with gradient rho

			/**
			 * by using the information from xcvar and batch var, we transform
			 * the derivatives with respect to gamma variable into the gradient
			 * of rho - so that to get rid of gamma permanently
			 * \param  TF  functional derivatives array in terms of gamma variable
			 * \reuturn funvals  the final functional values and derivatives array
			 */
			void transformGammaDerivIntoDRhoDeriv(const Double* TF,
					const XCFunc& xcfunc, const XCVar& xcvar, const BatchVar& bvar);

			///
			/// this function returns the fitting Gamma value for the modified 
			/// BR89 exchange hole model, see the following B13 Strong paper
			///
			Double b13StrongFreeAtomGammaBRVal(const UInt& Z) const;

			///
			/// for some functionals, it's definition may be based on the small tau and/or
			/// small exchange energy density (here in this program both tau and exrho are not
			/// with factor 1/2. We call the tau as small tau if it's with 1/2 factor). For example,
			/// TPSS functional, PSTS functional etc.
			///
			/// if the functional is defined by small tau etc. user needs to scale the tau/exrho
			/// on the entry of Fortran function of functional implementation (see psts.f). Therefore
			/// you will use the small tau/small exrho according to the definition. However, in 
			/// this case the result derivatives is usually with the small variable, then the result
			/// derivatives is two times larger. Therefore we need to scale it with 1/2.
			///
			/// the TPSS functional is an exception that we do not do it. For TPSS exchange, it's 
			/// converted from maple script and over there the functional derivatives is taken 
			/// for "big" tau rather than small tau. for TPSS correlation functional, the tau
			/// variable is defined in the functional wt used in TPSS correlation; the derivatives
			/// is also made on the big tau rather than the small one. Therefore in total we do not
			/// need to scale it anyway.
			///
			void scalFuncDeriv(const string& name, const XCVar& xcvar, Double* TD1F) const;

		public:

			/**
			 * constructor based on batch variable, xc functional information
			 * as well as variable infor
			 */
			BatchFunc(const BatchVar& bvar, const DenMtrx& den, 
					const XCFunc& xcfunc, const XCVar& xcvar,
					const XCIntJobInfor& infor);

			/**
			 * destructor
			 */
			~BatchFunc() { };

			/**
			 * main body to evaluate the functional derivatives etc.
			 */
			void doFuncDeriv(const BatchVar& bvar, const XCFunc& xcfunc, const XCVar& xcvar, 
					XCEnergyInfor& xcEList);

			/**
			 * get the pointer for the first order functional derivative
			 * \param varPos  the variable position
			 */
			const Double* getD1F(const UInt& varPos) const {
				return d1F.getPtr(0,varPos);
			};

			/**
			 * get the pointer for the exchange functional value
			 */
			const Double* getFunVals() const {
				return &funvals[0];
			};

			/**
			 * get the exchange correlation energy for current batch
			 */
			void getBatchExc(const BatchGrid& grid, Double& exc) const;

			/**
			 * for B13 strong functional, it uses hirshfeld weights for a modified BR89 
			 * exchange hole model so that to make the Neff normalized(Neff is used in B05 functional). 
			 *
			 * See the paper of 
			 * J. Chem. Phys. 138, 074109 (2013)
			 * Density functionals for static, dynamical, and strong correlation
			 * Axel D. Becke
			 *
			 * in appendix C
			 */
			void formBR89EXHoleWeights(const MolShell& ms, const SigAtomBasis& sigAtomBasis, 
					const Hirshfeld& hirshfeld);

			/**
			 * debug printing 
			 */
			void print();

	};
}

#endif
