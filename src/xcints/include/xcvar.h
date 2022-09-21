/**
 * \file    xcvar.h
 * \brief   control center for the variable information
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef XCVAR_H
#define XCVAR_H
#include <string>
#include "libgen.h"
#include "scalablevec.h"
#include "xcvarinfor.h"

namespace xcfunc {
	class XCFunc;
}

namespace xcvar {

	using namespace xcfunc;
	using namespace xcvarinfor;
	using namespace std;

	class XCVar {

		private:

			//
			// basic variable information
			//
			UInt maxVarDerivOrder;             ///< maximum variable derivatives order
			UInt maxFunDerivOrder;             ///< maximum functional derivatives order
			UInt maxBasisDerivOrder;           ///< maximum basis set derivative order
			UInt totalAlphaVar;                ///< how many alpha variables
			UInt totalVarNum;                  ///< how many variables actually used
			UInt totalNumVarWithGamma;         ///< how many variables - with gamma 
			UIntVec varTypes;                  ///< stores the types of variables

			//
			// functional derivatives information for variables without gamma
			//
			UIntVec d1Vars;                   ///< stores the variable positions in D1F
			UIntVec d2Vars;                   ///< stores the variable positions in D2F
			UIntVec d3Vars;                   ///< stores the variable positions in D3F

			//
			// functional derivatives for variables with gamma
			//
			UIntVec d1GVars;                  ///< stores the gamma variable positions in D1F
			UIntVec d2GVars;                  ///< stores the gamma variable positions in D2F
			UIntVec d3GVars;                  ///< stores the gamma variable positions in D3F

			/**
			 * setting the functional derivatives information
			 * \param order      functional derivatives order
			 * \param withGamma  whether we set the derivatives with Gamma variable
			 */
			void setVarDeriv(const UInt& order, bool withGamma);

			/**
			 * setting variable derivatives order and functional derivatives
			 * order
			 * we note that this is related to job and job order
			 */
			void setDerivOrder(const XCFunc& xcfunc, const UInt& job, const UInt& jobOrder);

		public:

			/**
			 * do we have rho? 
			 */
			bool hasRho() const;

			/**
			 * do we have gradient rho? - so as gamma var?
			 */
			bool hasGRho() const;

			/**
			 * do we have tau?
			 */
			bool hasTau() const;

			/**
			 * do we have lap?
			 */
			bool hasLap() const;

			/**
			 * do we have exrho?
			 */
			bool hasExRho() const;

			/**
			 * do we have beta variable?
			 */
			bool hasBetaVar() const {
				if(totalAlphaVar == totalVarNum) return false;
				return true;
			};

			/**
			 * return the number of variables actually used
			 */
			UInt getVarNum() const {
				return totalVarNum;
			};

			/**
			 * get variable number for gamma series
			 */
			UInt getNumVarWithGamma() const {
				return totalNumVarWithGamma;
			};

			/**
			 * get the number of functional derivatives for a 
			 * given order
			 */
			UInt getNumFuncDeriv(const UInt& order, bool withGamma = false) const;

			/**
			 * get the variable position from the given var - D1F
			 */
			UInt getVarPos(const UInt& var, bool withGamma = false) const; 

			/**
			 * get the variable position from the given var - D2F
			 */
			UInt getVarPos(const UInt& var1, const UInt& var2, bool withGamma = false) const; 

			/**
			 * from the index get the var 
			 */
			UInt getVar(const UInt& i, bool withGamma = false) const; 

			/**
			 * In the variable generation stage, we need to know
			 * what kind of "half Var" we need to generate(see batchvar.h)
			 * This is for the combination between density matrix
			 * and the gradient basis sets
			 */
			bool do1stGradHalfVar() const;

			/**
			 * In the variable generation stage, we need to know
			 * what kind of "half Var" we need to generate(see batchvar.h)
			 * This is for the combination between density matrix
			 * and the gradient basis sets
			 */
			bool do2edGradHalfVar() const;

			/**
			 * return the information of the maximum variable derivatives
			 */
			UInt getMaxVarDerivOrder() const {
				return maxVarDerivOrder;
			};

			/**
			 * return the information of the maximum functional derivatives
			 */
			UInt getMaxFunDerivOrder() const {
				return maxFunDerivOrder;
			};

			/**
			 * return the information of the maximum basis set derivatives
			 */
			UInt getMaxBasisDerivOrder() const {
				return maxBasisDerivOrder;
			};

			/**
			 * build the variable information for the functional
			 * calculation use - this array will be used in fortran section
			 */
			void setupVarInfor(Int* infor) const;

			/**
			 * in transformation work(here it means that we transformed the functional
			 * derivatives result from Gamma variable into the gradient rho variable),
			 * in terms of the close shell situation we may need to omit the 
			 * beta part variable(in close shell beta part is always the same with 
			 * alpha part, so we do not need to process it). Therefore, we do the 
			 * identify job here
			 */
			bool neglectVar(const UInt& var) const;

			/**
			 * constructor 
			 * \param  xcfunc   provide the xc functional information
			 * \param  nSpin    number of spin states
			 * \param  jobOrder the current job (ground state DFT, TDDFT etc.)
			 * \param  jobOrder the current job order (energy, gradient etc.)
			 */
			XCVar(const XCFunc& xcfunc, const UInt& nSpin, 
					const UInt& job, const UInt& jobOrder);

			///
			/// destructor
			///
			~XCVar() { };

			/**
			 * print function used in debug purpose
			 */
			void print() const;
	};

}

#endif

