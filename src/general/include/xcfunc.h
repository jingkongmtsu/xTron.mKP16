/**
 * \file    xcfunc.h
 * \brief   to gather functional information
 * \author  Fenglai Liu 
 */
#ifndef XCFUNC_H
#define XCFUNC_H
#include "libgen.h"
#include<string>
#include<vector>
#include "scalablevec.h"

namespace xcfunc {

	using namespace std;

	///
	/// here below let's first define the constant values used for the parameters in the XCFUNC
	///
	/// for example, the becke05 functional P and Q parameter values etc.
	///
	
	///
	/// default value for B13 Correlation functional method switch, see the variable of B13CoorMethod
	///
	const UInt     B13COORMETHOD_DEFAULT  = 0;

	/// B05_NDPAR_method let you pick modified nondynamic parallel spin correction.
	//  0 - original, 1 - min, 2 - jumping.
	const UInt     B05NDPARMETHOD_DEFAULT = 0;
	///
	/// default value for becke05_p
	///
	const Double   BECKE05_P_DEFAULT      = 115.0E0;

	///
	/// default value for becke05_q
	///
	const Double   BECKE05_Q_DEFAULT      = 120.0E0;

	///
	/// default value for kp14_alpha
	///
	const Double   KP14_ALPHA_DEFAULT     = 0.038E0;
	
	///
	/// default value for kp14_b
	///
	const Double   KP14_B_DEFAULT         = 1.355E0; 
	
	///
	/// default value for kp14_c_ndpar
	///
	const Double   KP14_C_NDPAR_DEFAULT   = 1.128E0;

	///
	/// default value for kp14_c_ndpar_cap
	///
	const Double   KP14_C_NDPAR_CAP_DEFAULT   = 0.5E0;

	/**
	 * \class   XCFunc
	 * \brief   describing the functional information
	 *
	 * For the xcfunc class, it accept one/two exchange-correlation name
	 * which should be already defined in the configuration file 
	 *
	 * Here we note that we only need one name or two names connected with space. 
	 * for example, Becke88 LYP means we want to do Becke 88 exchange together 
	 * with the LYP correlation. It's always exchange is first then the 
	 * correlation functional. 
	 *
	 * What's more, we need to say more about orbital functional. orbital functional
	 * means that this functional depends on the orbital. HF method, and all of 
	 * traditional post-HF methods are in this type
	 *
	 * more note about elementary functional. For a lot of cases that the functional
	 * itself is composed by some other functionals. These functional acts as "components"
	 * as the new functional. We call these functional as elementary functional.
	 *
	 * for maxFunDerivOrder, the default value is 1. Since the inplementation of the 
	 * functional always comes with functional derivatives in DFT part (this is usually
	 * used in the DFT code part)
	 */
	class XCFunc {

		private:

			////////////////////////////////////////////////////////////////////////////////////
			/// parameters information
			///
			/// sometimes the functional part needs the parameters. for example,
			/// in the Becke05 functional there's two parameters defined, called becke_p
			/// and becke_q. Both of the two parameters usually constant, but user can
			/// vary their values to debug the functional. We call these parameters as 
			/// functional parameters.
			///
			/// the parameters defined here are divided into two groups.
			///
			/// the first group is called method switch. Currently it's used to identify
			/// which form of B13 correlation functional we use. For the B13 coorelational,
			/// we have two forms, one is for the lambda averaged form, the other one 
			/// is lambda = 1 form. Different functional use different form. So you need to 
			/// specify it when use it.
			///
			/// the second group is the double parameters. For example the p, q values in
			/// the becke05 functional, or the three parameters used in KP14 functional etc.
			////////////////////////////////////////////////////////////////////////////////////

			///
			/// B13 correlational functional method switch value, only 1 or 2 accepted;
			///
			UInt B13CoorMethod;   

			///
			/// Choice for the nondynamic correlation parallel spin correction.
 			UInt B05NDParMethod;       

			///
			/// parameters for Becke05 functional, p and q
			///
			Double becke05_p;                 ///< the P value in the becke05 functional
			Double becke05_q;                 ///< the Q value in the becke05 functional

			///
			/// parameters for the KP14 functional 
			///
			Double kp14_alpha;                ///< alpha value in the KP14
			Double kp14_b;                    ///< b value in the KP14
			Double kp14_c_ndpar;              ///< c_ndpar value in the KP14
			Double kp14_c_ndpar_cap;          ///< c_ndpar_cap value in the KP14

			// functional information
			string exName;                    ///< the exchange functional name
			string ecName;                    ///< the correlation functional name
			UInt lenEx;                       ///< number of exchange components
			UInt lenEc;                       ///< number of correlation components
			UInt maxFunDerivOrder;            ///< maximum functional derivative order we can do
			bool islinear;                    ///< whether the components are linear combined
			bool hasRho;                      ///< whether has density
			bool hasGRho;                     ///< whether has gradient of density
			bool hasTau;                      ///< whether has kinetic density
			bool hasLap;                      ///< whether has Laplacian term
			bool hasExchRho;                  ///< whether has exchange energy density
			bool orbitalFunc;                 ///< whether this is orbital functional?
			vector<string> exchComponents;    ///< the components for exchange functional
			vector<string> corrComponents;    ///< the components for correlation functional
			DoubleVec  exchFactors;           ///< the factors for exchange functional
			DoubleVec  corrFactors;           ///< the factors for correlation functional

			/**
			 * after we get the name of the functional, we need to proceed its 
			 * details by exploring from the configuration file 
			 */
			void setupXCFunc();

		public:

			/**
			 * this is actually to read in functional name from input file
			 * for given section
			 */
			XCFunc(const string& input, const UInt& sec);

			/**
			 * this constructor is more clean, we just need a name, which is ex name
			 */
			XCFunc(const string& name0);

			/**
			 * this constructor is more clean, we just need a ex name and ec name
			 */
			XCFunc(const string& ex, const string& ec);

			/**
			 * destructor
			 */
			~XCFunc() { };

			/**
			 * debugging print
			 */
			void print() const;

			/**
			 * whether the given functional is linear?
			 */
			bool isLinear() const { return islinear; };

			/**
			 * whether we have correlation part?
			 */
			bool hasCorrelation() const {
				if (lenEc > 0) return true;
				return false;
			}

			/**
			 * whether the given functional has HF exchange 
			 * energy(so called "hybrid functional")
			 */
			bool isHybrid() const;

			/**
			 * if it has HF exchange in the functional, we are able to
			 * return the corresponding coefficients here
			 * if it does not, just return 0
			 */
			Double HFCoeff() const;

			/**
			 * return the number of elementary exchange functionals
			 */
			UInt getLenEx() const { return lenEx; };

			/**
			 * return the number of elementary correlation functionals
			 */
			UInt getLenEc() const { return lenEc; };

			/**
			 * get the total number of elementary functionals 
			 */
			UInt getElementaryFuncNumber() const { return lenEx + lenEc; };

			/**
			 * return exchange functionals
			 */
			const string& getEXName(UInt i) const{ return exchComponents[i];   };

			/**
			 * return correlation functionals
			 */
			const string& getECName(UInt i) const{ return corrComponents[i];  };

			/**
			 * return exchange functional coeffcients
			 */
			Double getEXFactor(UInt i) const{ return exchFactors[i];    };

			/**
			 * return correlation functional coeffcients
			 */
			Double getECFactor(UInt i) const{ return corrFactors[i];    };

			/**
			 * whether we have the rho variable?
			 */
			bool hasVarRho() const { return hasRho; };

			/**
			 * whether we have the gradient rho variable?
			 */
			bool hasVarGRho() const { return hasGRho; };

			/**
			 * whether we have the tau variable?
			 */
			bool hasVarTau() const { return hasTau; };

			/**
			 * whether we have the Laplacian variable?
			 */
			bool hasVarLap() const { return hasLap; };

			/**
			 * whether we have the exchange energy density?
			 */
			bool hasVarExRho() const { return hasExchRho; };

			/**
			 * whether it's an orbital functional?
			 */
			bool isOrbitalFunc() const { return orbitalFunc; };

			/**
			 * whether the given exchange functional 
			 * has hole function?
			 */
			bool hasExchangeHole() const; 

			/**
			 * whether the functional use B05 form of static correlation?
			 */
			bool useB05FormStaticCoor() const;

			/**
			 * whether this function use Hirshfeld weights?
			 * see the hrshfeld.cpp in xcints folder for more information
			 *
			 * right now it's only B13 strong functional uses the Hirshfeld weights
			 */
			bool useHirshfeldWeights() const {
				if (exName == "B13STRONG") {
					return true;
				}
				return false;
			};

			/**
			 * whether this functional use atom density matrix data?
			 * right now if functional needs hirshfeld weights, we need the atom density matrix
			 * data
			 */
			bool useAtomDenMtrxData() const {
				if (useHirshfeldWeights()) return true;
				return false;
			};

			/**
			 * whether this function is B13-Strong functional?
			 * 
			 * the B13-Strong functional is a kind of different from other functionals,
			 * it needs additional treatment in the batchfunc class. See the class of 
			 * batchfunc in xcints folder, for more details
			 */
			bool isB13Strong() const {
				if (exName == "B13STRONG") {
					return true;
				}
				return false;
			};

			/**
			 * whether it uses the B13 correlation functional
			 */
			bool useB13Coor() const;

			/**
			 * how many functional values we may have?
			 * this is helping batch functional class
			 * to identify the functional value number
			 */
			UInt getNFuncValues() const {

				// at least we have one defined
				UInt nFunc = 1;

				// we may or may not have correlation pieces
				if (hasCorrelation()) nFunc++;

				// if the exchange is actually the exchange hole 
				// function, then the functional value is spin-dependent
				// that means it has two values
				if (hasExchangeHole()) nFunc++;

				return nFunc;
			};

			///
			/// checking the input order with this maximum order
			/// the input order can not be higher than the maximum order
			///
			bool checkFuncDerivOrder(UInt order) const {
				if (order > maxFunDerivOrder) return false;
				return true;
			}

			/**
			 * return the name of exchange functional
			 */
			const string& getEXName() const { return exName; };

			/**
			 * return the name of correlation functional
			 */
			const string& getECName() const { return ecName; };

			/**
			 * return the functional name
			 */
			string getFuncName() const {
				string name = exName;
				if (ecName != "NONE") name = name + ecName;
				return name;
			};

			/**
			 * whether this is calculation without DFT involved?
			 */
			bool withoutDFTFunc() const;

			/**
			 * return the becke05 functional defined P value
			 */
			Double getBecke05PVal() const { return becke05_p; };

			/**
			 * return the becke05 functional defined Q value
			 */
			Double getBecke05QVal() const { return becke05_q; };

			/**
			 * return the KP14 functional defined alpha value
			 */
			Double getKP14AlphaVal() const { return kp14_alpha; };

			/**
			 * return the KP14 functional defined b value
			 */
			Double getKP14BVal() const { return kp14_b; };

			/**
			 * return the KP14 functional defined c_ndpar value
			 */
			Double getKP14CNDPARVal() const { return kp14_c_ndpar; };

			/**
			 * return the KP14 functional defined c_ndpar_cap value
			 */
			Double getKP14CNDPARCapVal() const { return kp14_c_ndpar_cap; };

			/**
			 * return the B13 correlation functional method switch
			 */
			UInt getB13CorrMethod() const { return B13CoorMethod; };

			/**
			 * return the B05 ND par method switch
			 */
			UInt getB05NDParMethod() const { return B05NDParMethod; };
	};

}

#endif
