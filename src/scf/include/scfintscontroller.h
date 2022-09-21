/**
 * \file    scfintscontroller.h
 * \brief   mornitoring the threshold etc. for integral calculation for SCF calculation
 * \author  Fenglai Liu 
 */
#ifndef SCFINTSCONTROLLER_H
#define SCFINTSCONTROLLER_H
#include "libgen.h"
#include "gintsinfor.h"
#include "xcintsinfor.h"
using namespace gintsinfor;
using namespace xcintsinfor;

namespace scfparam {
	class SCFParam;
}
namespace scfconv {
	class SCFConv;
}

namespace spinmatrix {
	class SpinMatrix;
}

namespace scfintscontroller {

	using namespace scfparam;
	using namespace scfconv;
	using namespace spinmatrix;

	//
	// here let's define the options for the variable
	// of scfOption
	//
	
	///
	/// this option refers that after initial criteria
	/// is met, then the user's choice(or default threshold)
	/// is switched on for the rest of SCF calculation
	///
	const UInt AFTER_CONVERGE_WITH_USER_CHOICE       = 1; 

	///
	/// this option refers that after initial criteria
	/// is met, then the threshold will be gradually
	/// increased until the target SCF convergence is met
	///
	const UInt AFTER_CONVERGE_WITH_STEPWIZE_PROGRESS = 2; 

	///
	/// this option refers that we do not use coarse
	/// initial criteria, from the beginning of SCF
	/// we just use user's (or default threshold) 
	/// for all of SCF calculation. That is to say,
	/// no scf controller is used
	///
	const UInt NO_CONTROLLER                         = 3; 


	/**
	 * this class manipulates the threshold setting etc. for GInts module
	 * and XCInts module during SCF calculation.
	 *
	 * The most expensive parts in SCF calculation, is the GInts4D and 
	 * XCInts parts calculation. However, when the SCF calculation is far away
	 * from convergence, the very high accuracy threshold demanded by the user
	 * is actually not quite practical for use. It takes more time, but it
	 * provides the same reliability as low threshold during initial convergence
	 * progress (for example, user may want the gints4d running under threshold
	 * of 1.0E-14, however; for the first several SCF iterations even 1.0E-7 is 
	 * good in general sense).
	 *
	 * Therefore, how we manipulate the integral threshold etc. for maintaining  
	 * accurate/robust SCF and at the same time keeping SCF fast?
	 *
	 * The first point is how we choose the criteria. We now use the rms difference
	 * between the current Fock matrix and the previous Fock matrix. You can see
	 * the updateCriteria function for more details.
	 *
	 * After setting the initial criteria, we set up the ratios between initial 
	 * criteria and the threshold value for doing xcints, gints4d jobs etc. 
	 * Once the SCF calculation meets the initial criteria, we may do further 
	 * refining afterwards.
	 *
	 * we note that we do not have threshold control for gints2D. We just use
	 * the most accurate one. This is because kinetic and NAI etc. forms the 
	 * backbone of the fock matrix, and it's not expensive at all comparing with
	 * gints4d etc. Therefore we hope it can be as accurate as possible.
	 *
	 * If the initial criteria is achieved, next step we have two choices:
	 * -  switch on the user demanded threshold for the rest SCF calculation,
	 *    if it's still in convergence; until the target SCF convergence criteria
	 *    is met;
	 * -  gradually raising up the criteria here, for example; from 1.0E-2 to
	 *    1.0E-3, 1.0E-4 etc. so that to revise the threshold value accordingly
	 *    until the target SCF convergence criteria is met
	 * 
	 * how to choose it will be decided by the user. However, in default we choose
	 * the second choice.
	 *
	 * there's only one case that the no controller situation is used, that when
	 * the SCF cycle is only 1; that means we only have one SCF cycle. In this 
	 * case we just use the original setting.
	 */
	class SCFIntsController {

		private:

			//
			// options to control the class here
			//
			UInt scfOption;               ///< do we switch on user's demand etc. see above macro?
			Double initialCriteria;       ///< initial value of criteria

			//
			// ratio values for GInts2D and GInts4D
			//
			Double csThreshRatio;         ///< the ratio for performing cauchy-schwarz check
			Double gints4DSPThreshRatio;  ///< the ratio for forming shell pair threshold for GInts4D
			Double gints4DThreshRatio;    ///< the ratio for forming GInts4D threshold

			//
			// ratio values for XCInts
			//
			Double xcintsThreshRatio;     ///< the ratio for forming xcints threshold

			//
			// criteria setting
			//
			Double criteria;              ///< this is current criteria

			//
			// original gints infor and the working gints infor
			// the original one is set by default or from user's choice
			//
			GIntsInfor  inputGIntsInfor;  ///< input analytical Gaussian integrals job information
			GIntsInfor  workGIntsInfor;   ///< working analytical Gaussian integrals job information

			//
			// original xcints infor and the working xcints infor
			// the original one is set by default or from user's choice
			//
			XCIntsInfor  inputXCIntsInfor;///< input XC numerical integrals job information
			XCIntsInfor  workXCIntsInfor; ///< working XC numerical integrals job information

			/**
			 * this is the working function to reset the integral information
			 * if withInitCriteria is true, then we use the initial criteria;
			 * else we will use the calculated criteria value
			 */
			void setIntsInfor(bool withInitCriteria);

		public:

			/**
			 * constructor
			 * \param  param the SCF parameters
			 */
			SCFIntsController(const SCFParam& param);

			/**
			 * default destructor
			 */
			~SCFIntsController() { };

			/**
			 * do we really use the SCF integral controller?
			 */
			bool withSCFIntsController() const { 
				if (scfOption == NO_CONTROLLER) return false;
				return true;
			};

			/**
			 * get the gints infor
			 */
			const GIntsInfor& getGIntsInfor() const { return workGIntsInfor; };

			/**
			 * get the xcints infor
			 */
			const XCIntsInfor& getXCIntsInfor() const { return workXCIntsInfor; };
	};
}

#endif

