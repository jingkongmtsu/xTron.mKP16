/**
 * \file   gintsinfor.h
 * \brief  collecting information for analytical integrals calculation
 * \author Fenglai Liu 
 */
#ifndef GINTSINFOR_H
#define GINTSINFOR_H
#include "libgen.h"
#include "hgp_os_ints.h"
#include "globalinfor.h"
#include "integralinfor.h"
#include "integraljobs.h"

namespace molecule {
	class Molecule;
};

namespace xcfunc {
	class XCFunc;
};

namespace gintsinfor {

	using namespace globalinfor;
	using namespace integralinfor;
	using namespace integraljobs;
	using namespace molecule;
	using namespace xcfunc;

	//
	// here we provide some method to change the threshold 
	// value according to the calculation requirement
	//
	const UInt  GINTS2D_SHELL_PAIR_THRESH   = 0;  ///< threshold index for forming shell pair (2d)
	const UInt  GINTS4D_SHELL_PAIR_THRESH   = 1;  ///< threshold index for forming shell pair (4d)
	const UInt  GINTS2D_THRESH              = 2;  ///< threshold index for gints2d
	const UInt  GINTS4D_THRESH              = 3;  ///< threshold index for gints4d

	/**
	 * \class GIntsInfor
	 * \brief collecting general integral calculation information
	 */
	class GIntsInfor {

		protected:

			//
			// for general integral calculation
			//
			GlobalInfor ginfor;          ///< a copy of global setting 
			bool closeShell;             ///< whether this is close shell system?
			UInt nSpin;                  ///< number of spin states in the calculation

			//
			// data for the exchange part
			// do we need to scale it?
			// and what is the scale factor?
			// this part of information comes from xcfunc
			//
			bool doScaleforK;            ///< do we really do scaling?
			Double kCoeff;               ///< scale coefficients

			//
			// threshold section
			//
			// because 2D calculation is usually very fast, therefore
			// GInts2D threshold value is in default very tight. So as
			// the Cauchy-Schwarz integral boundary calculation.
			//
			Double spGInts2DThresh;      ///< shell pair forming threshold for GInts2D       
			Double spGInts4DThresh;      ///< shell pair forming threshold for GInts4D       
			Double gints2DThresh;        ///< threshold for generating 2D analytical integrals for energy
			Double gints2DThreshDeriv;   ///< threshold for generating 2D analytical integrals for derivatives
			Double gints4DThresh;        ///< threshold for generating 4D analytical integrals for energy  
			Double gints4DThreshDeriv;   ///< threshold for generating 4D analytical integrals for derivatives

			//
			// parallel setting
			//
			// we note, that the following setting is largely depending on the parallel condition.
			// for example, the setting of nAtomBlocksInList should be adjusted according to the cache size
			// of CPU
			//
			UInt nAtomBlocksInList;      ///< if we use block matrix list to store the tmp result, how many blocks inside?

			//
			// MISC
			// 
			Double linearDependenceThresh;  ///< threshold detecting linear dependency of basis sets 

		public:

			/**
			 * constructor from the global infor and molecule data
			 * \param infor        global information for the job
			 * \param mol          molecular information
			 */
			GIntsInfor(const GlobalInfor& infor, const Molecule& mol);

			/**
			 * destructor
			 */
			~GIntsInfor() { };

			/**
			 * get the number of spin states
			 */
			UInt getNSpin() const {
				return nSpin;
			};

			/**
			 * whether it's close shell?
			 */
			bool isCloseShell() const { return closeShell; };

			/**
			 * return the inner copy of global information
			 */
			const GlobalInfor& getGlobalInfor() const {
				return ginfor;
			};

			/**
			 * pick up the threshold value for forming shell pair
			 * for the modules of GInts4D
			 */
			Double spGInts4DThreshVal() const {
				return spGInts4DThresh;
			};

			/**
			 * pick up the threshold value for forming shell pair
			 * for GInts2D module
			 */
			Double spGInts2DThreshVal() const {
				return spGInts2DThresh;
			};

			/**
			 * pick up the threshold value for forming 4D analytical integrals
			 */
			Double int4DThreshVal() const {
				return gints4DThresh;
			};

			/**
			 * pick up the threshold value for forming 2D analytical integrals
			 */
			Double int2DThreshVal() const {
				return gints2DThresh;
			};

			/**
			 * pick up the threshold value for forming 4D analytical integrals derivatives
			 */
			Double intDeriv4DThreshVal() const {
				return gints4DThreshDeriv;
			};

			/**
			 * pick up the threshold value for forming 2D analytical integrals derivatives
			 */
			Double intDeriv2DThreshVal() const {
				return gints2DThreshDeriv;
			};

			/**
			 * return the linear dependence thresh
			 */
			Double getLinearDepThresh() const {
				return linearDependenceThresh;
			};

			/**
			 * update the information from xcfunc 
			 */
			void updateFromXCFunc(const XCFunc& xcfunc);

			/**
			 * do we do scaling for exchange?
			 */
			bool doKScaling() const {
				return doScaleforK;   
			}
			
			/**
			 * return the exchange scaling coefficient
			 */
			Double getKCoeff() const {
				return kCoeff;
			};

			/**
			 * here we provide a method to change the 
			 * threshold value inside the gintsinfor
			 *
			 * i is the mode name of the threshold, like 
			 * GINTS4D_SHELL_PAIR_THRESH
			 *
			 * thresh is the given threshold value provided
			 */
			void changeThresh(UInt i, Double thresh); 

			/**
			 * whether the gints4d will use block form of result, 
			 * rather than the whole fock matrix?
			 *
			 * currently if we do parallel, we use the block form result.
			 * Since the block form result will then update to atomFockMtrx;
			 * then update to blockMtrxList (see gints4d.cpp). For parallel
			 * part each thread has a reference to the global result and 
			 * we update the global result in race condition.
			 *
			 * for serial part, we will do the update directly
			 *
			 */
			bool useBlockResult() const { 
				if (ginfor.useMultiThreads()) return true;
				return false; 
			};

			/**
			 * return the number of atom block results used in block matrix list
			 * see the gints4d.cpp
			 */
			UInt getNAtomBlockResults() const { return nAtomBlocksInList; };

	};


	/**
	 * \class GIntJobInfor
	 * \brief by given the job ID, we will form the practical information
	 * to guide the real calculation
	 */
	class GIntJobInfor : public GIntsInfor {

		private:

			//
			// basic integral job information
			//
			UInt intJob;               ///< what kind of integral job we do here 
			UInt jobOrder;             ///< the order for the job

			///
			/// integral engine information (memory usage + possible 
			/// derivative information)
			///
			IntegralInfor engineInfor;

		public:

			/**
			 * constructor 
			 * \param infor        general integral information
			 * \param intjob       the job ID for electron integrals
			 * \param jobOrder     energy, gradient or second derivatives calculation?
			 */
			GIntJobInfor(const GIntsInfor& infor, const UInt& intJob0, UInt jobOrder0 = 0);

			/**
			 * destructor
			 */
			~GIntJobInfor() { };

			/**
			 * return the integral job name
			 */
			UInt getIntJob() const { return intJob; };

			/**
			 * return the integral job order
			 */
			UInt getJobOrder() const { return jobOrder; };

			/**
			 * pick up the shell pair threshold value according to the job type
			 * currently we only distinguish the 4D case and the others
			 */
			Double pickUpSPThresh() const {
				if (fourBodyInts(intJob)) return spGInts4DThreshVal();
				return spGInts2DThreshVal();
			};

			/**
			 * pick up the threshold value according to the job type
			 * currently we only distinguish the 4D case and the others
			 */
			Double pickUpThresh() const {
				if (jobOrder == 0) {
					if (fourBodyInts(intJob)) return int4DThreshVal();
					return int2DThreshVal();
				}else{
					if (fourBodyInts(intJob)) return intDeriv4DThreshVal();
					return intDeriv2DThreshVal();
				}
			};

			/**
			 * whether the integral module will use heap memory?
			 * we will check the engine infor to get it
			 */
			bool useHeapMem(UInt maxL) const {
				if (engineInfor.getMaxMemLen(maxL) > 0) return true;
				return false;
			}; 

			/**
			 * return the number of elementary derivatives
			 */
			UInt getNElemDeriv() const { return engineInfor.getNElemDerivs(); };

			/**
			 * return the number of redundant derivatives
			 */
			UInt getNRedDeriv() const { return engineInfor.getNRedDerivs(); };

			/**
			 * set the memory size according to the maxL given
			 *
			 * because the memory length get from engine infor is just
			 * what the integral needs, therefore we increase a little bit;
			 * like 10 here below so that to avoid boundary check error 
			 * in class localmemscr in general folder
			 */
			UInt setMemSize(UInt maxL) const {
				return engineInfor.getMaxMemLen(maxL) + 10;
			}; 

			/**
			 * return the number of additional atom block results used in block matrix list
			 * see the block matrix list.cpp, this is to define the parameter of incLengthFactor
			 * in the constructor 
			 *
			 * remember, each digestion for K we have four matrix; for J we have two
			 */
			UInt getNAdditionalAtomBlockResults() const { 
				UInt moreMtrx = 6;
				if (onlyDoK(intJob)) moreMtrx = 4;
				if (onlyDoJ(intJob)) moreMtrx = 2;
				return moreMtrx;
			};

			/**
			 * return the single integral infor for the given LCode
			 */
			const SingleIntegralInfor& getIntegralInfor(const LInt& code) const{
				return engineInfor.getIntegralInfor(code);
			};
	};
}

#endif

