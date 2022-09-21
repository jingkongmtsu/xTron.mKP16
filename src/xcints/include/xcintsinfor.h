/**
 * \file    xcintsinfor.h
 * \brief   guiding class for calculations in xcints
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef XCINTSINFOR_H
#define XCINTSINFOR_H
#include "libgen.h"
#include "gintsinfor.h"
#include "integraljobs.h"
#include "globalinfor.h"
#include<string>

namespace molecule {
	class Molecule;
}

namespace xcintsinfor {

	using namespace molecule;
	using namespace globalinfor;
	using namespace gintsinfor;
	using namespace integraljobs;
	using namespace std;

	//
	// the step function type in doing Becke partition weights
	//
	const UInt BECKE_ORI_STEP_FUNCTION  = 1;   ///< original Becke step function
	const UInt STRATMANN_STEP_FUNCTION  = 2;   ///< stratmann's step function

	/**
	 * this class is used to collect/parse information about the foundamental DFT jobs in xcints
	 */
	class XCIntsInfor {

		protected:

			//
			// general settings for the XCInts job
			// 
			bool doPartitionWeights;     ///< do we do partition weights?
			bool insigGridsCheck;        ///< check that whether the grids are negligible
			bool doOddElectron;          ///< do we do odd electron population calculation?
			bool doXCEnergyProfing;      ///< do we do energy profiling work using the xcenergyinfor class? 
			UInt pruneMethod;            ///< if we use prune grid, this is the method name
			UInt stepFunctionType;       ///< which step function we take in doing partition weights? BECKE_ORI_STEP_FUNCTION?
			UInt gridChoice;             ///< what kind of grid we use, standard, coarse or fine? only for Baker pruned grids
			UInt nRad;                   ///< number of radial points for un-pruned grids
			UInt nAng;                   ///< number of angular points for un-pruned grids
			UInt nPtPerCPUThread;        ///< how many points per each CPU thread you expect?
			Double sigGridCheckThresh;   ///< threshold value to do significant grid check
			Double threshold;            ///< threshold value to determine the atom size etc.
			Double funcTol;              ///< tolerance value for determing small dft vars in functional derivatives calculation

			/*
			 * keep a copy of the spin information etc.
			 */
			UInt nSpin;                  ///< number of spin statess
			bool isCloseShell;           ///< whether this is close shell

			/*
			 * odd electron population parameters
			 */
			Double oddElecPar1;          ///< odd electron population parameter 1
			Double oddElecPar2;          ///< odd electron population parameter 2

			//JK. nDen here is the same as nSpin plus when nBeta==0 (the old single electron case.)
			UInt nDen;



			// debug printing choice
			UInt xcvar_debug;            ///< doing debug choice for xcvar
			UInt molegrids_debug;        ///< doing debug choice for atom grids as well as molegrids
			UInt sigatombasis_debug;     ///< doing debug choice for sigatombasis
			UInt batchgrid_debug;        ///< doing debug choice for batch grid
			UInt batchbasis_debug;       ///< doing debug choice for batch basis
			UInt batchvar_debug;         ///< doing debug choice for batch var
			UInt batchfunc_debug;        ///< doing debug choice for batch functional
			UInt batchxcmtrx_debug;      ///< doing debug choice for batch xc matrix
			UInt espints_debug;          ///< doing debug choice for espints
			UInt exrho_debug;            ///< doing debug choice for exrho

			// time printing choice
			bool printTimingForESP;      ///< whether to print timing for ESP integrals

		public:

			/**
			 * constructor
			 * \param  infor     global information
			 * \param  mol       molecule data
			 */
			XCIntsInfor(const GlobalInfor& ginfor, const Molecule& mol);

			///
			/// destrucotr
			///
			~XCIntsInfor() { };

			///
			/// debug print
			///
			void debugPrint() const;

			/**
			 * return the spin information
			 */
			UInt getNSpin() const { return nSpin; };

			/**
			 * is it close shell?
			 */
			bool closeShell() const { return isCloseShell; };

			///
			/// get the debug level of molgrids
			///
			UInt debugLevelMoleGrids() const { return molegrids_debug; };

			///
			/// do we do partitional weights? 
			///
			bool doParWeights() const { return doPartitionWeights; };

			///
			/// do we do odd electron population?
			///
			bool doOddElec() const { return doOddElectron; };

			///
			/// do we do xc energy profile/decomposition?
			///
			bool doXCEnergyProfile() const { return doXCEnergyProfing; };   

			///
			/// which step function for partition weights?
			///
			UInt getStepFunc() const { return stepFunctionType; };

			///
			/// get the number of points for each cpu thread
			///
			UInt getNPtPerCPUThread() const { return nPtPerCPUThread; };

			///
			/// do insig grid check?
			///
			bool doSigGridCheck() const { return insigGridsCheck; };

			///
			/// return prune grid method
			///
			UInt getPruneGridMethod() const { return pruneMethod; };

			///
			/// return grid choice
			///
			UInt getGridChoice() const { return gridChoice; };

			///
			/// return number of radial points
			///
			UInt getNAng() const { return nAng; };

			///
			/// return number of radial points
			///
			UInt getNRad() const { return nRad; };

			///
			/// return threshold value for determing atom size etc.
			///
			Double getThresh() const { return threshold; };

			///
			/// return threshold value for insignificant grid check
			///
			Double getInsigGridCheckThresh() const { return sigGridCheckThresh; };

			///
			/// get the tolerance of functional derivatives
			///
			Double getFuncTol() const { return funcTol; };

			///
			/// return the odd electron population parameters
			///
			Double getOddElecPar(UInt index) const { 
				if (index == 0) return oddElecPar1;; 
				return oddElecPar2;
			};

			///
			/// whether print out the esp part timing?
			///
			bool printESPIntsTiming() const { return printTimingForESP; };

			UInt getNDen() const { return nDen; };
	};

	/**
	 * based on the XCIntsInfor, this class will have the specific job information
	 */
	class XCIntJobInfor : public XCIntsInfor {

		private:

			///
			/// we save a copy of global infor
			///
			GlobalInfor globInfor;

			//
			// we need a copy of gintsinfor, so that to guide the espints calculation
			//
			GIntsInfor ginfor;           ///< a copy of gintsinfor

			//
			// job section
			//
			UInt jobName;                ///< the job for numerical integrals
			UInt order;                  ///< the job order, 0 is energy; 1 is for gradient; 2 for second derivatives 

			//
			// parallel related setting
			//
			bool useMultThreads;         ///< whether we do multiple threads for xcints (either use co-processor or CPU)
			UInt nCPUThreads;            ///< number of CPU threads for the xcints part

			//
			// batch size determination
			//
			UInt initBatchSize;          ///< the initial batch size, this is the guideline for the final batch size

		public:

			/**
			 * constructor
			 *
			 * this constructor will take all of XC part information from the input XCIntsInfor,
			 * so it does not do the parsing job anymore
			 */
			XCIntJobInfor(const GlobalInfor& globInfor0, const GIntsInfor& ginfor0, 
					const XCIntsInfor& xcinfor, const UInt& job, const UInt& jobOrder);

			///
			/// destrucotr
			///
			~XCIntJobInfor() { };

			///
			/// debug printing function
			///
			void print() const;

			/**
			 * whether we do debug for the given module
			 *
			 * for the multi-threading module, we can not print out the debugging infor
			 * in that case it will be a mess of the print out
			 */
			bool doDebug(const string& name) const;

			///
			/// return the gintsinfor information
			///
			const GIntsInfor& getGIntsInfor() const { return ginfor; };

			///
			/// return the copy of global infor
			///
			const GlobalInfor& getGlobalInfor() const { return globInfor; };

			///
			/// whether we use multi-threading mode?
			///
			bool useMultiThreads() const { return useMultThreads; };

			///
			/// get the number of threads for cpu
			///
			UInt getNCPUThreads() const { return nCPUThreads;  }; 

			///
			/// get the initial batch size
			///
			UInt getInitBatchSize() const { return initBatchSize;  }; 

			///
			/// return the integral job
			///
			UInt getIntJob() const { return jobName; };

			///
			/// return the job order
			///
			UInt getJobOrder() const { return order; };

			///
			/// whether we do numerical K job?
			///
			bool doNumericalK() const {
				if (doK(jobName)) return true;
				return false;
			};

			///
			/// whether we do numerical J job?
			///
			bool doNumericalJ() const {
				if (doJ(jobName)) return true;
				return false;
			};

			///
			/// get the number of energy density for numerical J/K
			///
			UInt getJKVarNum() const {
				UInt nVar = 0;
				if (doK(jobName)) {
					nVar += 1;
					if (getNSpin() == 2) nVar += 1; 
				}
				if (doJ(jobName)) nVar += 1;
				return nVar;
			};
	};
}

#endif

