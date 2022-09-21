/**
 * \file    diiscontroller.h
 * \brief   mornitoring the DIIS like procedure convergence progress
 * \author  Fenglai Liu 
 */
#ifndef DIISCONTROLLER_H
#define DIISCONTROLLER_H
#include "libgen.h"
#include "scfmacro.h"
#include "scfdiis.h"
#include "scalablevec.h"

namespace scfparam {
	class SCFParam;
}

namespace scfconv {
	class SCFConv;
}

namespace molecule {
	class Molecule;
}

namespace globalinfor {
	class GlobalInfor;
}

namespace denmtrx {
	class DenMtrx;
}

namespace spinmatrix{
	class SpinMatrix;
}

namespace oneemtrx {
	class OneEMtrx;
}

namespace diiscontroller {

	using namespace scfdiis;
	using namespace scfmacro;
	using namespace scfparam;
	using namespace scfconv;
	using namespace globalinfor;
	using namespace denmtrx;
	using namespace spinmatrix;
	using namespace oneemtrx;

	/**
	 * this class is used to mornitor the SCF convergence proress
	 * for DIIS type procedure, including DIIS, ADIIS/EDIIS 
	 *
	 * in general, we are able to perform single SCF algorithm like DIIS, EDIIS etc.
	 * alone, or to combine two algorithms together like EDIIS+DIIS, ADIIS+DIIS.
	 */
	class DIISController {

		private:

			//
			// debugging options
			//
			// for controller itself, 0 is normal print, 1 is for debugging print
			//
			bool printSCFEADIIS;              ///< whether print out scf EDIIS/ADIIS content if in use?
			bool printSCFDIIS;                ///< whether print out scf DIIS content if in use?
			UInt printController;             ///< level to print the given controller information?

			//
			// problem-solution checking choices
			//
			UInt nOscillationCycles;          ///< how many SCF cycles we used to check for oscillation status?
			Double oscillationThresh;         ///< threshold of energy for checking oscillation status
			///
			/// threshold of error for checking SCF is not converged 
			/// but in oscillation status
			///
			Double oscillationErrorThresh;    

			//
			// DIIS/EDIIS/ADIIS procedure information
			//
			UInt maxSCFCycles;                ///< we keep a copy of maximum SCF cycle data
			UInt currentJob;                  ///< current DIIS job
			UInt maxPrecEADIIS;               ///< precision digit for deriving result mo in EDIIS/ADIIS
			bool useExpensiveEADIISSolver;    ///< whether we use the expensive EADIIS solver?
			UInt maxEADIISSpace;              ///< the maximum dimension of vector/matrix in EDIIS/ADIIS
			UInt maxDIISSpace;                ///< the maximum dimension of vector/matrix in normal DIIS
			UInt DIISSpaceIncNumber;          ///< increment number on DIIS space if to enlarge DIIS space 
			UInt EADIISSpaceIncNumber;        ///< increment number on EADIIS space if to enlarge EADIIS space
			UInt scfVecSpaceChoice;           ///< the choice to use energy or diis err. forming scf index space 
			Double diisError;                 ///< current DIIS error
			UIntVec indexArray;               ///< scf index array used for DIIS/EADIIS for extropolation
			UIntVec historyIndexArray;        ///< a history record of indexArray for the given type of job
			
			///
			/// record of each history index array, each index array has
			/// two data fields. first one is beginning index in historyIndexArray,
			/// second one is the length of that given indexArray
			///
			UIntVec historyIndexRecord;       

			///
			/// this is the result of coefficients for extrapolation/interpolation
			///
			DoubleVec coefs;

			///
			/// available scf index space for picking up the result scf index 
			/// this is used by all of DIIS like methods
			///
			UIntVec scfIndexSpace;

			///
			/// normal Pulay's diis procedure
			/// we keep it inside because the error matrix is update per each
			/// SCF cycle.
			///
			SCFDIIS  scfDIIS;               

			///
			/// update the DIIS or ADIIS/EDIIS scf index
			/// according to the input energy/diis_error vector
			///
			/// the algorithm is to pick up the lowest energy or diis error 
			/// for the scf cycles as appeared in scfIndexSpace.
			///
			/// the result is scfIndexVec, the old indexArray is not changed yet
			///
			void updateIndexSpace(const SCFConv& scfconv, UIntVec& scfIndexVec) const;

			///
			/// after the updateIndexSpace, the newly generated scf vector space needed
			/// to be check to see what's kind of problem we may have
			///
			/// it returns the status defined in scfconvjob.h
			/// 
			UInt checkSCFIndexSpace(const SCFConv& scfconv, const UIntVec& newSCFVecSpace) const;

			///
			/// push the result new scf vector space into archive
			///
			void updateSpaceArchive(const UIntVec& newSCFVecSpace);

			///
			/// for restarting the algorithm (or switch to a new algorithm) we need to 
			/// reset the scf index space. In this function we will pick up the one 
			/// corresponding to lowest energy or diis error(the first SCF iteration is 
			/// not considered, see the comments inside the function) as first element
			/// after we restart the scf index space.
			///
			void resetSCFIndexSpace(const SCFConv& scfconv);

		public:

			/**
			 * constructor
			 */
			DIISController(const GlobalInfor& infor, const SCFParam& par);

			/**
			 * destructor
			 */
			~DIISController() { };

			/**
			 * print out the error matrix as well as coefficients
			 * for debug purpose
			 */
			void print() const;

			/**
			 * get current SCF job for DIIS like procedure
			 */
			UInt getJob() const { return currentJob; }; 

			/**
			 * set the job - from the upper class of SCFConv
			 * also all of index array etc. information will be cleared
			 */
			void setJob(UInt job) { 
				currentJob = job; 
				indexArray.clear();
				historyIndexArray.clear();
				historyIndexRecord.clear();
			};

			/**
			 * get the current scf index array information for the given job 
			 */
			const UIntVec& getIndexArray() const { return indexArray; }; 

			/**
			 * get the result coefficients array
			 */
			const DoubleVec& getCoefsArray() const { return coefs; }; 

			/**
			 * return the precision digit used for EDIIS/ADIIS procedure
			 */
			UInt getMaxEADIISPrec() const { return maxPrecEADIIS; };  

			/**
			 * get the error
			 */
			Double DIISError() const { return diisError; };

			/**
			 * whether you want more precise EADIIS result coefficients?
			 */
			bool withMorePreciseCoef() const { return useExpensiveEADIISSolver;  }; 

			/**
			 * get the length of terms for DIIS/EDIIS/ADIIS
			 */
			UInt getNTerms() const { return indexArray.size(); };

			/**
			 * as required by the top class of conv, we are going to reset 
			 * the index space array inside. this is a solution for DIIS
			 * like procesure to solve the scf convergence problem.
			 *
			 * this is used for switching to new SCF convergence method,
			 * or restarting the scf index array
			 *
			 * we note that this function does not provide a new scf index vector space,
			 * it just clean all of old ones
			 */
			void resetIndexSpace(const SCFConv& conv, UInt solution) {
				indexArray.clear();
				resetSCFIndexSpace(conv);

				// now depending on the solution required, let's see 
				// whether we need to clean the history
				if (solution == SWITCH_TO_NEW_SCF_CONV_METHOD) {
					historyIndexArray.clear();
					historyIndexRecord.clear();
				}
			};

			/**
			 * as required by the top class of conv, we are going to enlarge
			 * the index space array inside. this is a solution for DIIS
			 * like procesure to solve the scf convergence problem.
			 *
			 * we note that this function does not provide a new scf index 
			 * vector space, it just change the setting of the space limit
			 */
			void enlargeSCFIndexSpace() {
				if (useADIIS(currentJob) || useEDIIS(currentJob)) maxEADIISSpace += EADIISSpaceIncNumber;
				if (useDIIS(currentJob)) maxDIISSpace += DIISSpaceIncNumber;
			};

			/**
			 * this function just derive the scf index space 
			 * 
			 * firstly it derive the result;
			 * 
			 * then it will do the checking to see whether we have any problems(return value);
			 * 
			 * if everything going good then it will updat ethe archive then return good status
			 *
			 * if the saveResult is true, then even if the result status is not good, we will
			 * save the result index array and use it in later coefficients derivation
			 */
			UInt derivIndexSpace(const SCFConv& conv, bool saveResult); 

			/**
			 * this function will push the current scf index into the scf index space
			 * it should be called before deriving the scf index array in scfconv
			 */
			void incrementSCFIndexSpace(const SCFConv& conv);

			/**
			 * update DIIS information, this is done for every scf cycle
			 */
			void updateDIISInfor(const DenMtrx& denMtrx, const OneEMtrx& oneEMtrx, 
					const SpinMatrix& fock);

			/**
			 * this is the driver function to derive the new fock matrix
			 * through extrapolation/interpolation with DIIS like algorithm
			 *
			 * based on the derived scf index array, we will compute the 
			 * coefficients of the history fock matrix. Then finally we will
			 * compute the new Fock matrix
			 */
			void convSCF(const DenMtrx& denMtrx, const SCFConv& conv, SpinMatrix& fock);

	};
}

#endif

