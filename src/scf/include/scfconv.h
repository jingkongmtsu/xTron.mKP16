/**
 * \file    scfconv.h
 * \brief   monitoring the SCF convergence progress
 * \author  Fenglai Liu 
 */
#ifndef SCFCONV_H
#define SCFCONV_H
#include "libgen.h"
#include "globalinfor.h"
#include "histdataman.h"
#include "scfmacro.h"
#include "scfconvsol.h"
#include "diiscontroller.h"
#include "scalablevec.h"

namespace scfparam {
	class SCFParam;
}

namespace spinmatrix {
	class SpinMatrix;
}

namespace denmtrx {
	class DenMtrx;
}

namespace oneemtrx {
	class OneEMtrx;
}

namespace scfconv {

	using namespace scfmacro;
	using namespace globalinfor;
	using namespace histdataman;
	using namespace diiscontroller;
	using namespace scfparam;
	using namespace scfconvsol;
	using namespace spinmatrix;
	using namespace denmtrx;
	using namespace oneemtrx;

	/**
	 * this class is used to mornitor the SCF convergence proress
	 * it will call the DIIS, ADIIS/EDIIS etc. procedure to 
	 * perform the real SCF convergence process
	 */
	class SCFConv {

		private:

			//
			// general input data 
			//
			GlobalInfor infor;                ///< save a copy of global information
			Double  convCriteria;             ///< input convergence criteria 
			UInt maxSCFCycles;                ///< number of maximum scf cycles
			UInt nSpin;                       ///< keep a copy of spin information 
			UInt scfIndex;                    ///< current SCF index

			//
			// SCF convergence control information
			//
			UInt jobType;                     ///< current scf algorithm
			UInt firstJob;                    ///< user defined first scf algorithm
			UInt secondJob;                   ///< user defined second scf algorithm
			Double diisError;                 ///< current diis error
			Double secondJobSwitchOnLimit;    ///< value to determine when to switch second job
			UInt algorithmReplaceDIIS;        ///< if switch to new algorithm, what is to replace diis?
			UInt algorithmReplaceEDIIS;       ///< if switch to new algorithm, what is to replace ediis?
			UInt algorithmReplaceADIIS;       ///< if switch to new algorithm, what is to replace adiis?

			///
			/// sometimes during SCF we may have linear dependence problem,
			/// for example; in the DIIS the error matrix could be singular due to
			/// linear dependence. Therefore we set up this threshold value used
			/// to detect the linear dependence of these situations
			///
			Double scfLinearThresh; 

			//
			// result section of scf convergence 
			//
			DoubleVec energyList;             ///< storing each SCF cycle energy data
			DoubleVec diisErrors;             ///< diis error for the whole SCF

			//
			// SCF iteration status and possible problem and solution
			// all of status and solution information defined in scfmacro.h
			// we also keep a record of job information
			//
			UIntVec scfStatus;                ///< status for each SCF cycle 
			UIntVec scfSolution;              ///< solution for the problem status in each SCF cycle
			UIntVec jobInfor;                 ///< each SCF interation the SCF convergence algorithm

			//
			// history data section, will be used by the following 
			// SCF convergence algorithm
			//
			HistDataMan  FockMtrxList;        ///< Fock matrix history list
			HistDataMan  denMtrxList;         ///< density matrix history list

			///
			/// controller for algorithms like DIIS/ADIIS/EDIIS
			///
			DIISController  diisController; 

			///
			/// class to store the scf convergence problem-solution 
			///
			SCFConvSol scfConvSol;

			///
			/// this function is the driver the derive the scf index vector 
			/// for interpolating/extropolating the result fock matrix
			///
			void deriveIndexSpace(); 

			///
			/// reset the current job
			/// because of difficulty in current SCF algorithm, we 
			/// need to switch to anotherjob
			///
			/// the new algorithm will be derived according to the setting,
			/// for example; the setting of algorithmReplaceADIIS will set
			/// the algorithm to replace ADIIS
			///
			void resetJob();

			///
			/// check the available SCF convergence solutions
			/// if it's good then we use it, return true; else return false
			/// and we will go to next possible solution
			///
			bool checkUpSCFSolutions(UInt status, UInt solution) const;


		public:

			/**
			 * constructor
			 */
			SCFConv(const GlobalInfor& infor, const SCFParam& par);
			SCFConv(const GlobalInfor& infor, const SCFParam& par, const string& scfalgorithm);

			/**
			 * destructor
			 */
			~SCFConv() { };

			/**
			 * print out the error matrix as well as coefficients
			 * for debug purpose
			 */
			void print() const;

			////////////////////////////////////////////////////////////////
			//             member data fetching functions                 //
			////////////////////////////////////////////////////////////////

			/**
			 * return the number of spin states
			 */
			UInt getNSpin() const { return nSpin; };

			/**
			 * get current SCF job
			 */
			UInt getJob() const { return jobType; }; 

			/**
			 * get the current scf index of the job
			 */
			UInt getCurrentSCFIndex() const {
				return scfIndex;
			};

			/**
			 * whether this is first cycle of SCF?
			 */
			bool firstSCFCycle() const {
				return (scfIndex==0);
			};

			/**
			 * get current iteration energy
			 */
			Double getCurrentEnergy() const { 
				UInt scfIndex = getCurrentSCFIndex();
				return energyList[scfIndex]; 
			};

			/**
			 * get the history of density matrix
			 */
			const HistDataMan& getHistDenData() const {
				return denMtrxList;    
			};

			/**
			 * get the history of fock matrix
			 */
			const HistDataMan& getHistFockData() const {
				return FockMtrxList;    
			};

			/**
			 * get the history of energy data
			 */
			const DoubleVec& getHistEnergyData() const {
				return energyList;
			};

			/**
			 * get the history of diis error data
			 */
			const DoubleVec& getHistDIISErrorData() const {
				return diisErrors;    
			};

			/**
			 * get the history of energy data, just for a single index
			 */
			Double getHistEnergy(UInt i) const { return energyList[i]; };

			/**
			 * get the history of diis error data, just for a single index
			 */
			Double getHistDIISError(UInt i) const { return diisErrors[i]; };

			/**
			 * get the error
			 */
			Double DIISError() const { return diisError; };

			/**
			 * return the SCF linear dependence threshold value
			 */
			Double getSCFLinearThresh() const { return scfLinearThresh; };

			/**
			 * return the current lowest energy and the corresponding scf index
			 */
			void getLowestEnergyInfor(Double& E, UInt& index) const;

			////////////////////////////////////////////////////////////////
			//                 main functions of SCFConv                  //
			////////////////////////////////////////////////////////////////

			/**
			 * driver function to converge SCF
			 * \param currentIndex  current SCF iteration index
			 * \param E             total energy this iteration
			 * \param denMtrx       density matrix data
			 * \param oneEMtrx      archive for one electron matrix data
			 * \return fock         the interpolation/extrapolation will be performed for fock 
			 */
			void convSCF(const UInt& currentIndex, const Double& E, 
					const DenMtrx& denMtrx, const OneEMtrx& oneEMtrx, SpinMatrix& fock);

			/**
			 * is it no need to perform any SCF convergence algorithm?
			 */
			bool noNeedSCFConv() const {

				// this case we only do one cycle SCF
				if (maxSCFCycles == 1) return true;

				// the rest of other cases
				return false;
			};

			/**
			 * whether we are reaching the end of SCF cycle?
			 */
			bool scfEnds() const {

				// this case we do not do scf convergence 
				if (noNeedSCFConv()) return true;

				// this case we already get a converged result
				if (diisError<convCriteria) return true;

				// this case we meet the maximum SCF cycles requirement
				if (scfIndex>=maxSCFCycles-1) return true;

				// the rest of other cases
				return false;
			};

	};
}

#endif

