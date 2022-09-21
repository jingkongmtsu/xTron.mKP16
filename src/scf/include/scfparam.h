/**
 * \file    scfparam.h
 * \brief   collecting SCF parameters 
 * \author  Fenglai Liu 
 */
#ifndef SCFPARAM_H
#define SCFPARAM_H
#include "libgen.h"
#include "xcfunc.h"
#include "gintsinfor.h"
#include "xcintsinfor.h"
#include "globalinfor.h"
using namespace globalinfor;
using namespace xcfunc;
using namespace gintsinfor;
using namespace xcintsinfor;

namespace molecule {
	class Molecule;
}

namespace scfparam {

	using namespace molecule;

	//
	// here we name the way to generate scf guess
	//
	const UInt SCF_GUESS_READ   = 100;  ///< reading mo from file then form density matrix guesss
	const UInt SCF_GUESS_CORE   = 101;  ///< form the core guess
	const UInt SCF_GUESS_ATOMIC = 102;  ///< form the guess with atomic density matrix


	/**
	 * this class is used to collect/parse information about the SCF process
	 */
	class SCFParam {

		private:

			/*
			 * A COPY OF GLOBAL INFORMATION
			 */
			GlobalInfor gInfor;     ///< global settings

			/*
			 * SCF ANALYTICAL INTEGRAL SECTION
			 */
			GIntsInfor  gIntInfor;  ///< analytical Gaussian integrals job information

			/*
			 * SCF XC INTEGRAL SECTION
			 */
			XCIntsInfor  xcIntInfor;  ///< xc integrals job information

			/*
			 * SCF FUNCTIONAL SECTION
			 */
			XCFunc  xcFunc;         ///< xc functional information

			/*
			 * SCF GUESS SECTION
			 */
			UInt scfGuessMethod;    ///< the guess method, see above definition
			string moPath;          ///< if reading from the file, we need the definition

			/**
			 * density matrix data write/read
			 */
			string denMtrxPath;     ///< the path where to save the result density matrix data
			string initDenMtrxPath; ///< the path where to save the initial density matrix data
			string guessDenMtrxPath;///< the path where to read in the guess 

			/*
			 * SCF PROPERTY SECTION
			 */
			bool  withFile;         ///< do we use file to keep result etc.?
			bool  withPostSCF;      ///< do we do post SCF process? sometimes we need to disable the post SCF
			bool isCloseShell;      ///< whether the alpha and beta will share same MO?
			UInt nSpin;             ///< number of spin states in SCF process
			Double convCriteria;    ///< SCF convergence criteria
			UInt maxSCFCycles;      ///< maximum scf cycles
			UInt section;           ///< which geometry section it is

			/*
			 * SCF CONTROL SECTION
			 */
			bool doDipoleMOMPostSCF;   ///< whether you do dipole moment calculation in the post SCF process?
			bool doVDWPostSCF;         ///< whether we do vdw calculation in the post SCF process?

			/*
			 * SCF DEBUG SECTIOn
			 *
			 * we are trying to provide a detailed printing list for SCF progress
			 */
			bool keepSCFSilent;        ///< whether the SCF will be silent, do not print out anything?
			bool doJKSeparately;       ///< do you want to manually do exchange/Coulomb separately?
			bool printSCFTiming;       ///< do we print out the timing of SCF part?
			bool saveJKMtrx;           ///< whether you want to save JK matrix?
			bool saveXCMtrx;           ///< whether you want to save XC matrix?
			bool printInitDenMtrx;     ///< whether to print out the initial density matrix before SCF?
			bool printDenMtrx;         ///< whether to print out density matrix during SCF(not converged yet)?
			bool printResultDenMtrx;   ///< whether to print out the result density matrix at end of SCF?
			bool printCoreMtrx;        ///< whether to print out core matrix?
			bool printJKMtrx;          ///< whether to print out JK matrix?
			bool printXCMtrx;          ///< whether to print out XC matrix?
			bool printCalFockMtrx;     ///< whether to print out the calculated Fock matrix?
			bool printSCFFockMtrx;     ///< whether to print out the Fock matrix for deriving new density matrix?
			bool printSCFConv;         ///< whether to print out the scf conv class content?
			UInt printMO;              ///< the level to print out MO data: 1-only MO energy; 2-mo vectors; 3 all

		public:

			/**
			 * constructor
			 * \param  infor     global information
			 * \param  mol       molecule data
			 */
			SCFParam(const GlobalInfor& infor, const Molecule& mol);

			/**
			 * default destructor
			 */
			~SCFParam() { };

			/**
			 * an convinient member function to return the number 
			 * of spin states
			 */
			UInt getNSpin() const { return nSpin; };

			/**
			 * whether we perform the calculation in close shell way?
			 * Close shell means that whether the alpha electrons 
			 * and beta electrons are sharing the same part of MO
			 * if it is, then it's close shell; else it's open 
			 * shell
			 *
			 * this is different from the number of spin states
			 */
			bool closeShell() const { return isCloseShell; };

			/**
			 * get the section number
			 */
			UInt getSec() const { return section; };

			/**
			 * do we use File in SCF?
			 */
			bool useFile() const { return withFile; };

			/**
			 * overwrite the option of use file
			 * force to store everything into memory
			 *
			 * Use it with caution!
			 */
			void keepInMem() { withFile = false; };

			/**
			 * sometimes we need to disable the post scf to avoid trouble
			 * for example, if we want to run vdw calculation with xdm; since
			 * in xdm it runs the atom SCF calculation again then if it runs
			 * post SCF again then it begins to form an endless recursive loop;
			 * so SCF never ends
			 *
			 * for such cases we disable the post SCF process
			 */
			void disablePostSCF() { withPostSCF = false; };

			/**
			 * do we do post scf work?
			 */
			bool doPostSCFWork() const { return withPostSCF; };

			/**
			 * get the maximum SCF cycles
			 */
			UInt getMaxSCFCycles() const { return maxSCFCycles; };

			/**
			 * some times you may need to change the maximum SCF cycles by hard assignment
			 */
			void changeMaxSCFCycles(UInt newVal) { maxSCFCycles = newVal; };

			/**
			 * in default we always like to have SCF print out
			 *
			 * however, for some cases like atom density matrix data generation,
			 * you may like it to be in silence. We change it here
			 */
			void keepSCFSilence() { keepSCFSilent = true; };

			/**
			 * get the SCF convergence criteria
			 */
			Double getConvCriteria() const { return convCriteria; };

			/**
			 * return scf guess method
			 */
			UInt getGuess() const { return scfGuessMethod; };

			/**
			 * do we print out mo?
			 */
			bool printMOData() const { 
				if (! keepSCFSilent && printMO > 0) return true;
				return false;
			};

			///
			/// return the MO print option 
			///
			UInt getMOPrintOption() const { return printMO; };

			/**
			 * do you have the guess MO path defined?
			 */
			bool hasGuessMO() const {
				if (moPath == "NONE") return false;
				return true;
			};

			/**
			 * return scf guess mo path
			 */
			const string& getGuessMOPath() const { return moPath; };

			/**
			 * get the xcfunc infor
			 */
			const XCFunc& getXCFunc() const { return xcFunc; };

			/**
			 * some times you may also want to change the xc functional information
			 */
			void updateXCFunc(const XCFunc& xcfunc) {
				xcFunc = xcfunc;
			};

			/**
			 * get the gints infor
			 */
			const GIntsInfor& getGIntsInfor() const { return gIntInfor; };

			/**
			 * get the xcints infor
			 */
			const XCIntsInfor& getXCIntsInfor() const { return xcIntInfor; };

			/**
			 * get the copy of global infor
			 */
			const GlobalInfor& getGlobalInfor() const { return gInfor; };

			/**
			 * whether you do dipole moment calculation in the post scf process?
			 */
			bool doDipole() const { return doDipoleMOMPostSCF; };

			/**
			 * print out the jk matrix?
			 */
			bool printJKMatrix() const { 
				if (! keepSCFSilent) return printJKMtrx; 
				return false;
			};

			/**
			 * print out the core matrix?
			 */
			bool printCoreMatrix() const { 
				if (! keepSCFSilent) return printCoreMtrx; 
				return false;
			};

			/**
			 * print out the XC matrix?
			 */
			bool printXCMatrix() const { 
				if (! keepSCFSilent) return printXCMtrx; 
				return false;
			};

			/**
			 * print out the density matrix during SCF process?
			 */
			bool printDensityMatrixInSCF() const { 
				if (! keepSCFSilent) return printDenMtrx; 
				return false;
			};

			/**
			 * print out the initial density matrix?
			 */
			bool printInitialDensityMatrix() const { 
				if (! keepSCFSilent) return printInitDenMtrx; 
				return false;
			};

			/**
			 * print out the result density matrix?
			 */
			bool printResultDensityMatrix() const { 
				if (! keepSCFSilent) return printResultDenMtrx; 
				return false;
			};

			/**
			 * print out the calcualted fock matrix?
			 */
			bool printFockMatrix() const { 
				if (! keepSCFSilent) return printCalFockMtrx; 
				return false;
			};

			/**
			 * print out the derived fock matrix used for generating new MO?
			 */
			bool printSCFFockMatrix() const { 
				if (! keepSCFSilent) return printSCFFockMtrx; 
				return false;
			};

			/**
			 * print out scfconv class content?
			 */
			bool printSCFConvClass() const { 
				if (! keepSCFSilent) return printSCFConv; 
				return false;
			};

			/**
			 * whether SCF print out in silence?
			 */
			bool isSCFSilent() const { return keepSCFSilent; };

			/**
			 * do you want to do JK apart?
			 */
			bool doJKApart() const { return doJKSeparately; }    

			/**
			 * do you want to save JK matrix on disk?
			 */
			bool saveJKOnDisk() const { return saveJKMtrx; };

			/**
			 * do you want to save XC matrix on disk?
			 */
			bool saveXCOnDisk() const { return saveXCMtrx; };

			/**
			 * do you want to save the result density matrix on disk?
			 */
			bool saveResultDenMtrxOnDisk() const {
				if (denMtrxPath == "NONE") return false;
				return true;
			};

			/**
			 * do you want to save the initial density matrix on disk?
			 */
			bool saveInitDenMtrxOnDisk() const {
				if (initDenMtrxPath == "NONE") return false;
				return true;
			};

			/**
			 * return the path to save result density matrix
			 */
			const string& getResultDenMtrxPath() const { return denMtrxPath; };

			/**
			 * do you have the guess density matrix path defined?
			 */
			bool hasGuessDenMtrx() const {
				if (guessDenMtrxPath == "NONE") return false;
				return true;
			};

			/**
			 * return the path to read in the guess density matrix
			 */
			const string& getGuessDenMtrxPath() const { return guessDenMtrxPath; };

			/**
			 * return the path to save initial density matrix
			 */
			const string& getInitDenMtrxPath() const { return initDenMtrxPath; };

			/**
			 * do you want to print out SCF timing data?
			 */
			bool printSCFTimingData() const { return printSCFTiming; };

			/**
			 * do you want to do vdw in the post SCF?
			 */
			bool doVDWInPostSCF() const { return doVDWPostSCF; };

	};
}

#endif

