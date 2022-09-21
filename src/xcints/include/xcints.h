/**
 * \file    xcints.h
 * \brief   top interface class for doing DFT work
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef XCINTS_H
#define XCINTS_H
#include "tbb/tbb.h"
#include "libgen.h"
#include "scalablevec.h"
#include "xcintsinfor.h"
#include "xcfunc.h"
#include "xcvar.h"
#include "molegrids.h"
#include "shellsize.h"
#include "spinmatrix.h"

using namespace xcintsinfor;
using namespace xcfunc;
using namespace xcvar;
using namespace molegrids;
using namespace shellsize;
using namespace spinmatrix;

namespace xcenergyinfor {
	class XCEnergyInfor; 
}

namespace oddelec {
	class ODDElec;
}

namespace atomdenmtrx {
	class AtomDenMtrx;
}

namespace denmtrx {
	class DenMtrx;
}

namespace globalinfor {
	class GlobalInfor;
}

namespace gintsinfor {
	class GIntsInfor;
}

namespace shell {
	class MolShell;
}

namespace sigshellpairinfor {
	class SigMolShellPairInfor;
};

namespace batchgrid {
  class BatchGrid;
};

namespace batchbasis {
  class BatchBasis;
};

namespace sigatombasis {
  class SigAtomBasis;
};

namespace batchvar {
  class BatchVar;
};

namespace xcints {

	using namespace xcenergyinfor; 
	using namespace sigshellpairinfor; 
	using namespace atomdenmtrx;
	using namespace globalinfor; 
	using namespace gintsinfor;
	using namespace denmtrx;
	using namespace oddelec; 
	using namespace shell;
  using namespace batchgrid;
  using namespace batchbasis;
  using namespace sigatombasis;
  using namespace batchvar;

	/**
	 * \class TBB_XCIntsMtrx
	 *
	 * Here in this class, each threads will compute a portion of xc matrix etc.
	 * then in using the mutex it will add the result into the output global results
	 */
	class TBB_XCIntsMatrix {

		private:

			//
			// reference for the input 
			//
			const MolShell& ms;                    ///< input shell
			const SigMolShellPairInfor& spData;    ///< shell pair information data
			const MoleGrids& molGrids;             ///< molecular grids
			const MolShellSize& molShellSize;      ///< molecular shell size data
			const XCFunc& xcfunc;                  ///< xc functional information
			const XCVar& xcvar;                    ///< xc variable information
			const XCIntJobInfor& infor;            ///< information center for real calculation
			const DenMtrx& denMtrx;                ///< density matrix
			const DoubleVec& cartDen;              ///< Cartesian form density matrix in significant shell pair form
			const AtomDenMtrx& atomDenMtrx;        ///< density matrix data for free atom

			//
			// the output data reference
			//
			Double& exc;                           ///< xc energy
			Double& alphaEx;                       ///< alpha exchange energy 
			Double& betaEx;                        ///< beta exchange energy 
			SpinMatrix& xcMtrx;                    ///< xc result matrix/ K or J+K matrix
			Mtrx& jMtrx;                           ///< Coulomb matrix
			DoubleVec& s2VecJMtrx;                 ///< the s2 form of J matrix part

			//
			// possible odd electron data
			//
			ODDElec& oddElec;                      ///< odd electron population results

			//
			// possible xc energy infor (XC energy decomposition data)
			//
			XCEnergyInfor& xcEnergyProfile;         ///< XC energy decomposition data

		public:

			/**
			 * constructor 
			 *
			 * actually we pass all of input/output reference inside
			 */
			TBB_XCIntsMatrix(const MolShell& ms0, const SigMolShellPairInfor& spData0, 
					const MoleGrids& molGrids0, const MolShellSize& molShellSize0, const XCFunc& xcfunc0, 
					const XCVar& xcvar0, const XCIntJobInfor& infor0, 
					const DenMtrx& den, const DoubleVec& cartDen0, 
					const AtomDenMtrx& atomDenMtrx0,
					Double& exc0, Double& alphaEx0, Double& betaEx0,
					SpinMatrix& xcMtrx0, Mtrx& jMtrx0, DoubleVec& s2VecJMtrx0, 
					ODDElec& oddElec0, XCEnergyInfor& xcEnergyProfile0):ms(ms0),spData(spData0),
			molGrids(molGrids0),molShellSize(molShellSize0),
			xcfunc(xcfunc0),xcvar(xcvar0),infor(infor0),denMtrx(den),cartDen(cartDen0),
			atomDenMtrx(atomDenMtrx0),exc(exc0),alphaEx(alphaEx0),betaEx(betaEx0),
			xcMtrx(xcMtrx0),jMtrx(jMtrx0),s2VecJMtrx(s2VecJMtrx0),oddElec(oddElec0),
			xcEnergyProfile(xcEnergyProfile0) { };

			/**
			 * destructor
			 */
			~TBB_XCIntsMatrix() { };

			/**
			 * functional operator to generate the results
			 * this is used as parallel_for together with the mutex
			 */
			void operator()(const blocked_range<UInt>& r) const; 

			/**
			 * splitting constructor used by TBB
			 * we copy all of reference data 
			 */
			TBB_XCIntsMatrix(const TBB_XCIntsMatrix& tbb_xcints, split):ms(tbb_xcints.ms),
			spData(tbb_xcints.spData),molGrids(tbb_xcints.molGrids),molShellSize(tbb_xcints.molShellSize),
			xcfunc(tbb_xcints.xcfunc),xcvar(tbb_xcints.xcvar),infor(tbb_xcints.infor),
			denMtrx(tbb_xcints.denMtrx),cartDen(tbb_xcints.cartDen),atomDenMtrx(tbb_xcints.atomDenMtrx),
			exc(tbb_xcints.exc),alphaEx(tbb_xcints.alphaEx),betaEx(tbb_xcints.betaEx),xcMtrx(tbb_xcints.xcMtrx),
			jMtrx(tbb_xcints.jMtrx),s2VecJMtrx(tbb_xcints.s2VecJMtrx),oddElec(tbb_xcints.oddElec), 
			xcEnergyProfile(tbb_xcints.xcEnergyProfile) { };
	};

	/**
	 * top interface class for doing DFT work 
	 */
	class XCInts {

		protected:

			//
			// energy result
			//
			Double exc;             ///< energy of exchange + correlation

			//
			// class members
			//
			XCIntJobInfor infor;    ///< information of xcints
			XCFunc      xcfunc;     ///< xcfunctional information
			XCVar       xcvar;      ///< xc variable information
			MolShellSize size;      ///< molecular shell size data
			MoleGrids   molGrids;   ///< atom and molecule grids setting

			//
			// result matrix of XCInts
			//
			SpinMatrix  xcMtrx;     ///< it may be used to store the xc matrix, K matrix or J+K matrix
			Mtrx        jMtrx;      ///< if numerical J is required and doing seperately, we do it here

		public:

			/**
			 * construtor for the XCInts - only for information processing
			 * \param gInfor      global setting information
			 * \param xcfunc0     functional information
			 * \param ms          shell data
			 * \param xcinfor     DFT information center
			 * \param gintsinfor  analytical integrals(gints) information center
			 * \param job         job name
			 * \param joborder    the order to do energy/gradient calculation
			 */
			XCInts(const GlobalInfor& gInfor, const XCFunc& xcfunc0, const MolShell& ms, 
					const XCIntsInfor& xcinfor, const GIntsInfor& gintsInfor, const UInt& job, 
					const UInt& jobOrder);

			/**
			 * destructor
			 */
			~XCInts() { };

			/**
			 * top interface to form the XC Matrix in normal way
			 * \param   molShell  shell data
			 * \param   denMtrx   density matrix
			 * \return  result    result matrix waiting to be updated
			 */
			void doXCMtrx(const MolShell& ms, const DenMtrx& denMtrx, SpinMatrix& result, 
					bool printTiming, bool printMatrix);

			/**
			 * top interface to form the JK Matrix in numerical way
			 * \param   molShell  shell data
			 * \param   denMtrx   density matrix
			 * \return  result    result matrix waiting to be updated
			 */
			void doJKMtrx(const MolShell& ms, DenMtrx& denMtrx, SpinMatrix& result,
					bool printTiming, bool printMatrix);

			/**
			 * whether you are saving the result to disk?
			 */
			void saveToDisk(const string& dataPath, const SpinMatrix& result) const;

			/**
			 * return the exc
			 */
			Double getExc() const { return exc; };
	};

}

#endif

