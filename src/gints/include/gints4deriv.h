/**
 * \file   gints4deriv.h
 * \brief  working class to handle the 4D analytical integrals deriv calculation
 * \author Fenglai Liu 
 */
#ifndef GINTS4DERIV_H
#define GINTS4DERIV_H
#include "tbb/tbb.h"
#include "shell.h"
#include "libgen.h"
#include "matrix.h"
#include "gintsinfor.h"
#include "scalablevec.h"
#include "sigshellpairinfor.h"

namespace denmtrx {
	class DenMtrx;
};

namespace molecule {
	class Molecule;
};

namespace shellpair {
	class ShellPair;
};

namespace denmtrxinfor {
	class DenMtrxInfor;
};

namespace symmjkdigest {
	class SymmJKDerivDigest; 
}

namespace gints4deriv {

	using namespace tbb;
	using namespace shell;
	using namespace matrix;
	using namespace denmtrx;
	using namespace molecule;
	using namespace shellpair;
	using namespace gintsinfor;
	using namespace denmtrxinfor;
	using namespace symmjkdigest;
	using namespace sigshellpairinfor;

	/**
	 * \class TBB_GInts4Deriv
	 *
	 * working class to do the derivatives calculation for derivatives on (ab|cd)
	 */
	class TBB_GInts4Deriv {

		private:

			//
			// basic data needed for calculation
			//
			bool sameShellPair;                   ///< whether the bra and ket are coming from same source?
			UInt nAtoms;                          ///< total number of atoms, get from input molecule information
			UInt nSpin;                           ///< number of spin states
			Double thresh;                        ///< threshold value for gints4d derivatives work
			Double spThresh;                      ///< threshold value for generating shell pairs
			UInt braBatchIndex;                   ///< the bra side batch index
			UInt ketBatchIndex;                   ///< the ket side batch index

			//
			// reference for the input 
			//
			const MolShell&     bra1;             ///< bra1 shell
			const MolShell&     bra2;             ///< bra2 shell
			const MolShell&     ket1;             ///< ket1 shell
			const MolShell&     ket2;             ///< ket2 shell
			const SigMolShellPairInfor& braInfo;  ///< shell pair information for bra side
			const SigMolShellPairInfor& ketInfo;  ///< shell pair information for ket side
			const GIntJobInfor& ginfor;           ///< information center for real calculation
			const Mtrx&         alphaDen;         ///< alpha density matrix
			const Mtrx&         betaDen;          ///< beta density matrix
			const DenMtrxInfor& denMtrxInfor;     ///< density matrix infor

			//
			// the output data
			//
			Mtrx derivData;                       ///< for deriv order = 1, (natoms,3)
			                                      ///< for deriv order = 2, (3*natoms,3*natoms)

		public:

			/**
			 * constructor 
			 * \param spStatus       whether bra and ket are same shell pair infor?
			 * \param braBatchIndex  the bra side sigAtomShellPair index
			 * \param ketBatchIndex  the ket side sigAtomShellPair index
			 * \param mol            molecule information
			 * \param ginfor         information center for gints job
			 * \param bra,ket        shell pair infor class for bra and ket side
			 * \param b1,b2,k1,k2    reference for shell data
			 * \param alpha,beta     input density matrices
			 * \param denMtrxInfor   the density matrix information
			 */
			TBB_GInts4Deriv(bool spStatus, const UInt& braBatchIndex0, const UInt& ketBatchIndex0,
					const Molecule& mol, const GIntJobInfor& ginfor0, const SigMolShellPairInfor& bra,
				  const SigMolShellPairInfor& ket, const MolShell& b1, const MolShell& b2, 
				  const MolShell& k1, const MolShell& k2, const Mtrx& alpha, const Mtrx& beta,
				  const DenMtrxInfor& denMtrxInfor0);
			
			/**
			 * destructor
			 */
			~TBB_GInts4Deriv() { };

			/**
			 * functional operator to generate the integral results
			 * this is used as parallel_reduce
			 */
			void operator()(const blocked_range<UInt>& r); 

			/**
			 * splitting constructor used by TBB
			 * we copy all of reference data except result
			 * which should be always empty initially
			 */
			TBB_GInts4Deriv(const TBB_GInts4Deriv& tbb_gints4deriv, split);

			/**
			 * for assemble results in threads together
			 */
			void join(TBB_GInts4Deriv& tbb_gints4deriv) {  
				derivData.add(tbb_gints4deriv.derivData);
			};

			/**
			 * get the result matrix 
			 */
			const Mtrx& getResultMatrix() const {
				return derivData;
			};
	};

	/**
	 * \class GInts4Deriv
	 * \brief top level class to mornitor 4D analytical integral derivatives calculation
	 */
	class GInts4Deriv {

		private:

			GIntJobInfor ginfor;               ///< information center for real calculation

			void do4DAODeriv(const SigMolShellPairInfor& braInfor,
					const SigMolShellPairInfor& ketInfor,
					const MolShell& bra1, const MolShell& bra2, 
					const MolShell& ket1, const MolShell& ket2, 
					const Mtrx& alphaDen, const Mtrx& betaDen, 
					const DenMtrxInfor& denMtrxInfor, Mtrx& result) const;

		public:

			////////////////////////////////////////////////////////////
			//                  constructor etc.                      //
			////////////////////////////////////////////////////////////

			/**
			 * constructor - build shell pair information based on one shell
			 */
			GInts4Deriv(const GIntsInfor& infor0, const UInt& job, 
					const UInt& order):ginfor(infor0,job,order) { };

			/**
			 * destructor
			 */
			~GInts4Deriv() { };

			///
			/// form JK matrix for the given one shell
			/// so all of shell pair information are same
			/// the input density matrix will be transformed into Cartesian
			/// inside, but finally they will be transformed back
			///
			void doJKDeriv(const MolShell& s, const Molecule& mol, DenMtrx& den, 
					Mtrx& result, bool printTiming, bool printResult); 
	};
}

#endif

