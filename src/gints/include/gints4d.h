/**
 * \file   gints4d.h
 * \brief  working class to handle the 4D analytical integrals calculation
 * \author Fenglai Liu 
 */
#ifndef GINTS4D_H
#define GINTS4D_H
#include "tbb/tbb.h"
#include "shell.h"
#include "libgen.h"
#include "matrix.h"
#include "gintsinfor.h"
#include "spinmatrix.h"
#include "scalablevec.h"
#include "sigshellpairinfor.h"

namespace denmtrx {
	class DenMtrx;
};

namespace shellpair {
	class ShellPair;
};

namespace denmtrxinfor {
	class DenMtrxInfor;
};

namespace symmjkdigest {
	class SymmJKDigest; 
}

namespace blockmatrixlist {
	class BlockMtrxList; 
}

namespace gints4d {

	using namespace tbb;
	using namespace shell;
	using namespace matrix;
	using namespace denmtrx;
	using namespace shellpair;
	using namespace spinmatrix;
	using namespace gintsinfor;
	using namespace denmtrxinfor;
	using namespace symmjkdigest;
	using namespace blockmatrixlist;
	using namespace sigshellpairinfor;

	/**
	 * \class TBB_GInts4D
	 *
	 * the result is in form of matrix
	 */
	class TBB_GInts4D {

		private:

			//
			// basic data needed for calculation
			//
			bool sameShellPair;                   ///< whether the bra and ket are coming from same source?
			UInt nSpin;                           ///< number of spin states
			Double thresh;                        ///< threshold value for gints4d work
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
			Mtrx& alphaMatrix;                    ///< alpha  matrix
			Mtrx& betaMatrix;                     ///< beta  matrix

		public:

			/**
			 * constructor 
			 * \param spStatus       whether bra and ket are same shell pair infor?
			 * \param braBatchIndex  the bra side sigAtomShellPair index
			 * \param ketBatchIndex  the ket side sigAtomShellPair index
			 * \param ginfor         information center for gints job
			 * \param bra,ket        shell pair infor class for bra and ket side
			 * \param b1,b2,k1,k2    reference for shell data
			 * \param alpha,beta     input density matrices
			 * \param denMtrxInfor   the density matrix information
			 * \return alphaMtrx     alpha result matrix
			 * \return betaMtrx      beta result matrix
			 */
			TBB_GInts4D(bool spStatus, const UInt& braBatchIndex0, const UInt& ketBatchIndex0, 
					const GIntJobInfor& ginfor0, const SigMolShellPairInfor& bra, 
					const SigMolShellPairInfor& ket, const MolShell& b1, const MolShell& b2, 
					const MolShell& k1, const MolShell& k2, const Mtrx& alpha, const Mtrx& beta,
					const DenMtrxInfor& denMtrxInfor0, Mtrx& alphaMtrx, Mtrx& betaMtrx):sameShellPair(spStatus),
			nSpin(ginfor0.getNSpin()),thresh(ginfor0.pickUpThresh()),spThresh(ginfor0.pickUpSPThresh()),
			braBatchIndex(braBatchIndex0),ketBatchIndex(ketBatchIndex0),bra1(b1),bra2(b2),ket1(k1),ket2(k2),
			braInfo(bra),ketInfo(ket),ginfor(ginfor0),alphaDen(alpha),betaDen(beta),denMtrxInfor(denMtrxInfor0),
			alphaMatrix(alphaMtrx),betaMatrix(betaMtrx) { };

			/**
			 * splitting constructor used by TBB
			 * we copy all of reference data except result
			 * which should be always empty initially
			 */
			TBB_GInts4D(const TBB_GInts4D& tbb_gints4d, split):sameShellPair(tbb_gints4d.sameShellPair),
			nSpin(tbb_gints4d.nSpin),thresh(tbb_gints4d.thresh),spThresh(tbb_gints4d.spThresh),
			braBatchIndex(tbb_gints4d.braBatchIndex),ketBatchIndex(tbb_gints4d.ketBatchIndex),
			bra1(tbb_gints4d.bra1),bra2(tbb_gints4d.bra2),ket1(tbb_gints4d.ket1),
			ket2(tbb_gints4d.ket2),braInfo(tbb_gints4d.braInfo),ketInfo(tbb_gints4d.ketInfo),
			ginfor(tbb_gints4d.ginfor),alphaDen(tbb_gints4d.alphaDen),betaDen(tbb_gints4d.betaDen), 
			denMtrxInfor(tbb_gints4d.denMtrxInfor),alphaMatrix(tbb_gints4d.alphaMatrix), 
			betaMatrix(tbb_gints4d.betaMatrix) { };

			/**
			 * destructor
			 */
			~TBB_GInts4D() { };

			/**
			 * functional operator to generate the integral results
			 */
			void operator()(const blocked_range<UInt>& r) const; 
	};

	/**
	 * \class GInts4D
	 * \brief top level class to mornitor 4D analytical integral calculation
	 */
	class GInts4D {

		private:

			GIntJobInfor ginfor;               ///< information center for real calculation

			////////////////////////////////////////////////////////////
			//                                                        //
			// for serial calculation, the above tbb structure is not //
			// efficient. Therefore, by re-structure the code of      //
			// TBB_GInts4DMatrix, now we provide the function for     //
			// doing 4D AO matrix calculation in serial mode.         //
			// this function is used for testing purpose              //
			//                                                        //
			////////////////////////////////////////////////////////////
			void do4DAOMtrx(const SigMolShellPairInfor& braInfor,
					const SigMolShellPairInfor& ketInfor,
					const MolShell& bra1, const MolShell& bra2, 
					const MolShell& ket1, const MolShell& ket2, 
					const Mtrx& alphaDen, const Mtrx& betaDen, 
					const DenMtrxInfor& denMtrxInfor,
					Mtrx& alphaResult, Mtrx& betaResult) const;

			///
			/// if user want to save the JK matrix into disk, this is where to go
			///
			void saveJKMtrx(UInt section, UInt nSpin, const GlobalInfor& globInfor, 
					const Mtrx& A, const Mtrx& B) const;

		public:

			////////////////////////////////////////////////////////////
			//                  constructor etc.                      //
			////////////////////////////////////////////////////////////

			/**
			 * constructor - build shell pair information based on one shell
			 */
			GInts4D(const GIntsInfor& infor0, const UInt& job):ginfor(infor0,job,0) { };

			/**
			 * destructor
			 */
			~GInts4D() { };

			///
			/// form JK matrix for the given one shell
			/// so all of shell pair information are same
			/// the input density matrix will be transformed into Cartesian
			/// inside, but finally they will be transformed back
			///
			void doJKMtrx(const MolShell& s, DenMtrx& den, SpinMatrix& intMtrx, 
					Double& eJK, bool printTiming, bool printMatrix, bool saveMatrix); 

			///
			/// recover the given data matrix from hard disk
			/// \param  dataPath   where the data is stored
			///
			void recover(const string& dataPath, SpinMatrix& M) const;
	};
}

#endif

