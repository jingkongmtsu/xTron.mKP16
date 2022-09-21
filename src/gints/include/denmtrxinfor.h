/**
 * \file    denmtrxinfor.h
 * \brief   abstract the information from density matrices
 * \author  Fenglai Liu 
 *
 */
#ifndef DENMTRXINFOR_H
#define DENMTRXINFOR_H
#include "libgen.h"
#include "spinmatrix.h"
using namespace spinmatrix;

namespace denmtrx {
	class DenMtrx;
}

namespace shell {
	class MolShell;
}

namespace gintsinfor {
	class GIntJobInfor;
}

namespace shellpair {
	class ShellPair;
}

namespace sigshellpairinfor {
	class SigAtomShellPairInfor;
}

namespace globalinfor{
	class GlobalInfor;
}

namespace denmtrxinfor {

	using namespace denmtrx;
	using namespace shell;
	using namespace shellpair;
	using namespace gintsinfor;
	using namespace globalinfor;
	using namespace sigshellpairinfor;

	/**
	 * \class TBB_DenMtrxInfor
	 *
	 * This is a working class to form the density matrix infor data
	 * by using the threads
	 */
	class TBB_DenMtrxInfor {

		private:

			//
			// basis set type
			//
			UInt type;                            ///< what type of basis set it's?

			//
			// reference for the input 
			//
			const MolShell&     rs;               ///< row shell information
			const MolShell&     cs;               ///< column shell information
			const DenMtrx&      denMtrx;          ///< the original density matrix data

			//
			// reference for the output
			//
			SpinMatrix&     denMtrxAtomInfor;     ///< the infor data on atom-atom block
			SpinMatrix&     denMtrxShellInfor;    ///< the infor data on shell-shell block

		public:

			/**
			 * constructor 
			 */
			TBB_DenMtrxInfor(const UInt& type0, const MolShell& rowShell, 
					const MolShell& colShell, const DenMtrx& denMtrx0, SpinMatrix& atomInfor, 
					SpinMatrix& shellInfor):type(type0),rs(rowShell),cs(colShell),denMtrx(denMtrx0),
			denMtrxAtomInfor(atomInfor),denMtrxShellInfor(shellInfor) { };

			/**
			 * destructor
			 */
			~TBB_DenMtrxInfor() { };

			/**
			 * functional operator to generate the integral results
			 * this is used as parallel_for
			 */
			void operator()(const blocked_range<UInt>& r) const; 
	};

	/**
	 * this class is used to abstract the density matrix information  
	 * used for integrals forming. This class is based on the "dense"
	 * form of density matrix.
	 *
	 * the content inside the class focus on two data matrix.
	 *
	 * The density matrix atom infor is a matrix in (natoms,natoms)
	 * dimension. This matrix extract the maximum absolute value of 
	 * block density matrix per each atom shell pair.
	 *
	 * The density matrix shell infor is a matrix in (nshells,nshells)
	 * dimension. This matrix extract the maximum absolute value of 
	 * block density matrix per each shell pairs for all of shells.
	 *
	 * we note that generally we do not assume that the two input shell
	 * should be same.
	 */
	class DenMtrxInfor {

		private:
			
			//
			// information
			//
			UInt upperDimEnableThreads;     ///< dimension limit to enable threads mode in forming data
			UInt type;                      ///< what kind of type of the dimension of density matrix?

			//
			// data part
			//
			DoubleVec maxDenMtrxPerRow;     ///< maximum element for each row of denMtrxShellInfor
			DoubleVec maxDenMtrxPerCol;     ///< maximum element for each column of denMtrxShellInfor
			SpinMatrix  denMtrxAtomInfor;   ///< stores the abs max value per atom shell pair
			SpinMatrix  denMtrxShellInfor;  ///< stores the abs max value per shell pair

			/**
			 * this is the working function to form density matrix infor data in single threads
			 */
			void directFormPMaxInfor(const MolShell& rs, const MolShell& cs, 
					const DenMtrx& denMtrx);

			/**
			 * this is the driver function to threadly form the density matrix infor data
			 * per shell-shell block and atom-atom block
			 */
			void formPMaxInforInThreads(const GlobalInfor& infor, const MolShell& rs, 
					const MolShell& cs, const DenMtrx& denMtrx);

		public:

			///
			/// constructor to build all of data
			/// 
			DenMtrxInfor(const GlobalInfor& infor, const MolShell& rs, 
					const MolShell& cs, const DenMtrx& denMtrx);

			///
			/// this constructor is used to hold basic information of data
			/// the real building process is used in init function
			///
			DenMtrxInfor(UInt nSpin):upperDimEnableThreads(100),type(-1),
			denMtrxAtomInfor(nSpin),denMtrxShellInfor(nSpin) { };

			///
			/// to build the data through input shell and density matrices
			///
			void init(const GlobalInfor& infor, const MolShell& rs, 
					const MolShell& cs, const DenMtrx& denMtrx);

			///
			/// return the maximum PMax value for the row shell
			///
			Double getPmaxOnRowShell(UInt index, UInt iSpin) const { 
				const Mtrx& M = denMtrxShellInfor.getMtrx(0);
				UInt nRow = M.getRow();
				return maxDenMtrxPerRow[index+iSpin*nRow];
			};

			///
			/// return the maximum PMax value for the col shell
			///
			Double getPmaxOnColShell(UInt index, UInt iSpin) const { 
				const Mtrx& M = denMtrxShellInfor.getMtrx(0);
				UInt nCol = M.getCol();
				return maxDenMtrxPerCol[index+iSpin*nCol];
			};

			///
			/// return the number of spin states for the density matrix
			/// because it's same in terms of atom/shell, so we just
			/// check the 
			///
			UInt getNSpin() const { return denMtrxAtomInfor.getNSpin(); };

			///
			/// destructor
			///
			~DenMtrxInfor() { };

			///
			/// function to get the biggest PMax value for the J type digestion
			/// For J type of digestion, we will have P12*(12|34) and P34*(12|34)
			/// so we will return the maximum P value by checking both P12 and P34
			/// 
			///
			/// the value returns will be sum between alpha and beta, since the 
			/// coulomb type of digestion is spin sumed
			///
			/// we note, that here the switch status between braSP and ketSP does 
			/// not affect the result
			///
			/// \param braSP            : bra side shell pair data
			/// \param ketSP            : ket side shell pair data
			///
			Double getJPMax(const ShellPair& braSP, const ShellPair& ketSP) const;

			///
			/// function to get the biggest PMax value for the K type digestion
			/// For K type of digestion, we will have:
			/// - P13*(12|34) 
			/// - P14*(12|34)
			/// - P23*(12|34) 
			/// - P24*(12|34)
			/// 4 situations for each spin state
			///
			/// we will return the maximum P value by checking all of P13, P14 etc.
			/// for both alpha and beta spin states (this is not spin sum over case)
			//
			/// we note, that here the switch status between braSP and ketSP does 
			/// not affect the result, referring that  (12|34) <=> (34|12)
			/// you can see that the maximum p value does not alter
			///
			///
			Double getKPMax(const ShellPair& braSP, const ShellPair& ketSP) const;

			///
			/// function to get the biggest PMax Pair value for the J type digestion
			/// For J type of digestion, this is basically P12*P34
			///
			/// the value returns will be sum between alpha and beta, since the 
			/// coulomb type of digestion is spin sumed
			///
			/// \param braSP            : bra side shell pair data
			/// \param ketSP            : ket side shell pair data
			///
			Double getJPPairMax(const ShellPair& braSP, const ShellPair& ketSP) const;

			///
			/// function to get the biggest PMax value for the K type digestion
			/// For K type of digestion, we will have:
			/// - P13*(12|34)*P24 
			/// - P14*(12|34)*P23
			/// 2 situations for each spin state
			///
			/// we will return the maximum P value by checking both of two cases
			/// for both alpha and beta spin states (this is not spin sum over case)
			//
			/// we note, that here the switch status between braSP and ketSP does 
			/// not affect the result, referring that  (12|34) <=> (34|12)
			/// you can see that the maximum p value does not alter
			///
			///
			Double getKPPairMax(const ShellPair& braSP, const ShellPair& ketSP) const;

			///
			/// function to test that whether the corresponding coulomb digestion is significant for 
			/// the atom shell pair data
			///
			bool isJSig(const Double& thresh, const SigAtomShellPairInfor& sigBraAtomSPInfor, 
					const SigAtomShellPairInfor& sigKetAtomSPInfor) const;

			///
			/// function to test that whether the corresponding exchange digestion for the given spin state
			/// is significant for the atom shell pair data
			///
			bool isKSig(const UInt& iSpin, const Double& thresh, const SigAtomShellPairInfor& sigBraAtomSPInfor, 
					const SigAtomShellPairInfor& sigKetAtomSPInfor) const;

			///
			/// function to test that whether the corresponding coulomb integral derivatives 
			/// digestion is significant for the atom shell pair data
			///
			bool isJDerivSig(const Double& thresh, const SigAtomShellPairInfor& sigBraAtomSPInfor, 
					const SigAtomShellPairInfor& sigKetAtomSPInfor) const;

			///
			/// function to test that whether the corresponding exchange integral derivatives
			/// digestion for the given spin state is significant for the atom shell pair data
			///
			bool isKDerivSig(const UInt& iSpin, const Double& thresh, const SigAtomShellPairInfor& sigBraAtomSPInfor, 
					const SigAtomShellPairInfor& sigKetAtomSPInfor) const;

	};

}

#endif
