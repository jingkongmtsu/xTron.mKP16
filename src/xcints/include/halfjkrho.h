/**
 * \file    halfjkrho.h
 * \brief   produce half of exchange/Coulomb energy density (monitoring ESP integrals calculation)
 * \author  Fenglai Liu 
 */
#ifndef HALFJKRHO_H
#define HALFJKRHO_H
#include "libgen.h"
#include "matrix.h"
#include "spinmatrix.h"
#include "scalablevec.h"
#include "blas.h"
#include "blas1.h"

namespace shell {
	class MolShell;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace batchgrid {
	class BatchGrid;
}

namespace batchbasis {
	class BatchBasis;
}

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace denmtrx {
	class DenMtrx;
}

namespace xcvar {
	class XCVar;
}

namespace sigshellpairinfor {
	class SigMolShellPairInfor;
};

namespace halfjkrho {

	using namespace shell;
	using namespace sigatombasis;
	using namespace xcintsinfor;
	using namespace batchgrid;
	using namespace matrix;
	using namespace xcvar;
	using namespace denmtrx;
	using namespace spinmatrix;
	using namespace batchbasis;
	using namespace blas;
	using namespace sigshellpairinfor; 

	class HalfJKRho {

		private:

			///////////////////////////////////////////////
			//                data section               //
			///////////////////////////////////////////////

			// dimension information etc.
			UInt nGrids;                  ///< number of grid points for this batch
			UInt nSigBasis;               ///< number of significant basis sets

			/**
			 * value for producing exchange energy density
			 * we call it as halfExRho
			 *
			 * \f$ ex_{\nu}(r)  = \sum_{i}^{occ} \phi_{i}(r)
			 * \int \frac{\phi_{i}(r^{'})\phi_{\nu}(r^{'})}{|r-r^{'}|} dr^{'}\f$
			 *
			 * i is the occupied orbital, nu is the basis set index
			 *
			 * this result is in dimension of (nGrids,nBas)
			 * nBas is the whole shell dimension
			 *
			 * we hold this dimension is because after ESP integral calculation
			 * the basis dimension needs to be transformed into normal order
			 * so that to further combined with basis set value to form 
			 * the matrix.
			 */
			SpinMatrix halfExRho;

			/**
			 * value of half Coulomb energy density vector, we call it as halfCouRhoVec
			 *
			 * \f$ ec(r) = \sum_{\mu\nu}P_{\mu\nu}
			 * \int \frac{\phi_{\mu}(r^{'})\phi_{\nu}(r^{'})}{|r-r^{'}|} dr^{'}\f$
			 *
			 * the result is in nGrids
			 */
			DoubleVec halfCouRhoVec;

			/**
			 * values of half Coulomb energy density matrix, however it's in "shell pair vector"
			 * form (see the function of convertMatrixIntoSigSPVec in sigshellpairinfor.h)
			 *
			 * \f$ ec_{\mu\nu}(r) = \rho(r)\int \frac{\phi_{\mu}(r^{'})
			 * \phi_{\nu}(r^{'})}{|r-r^{'}|} dr^{'}\f$
			 *
			 * this part is needed in forming the numerical J matrix.
			 */
			DoubleVec halfCouRhoVecMatrix;

		public:

			/**
			 * initilize the data 
			 * \param infor  :   the XCInts job information center
			 * \param xcvar  :   whether we do exchange energy density in DFT hyper functional
			 * \param grid   :   batch grid information
			 * \param molShl :   the molecule shell information
			 * \param sigList:   significant atom and basis sets
			 * \param cartDen:   the shell pair vector form density matrix in Cartesian order
			 */
			HalfJKRho(const XCIntJobInfor& infor, const XCVar& xcvar,
					const BatchGrid& grid, const MolShell& molShl, 
					const SigAtomBasis& sigList, const DoubleVec& cartDen);

			/**
			 * destructor
			 */
			~HalfJKRho() { };

			/**
			 * print function used in debug purpose
			 */
			void print() const;

			/**
			 * return the number of significant basis sets
			 */
			UInt getNSigBasis() const{
				return nSigBasis;
			};

			/**
			 * return the number of grid points
			 */
			UInt getNGrids() const {
				return nGrids;
			};

			/**
			 * return the half exchange energy density
			 */
			const SpinMatrix& getHalfExRho() const { return halfExRho; };

			/**
			 * whether the half exrho data is really defined?
			 * 
			 * if the half exrho matrix has dimension information, then
			 * it's defined; else it's not defined
			 */
			bool hasHalfExRho() const;

			/**
			 * whether the half Coulomb energy density data is really defined?
			 * 
			 * if the array has dimension information, then
			 * it's defined; else it's not defined
			 */
			bool hasHalfCouRho() const {
				if (halfCouRhoVec.size() > 0) return true;
				return false;
			};

			/**
			 * return the half Coulomb energy density vector
			 */
			const DoubleVec& getHalfCouRhoVec() const { return halfCouRhoVec; };

			/**
			 * return the half Coulomb energy density matrix
			 */
			const DoubleVec& getHalfCouRhoVecMatrix() const { return halfCouRhoVecMatrix; };

			/**
			 * this is the driver function to produce the half exchange/Coulomb rho 
			 * \param infor  :   the XCInts job information center
			 * \param xcvar  :   whether we do exchange energy density in DFT hyper functional
			 * \param ms     :   the molecule shell information
			 * \param grid   :   batch grid information
			 * \param sigList:   significant atom and basis sets
			 * \param den    :   the input whole density matrix
			 * \param cartDen:   the shell pair vector form density matrix in Cartesian order
			 * \param spData :   the significant shell pair information
			 */
			void doHalfJKRho(const XCIntJobInfor& infor, const XCVar& xcvar, const MolShell& ms, 
					const BatchGrid& grid, const SigAtomBasis& sigList, const BatchBasis& basis,
					 const DenMtrx& den, const DoubleVec& cartDen, const SigMolShellPairInfor& spData);
			void doHalfJKRhoXHole(const XCIntJobInfor& infor, const XCVar& xcvar, const MolShell& ms,
					 const BatchGrid& grid, const SigAtomBasis& sigList, const BatchBasis& basis,
                                         const DenMtrx& den, const DoubleVec& cartDen, const SigMolShellPairInfor& spData, Double sValue); //yw

			///
			/// update this batch halfCouRhoVecMatrix to the global result
			///
			void updateHalfCouRhoVecMatrix(DoubleVec& s2VecJMtrx) const {
				vaxpy(&halfCouRhoVecMatrix[0],&s2VecJMtrx[0],ONE,halfCouRhoVecMatrix.size());
			};

	};

}


#endif
