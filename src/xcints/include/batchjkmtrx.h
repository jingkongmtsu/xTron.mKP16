/**
 * \file    batchjkmtrx.h
 * \brief   update the JK Matrix in batches
 * \author  Fenglai Liu 
 */
#ifndef BATCHJKMTRX_H
#define BATCHJKMTRX_H
#include "libgen.h"
#include "dftmatrix.h"
#include "scalablevec.h"

namespace batchbasis {
	class BatchBasis;
}

namespace batchgrid {
	class BatchGrid;
}

namespace xcintsinfor {
	class XCIntJobinfor;
}

namespace halfjkrho {
	class HalfJKRho;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace shell {
	class MolShell;
}

namespace batchjkmtrx {

	using namespace shell;
	using namespace batchbasis;
	using namespace dftmatrix;
	using namespace batchgrid;
	using namespace sigatombasis;
	using namespace xcintsinfor;
	using namespace halfjkrho;

	///
	/// this is to calculate batch JK matrix (exchange/Coulomb)
	///
	/// we note that if the JK are doing together, then because the K matrix's dimension contains
	/// J matrix, therefore we will update the J matrix into the K matrix and the result are not
	/// saperated
	///
	class BatchJKMtrx {

		private:

			UInt intJob;        ///< the integral job name
			DFTMatrix  JMtrx;   ///< this is the batch J matrix (Coulomb),  in dimension of (nSigBas,nSigBas)
			DFTMatrix  KMtrx;   ///< this is the batch K matrix (Exchange), in dimension of (nBas,nSigBas)

		public:

			/**
			 * this is the default constructor to form the Fock matrix
			 * for the exchange rho 
			 */
			BatchJKMtrx(const MolShell& ms, const SigAtomBasis& sigList, const HalfJKRho& halfJKRho, 
					const BatchBasis& basis, const BatchGrid& grid, const XCIntJobInfor& infor);

			///
			/// destructor
			///
			~BatchJKMtrx() { };

			/**
			 * update the JK matrix into the result Fock matrix
			 *
			 * if only Coulomb matrix is calculated, then we only update the jMtrx
			 */
			void updateJKMtrx(const XCIntJobInfor& infor, const SigAtomBasis& sigList, 
					SpinMatrix& M, Mtrx& jMtrx) const; 
	};
}

#endif
