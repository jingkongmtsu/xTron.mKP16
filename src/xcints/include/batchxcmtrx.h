/**
 * \file    batchxcmtrx.h
 * \brief   update the XC Matrix in batches
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef BATCHXCMTRX_H
#define BATCHXCMTRX_H
#include "libgen.h"
#include "dftmatrix.h"
#include "scalablevec.h"

namespace batchbasis {
	class BatchBasis;
}

namespace batchvar {
	class BatchVar;
}

namespace batchgrid {
	class BatchGrid;
}

namespace batchfunc {
	class BatchFunc;
}

namespace halfjkrho {
	class HalfJKRho;
}

namespace xcvar {
	class XCVar;
}

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace shell {
	class MolShell;
}

namespace batchxcmtrx {

	using namespace shell;
	using namespace batchbasis;
	using namespace batchvar;
	using namespace batchfunc;
	using namespace dftmatrix;
	using namespace batchgrid;
	using namespace xcvar;
	using namespace sigatombasis;
	using namespace halfjkrho;
	using namespace xcintsinfor;

	class BatchXCMtrx : public DFTMatrix {

		//
		// batch xc matrix does not have it's own member data
		// it's same with DFT matrix
		//

		private:

			/**
			 * build the Fock matrix for given spin state
			 * \f[  F_{\munu} = \sum w_{i}(r) \frac{\partial F(r)}
			 *                       {\partial \xi(r)}\frac{\partial \xi(r)}
			 *                       {\partial P_{\mu\nu}} \f]
			 * w is the weight, \f$\xi\f$ is the density variable. 
			 */
			void buildFxcMtrx(const BatchBasis& basis, 
					const XCVar& xcvar, const UInt& iSpin, const MtrxVec& halfFxc);

			/**
			 * build Half of the Fock matrix. Actually here we multiply half of the
			 * \f$\frac{\partial \xi(r)}{\partial P_{\mu\nu}}\f$ for each variable
			 * We label it as \f$\phi(r)\f$.
			 * \f[  halfFxc_{\mu} = \sum w_{i}(r) \frac{\partial F(r)}
			 *                       {\partial \xi(r)} \phi(r)\f]
			 */
			void buildHalfFxcMtrx(const BatchBasis& basis, const BatchFunc& func,
					const BatchGrid& grid, const XCVar& xcvar, const UInt& iSpin, MtrxVec& halfFxc) const;

		public:

			/**
			 * this is the default constructor to form the Fock matrix
			 * in normal way
			 */
			BatchXCMtrx(const SigAtomBasis& sigList, const BatchBasis& basis, 
					const BatchFunc& func, const BatchGrid& grid, 
					const XCVar& xcvar, const XCIntsInfor& infor);

			///
			/// destructor
			///
			~BatchXCMtrx() { };

			/**
			 * update the XC matrix(convert the result from significant order
			 * to normal basis set order)
			 */
			void updateXCMtrx(const SigAtomBasis& sigList, SpinMatrix& M) const {
				for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
					intoNormOrder(sigList,sigList,iSpin,M.getMtrx(iSpin));
				}
			};

	};

	///
	/// this is to perform batch xc matrix in terms of exchange energy density
	/// for this type of matrix, it always has a dimension which is in shell
	/// dimension
	///
	class BatchEXRhoMtrx : public DFTMatrix {

		public:

			/**
			 * this is the default constructor to form the Fock matrix
			 * for the exchange rho 
			 */
			BatchEXRhoMtrx(const MolShell& ms, const SigAtomBasis& sigList, 
					const BatchBasis& basis, const HalfJKRho& halfJKRho, 
					const BatchFunc& func, const BatchGrid& grid, 
					const XCVar& xcvar, const XCIntJobInfor& infor);

			///
			/// destructor
			///
			~BatchEXRhoMtrx() { };

			/**
			 * update the XC matrix(convert the result from significant order
			 * to normal basis set order)
			 */
			void updateXCMtrx(const SigAtomBasis& sigList, SpinMatrix& M) const {
				for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
					intoNormOrder(sigList,iSpin,M.getMtrx(iSpin));
				}
			};
	};

	class BatchXCVect : public DFTVect {

		private:

			/**
			 * build the Fock matrix in RI way
			 * \f$ F_{XC} = \sum w_{i}(r) f_{XC}^{\rho}(\rho(r)) \phi_{P}(r) \f$
			 * We note that currently we only know how to process rho and gradient
			 * rho. Here phi is the aux basis set/gradient basis set, the whole 
			 * result is in dimension of nSigAuxBasis.
			 */
			void buildXCVect(const BatchBasis& basis, const BatchFunc& func,
					const BatchGrid& grid, const XCVar& xcvar, const UInt& iSpin); 

		public:

			/**
			 * build the XC vector
			 */
			BatchXCVect(const SigAtomBasis& sigList, const BatchBasis& basis, 
					const BatchFunc& func, const BatchGrid& grid, 
					const XCVar& xcvar, const XCIntsInfor& infor);

			/**
			 * destructor
			 */
			~BatchXCVect() { };

			/**
			 * update the XC vector(convert the result from significant order
			 * to normal basis set order) for the given spin state. 
			 */
			void updateXCVect(const SigAtomBasis& sigList, DoubleVec& result) const {
				intoNormOrder(result,sigList);
			};
	};
}

#endif
