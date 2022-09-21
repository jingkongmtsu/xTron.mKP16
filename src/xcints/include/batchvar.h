/**
 * \file    Batchvar.h
 * \brief   produce density variables in the given batch
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef BATCHVAR_H
#define BATCHVAR_H
#include "libgen.h"
#include "scalablevec.h"
#include "matrix.h"

namespace batchbasis {
	class BatchBasis;
}

namespace halfjkrho {
	class HalfJKRho;
}

namespace batchgrid {
	class BatchGrid;
}

namespace xcvar {
	class XCVar;
}

namespace shell {
	class MolShell;
}

namespace dftmatrix {
	class DFTMatrix;
	class DFTVect;
}

namespace denmtrx {
	class DenMtrx;
}

namespace xcintsinfor {
	class XCIntJobInfor;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace batchvar {

	using namespace shell;
	using namespace batchbasis;
	using namespace batchgrid;
	using namespace xcvar;
	using namespace matrix;
	using namespace denmtrx;
	using namespace xcintsinfor;
	using namespace dftmatrix;
	using namespace sigatombasis;
	using namespace halfjkrho; 

	/**
	 * The basic idea to build variables in DFT module is quite simple.
	 * In the first step, according to the variable information we combine
	 * the data in batchbasis together with the density matrix. In this
	 * step dgemm is used. In the second step, the converted data(we called
	 * it as "halfVar") is then combined with the batchbasis data again, so
	 * to produce final varible. Such process applies to all of variables
	 * except the exchange energy density, which is more sophisticated. 
	 */
	class BatchVar {

		private:

			// dimension information
			UInt nGrids;                  ///< number of grid points for this batch
			UInt nSigBasis;               ///< number of significant basis sets

			// data section
			DoubleVec vars;               ///< variables 

			/**
			 * build half variable: \f$ result = P_{\mu\nu}\phi_{\nu}\f$
			 * namely to combine the density matrix with the basi set 
			 * value or its derivatives.
			 * \param  DenMtrx :  alpha/beta density matrix
			 * \param  basis   :  batch basis set data
			 * \param  xcvar   :  do density matrix + 2ed basis set derivatives?
			 * \return halfvar :  the result after combination
			 */
			void buildHalfVar(const BatchBasis& basis, const Mtrx& denMtrx, 
					const XCVar& xcvar, MtrxVec& halfVar) const; 

			/**
			 * build the variable for the given spin state
			 * \param  basis   :  batch basis set data
			 * \param halfvar  :  half variable with density matrix
			 * \param iSpin    :  0 is alpha and 1 is beta
			 * \param xcvar    :  the xc variable information
			 */
				void buildVar(const BatchBasis& basis, const MtrxVec& halfVar, 
						const UInt& iSpin, const XCVar& xcvar); 

			/**
			 * build the DFT variable in RI way for the given
			 * spin state
			 * \param  basis   :  batch aux basis set data
			 * \param xcvar    :  the xc variable information
			 * \param  DenMtrx :  alpha/beta density matrix
			 * \param iSpin    :  0 is alpha and 1 is beta
			 */
			void buildRIVar(const BatchBasis& basis, const XCVar& xcvar, 
					const Double* denMtrx, const UInt& iSpin);

			/**
			 * For the approximated DFT variables, there's potential 
			 * issue that the calculation result on grid could be 
			 * qualitatively wrong. For example, density will always
			 * between 0 and 1 on each grid, however; the approximated
			 * density could be less than 0. Therefore, we need to
			 * smooth it
			 */
			void smoothDFTVar(const XCVar& xcvar);

		public:
			const DenMtrx& denMtrx;
			/**
			 * constructor for the batch variable class in common way
			 * \param sigList:     significant list, which is used to convert density matrix
			 * \param basis  :     batch basis set 
			 * \param xcvar  :     providing variable information
			 * \param den    :     density matrix 
			 */
			BatchVar(const SigAtomBasis& siglist, const BatchBasis& basis, 
					const XCVar& xcvar, const DenMtrx& den);

			/**
			 * constructor for building the batch variables
			 * comparing with the above constructor, this one uses the 
			 * density matrix in significant order directly
			 */
			BatchVar(const BatchBasis& basis, const XCVar& xcvar, 
					const DFTMatrix& den, const DenMtrx&);

			/**
			 * constructor for the batch variable class in RI way
			 * in this way, the density matrix is only one dimension vector
			 * \param basis:     batch basis set 
			 * \param xcvar:     providing variable information
			 * \param vec  :     density matrix in vector form for alpha and beta
			 */
			BatchVar(const BatchBasis& basis, const XCVar& xcvar, const DFTVect& vec,
			         const DenMtrx&);

			/**
			 * constructor for the batch variable, only for numerical J/K
			 * \param sigList:     significant list, which is used to convert density matrix
			 * \param basis  :     batch basis set 
			 * \param infor  :     providing numerical JK information
			 * \param den    :     density matrix 
			 */
			BatchVar(const MolShell& ms, const SigAtomBasis& siglist, 
					const BatchBasis& basis, const HalfJKRho& halfJKRho,
					const XCIntJobInfor& infor, const DenMtrx& den);

			///
			/// destructor
			///
			~BatchVar() { };

			///
			/// build the exchange energy density
			///
			void buildExRho(const MolShell& ms, const SigAtomBasis& sigList, 
					const DenMtrx& denMtrx, const BatchBasis& basis, const HalfJKRho& halfJKRho, 
					const XCVar& xcvar);

			/**
			 * return the number of grids
			 */
			UInt getNGrids() const {
				return nGrids;
			};

			/**
			 * return the variable array
			 */
			const Double* getVar(UInt pos) const {
				if(pos == static_cast<UInt>(-1)) return NULL;  // we do not have this variable
				return &vars[pos*nGrids];
			};

			/**
			 * return the electron number estimation for alpha/beta density
			 * \param  iSpin  alpha or beta?
			 * \param  xcvar  get the variable storage information
			 * \param  bg     get the weights information
			 */
			Double getNEle(const UInt& iSpin, const XCVar& xcvar,
					const BatchGrid& bg) const;

			/**
			 * return the sum of exchange energy density(for DFT functional case)
			 * this is used in debug purpose
			 */
			Double getEX(const UInt& iSpin, const XCVar& xcvar,
					const BatchGrid& bg) const;

			/**
			 * return the sum of exchange energy density
			 * this is used in debug purpose for numerical JK
			 */
			Double getEX(const UInt& iSpin, const XCIntJobInfor& infor, const BatchGrid& bg) const;

			/**
			 * return the batch Coulomb energy density
			 */
			Double getEJ(const XCIntJobInfor& infor, const BatchGrid& bg) const;

      //compute exact exchange hole.
      void xhole(const XCVar& xcvar, vector<DoubleVec>& xholeVal) const;

			/**
			 * printing function used in debug purpose
			 */
			void print(const XCVar& xcvar) const;

			/**
			 * debug printing for the numerical JK variable
			 */
			void print(const XCIntJobInfor& infor) const;

			void varForFort(const Double*& rhoA, const Double*& rhoB,
			      const Double*& DRA, const Double*& DRB, const Double*& TA,
			      const Double*& TB, const Double*& LA, const Double*& LB,
			      const Double*& EXRA, const Double*& EXRB, const XCVar& xcvar) const;
	};
}

#endif
