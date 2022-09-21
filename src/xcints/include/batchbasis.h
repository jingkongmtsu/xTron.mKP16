/**
 * \file    Batchbasis.h
 * \brief   produce basis set value and derivatives etc. in the given batch
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef BATCHBASIS_H
#define BATCHBASIS_H
#include "libgen.h"
#include "matrix.h"
#include "spinmatrix.h"
#include "scalablevec.h"

namespace shell {
	class MolShell;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace batchgrid {
	class BatchGrid;
}

namespace denmtrx {
	class DenMtrx;
}

namespace xcvar {
	class XCVar;
}

namespace batchbasis {

	using namespace shell;
	using namespace sigatombasis;
	using namespace batchgrid;
	using namespace matrix;
	using namespace xcvar;
	using namespace denmtrx;
	using namespace spinmatrix;

	class BatchBasis {

		private:

			///////////////////////////////////////////////
			//                data section               //
			///////////////////////////////////////////////

			// dimension information etc.
			UInt maxDerivOrder;           ///< maximum order of basis set derivatives
			UInt nGrids;                  ///< number of grid points for this batch
			UInt nSigBasis;               ///< number of significant basis sets

			/**
			 * values and derivative values for the basis sets in this batch
			 * the data is arranged in the order below:
			 * basis set values->1st derivative->2ed derivative
			 * ->3rd derivative->4th derivative 
			 * (ngrids, nsigBasis, nderivs)
			 */
			MtrxVec phi;

			///////////////////////////////////////////////
			//   numerical functions used to create the  //
			//      basis set values and derivatives     //
			///////////////////////////////////////////////

			/**
			 * This function is used to create the radial part data for 
			 * a given shell(not a composite one).
			 * \param  np         number of Gaussian primitives in this shell
			 * \param  coe        coefficients for Gaussian primitives
			 * \param  e          exponent array for Gaussian primitives
			 * \param  pts        grid points
			 * \param  xyz        xyz for the center of the shell
			 * \return rad        radial data in dimension of (ngrids, maxOrder+1)
			 */
			void radDFTBasis(const UInt& np, const Double* coe, const Double* e, 
					const Double* pts, const Double* xyz, DoubleVec& rad) const;

			/**
			 * This function generate the angular data for the given atom(characterized 
			 * by c)
			 * \param  lmax   maximum of the anglar momentum for this atom
			 * \param  pts    grid point coordinate
			 * \param  c      coordinate for the atom center
			 * \return ang    angular data for the basis set (nBas,nGrids)
			 */
			void angDFTBasis(const UInt& lmax, const Double* pts, const Double* c, DoubleVec& ang) const;

			/**
			 * This function simply combine the angular and radial part data
			 * together to form the basis set value. 
			 * \return basis basis set value in dimension of (nCarBas,nGrids)
			 */
			void dftBasis0(const UInt& L, const UInt& nTotalBas, 
					const DoubleVec& ang, const DoubleVec& rad, Double* basis) const; 

			/**
			 * This function gives the first order derivative value for the basis set
			 * \param  ng        number of grid points
			 * \param  nTolBas   number of Cartesian basis set the ang array (first dimension)
			 * \param  L         angular momentum for this shell
			 * \return basis     basis set 1st deriatives in dimension of (nCarBas,nGrids,3)
			 */
			void dftbasisderiv1(const UInt& ng, const UInt& L, const UInt& nTolBas, 
					const Double* ang, const Double* rad, Double* basis) const;

			/**
			 * This function gives the second order derivative value for the basis set
			 * \param  ng        number of grid points
			 * \param  nTolBas   number of Cartesian basis set the ang array (first dimension)
			 * \param  L         angular momentum for this shell
			 * \return basis     basis set 2ed deriatives in dimension of (nCarBas,nGrids,6)
			 */
			void dftbasisderiv2(const UInt& ng, const UInt& L, const UInt& nTolBas, 
					const Double* ang, const Double* rad, Double* basis) const;

			/**
			 * This function gives the third order derivative value for the basis set
			 * \param  ng        number of grid points
			 * \param  nTolBas   number of Cartesian basis set the ang array (first dimension)
			 * \param  L         angular momentum for this shell
			 * \return basis     basis set 3rd deriatives in dimension of (nCarBas,nGrids,10)
			 */
			void dftbasisderiv3(const UInt& ng, const UInt& L, const UInt& nTolBas, 
					const Double* ang, const Double* rad, Double* basis) const;

			/**
			 * This function gives the fourth order derivative value for the basis set
			 * \param  ng        number of grid points
			 * \param  nTolBas   number of Cartesian basis set the ang array (first dimension)
			 * \param  L         angular momentum for this shell
			 * \return basis     basis set 4th deriatives in dimension of (nCarBas,nGrids,15)
			 */
			void dftbasisderiv4(const UInt& ng, const UInt& L, const UInt& nTolBas, 
					const Double* ang, const Double* rad, Double* basis) const;

		public:

			/**
			 * initilize the data in batch basis
			 */
			BatchBasis(const XCVar& xcvar, const BatchGrid& grid, 
					const MolShell& molShl, const SigAtomBasis& sigList);

			/**
			 * some times we need a empty batch basis set to recieve data
			 * \param maxOrder  maximum derivative order
			 * \param ng        nGrids
			 * \param nb        number of significant basis sets
			 * \param nSpin     number of spin states
			 *
			 * for this constructor, we can not initilize the exchange energy density part
			 * since there are currently no information pass in (we need the MolShell pass in)
			 */
			BatchBasis(const UInt& nSpin, const UInt& maxOrder, const UInt& ng, const UInt& nb);

			/**
			 * destructor
			 */
			~BatchBasis() { };

			/**
			 * print function used in debug purpose
			 */
			void print() const;

			/**
			 * return the number of significant basis sets
			 */
			UInt getNSigBasis() const{ return nSigBasis; };

			/**
			 * return the number of grid points
			 */
			UInt getNGrids() const { return nGrids; };

			/**
			 * get the maximum basis set deriv order
			 */
			UInt getMaxBasisDerivOrder() const { return maxDerivOrder; };

			/**
			 * this is the main driver function to calculation basis set values on grids
			 * as well as basis set derivative values
			 */
			void setupPhi(const BatchGrid& grid, const MolShell& ms, const SigAtomBasis& sigList);

			/**
			 * return the phi according to the order given
			 */
			const Mtrx& getPhi(UInt deriv_order) const { return phi[deriv_order]; };

			/**
			 * form the batch basis set based on another batch basis set
			 * data. Actually we do split work here, only copy the batch basis
			 * set data related to the given iSigAtom 
			 */
			void batchBasisSplit(const BatchBasis& bbasis, 
					const SigAtomBasis& sigList, const UInt& iSigAtom);
	};

}


#endif
