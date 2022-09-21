/**
 * \file    hirshfeld.h
 * \brief   classes used to describe the Hirshfeld partition scheme
 * \author  Fenglai Liu and Jing Kong
 * \note
 *
 * Hirshfeld references:
 * F. L. Hirshfeld, Theor. Chim. Acta Vol 44, Page 129, 1977
 */
#ifndef HIRSHFELD_H
#define HIRSHFELD_H
#include "libgen.h"
#include "matrix.h"
#include "xcfunc.h"
#include "scalablevec.h"

namespace shell {
	class MolShell;
}

namespace xcintsinfor {
	class XCIntsInfor;
}

namespace atomdenmtrx {
	class AtomDenMtrx;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace batchbasis {
	class BatchBasis;
}

namespace batchgrid {
	class BatchGrid;
}

namespace batchvar {
	class BatchCar;
}

namespace hirshfeld {

	using namespace matrix;
	using namespace shell;
	using namespace xcfunc;
	using namespace atomdenmtrx;
	using namespace xcintsinfor;
	using namespace batchgrid;
	using namespace sigatombasis;
	using namespace batchbasis;
	using namespace batchvar;

	class Hirshfeld {

		private:

			Double thresh;     ///< threshold value to oversee weights calculation

			/**
			 * the faked xcfunctional information
			 * actually hirshfeld will only need density in terms of 
			 * the free atom density matrix. For correctly construct
			 * the density, we need the xcfunc information
			 *
			 * here we just pretend to be slater functional
			 */
			XCFunc xcfunc;

			/** 
			 * Hirshfeld weights: 
			 * \f$w(r) = \frac{\rho_{atom}(r)}{\sum_{natoms} \rho_{atom}(r)}\f$
			 * dimension is (ngrids,natoms)
			 *
			 * \rho_{atom}(r) is spin resolved free atom density based on 
			 * the current batch basis set value on the given batch grids
			 * the sum go over all of significant atoms in the current batch. 
			 */
			Mtrx weights;  

			/**
			 * this is free atom volumn for the given batch
			 * \f$freeAtomVol(i) = \int r^{3} \rho_{atom}(r) dr\f$
			 * this is used for polorizability calculation
			 *
			 * we need to calculate the free atom volumn here because
			 * the rho_{atom} here is for free atom density
			 */
			DoubleVec freeAtomVol;

		public:

			/**
			 * initilize the Hirshfeld weights data etc.
			 */
			Hirshfeld(const XCIntJobInfor& infor, const SigAtomBasis& sigList, 
					const BatchGrid& bg);

			/**
			 * destructor
			 */
			~Hirshfeld() { };

			/**
			 * form the weights as well as the free atom volumn
			 * \param ms      the whole shell information - we need it to form new batch basis
			 * \param sigList significant atom, basis set information
			 * \param basis   the original batch basis set is used to form new one for free
			 *                atom density
			 * \param bg      provides batch points and weights information 
			 * \param atomDenMatrices density matrix data for free atom
			 */
			void formWtsFreeAtomVol(const MolShell& ms, const SigAtomBasis& sigList, 
					const BatchBasis& basis, const BatchGrid& bg, 
					const AtomDenMtrx& atomDenMatrices);

			/**
			 * for given significant atom, we return its corresponding weights
			 * vector
			 */
			const Double* getWts(const UInt& iSigAtom) const {
				return weights.getPtr(0,iSigAtom);
			};

			/**
			 * for given significant atom, return it's free vol
			 */
			Double getFreeAtomVol(const UInt& iSigAtom) const {
				return freeAtomVol[iSigAtom];
			};
	};
}

#endif

