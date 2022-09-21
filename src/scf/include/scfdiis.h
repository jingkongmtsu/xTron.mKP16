/**
 * \file    scfdiis.h
 * \brief   produce the SCF DIIS result
 * \author  Fenglai Liu 
 * \note
 *
 * The implementation here is following the references below:
 *
 * (oringinal DIIS idea)
 * Peter Pulay,
 * Convergence acceleration of iterative sequences. the case of scf iteration
 * Chemical Physics Letters, 1980, volume 73, pages 393 - 398
 *
 * (this is for close shell case)
 * Peter Pulay,
 * Improved SCF convergence acceleration
 * Journal of Computational Chemistry, 1982, volume 3, pages 556--560
 *
 * (this is for open shell case)
 * Tracy P. Hamilton and Peter Pulay
 * Direct inversion in the iterative subspace (DIIS) optimization of open‐shell, 
 * excited‐state, and small multiconfiguration SCF wave functions
 * J. Chem. Phys. vol 84, page 5728 (1986) 
 *
 */
#ifndef SCFDIIS_H
#define SCFDIIS_H
#include "libgen.h"
#include "matrix.h"
#include "globalinfor.h"
#include "histdataman.h"
#include "scalablevec.h"
using namespace matrix;
using namespace globalinfor;
using namespace histdataman;

namespace denmtrx {
	class DenMtrx;
}

namespace spinmatrix {
	class SpinMatrix;
}

namespace oneemtrx {
	class OneEMtrx;
}

namespace scfparam {
	class SCFParam;
}

namespace scfdiis {

	using namespace spinmatrix; 
	using namespace denmtrx;
	using namespace oneemtrx;
	using namespace scfparam;

	/**
	 * this class is used to obtain the SCF DIIS result
	 */
	class SCFDIIS {

		private:

			GlobalInfor infor;           ///< save a copy of global information
			Double error;                ///< DIIS error
			Mtrx errMtrx;                ///< error matrix formed based on error vector
			HistDataMan errVecList;      ///< entire error vectors list

			/**
			 * this function is used to generate the error vector for each 
			 * DIIS process based on the calculated Fock matrix etc. 
			 *
			 * In this function, error vector is defined as:
			 * \f$ E = FPS - SPF\f$
			 * F is the Fock matrix for given spin state, P is 
			 * corresponding density matrix, and S is the overlap
			 * matrix. finally, we will transform E into orthogonal
			 * basis set by using orthogonal matrix(O)
			 */
			void errorVec(const Mtrx& S, const Mtrx& O, const DenMtrx& den, 
					const SpinMatrix& fock, SpinMatrix& errVec) const;

			/**
			 * update the error based on the errVec
			 */
			void updateError(const SpinMatrix& errVec); 

			/**
			 * this is to form the error matrix for DIIS
			 */
			void formErrMtrx(const SpinMatrix& errVec);

		public:

			/**
			 * constructor
			 */
			SCFDIIS(const SCFParam& par);

			/**
			 * destructor
			 */
			~SCFDIIS() { };

			/**
			 * get the DIIS error
			 */
			Double diisError() const { return error; };

			/**
			 * main function to update DIIS error and error matrix etc.
			 * we will update the error matrix for all SCF cycles
			 */
			void updateDIISInfor(const OneEMtrx& oneEMtrx, const DenMtrx& den, const SpinMatrix& fock);

			/**
			 * update the DIIS coefficients
			 *
			 * we will form the error matrix given by the input scfIndexArray
			 * that means, this matrix is only a part of errMtrx
			 *
			 * since the DIIS coefficients array dimension will be nDim+1
			 * therefore we only copy the nDim elements
			 */
			void updateCoefs(const UIntVec& scfIndexArray, DoubleVec& coe) const;

			/**
			 * print out the error matrix as well as coefficients
			 * for debug purpose
			 */
			void print() const;
	};
}

#endif

