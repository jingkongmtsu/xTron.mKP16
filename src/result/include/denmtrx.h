/**
 * \file    denmtrx.h
 * \brief   describe the density matrices
 * \author  Fenglai Liu 
 *
 */
#ifndef DENMTRX_H
#define DENMTRX_H
#include<string>
#include "libgen.h"
#include "spinmatrix.h"
#include "globalinfor.h"
using namespace spinmatrix;
using namespace globalinfor;
using namespace std;

namespace mo {
	class MO;
}

namespace shell {
	class MolShell;
}

namespace molecule {
	class Molecule;
}

namespace oneemtrx {
	class OneEMtrx;
}

namespace scf {
	class SCF;
}

namespace denmtrx {

	using namespace mo;
	using namespace scf;
	using namespace shell;
	using namespace molecule;
	using namespace oneemtrx;

	// declare the denMtrx class
	class DenMtrx;

	///
	/// here we define the vector of density matrices
	///
	typedef std::vector<DenMtrx,tbb::scalable_allocator<DenMtrx> > DenMtrxVec;

	/**
	 * \brief describing the AO based density matrix
	 */
	class DenMtrx : public SpinMatrix {

		private:

			GlobalInfor infor;        ///< a copy of global information
			UInt section;             ///< geometry section number
			UInt nAlpha;              ///< number of alpha electrons
			UInt nBeta;               ///< number of beta electrons

		public:

			///
			/// constructor to build the blank density matrix
			/// 
			DenMtrx(const GlobalInfor& infor0, const Molecule& mol, 
					const MolShell& rs, const MolShell& cs, const UInt& nSpin);

			///
			/// destructor
			///
			~DenMtrx() { };

			///
			/// get number of electrons
			///
			UInt getNEle(UInt iSpin) const {
				if (iSpin == 0) return nAlpha;
				return nBeta;
			};

			///
			/// get the section number
			///
			UInt getSec() const { return section; };

			/**
			 * forming the density matrix from the given MO space
			 * P = C*C^{T}
			 */
			void formDenMtrx(const MO& mo);

			/**
			 * forming the core guess
			 *
			 * in this version the Fock matrix (core) and the 
			 * orthogonal matrix is formed inside
			 */
			void coreGuess(const MolShell& ms, const Molecule& mol);

			/**
			 * forming the core guess with input core matrix etc.
			 */
			void coreGuess(const MolShell& ms, const Molecule& mol, const OneEMtrx& oneEMtrx);

			/**
			 * forming the guess with atomic density matrix method
			 */
			void atomicDenMtrxGuess(const MolShell& ms, const Molecule& mol);

			/**
			 * the driver function to form guess for scf
			 */
			void formSCFGuess(const MolShell& ms, const Molecule& mol, const SCF& scf);

			/**
			 * write the current density matrix data into binary file in standard dir structure
			 */
			void writeToDisk() const;

			/**
			 * write the current density matrix data into binary file according to the given 
			 * folder
			 */
			void writeToDisk(const string& dir) const;

			/**
			 * recover the density matrix data from the given folder name
			 */
			void recover(const string& dir);

      ///
      /// this is a tool function to see that whether this is single
      /// electron system?
      ///
      bool isSingleEleSystem() const {
        return (nAlpha == 0 || nBeta == 0);
      };
	};

}

#endif
