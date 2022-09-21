/**
 * \file    mo.h
 * \brief   describing the MO generated in the SCF procedure
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef MO_H
#define MO_H
#include "libgen.h"
#include "globalinfor.h"
#include "spinmatrix.h"
using namespace globalinfor;
using namespace spinmatrix;

namespace shell {
	class MolShell;
}

namespace molecule {
	class Molecule;
}

namespace oneemtrx {
	class OneEMtrx;
}

namespace fracspininfor {
	class  FracSpinInfor;
}

namespace mo {

	using namespace shell;
	using namespace molecule;
	using namespace oneemtrx;
	using namespace fracspininfor;

	/**
	 * \class   MO
	 * \brief   molecuar orbitals generated in the SCF 
	 */
	class MO : public SpinMatrix {

		protected:

			GlobalInfor infor;        ///< a copy of global information
			UInt section;             ///< geometry section number
			UInt shellCode;           ///< the shell code for the basis set
			UInt nAlpha;              ///< number of occupied orbitals for alpha state
			UInt nBeta;               ///< number of occupied orbitals for beta state
			DoubleVec alphaEnergy;    ///< MO energy for alpha orbitals
			DoubleVec betaEnergy;     ///< MO energy for alpha orbitals

			/**
			 * constructing the MO for one spin state. Basically we will solve the 
			 * Fock matrix inside and get MO directly
			 * We note that the Fock matrix will be de stroyed inside
			 * \param iSpin   which spin state we refer to
			 * \param ortho   orthogonal matrix in solving the MO
			 * \param Fock    the result Fock matrix in AO form
			 * \param withCSCEqI whether you want the result mo satisfy C^{+}SC = I
			 *
			 * the orthogonal matrix itself must be a full matrix
			 * this is required by the multiplication inside this function
			 * Please make sure of this
			 *
			 * the solving MO process is to solve a symmetrical eigen-vector problem.
			 * For DFT and Hartree-Fock equation, the function is in form of \f$ FC = SCE \f$.
			 * where F is the result Fock matrix, S is two body overlap matrix, C is the
			 * MO and E is the orbital energy.
			 *
			 * To solve this equation, firstly we need some matrix satisfy \f$T^{+}ST = I\f$ so that
			 * to drop the S in the equation. \f$T^{+}\f$ means the transpose of the matrix T. 
			 * Also we note that in this case $T^{-} = T^{+}S$. 
			 * Suggesting we have the T, so we have \f$ (T^{+}FT)(T^{-}C) = T^{+}SCE\f$ 
			 * which leads to \f$ (T^{+}FT)(T^{-}C) = T^{-}CE\f$ and it's actually
			 * \f$ F^{'}C^{'} = C^{'}E\f$. So solving this eigen function it's able to
			 * derive the \f$ C^{'} \f$, and because F is symmetrical matrix the result
			 * MO satisfy \f$ (C^{'})^{+}C^{'} = I \f$.
			 *
			 * however, in practise we use the mo satisfy the condition that \f$ C^{+}SC = I\f$.
			 * In this case, it's easy to see that \f$ C^{'} = T^{-}C\f$ so that \f$ C = TC^{'}\f$.
			 * this is the reason why in the end we multiply the mo with orthogonal matrix T
			 * if withCSCEqI is true.
			 *
			 * For both \f$ C^{'} \f$ and C they are proved to be equivalnt on
			 * physics. If the Fock matrix is in full rank (not singular), the
			 * best matrix for T is the lowdin orthogonalation matrix of
			 * \f$ S^{-1/2} \f$. For more information about the physical meaning of T,
			 * please refer to the paper "On Löwdin's method of symmetric orthogonalization"
			 * by I. Mayer in International journal of quantum chemistry, vol 90 pages 63-65, 2002.
			 */
			// tvoid solveFock(bool withCSCEqI, const UInt& iSpin, const Mtrx& ortho, Mtrx& Fock);
			void solveFock(const UInt& iSpin, const Mtrx& ortho, Mtrx& Fock);

		public:

			/**
			 * constructor for the MO - build the dimension etc. for the MO
			 * in default the nOrb = nBas 
			 */
			MO(const GlobalInfor& infor0, const MolShell& ms, 
					const Molecule& mol, const UInt& nSpin);

			/**
			 * destructor
			 */
			~MO() { };

			/**
			 * whether there are redundancy in the basis sets so that
			 * the number of MO does not equal to the number of basis sets?
			 */
			bool isRedundantBasis() const {
				// we only test the alpha state
				// the beta one should be same with alpha
				if (mtrxList[0].getCol() != mtrxList[0].getRow()) return true;
				return false;
			};

			/**
			 * return the number of MO
			 */
			UInt getNOrb(UInt iSpin) const { 
				if (iSpin == 0) return mtrxList[0].getCol(); 
				return mtrxList[1].getCol(); 
			};

			/**
			 * return the number of basis set
			 */
			UInt getNBas(UInt iSpin) const { 
				if (iSpin == 0) return mtrxList[0].getRow(); 
				return mtrxList[1].getRow(); 
			};

			/**
			 * return the number occupied orbitals
			 */
			UInt getNOcc(UInt iSpin) const { 
				if (iSpin == 0) return nAlpha; 
				return nBeta;
			};

			/**
			 * return the LUMO index 
			 * the energy is in ascending order
			 * for the MO 
			 */
			UInt getLUMOIndex(UInt iSpin) const {
				if (iSpin == 0) return nAlpha; 
				return nBeta;
			};

			/**
			 * get the E(LUMO) - E(HOMO) for the given spin state
			 */
			Double getLUMOHOMODiff(UInt iSpin) const; 				

			/**
			 * write the result MO data into Disk files
			 */
			void writeToDisk() const;

			/**
			 * recover the mo data from disk files
			 * dir is the given mo file path
			 */
			void recover(const string& moPath); 

			///
			/// get the energy vector
			///
			const DoubleVec& getEnergyVec(UInt iSpin) const {
				if (iSpin == 0) return alphaEnergy;
				return betaEnergy;
			};

			///
			/// get the energy vector in non-const way
			///
			DoubleVec& getEnergyVec(UInt iSpin) {
				if (iSpin == 0) return alphaEnergy;
				return betaEnergy;
			};

			///
			/// form MO data for given Fock matrix
			/// the orthogonal matrix will be given through oneEMtrx
			///
			/// we note, the ortho matrix must be full matrix
			///
			void formMO(const MolShell& ms, const Molecule& mol, const OneEMtrx& oneEMtrx, SpinMatrix& Fock);

			///
			/// form MO data for given Fock matrix
			///
			/// this function is used for generating core guess 
			///
			/// the beta part is sharing the alpha mo data
			/// 
			/// we note, the ortho matrix must be full matrix
			///
			void formMO(const Mtrx& ortho, Mtrx& Fock);

			///
			/// form the fractional spined MO data based on the input infor
			///
			/// here one thing we need to note here. There are two source of molecule information for 
			/// this mo to refer, one is the normal molecule, the other one is the "newMol" contained
			/// in the FracSpinInfor class. Here if you want to form the fractional spin mo by calling
			/// this function, you must pass the "normal molecule" information first to the contructor
			/// so that the mo knows the normal alpha and beta electron numbers. This will help to 
			/// determine the HOMO index used in this function.
			///
			/// after the mo is scaled according to the fractional spin information, we will modify
			/// the alpha and beta electron numbers according to the newMol contained in the FracSpinInfor
			/// class. So be careful about this case!!!
			///
			void formFracSpinMO(const FracSpinInfor& infor);

			///
			/// printing function 
			///
			/// the level is set to 3, only level = 1, 2, 3 accept. Other numbers do not
			/// recognized
			/// 
			/// level == 1: print out the MO energy 
			/// 
			/// level == 2: print out the MO matrix data
			/// 
			/// level == 3: print out the MO energy and matrix data
			///
			///
			void print(UInt level) const;
	};

}

#endif
