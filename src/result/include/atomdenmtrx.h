/**
 * \file    atomdenmtrx.h
 * \brief   classes used to generate the density matrix for free atom
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef ATOMDENMTRX_H
#define ATOMDENMTRX_H
#include "libgen.h"
#include "globalinfor.h"
#include "denmtrx.h"

namespace molecule {
	class Molecule;
}

namespace shell {
	class MolShell;
}

namespace atomdenmtrx {

	using namespace molecule;
	using namespace shell;
	using namespace denmtrx;
	using namespace globalinfor;

	class AtomDenMtrx {

		private:

			//
			// information center
			//
			GlobalInfor  infor;          ///< save a copy of global information center
			UInt section;                ///< the corresponding section of molecule

			//
			// these are the settings to generate the atom density matrix data
			//
			string scfEXMethod;          ///< the SCF exchange method to generate the atom density matrix
			string scfECMethod;          ///< the SCF exchange method to generate the atom density matrix
			UInt maxSCFCycles;           ///< maximum SCF cycles for SCF running
			bool scfRunInSilence;        ///< whether the atom SCF running in silence, for debug purpose

			//
			// result data
			//
			UIntVec atomTypes;           ///< atomic number used to represent atom types
			DenMtrxVec atomDenMtrices;   ///< atom density matrix data   

		public:

			/**
			 * constructing the atom density matrix
			 */
			AtomDenMtrx(const GlobalInfor& infor, const Molecule& mol);

			/**
			 * this is the build the atomic density matrix basic information
			 * this contructor is used later together with then recover function
			 * so that to recover the density matrix data inside
			 *
			 * right now this constructor only build an empty atom density matrix class
			 */
			AtomDenMtrx(const GlobalInfor& infor0, const UInt& sec):infor(infor0),
			section(sec),scfEXMethod("HF"),scfECMethod("NONE"),maxSCFCycles(100),scfRunInSilence(true) { };

			/**
			 * destructor
			 */
			~AtomDenMtrx() { };

			/**
			 * for each atom in the molecule, we will form it's SCF converged 
			 * density matrix 
			 */
			void formAtomDenMtrxInSCF(const Molecule& mol, const MolShell& ms);

			/**
			 * write the atom density matrix into hard disk
			 */
			void writeToDisk() const;

			/**
			 * recover the density matrix data from the binary files from hard disk
			 *
			 * we need the mol shell data to fetch the basis set information for each
			 * atom density matrix
			 */
			void recover(const MolShell& ms); 

			/**
			 * to average the alpha and beta density matrix for all of atoms
			 * in the vector
			 *
			 * Basically, it's new alpha = (alpha_den+beta_den)/2
			 *
			 * so alpha part of density matrix will be destroyed, beta is still 
			 * the same; We do it only for open shell atom 
			 */
			void averageSpinDenMtrx();

			/**
			 * get the atom density matrix data according to 
			 * given atomic number
			 */
			const DenMtrx& getAtomDenMtrx(const UInt& atomic) const;

			/**
			 * for debug printing
			 */
			void print() const;

	};
}

#endif

