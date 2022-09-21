/**
 * \file  molecule.h
 * \brief This class contains the molecular data for the whole program
 * \author Fenglai Liu and Jing Kong
 */

#ifndef MOLECULE_H
#define MOLECULE_H
#include <vector>
#include <string>
#include "libgen.h"
#include "atom.h"
#include "tbb/scalable_allocator.h"
#include "scalablevec.h"
using namespace std;
using namespace atom;
using namespace tbb;

// define the vector of atoms
typedef std::vector<Atom,tbb::scalable_allocator<Atom> > AtomVec;

namespace molecule {

	/**
	 * \class molecule
	 *
	 * Here we need to explain the meaning of "section" in molecule class.
	 * Since in one job we may have many molecule data defined, therefore it's
	 * necessary to point out which section this data points to.  Section, is
	 * like the "ID" to each given molecule from the input file.  We note, that
	 * for cluster the section is 0, for molecule it's >= 1. 
	 */
	class Molecule {

		protected:

			///////////////////////////////////////////
			//           data definition             //
			///////////////////////////////////////////
			bool unrestricted;    ///< whether the molecule is force to be unrestricted?
			bool unitInAngstrom;  ///< whether the geometry data is in angstrom
			UInt section;         ///< section number for molecule
			UInt multiplicity;    ///< spin multiplicity
			UInt nalpha;          ///< number of alpha electrons
			UInt nbeta;           ///< number of beta electrons
			Int charge;           ///< charge
			AtomVec atoms;        ///< all of atom's data

			/**
			 * get the location of given section
			 * for section label > 0; it's molecule section
			 * label = 0 it's cluster section
			 *
			 * also here we will determine whether the geometry data unit is in Bohr or Angstrom
			 */
			UInt getSecLoc(ifstream& inf);

			/**
			 * reading geometry data from the given input stream. loc refers
			 * to the location of the molecule piece. We note that the 
			 * atom data is directly pushing into the atoms array 
			 */
			void readAtomDataInXYZ(ifstream& inf, const UInt& loc);

			/**
			 * reading geometry data from the given input file. loc refers
			 * to the location of the molecule piece. Z-matrix format of 
			 * data would be interally read in and parsed into XYZ format,
			 * and finally store the data into atoms array
			 */
			void readAtomDataInZMatrix(const string& input);

			/**
			 * general interface to read in geometry data 
			 */
			void readAtomData(const string& input, const UInt& loc);

			/**
			 * reading charge and multiplicity from the given input stream. 
			 */
			void readChargeMult(ifstream& inf, const UInt& loc);

			/**
			 * set nalpha and nbeta according to spin multiplicity and the charge value
			 */
			void setNAlphaBeta(); 

			/**
			 * checking the geometry data
			 * in case that we have two atoms that are too close with each other
			 * which must be wrong
			 */
			bool checkDistance() const;

			/**
			 * update the atom index information for each atom in the array
			 * this is done after atom data is complete
			 */
			void updateAtomIndexInfor();

		public:

			/**
			 * constructor to build the molecule data from input file
			 * \param   whichGeom which geometry we are referring
			 */
			Molecule(const string& input, UInt whichGeom = 1);

			/**
			 * constructor to build one atom molecule
			 * the atom is free atom, therefore a lot of its 
			 * data could be determined from element part of codes
			 *
			 * the only thing needs to input is atomic number
			 */
			Molecule(const UInt& atomic);

			/**
			 * destructor
			 */
			~Molecule() { };

			///
			/// for debug purpose
			///
			void print() const;

			/**
			 * get the corresponding atom
			 */
			const Atom& getAtom(UInt atom) const {
				return atoms[atom];
			};

			/**
			 * get the given atom xyz
			 */
			const Double* getXYZ(UInt atom) const;

			/**
			 * get the given atom atomic number
			 */
			UInt getAtomic(UInt atom) const;

			/**
			 * get the given atom atomic symbol
			 */
			string getSymbol(UInt atom) const; 

			/**
			 * we sort out the atom types and store it inside 
			 * the the array of atom types
			 */
			void getAtomTypes(UIntVec& atomTypes) const;

			/**
			 * return the number of atom types
			 */
			UInt getNAtomTypes() const;

			/**
			 * return the number of atoms 
			 */
			UInt getNAtoms() const {
				return atoms.size();
			};

			/**
			 * return the number of alpha electrons
			 */
			UInt getNAlpha() const {
				return nalpha;
			};

			/**
			 * return the number of beta electrons
			 */
			UInt getNBeta() const {
				return nbeta;
			};

			/**
			 * return the number of electrons depending 
			 * on the spin state given
			 */
			UInt getNEle(const UInt& iSpin) const {
				if (iSpin == 0) return nalpha;
				return nbeta;
			};

			/**
			 * return the total number of electrons 
			 */
			UInt getTotalEle() const {
				return nalpha+nbeta;
			};

			///
			/// reset the number of occupied orbitals for the given spin state
			/// this is used for special purpose such as fractional spin calculation
			/// Please use it carefully!!!
			///
			void changeEleNum(UInt iSpin, UInt n) {
				if (iSpin == 0) {
					nalpha = n;
				}else{
					nbeta  = n;
				}
			};

			/**
			 * return the charge
			 */
			Int getCharge() const {
				return charge;
			};

			/**
			 * return the multiplicity
			 */
			UInt getMult() const {
				return multiplicity;
			};

			/**
			 * return section number
			 */
			UInt getSec() const {
				return section;
			};

			/**
			 * sometimes you may want to update the section number
			 * and if this is want you want, here we provide the function
			 *
			 * use it with caution!
			 */
			void updateSec(UInt newSec) {
				section = newSec;
			};

			/**
			 * judging by the alpha and beta electrons
			 * number, we will decide the number of 
			 * spin states
			 */
			UInt getNSpin() const {
				if (unrestricted) return 2;
				if (nalpha != nbeta) {
					if (nalpha == 0 || nbeta == 0) return 1;
					return 2;
				}else{
					return 1;
				}
			};

			/**
			 * get the close shell status by counting the 
			 * nalpha and nbeta
			 */
			bool isCloseShell() const {
				if (unrestricted) return false;
				if (nalpha != nbeta) {
					return false;
				}else{
					return true;
				}
			};

			/**
			 * is it a single electron system?
			 */
			bool isSingEleSystem() const {
				return (nalpha + nbeta == 1);
			};

			/**
			 * what is the input geometry unit?
			 */
			bool isInputGeomDataInAng() const {
				return unitInAngstrom; 
			};

			///
			/// force unrestricted calculation
			///
			void forceUNRestricted();

			/**
			 * get the distance between two atoms
			 */
			Double getDistance(const UInt& atom1, const UInt& atom2) const;

			/**
			 * get the nuclear replusion energy
			 */
			Double getNuclearRepulsion() const;

			/**
			 * calculate the distance matrix between the atoms
			 * order 0: only distance between atoms is given
			 * order 1: derivatives for distance with respect to the x,y and z
			 */
			void getAtomDistanceMatrix(Double* dmatrix, const UInt& order) const; 
	};


}

#endif
