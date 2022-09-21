/**
 * \file  zmatrix.h
 * \brief This class handles the zmatrix geometry data
 * \author Fenglai Liu 
 */
#ifndef __ZMATRIX_H__
#define __ZMATRIX_H__
#include<vector>
#include "libgen.h"
#include "tbb/tbb.h"
#include "scalablevec.h"

namespace zmatrix {

	///
	/// declare the ZMatrixItem class
	///
	class ZMatrixItem;

	///
	/// here we define the vector of density matrices
	///
	typedef std::vector<ZMatrixItem,tbb::scalable_allocator<ZMatrixItem> > ZMatrixVec;

	///
	///  this class describe one item of z matrix information (one item means for one atom)
	///
	class ZMatrixItem {

		private:

			UInt atomic;///< the atomic value for this atom (atom1)
			UInt atom2; ///< the atom index that forms bond: bond distance is (atom1,atom2)
			UInt atom3; ///< the atom index that forms angle: angle is (atom1,atom2,atom3)
			UInt atom4; ///< the atom index that forms dihedral: dihedral defined as (atom1,atom2,atom3,atom4)
			Double bond;     ///< the bond distance value
			Double angle;    ///< the angle value formed by (atom1,atom2,atom3)
			Double dihedral; ///< the dihedral value formed by (atom1,atom2,atom3,atom4)

		public:

			///
			/// contructor, we directly fill in all of values
			///
			ZMatrixItem(UInt atomic0, UInt atom2Index, UInt atom3Index, UInt atom4Index,
					Double bondVal, Double angleVal, Double dihedralVal):atomic(atomic0),
			atom2(atom2Index),atom3(atom3Index),atom4(atom4Index),bond(bondVal),
			angle(angleVal),dihedral(dihedralVal) { };

			///
			/// contructor, this is the first atom so that only atom1 Index
			/// is provided
			///
			ZMatrixItem(UInt atomic0):atomic(atomic0),atom2(-1),atom3(-1),atom4(-1),bond(ZERO),
			angle(ZERO),dihedral(ZERO) { };

			///
			/// contructor, this is the second atom so that only atom1 and atom2 present
			///
			ZMatrixItem(UInt atomic0, UInt atom2Index, Double bondVal):atomic(atomic0),
			atom2(atom2Index),atom3(-1),atom4(-1),bond(bondVal),angle(ZERO),dihedral(ZERO) { };

			///
			/// contructor, this is the third atom so that we only have bond and angle
			/// information provided
			///
			ZMatrixItem(UInt atomic0, UInt atom2Index, UInt atom3Index, 
					Double bondVal, Double angleVal):atomic(atomic0),
			atom2(atom2Index),atom3(atom3Index),atom4(-1),bond(bondVal),
			angle(angleVal),dihedral(ZERO) { };

			///
			/// desctructor
			///
			~ZMatrixItem() { };

			///
			/// whether this is the first atom?
			///
			bool is1stAtom() const {
				UInt nullVal = static_cast<UInt>(-1);
				if (atom2 == nullVal && atom3 == nullVal && atom4 == nullVal) return true;
				return false;
			};

			///
			/// whether this is the second atom?
			///
			bool is2edAtom() const {
				UInt nullVal = static_cast<UInt>(-1);
				if (atom2 != nullVal && atom3 == nullVal && atom4 == nullVal) return true;
				return false;
			};

			///
			/// whether this is the third atom?
			///
			bool is3rdAtom() const {
				UInt nullVal = static_cast<UInt>(-1);
				if (atom2 != nullVal && atom3 != nullVal && atom4 == nullVal) return true;
				return false;
			};

			///
			/// whether this is the normal atom index, with all of data field defined?
			///
			bool isNormalAtom() const {
				UInt nullVal = static_cast<UInt>(-1);
				if (atom2 != nullVal && atom3 != nullVal && atom4 != nullVal) return true;
				return false;
			};

			///
			/// return this atom index
			///
			UInt getAtomic() const { return atomic; };

			///
			/// return second atom index
			///
			UInt getAtom2Index() const { return atom2; };

			///
			/// return third atom index
			///
			UInt getAtom3Index() const { return atom3; };
			
			///
			/// return fourth atom index
			///
			UInt getAtom4Index() const { return atom4; };

			///
			/// return bond distance
			///
			Double getBondDistance() const { return bond; };

			///
			/// return angle value
			///
			Double getAngleVal() const { return angle; };

			///
			/// return dihedral value
			///
			Double getDihedralVal() const { return dihedral; };

			///
			/// debug print
			///
			void print() const;
	};

	///
	///  this class describe the whole zmatrix set for the given molecule
	///
	class ZMatrix {

		private:

			ZMatrixVec zmatrixList;  ///< the zmatrix coordinate system

			///
			/// for an input stream, return the location that the given 
			/// geometry section is defined
			///
			UInt getSecLoc(ifstream& inf, UInt section) const;

			///
			/// read in variable values if in zmatrix, there are variables
			/// contained
			///
			void readInVarValues(ifstream& inf, UInt loc, 
					vector<string>& vars, DoubleVec& vals) const;

			///
			/// read in the z matrix geometry data
			///
			void readInGeom(ifstream& inf, UInt loc,
					const vector<string>& vars, const DoubleVec& vals);

		public:

			///
			/// contructor, we read in the data from the given file
			///
			ZMatrix(const string& file, UInt section);

			///
			/// desctructor
			///
			~ZMatrix() { };

			///
			/// debug print
			///
			void print() const;

			///
			/// return the number of atoms
			///
			UInt getNAtoms() const { return zmatrixList.size(); };

			///
			/// return the z matrix item
			///
			const ZMatrixItem& getZMatrix(UInt i) const { return zmatrixList[i]; };

			///
			/// change the z matrix data into input XYZ coordinates
			///
			void loadInToXYZ(DoubleVec& coord) const;
	};

}

#endif
