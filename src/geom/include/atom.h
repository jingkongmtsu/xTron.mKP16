/**
 * \file  atom.h
 * \brief This class contains the data to define atom
 *
 * Natrually, for each atom's data it is represented by items below:
 *   - atomic number;
 *   - atomic symbol(equivalent to atomic number);
 *   - atomic coordinates
 *   - atom index in the molecule
 *   - etc.
 *
 * \author Fenglai Liu and Jing Kong
 */

#ifndef ATOM_H
#define ATOM_H
#include "libgen.h"
#include<cmath>

namespace atom {

	/**
	 * \class Atom
	 * \brief Atoms in chemistry
	 */
	class Atom {

		protected:

			Double xyz[3];   ///< coordinates 
			UInt atomic;     ///< atomic number
			UInt atomIndex;  ///< atom index ine molecule

		public:

			/**
			 * constructor 
			 * in default the atom index is set to be -1
			 * which means unreachable
			 * this is logical because in forming atom object
			 * we do not need to know it's position in molecule
			 * this will be formed later
			 */
			Atom(const Double* xyz0, const UInt& atomic0) {
				xyz[0] = xyz0[0];
				xyz[1] = xyz0[1];
				xyz[2] = xyz0[2];
				atomic = atomic0;
				atomIndex = -1;
			};

			/**
			 * default constructor
			 * in default atomic number is set to something unreachable
			 * atom index is also unreachable
			 */
			Atom() {
				xyz[0] = ZERO;
				xyz[1] = ZERO;
				xyz[2] = ZERO;
				atomic = -1;
				atomIndex = -1;
			};

			/**
			 * destructor
			 */
			~Atom(){ };

			/**
			 * define the operator <
			 */
			bool operator<(const Atom& a) const {
				return (atomic<a.atomic);
			};

			/**
			 * update the atom index information
			 * this should be only called from molecule or 
			 * it's derived class!
			 */
			void updateAtomIndex(const UInt& index) {
				atomIndex = index;
			};

			/**
			 * get the xyz of this atom
			 */
			const Double* getXYZ() const {
				return xyz;
			};

			/**
			 * get atomic number
			 */
			UInt getAtomic() const {
				return atomic;
			};

			/**
			 * get atom index
			 */
			UInt getAtomIndex() const {
				return atomIndex;
			};

			///
			/// get distance between an arbitrary point with this atom
			///
			Double getDistance(const Double* P) const {
				Double d2 = (xyz[0]-P[0])*(xyz[0]-P[0]) + (xyz[1]-P[1])*(xyz[1]-P[1])
					+ (xyz[2]-P[2])*(xyz[2]-P[2]);
				return sqrt(d2);
			};

			/**
			 * sometimes we may set the xyz to a new one
			 */
			void updateXYZ(const Double& X, const Double& Y, const Double& Z) {
				xyz[0] = X;
				xyz[1] = Y;
				xyz[2] = Z;
			};

	};
}

#endif
