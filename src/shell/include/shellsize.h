/**
 * \file    shellsize.h
 * \brief   collecting shell size data 
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef SHELLSIZE_H
#define SHELLSIZE_H
#include "libgen.h"
#include <vector>
#include "scalablevec.h"
#include "tbb/scalable_allocator.h"

namespace shell {
	class Shell;
	class AtomShell;
	class MolShell;
}

namespace shellsize {

	using namespace shell;
	using namespace std;

	///
	/// declare
	///
	class AtomShellSize;

	///
	/// typedef the vector of atom shell size
	///
	typedef std::vector<AtomShellSize,tbb::scalable_allocator<AtomShellSize> >   AtomShellSizeVec;

	/**
	 * collect shell radius as well as determine atom size
	 * for a given atom type
	 */
	class AtomShellSize {

		private:

			UInt atomic;                   ///< atom type
			Double atomRadius;             ///< atom radius
			DoubleVec shellRadius;         ///< shell radius for all shells in this atom

			///
			/// creating shell radius of r by determining that:
			/// \f$\psi_{L}(r) = 0\f$
			/// \psi is the shell which has angular momentum of L. For the
			/// composite shell, we just return the maximum value we can
			/// find among the sub-shells.
			///
			/// the \f$\psi_{L}(r)\f$ function actually we take the form of :
			/// \f$\psi_{L}(r) = \sum_{i}c_{i}r^{L}e^{-\alpha_{i}r^{2}}\f$
			/// it's obvious that this function will have largerr shell radius (
			/// so as atom radius) than the normal Cartesian type of basis set functions.
			/// This is because that for r, it's r>=x, r>=y and r>=z.
			/// 
			/// the function uses Newton-Rhaphson process to compute shell radius. The real
			/// function used is actually the log format of psi, which is 
			/// \f$\log(\psi_{L}(r))\f$.
			///
			Double setShellRadius(const Shell& s, const Double& thresh) const;

		public:

			///
			/// constructor - form the shell size inside
			///
			AtomShellSize(const AtomShell& s, const Double& thresh);

			///
			/// default destructor
			///
			~AtomShellSize(){ };

			/**
			 * print for debug choice
			 */
			void print() const;

			/**
			 * get the shell radius
			 * \param  iShell   local shell index within the atom - starting from 0
			 */
			Double getShellRadius(const UInt& iShell) const {
				return shellRadius[iShell];
			}; 

			/**
			 * get the atom size
			 */
			Double getAtomSize() const { return atomRadius; };

			/**
			 * get the atomic number
			 */
			UInt getAtomic() const { return atomic; };

			/**
			 * for given distance and another shell size data,
			 * let's determine that whether the two shells are out of 
			 * shell radius?
			 * \param iShell  the shell index (local) for this atom shell size object
			 * \param jShell  the shell index (local) for another atom shell size object
			 * \param size2   another atom shell size object
			 * \param dis     distance that we want to compare
			 */
			bool outofShellRadius(const UInt& iShell, const UInt& jShell, 
					const AtomShellSize& size2,const Double& dis) const {
				Double thisShellRadius = getShellRadius(iShell);
				Double thatShellRadius = size2.getShellRadius(jShell);
				if (dis>=thisShellRadius+thatShellRadius) return true;
				return false;
			};

			/**
			 * for given distance and another shell size data,
			 * let's determine that whether we are out of atom radius?
			 * \param size2   another atom shell size object
			 * \param dis     distance that we want to compare
			 */
			bool outofAtomRadius(const AtomShellSize& size2, const Double& dis) const {
				Double thisRadius = getAtomSize();
				Double thatRadius = size2.getAtomSize();
				if (dis>=thisRadius+thatRadius) return true;
				return false;
			};
	};

	/**
	 * collect the shell size for a given molecular shell
	 */
	class MolShellSize {

		private:

			AtomShellSizeVec atomShellSizes;   ///< record atom distance in terms of element type

		public:

			///
			/// constructor
			///
			MolShellSize(const MolShell& ms,const Double& thresh);

			///
			/// default destructor
			///
			~MolShellSize(){ };

			/**
			 * print for debug choice
			 */
			void print() const;

			/**
			 * get the atom shell size information by given atomic type
			 */
			const AtomShellSize& getAtomShellSizeInfor(const UInt& atomic) const; 

			/**
			 * get the number of atom shell data
			 */
			UInt getNAtomShellData() const { return atomShellSizes.size(); };

			/**
			 * here we provide a easy function to reach the atom size
			 * just provide the atomic number
			 */
			Double getAtomSize(const UInt& atomic) const {
				const AtomShellSize& infor = getAtomShellSizeInfor(atomic);
				return infor.getAtomSize();
			};

			/**
			 * extract the atom type information (by atomic number)
			 * from the vector of atom shell data
			 */
			void getAtomTypeInfor(UIntVec& infor) const;

			/**
			 * get the number of atom types
			 */
			UInt getNAtomTypes() const;
	};

}


#endif
