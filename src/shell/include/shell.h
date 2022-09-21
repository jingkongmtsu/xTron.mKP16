/**
 * \file    shell.h
 * \brief   describing the shell structure used in the program
 * \author  Fenglai Liu and Jing Kong
 *
 * For the shell data, we follow the natrual way to do it:
 * shell -> atom shell -> molecule shell
 *
 * \note
 *
 * currently, the spliting mode for shells (that is, we are 
 * able to form a subshell from a shell) is not set up yet.
 * will do it in the future.
 *
 * keep an eye on the operator == for MolShell class. In the 
 * future, we may need to revise it.
 */
#ifndef SHELL_H
#define SHELL_H
#include "libgen.h"
#include "atom.h"
#include <string>
#include "scalablevec.h"
#ifdef SHELL_PHI_DEBUG
#include "atomshelltypedatainfor.h" 
#endif

namespace parseshell {
	class RawShell;
	class RawAtomShell;
	class RawMolShell;
}

namespace molecule {
	class Molecule;
}

namespace shell {

	using namespace parseshell;
	using namespace atom;
	using namespace molecule;
	using namespace std;
#ifdef SHELL_PHI_DEBUG
	using namespace atomshelltypedatainfor; 
#endif

	///
	/// declare class
	///
	class Shell;
	class AtomShell;

	///
	/// typedef the vector of shell and atom shells
	///
	typedef std::vector<Shell,tbb::scalable_allocator<Shell> > ShellVec;
	typedef std::vector<AtomShell,tbb::scalable_allocator<AtomShell> > AtomShellVec;

	/**
	 * \class   Shell
	 * \brief   data sctructure for unit shell
	 */
	class Shell {

		private:

			bool isDummy;          ///< whether the shell is dummy shell
			UInt localShellIndex;  ///< the local shell index
			UInt globalShellIndex; ///< the global shell index across the whole molshell
			UInt normBasStart;     ///< The starting normal basis set index for this shell
			UInt cartBasStart;     ///< The starting Cartesian basis set index for this shell
			UInt localNormBasStart;///< The local version (only in this atom shell) of normBasStart
			UInt localCartBasStart;///< The local version (only in this atom shell) of cartBasStart
			UInt nBas;             ///< number of basis sets
			UInt nCarBas;          ///< number of cartesian format basis sets
			UInt lmin;             ///< minimum type of angular momentum in the shell
			UInt lmax;             ///< maximum type of angular momentum in the shell
			UInt nPrim;            ///< number of primitive functions
			DoubleVec c;           ///< coefficients array (number of primitives*(lmax-lmin+1))
			DoubleVec e;           ///< exponents array (number of primitives)

		public:

			/**
			 * form the shell data from  the raw shell
			 * \param sr data from reading in shell
			 */
			Shell(const RawShell& rs); 

			/**
			 * \brief form an empty shell
			 * index are set to be -1 for avoiding reaching it
			 */
			Shell():isDummy(false),localShellIndex(-1),globalShellIndex(-1),normBasStart(-1),
			cartBasStart(-1),localNormBasStart(-1),localCartBasStart(-1),
			nBas(0),nCarBas(0),lmin(0),lmax(0),nPrim(0) { };

			/**
			 * form a shell data used to hold the memory
			 * \param maxNP  maximum number of primitive functions
			 * \param maxNL  maximum number of L contraction, like D is 1, SP = 2, SPD = 3 etc.
			 */
			Shell(UInt maxNP, UInt maxNL):isDummy(false),localShellIndex(-1),globalShellIndex(-1),
			normBasStart(-1),cartBasStart(-1),localNormBasStart(-1),localCartBasStart(-1),
			nBas(0),nCarBas(0),lmin(0),lmax(0),nPrim(0),
			c(maxNP*maxNL,ZERO),e(maxNP,ZERO) { };

			/**
			 * automatic destructor
			 */
			~Shell() { };

			/**
			 * form a dummy shell
			 */
			void formDummyShell();

			/**
			 * forming the basis set index for this shell
			 */
			void formShellOffset(UInt basOffset, UInt carBasOffset) {
				normBasStart = basOffset;     
				cartBasStart = carBasOffset;    
			};

			/**
			 * forming shell index 
			 * this call must be made from its parent atom shell
			 */
			void formShellIndex(UInt offset, UInt index) {
				localShellIndex = index;
				globalShellIndex = offset+index;
			};

			/**
			 * forming the local basis set index for this shell
			 * this should be only called in atom shell class
			 */
			void formLocalShellOffset(UInt basOffset, UInt carBasOffset) {
				localNormBasStart = basOffset;     
				localCartBasStart = carBasOffset;    
			};

			/**
			 * this function will reset the shell data from the compact form
			 * \param lmin0, lmax0: the shell Lmin and Lmax
			 * \param shellStatus : whether the shell is Cartesian or spherical(pure)?
			 * \param np          : number of primitives
			 * \param coef        : the coefficients array for this shell
			 * \param primExp     : exponential factors for the primitive functions
			 */
			void resetShellData(const UInt& lmin0, const UInt& lmax0, const UInt& shellStatus, 
					const UInt& np, const Double* coef, const Double* primExp);

			/**
			 * debugging printing 
			 */
			void print(UInt iprint=0) const;

			/**
			 * get the min L
			 */
			UInt getLmin() const {
				return lmin;
			};

			/**
			 * get the max L
			 */
			UInt getLmax() const {
				return lmax;
			};

			/**
			 * get the number of primitive Gaussian functions
			 */
			UInt getNPrim() const {
				return nPrim;
			};

			/**
			 * whether the shell is in spherical form, 5D,7F etc.
			 * we note, that S and P shell should be treated as Cartesian
			 */
			bool isPure() const {
				if (lmax<=1) return false;
				if (nBas == nCarBas) {
					return false;
				}else{
					return true;
				}
			};

			/**
			 * whether the shell is dummy or not
			 */
			bool isDummyShell() const {
				return isDummy;
			};

			/**
			 * whether the shell is empty?
			 */
			bool empty() const {
				return (nBas == 0);
			};

			/**
			 * get the exponent array 
			 */
			const Double* getExp() const {
				return &e.front();
			};

			/**
			 * get the coefficient array 
			 */
			const Double* getCoe() const {
				return &c.front();
			};

			/**
			 * get the coefficient array by giving the angular momentum L
			 */
			const Double* getCoe(const UInt& L) const {
				return &c[(L-lmin)*getNPrim()];
			};

			/**
			 * get the number of normal basis set
			 */
			UInt getNBas() const {
				return nBas;
			};

			/**
			 * get the number of basis set according to the type given
			 */
			UInt getNBas(const UInt& type) const {
				if (type == TYPE_CART) return nCarBas;
				return nBas;
			};

			/**
			 * get the number of Cartesian basis set
			 */
			UInt getNCarBas() const {
				return nCarBas;
			};

			/**
			 * whether it's compound shell
			 */
			bool isCompoundShell() const {
				if (lmin == lmax) {
					return false;
				}else{
					return true;
				}
			};

			/**
			 * whether the shell is in high L (L>=D)?
			 */
			bool isHighL() const { return (lmax>=2); };

			/**
			 * This function return the basis set index in global sense 
			 * (that is to say, in the whole molecule shell order)
			 * iBas is the basis set index in this shell
			 */
			UInt getBasisIndex(const UInt& iBas, const UInt& type) const;

			/**
			 * This function return the basis set index in local sense 
			 * (that is to say, in it's parent atom shell order)
			 * iBas is the basis set index in this shell
			 */
			UInt getLocalBasisIndex(const UInt& iBas, const UInt& type) const;

			/**
			 * This function return the local shell index within an atom shell
			 */
			UInt getLocalShellIndex() const { return localShellIndex; };

			/**
			 * This function return the global shell index in an molecule shell
			 */
			UInt getGlobalShellIndex() const { return globalShellIndex; };

			/**
			 * operator for comparing two shells in angular momentum
			 */
			bool operator>=(const Shell& s) const {
				if(lmax == s.getLmax()) {
					return (lmin >= s.getLmin()); 
				}else{ 
					return (lmax >= s.getLmax());
				}
			};

			/**
			 * operator for comparing two shells in angular momentum
			 * will be used in the stable sort
			 */
			bool operator<(const Shell& s) const {
				if(lmax == s.getLmax()) {
					return (lmin < s.getLmin());
				}else {
					return (lmax < s.getLmax());
				}
			};

			/**
			 * this is a debug function, to see whether two shell data
			 * are exactly same
			 */
			bool sameShellData(const Shell& s) const;
	};

	/**
	 * \class   AtomShell
	 * \brief   Shell data for the given ATOM
	 */
	class AtomShell : public Atom {

		private:

			bool allCartesian;  ///< whether shells with L>=D are in Cartesian form
			bool allSpherical;  ///< whether shells with L>=D are in pure form
			UInt nBas;          ///< number of basis sets
			UInt nCarBas;       ///< number of basis sets in cartesian format
			UInt normBasStart;  ///< The starting normal basis set index for this atom shell
			UInt cartBasStart;  ///< The starting Cartesian basis set index for this atom shell
			UInt maxL;          ///< maximum type of angular momentum in this shell
			UInt maxNP;         ///< maximum contraction degree
			UInt shlIndexStart; ///< the beginning index of shell in terms of whole molecule shell
			UInt nShells;       ///< number of shells for this atom shell
			ShellVec shells;    ///< shell data

			/**
			 * from the local shell/basis set offset for this atom
			 * this is done right after the whole atom shell is finished
			 * should be called each time you want to swap shells inside
			 * the atom shells
			 */
			void formLocalShellOffset();

		public:

			/**
			 * constructor for atom shell
			 * \param  as:     raw shell information
			 * \param  atom:   atom information
			 * \param  splitSPShell: do we need to split the sp shell into a S shell and a P shell
			 */
			AtomShell(const RawAtomShell& as, const Atom& atom, bool splitSPShell);

			/**
			 * construct to build an empty atom shell
			 * here index is set to -1 so that to avoid accidentally reach it
			 * so as the basis set starting index
			 */
			AtomShell():allCartesian(false),allSpherical(false),
			nBas(0),nCarBas(0),normBasStart(-1),cartBasStart(-1),
			maxL(0), maxNP(0),shlIndexStart(-1),nShells(0) { };

			/**
			 * constructor to build atom shell data to hold the memory
			 * \param maxNP  maximum number of primitive functions
			 * \param maxNL  maximum number of L contraction, like D is 1, SP = 2, SPD = 3 etc.
			 * \param maxNShell  maximum number of shells
			 */
			AtomShell(UInt maxNP, UInt maxNL, UInt maxNShell):allCartesian(false),allSpherical(false),
			nBas(0),nCarBas(0),normBasStart(-1),cartBasStart(-1),
			maxL(0),maxNP(0),shlIndexStart(-1),nShells(0) { 
				shells.reserve(maxNShell);
				Shell s(maxNP,maxNL);
				for(UInt i=0; i<maxNShell; i++) shells.push_back(s);
			};

			///
			/// whether it's empty?
			///
			bool empty() const { return (getNShell() == 0); };

			///
			/// destructor
			///
			~AtomShell() { };

			/**
			 * debug printing 
			 */
			void print(UInt iprint=0) const;

			///
			/// is the given atom shell is dummy shell?
			///
			bool isDummyAtomShell() const {
				if (getNShell() == 1 && shells[0].isDummyShell()) return true;
				return false;
			};

			///
			/// form a dummy shell from empty shell
			///
			void formDummyAtomShell();

			/**
			 * from the shell offset for this atom
			 * this is only for global basis set offset 
			 * which is determined by molShell
			 * therefore this function should be only called
			 * from molshell object
			 */
			void formShellOffset(UInt basOffset, UInt carBasOffset);

			/**
			 * this function is forming the shell index for all of shell data
			 * in the shell vector. for using this function, we assume that
			 * shelIndexStart has been properly initilized. we will check it 
			 * before forming shell index
			 */
			void formShellIndex();

			/**
			 * this function will restore the atom shell data from the input data 
			 *
			 * this is to convert the shell data from "compact" form into normal form
			 * with this function, the user should build the atom shell object with the 
			 * "memory holder" constructor. Inside we do not change the memory, we just
			 * use the one already being allocated
			 *
			 * \param atomShellIndex  the index for this atom shell in the whole molecule shell
			 * \param coord           the geometry data, we use atomShellIndex to retrieve the xyz of the atom
			 * \param compactInfor    the compact data of this atom shell
			 * \param cartBasOffset   the cartesian basis set offset data for this atom shell
			 * \param normBasOffset   the normal    basis set offset data for this atom shell
			 * \param shellOffset     the global shell offset for this atom shell
			 */
#ifdef SHELL_PHI_DEBUG
			void formAtomShellFromCompactForm(const UInt& atomShellIndex, const Double* coord, 
					const AtomShellTypeDataInfor& compactInfor, 
					UInt cartBasOffset, UInt normBasOffset, UInt shellOffset);
#endif

			/**
			 * get shell data for the corresponding ishell
			 * ishell is the index for the shell array
			 */
			const Shell& getShell(const UInt& iShell) const {
				return shells[iShell];
			};

			/**
			 * return the atom index
			 */
			UInt getAtomShellIndex() const { return atomIndex; };

			/**
			 * get the number of shells
			 */
			UInt getNShell () const {
				return nShells;
			};

			/**
			 * get the number of basis sets
			 */
			UInt getNBas () const {
				return nBas;
			};

			/**
			 * get the number of basis sets according to the type given
			 */
			UInt getNBas (const UInt& type) const {
				if (type == TYPE_NORM) return nBas;
				return nCarBas;
			};

			/**
			 * get the number of basis sets in cartesian format
			 */
			UInt getNCarBas() const {
				return nCarBas;
			};

			/**
			 * get basis set starting index according to the type
			 */
			UInt getBasisStartIndex(const UInt& type) const;

			/**
			 * get shell starting index for this atom
			 */
			UInt getShellStartIndex() const { return shlIndexStart; };

			/**
			 * update shell starting index for this atom
			 * this function can be only called from MolShell class!
			 */
			void updateShellStartIndex(UInt index) { shlIndexStart = index; };

			/**
			 * get the maximum angular momentum
			 */
			UInt getMaxL() const {
				return maxL;
			};

			/**
			 * get the maximum contraction degree
			 */
			UInt getMaxNP() const {
				return maxNP;
			};

			/**
			 * do we have high L (L>=D)?
			 */
			bool hasHighL() const;

			/**
			 * whether the given atom shell is all Cartesian?
			 */
			bool allCart() const {
				return allCartesian;
			};

			/**
			 * whether the given atom shell is all pure?
			 */
			bool allPure() const {
				return allSpherical;
			};

			/**
			 * return the composite shell degree
			 * for example,
			 * S   shell: 1
			 * SP  shell: 2
			 * SPD shell: 3
			 * etc.
			 */
			UInt getMaxCompDegree() const;

			/**
			 * this is a debug function
			 * it's used to test that whether the two atom shells are same
			 */
			bool sameAtomShell(const AtomShell& as) const;
	};

	/**
	 * \class   MolShell
	 * \brief   Shell data for the given MOLECULE
	 *
	 * In the mol shell class, we use a "code" to represent the basis set.
	 * the code is generated in the shellprop.h
	 */
	class MolShell {

		private:

			bool allCartesian;             ///< whether shells with L>=D are in Cartesian form
			bool allSpherical;             ///< whether shells with L>=D are in pure form
			UInt nShell;                   ///< number of shells
			UInt nBas;                     ///< number of basis sets
			UInt nCarBas;                  ///< number of basis sets in cartesian format
			UInt maxL;                     ///< maximum type of angular momentum in this shell
			UInt maxNP;                    ///< maximum contraction degree
			UInt code;                     ///< the basis set code 
			AtomShellVec atomShells;       ///< shell data in atom

			/**
			 * this is the main function to form the molShell data from the raw data
			 */
			void formMolShellData(const Molecule& mol, const RawMolShell& rms);

			/**
			 * form the basis set index for all of shells
			 */
			void formShellIndex();

			/**
			 * appending the atom shell data into the molshell data
			 * we do not provide function to resort shell index inside
			 * it's most likely that you still want to keep it (for example,
			 * when you construct an empty shell and form a subset of some
			 * given shell).
			 *
			 * Therefore, if you want to reform shell index, you need to do
			 * it in other way
			 */
			void appendAtomShell(const AtomShell& atomShl);

		public:

			/**
			 * constructor to form shell from input file
			 * \param  input  user input file
			 * \param  part   which part of basis set we are going to read in? 
			 *                part = 0, main basis set; part >= 1, aux basis set
			 * \param  mol0   molecular geometry
			 * \param  debugRawShell  if you want to switch on the debug of the raw shell data, turn this on
			 */
			MolShell(const string& input, const Molecule& mol0, UInt part = 0, bool debugRawShell = false);

			/**
			 * constructor for empty MolShell - shell code is unreachable
			 */
			MolShell():allCartesian(true),allSpherical(true),nShell(0),nBas(0),
			nCarBas(0),maxL(0),maxNP(0),code(-1) { };

			///
			/// construct on atom shell
			///
			MolShell(const AtomShell& as);

			///
			/// whether it's empty?
			///
			bool empty() const { return (getNAtomShells() == 0); };

			/**
			 * form a dummy shell based on the empty shell
			 */
			void formDummyShell();

			///
			/// destructor
			///
			~MolShell() { };

			///
			/// debug printing
			///
			void print(UInt iprint=0) const;

			/**
			 * compare that whether s is identical to this shells
			 * note that we only use the code of basis set!
			 */
			bool operator==(const MolShell& s) const {
				return (code == s.code); 
			};

			/**
			 * whether this is dummy shell?
			 */
			bool dummy() const {
				if (getNAtomShells() == 1 && atomShells[0].isDummyAtomShell()) return true;
				return false;
			};

			/**
			 * get the shell code
			 */
			UInt getCode() const { return code;};

			/**
			 * get atom shell data 
			 */
			const AtomShell& getAtomShell(const UInt& i) const {
				return atomShells[i];
			};

			/**
			 * get the length of atom shells
			 */
			UInt getNAtomShells() const {
				return atomShells.size();
			};

			/**
			 * get the number of shells
			 */
			UInt getNShell() const {
				return nShell;
			};

			/**
			 * get the number of basis sets
			 */
			UInt getNBas() const {
				return nBas;
			};

			/**
			 * get the number of basis sets according to the type given
			 */
			UInt getNBas (const UInt& type) const {
				if (type == TYPE_CART) return nCarBas;
				return nBas;
			};

			/**
			 * get the number of basis sets in cartesian format
			 */
			UInt getNCarBas() const {
				return nCarBas;
			};

			/**
			 * get max angular momentum 
			 */
			UInt getMaxL() const {
				return maxL;
			};

			/**
			 * get max contraction degree
			 */
			UInt getMaxNP() const {
				return maxNP;
			};

			/**
			 * whether it's all cartesian?
			 */
			bool allCart() const {
				return allCartesian;
			};

			/**
			 * whether it's all pure?
			 */
			bool allPure() const {
				return allSpherical;
			};

			/**
			 * whether we have the given angular momentum in molecule shell data?
			 */
			bool hasAng(const UInt& lmin, const UInt& lmax) const; 

			/**
			 * return the number of types in the atom shell list
			 */
			UInt getNAtomShellTypes() const;

			/**
			 * return the maximum number of basis set for the atom
			 * this is done by loop over atom types
			 * since each type of atom has same basis set functions
			 * \param type  cart, pure or norm?
			 */
			UInt getMaxBasis(const UInt& type) const;

			/**
			 * return the maximum number of shells for the atom
			 * this is done by loop over atom types
			 * since each type of atom has same basis set functions
			 */
			UInt getMaxNShell() const;

			/**
			 * get the maximum number of sub shell across the molecule shell
			 */
			UInt getMaxLContraction() const;
	};

}

#endif
