/**
 * \file    sigatombasis.h
 * \brief   produce significant atom and basis set list
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef SIGATOMBASIS_H
#define SIGATOMBASIS_H
#include "libgen.h"
#include "scalablevec.h"

namespace shell {
	class AtomShell;
	class MolShell;
}

namespace shellsize {
	class MolShellSize;
}

namespace batchgrid {
	class BatchGrid;
}

namespace xcintsinfor {
	class XCIntsInfor;
}

namespace sigatombasis {

	using namespace shell;
	using namespace shellsize;
	using namespace batchgrid;
	using namespace xcintsinfor;
	using namespace std;

	class SigAtomBasis {

		private:

			///////////////////////////////////////////////
			//                data section               //
			///////////////////////////////////////////////

			// job information
			bool withAtomApprox;         ///< basis set evalulated only for atom=batchatom 

			// sigfinicant array information
			UIntVec sigAtoms;            ///< significant atoms list
			UIntVec sigBas;              ///< significant basis set list(in global normal basis set index)
			UIntVec sigShell;            ///< significant shell list(in local shell index order per atom)

			/**
			 * this array set up mapping between sigAtoms and sigBas. it's total
			 * length is nSigAtoms+1. for each significant atom, it stores the 
			 * index of sigBas:
			 * significant atom i (starting from 0)-> 
			 * atomBasisMapping[i]    : starting index in sigBas for this sig atom
			 * atomBasisMapping[i+1]  : ending index in sigBas for this atom
			 *
			 * for example, atomBasisMapping[3] = 7; it means that for for significant
			 * atom, it's significant basis set starts from sigBas[7]; and 
			 * atomBasisMapping[4] = 20, it means that the significant basis set for
			 * this atom ends on sigBas[20]
			 */
			UIntVec atomBasisMapping;     

			/**
			 * this array set up mapping between sigAtoms and sigShell. it's total
			 * length is nSigAtoms+1. for each significant atom, it stores the 
			 * index of sigShell:
			 * significant atom i (starting from 0)-> 
			 * atomShellMapping[i]    : starting index in sigShell for this sig atom
			 * atomShellMapping[i+1]  : ending index in sigShell for this atom
			 *
			 * for example, atomShellMapping[3] = 7; it means that for for significant
			 * atom, it's significant shell starts from sigShell[7]; and 
			 * atomShellMapping[4] = 20, it means that the significant shell for
			 * this atom ends on sigShell[20]. From sigShell[7] to sigShell[20] 
			 * it will gives the significant shells for this atom.
			 */
			UIntVec atomShellMapping;     

			/**
			 * forming significant atom and shell data under atom center approximation
			 * in this approximation, there's only one significant atom which is 
			 * the batch atom itself (use with withAtomApprox)
			 */
			void formSigDataWithAtomCenterApprox(const BatchGrid& grid, const MolShell& ms);

			/**
			 * forming significant atom and shell data in normal way
			 */
			void formSigData(const BatchGrid& grid, const MolShell& ms, 
					const MolShellSize& shellSize); 

		public:

			/**
			 * constructor - directly form the significant atom and 
			 * basis set list
			 */
			SigAtomBasis(const BatchGrid& grid, const XCIntJobInfor& infor,
					const MolShell& ms, const MolShellSize& shellSize); 

			/**
			 * destructor
			 */
			~SigAtomBasis() { };

			/**
			 * print function used in debug purpose
			 */
			void print() const;

			/**
			 * judge whether the given atom is significant
			 * iAtom is in global index in molshell
			 */
			bool isSigAtoms(const UInt& iAtom) const;

			/**
			 * judge whether the given shell is significant
			 * We determine it by the given basis set offset
			 * the basoffset should be global basis set index
			 */
			bool isSigShell(const UInt& basOffset) const;

			/**
			 * this function returns the information about significant 
			 * basis set list
			 */
			const UIntVec& getGlobalSigBasisIndex() const { return sigBas; };

			/**
			 * this function returns the information about significant 
			 * shell list
			 */
			const UIntVec& getSigShellIndex() const { return sigShell; };

			/**
			 * return the global basis set index (in normal type)
			 * for the given significant atom and local significant basis set index
			 * \iSigAtom    local significant atom index
			 * \iSigbasis   local significant basis index
			 */
			UInt getGlobalBasisIndex(const UInt& iSigAtom, const UInt& iSigBasis) const {
				UInt basOffset = getSigBasisBeginIndex(iSigAtom);
				UInt index = basOffset + iSigBasis;
				return sigBas[index];
			};

			/**
			 * return the significant basis set index by given the global basis set index
			 */
			UInt getSigBasisIndex(const UInt& iBasis) const; 

			/**
			 * return the total number of significant basis sets
			 */
			UInt getNSigBasis() const{
				return sigBas.size();
			};

			/**
			 * return the number of significant basis sets for a given atom
			 * the iAtom index should be the local index of sigAtoms
			 */
			UInt getNSigBasis(const UInt& iAtom) const{
				UInt basBeginIndex = atomBasisMapping[iAtom];
				UInt basEndIndex   = atomBasisMapping[iAtom+1];
				return basEndIndex - basBeginIndex;
			};

			/**
			 * return the begining index of basis set in sigBas array
			 * for this significant atom
			 * the iAtom index should be the local index of sigAtoms
			 */
			UInt getSigBasisBeginIndex(const UInt& iAtom) const{
				return atomBasisMapping[iAtom];
			};

			/**
			 * return the number of significant shells for a given atom
			 * the iAtom index should be the local index of sigAtoms
			 */
			UInt getNSigShell(const UInt& iAtom) const{
				UInt shlBeginIndex = atomShellMapping[iAtom];
				UInt shlEndIndex   = atomShellMapping[iAtom+1];
				return shlEndIndex - shlBeginIndex;
			};

			/**
			 * return the begining local index of shell in sigShell array
			 * for this significant atom
			 * the iAtom index should be the local index of sigAtoms
			 */
			UInt getSigShellBeginIndex(const UInt& iAtom) const{
				return atomShellMapping[iAtom];
			};

			/**
			 * return the number of significant atoms
			 */
			UInt getNSigAtoms() const{
				return sigAtoms.size();
			};

			/**
			 * return the significant atom list
			 */
			const UIntVec& getSigAtoms() const{
				return sigAtoms;
			};

			/**
			 * whether no significant basis set found in this batch
			 */
			bool empty() const{
				if (sigBas.size() == 0) return true;
				return false;
			};

			/**
			 * obtain the local basis set index for the given atom
			 * suppose that the global basis set corresponding to the 
			 * iSigAtom and iSigBasis set is X, and the basis set 
			 * starting index for the atom shell is Y, then the local
			 * index is just X - Y
			 */
			UInt getLocalAtomBasisIndex(const AtomShell& as, const UInt& iSigAtom, 
					const UInt& iSigBasis) const;
	};
}


#endif
