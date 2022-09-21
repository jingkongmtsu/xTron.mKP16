/**
 * \file    shellpair.h
 * \brief   describing the shell pairs data for two given atom shells
 * \author  Fenglai Liu
 */
#ifndef SHELLPAIR_H
#define SHELLPAIR_H
#include "libgen.h"
#include<vector>
#include "tbb/scalable_allocator.h"
#include "scalablevec.h"
#ifdef GINTS_PHI_DEBUG
//#include "compactshellpairdata_phi.h"
#endif

namespace shell {
	class shell;
	class AtomShell;
}

namespace sigshellpairinfor {
	class SigAtomShellPairInfor;
}

namespace shellpair {

#ifdef GINTS_PHI_DEBUG
	// using namespace compactshellpairdata_phi; 
#endif
	using namespace sigshellpairinfor; 
	using namespace shell;
	using namespace std;

	//
	// forward declare
	//
	class ShellPair;
	class AtomShellPair;

	// redefine the vector for shell pair array 
	typedef std::vector<ShellPair,tbb::scalable_allocator<ShellPair> > ShellPairVec;

	/**
	 * \class ShellPair
	 *
	 * This class is basically set up shell pair data for integral calculation
	 * the result data listed below, must be in accordance with the integral
	 * file requirement (see the function of switchShellPair defined in 
	 * shellprop.h and it's comments above).
	 *
	 * Namely, the shell data below must have shell i >= shell j.
	 * which is related to rule 1 and rule 2 in shellprop.h
	 * in our code, the shell i could be also referred as "row shell",
	 * and the shell j is also referred as "col shell". This is because
	 * in integral calculation, shell i usually forms the row side,
	 * and shell j forms the column data side.
	 *
	 * the shell pair data is used in scratch way. That is to say, we set up
	 * it first(the most important thing, is to set up enough memory for future
	 * calculation), then reuse it each time we form the shell pair data.
	 *
	 */
	class ShellPair {

		private:

			//
			// angular momentum information 
			//
			bool shellInversed;         ///< whether we inverse the row and col shell
			bool pureShellI;            ///< shell i is pure 
			bool pureShellJ;            ///< shell j is pure 
			UInt iLmin;                 ///< Lmin for shell i
			UInt iLmax;                 ///< Lmax for shell i
			UInt jLmin;                 ///< Lmin for shell j
			UInt jLmax;                 ///< Lmax for shell j
			UInt LCode;                 ///< angular momentum code for this shell pair

			//
			// basis set information
			// We store the local basis set order, which is basis set index inside
			// the atom shell
			// also we store the global basis set order, which is indexed across 
			// the whole molecule shell
			//
			UInt iBasOffset;            ///< shell i basis offset - in normal basis set order
			UInt jBasOffset;            ///< shell j basis offset - in normal basis set order
			UInt iGBasOffset;           ///< shell i global basis offset - normal basis set order
			UInt jGBasOffset;           ///< shell j global basis offset - normal basis set order
			UInt iNBas;                 ///< number of normal basis set for shell i
			UInt jNBas;                 ///< number of normal basis set for shell j
			UInt iCarBasOffset;         ///< shell i basis offset - in Cartesian basis set order
			UInt jCarBasOffset;         ///< shell j basis offset - in Cartesian basis set order
			UInt iGCarBasOffset;        ///< shell i global basis offset - in Cartesian basis set order
			UInt jGCarBasOffset;        ///< shell j global basis offset - in Cartesian basis set order
			UInt iNCarBas;              ///< number of cartesian basis set for shell i
			UInt jNCarBas;              ///< number of cartesian basis set for shell j

			//
			// shell index information
			// including local index (index counted inside the atom shell)
			// and global index (index counted inside the molecule shell)
			//
			UInt iShellIndex;           ///< shell i local index
			UInt jShellIndex;           ///< shell j local index
			UInt iGShellIndex;          ///< shell i global index
			UInt jGShellIndex;          ///< shell j global index

			//
			// The primitive pairs data in compact form
			// most of time we do not need e2diff, we only evaluate it for kinetic
			// integral etc.
			// however, since to compute it does not take anything; therefore we compute
			// it every time
			//
			UInt np2;                   ///< the practical length of primitive pairs
			DoubleVec c2;               ///< coefficients product: ic*jc
			DoubleVec e2;               ///< exponent: 1/(alpha+beta) 
			DoubleVec fac;              ///< pre-factors for basic integrals - ovrelap integrals
			DoubleVec P;                ///< new centers created by two Gaussian primitives
			DoubleVec e2diff;           ///< difference between exponents (shell i - shell j)
			Double A[3];                ///< center for shell i
			Double B[3];                ///< center for shell j

			//
			// this is a scratch array used to record the significance 
			// of primitive pairs, should not be used in integral calculation
			//
			UIntVec sigPrimPairInfor;   ///< if it's significant, then it's 1, else it's 0

			//
			// check the sparsity of primitive pairs
			// also push in primitive data
			//
			void formPrimPairsData(const Shell& is, const Shell& js, const Double& thresh); 

		public:

			/**
			 * constructor
			 * set up the scratch data
			 *
			 * for set up the array, we need maxNP2; which is the maximum primitive 
			 * pairs number. and maxNL, which is the maximum length of subshells 
			 * in a given composite shell (for example, SP is maxNL=2, SPD is 3)
			 */
			ShellPair(UInt maxNP2, UInt maxNL):shellInversed(false),pureShellI(false),
			pureShellJ(false),iLmin(-1),iLmax(-1),jLmin(-1),jLmax(-1),               
			LCode(0),iBasOffset(-1),jBasOffset(-1),iGBasOffset(-1),jGBasOffset(-1),
			iNBas(0),jNBas(0),iCarBasOffset(-1),jCarBasOffset(-1),iGCarBasOffset(-1),
			jGCarBasOffset(-1),iNCarBas(0),jNCarBas(0),
			iShellIndex(-1),jShellIndex(-1),iGShellIndex(-1),jGShellIndex(-1),
			np2(0),c2(maxNL*maxNL*maxNP2),e2(maxNP2),fac(maxNP2),P(3*maxNP2),        
			e2diff(maxNP2),sigPrimPairInfor(maxNP2) {  };   

			/**
			 * default destructor
			 */
			~ShellPair() { };

			/**
			 * initilize the shell pair data 
			 * \param  is,js: two shells for forming shell pair data
			 * \param  xyz:   the coordinates for two shell's center
			 * \param  threshold: threshold value to create shell pairs 
			 */
			void init(const Shell& is, const Shell& js, 
					const Double* ixyz, const Double* jxyz, const Double& threshold);

			/**
			 * debug printing 
			 */
			void print(UInt iprint = 0) const;

			/**
			 * get contraction degree
			 * which is always the length of e2
			 */
			UInt getNP2() const {
				return np2;
			};

			/**
			 * get angular momentum infor
			 */
			UInt getLCode() const {
				return LCode;
			};

			/**
			 * get the L information
			 * pay attention that L may be reversed compared with input
			 * shell information
			 */
			void getL(UInt& ilmin, UInt& ilmax, UInt& jlmin, UInt& jlmax) const {
				ilmin = iLmin;
				ilmax = iLmax;
				jlmin = jLmin;
				jlmax = jLmax;
			};

			/**
			 * get the max L between two lmax value 
			 */
			UInt getLmax() const {
				return (iLmax>jLmax ? iLmax : jLmax);
			};

			/**
			 * get local basis set offset - normal order
			 */
			void getLocalBasOffSet(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iBasOffset;
				colIndex = jBasOffset;
			};

			/**
			 * get global basis set offset - normal order
			 */
			void getGlobalBasOffSet(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iGBasOffset;
				colIndex = jGBasOffset;
			};

			/**
			 * get local basis set offset - Cartesian order
			 */
			void getLocalCarBasOffSet(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iCarBasOffset;
				colIndex = jCarBasOffset;
			};

			/**
			 * get global basis set offset - Cartesian order
			 */
			void getGlobalCarBasOffSet(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iGCarBasOffset;
				colIndex = jGCarBasOffset;
			};

			/**
			 * get local shell index
			 */
			void getLocalShellIndex(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iShellIndex;
				colIndex = jShellIndex;
			};

			/**
			 * get global basis set offset - Cartesian order
			 */
			void getGlobalShellIndex(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iGShellIndex;
				colIndex = jGShellIndex;
			};

			/**
			 * get basis set numbers in Cartesian format
			 */
			void getNCarBas(UInt& rowNBas, UInt& colNBas) const {
				rowNBas = iNCarBas;
				colNBas = jNCarBas;
			};

			/**
			 * get basis set numbers in normal format
			 */
			void getNBas(UInt& rowNBas, UInt& colNBas) const {
				rowNBas = iNBas;
				colNBas = jNBas;
			};

			/**
			 * get shell i's center
			 */
			const Double* getA() const {
				return A;
			};

			/**
			 * get shell j's center
			 */
			const Double* getB() const {
				return B;
			};

			/**
			 * get new center
			 */
			const Double* getP() const {
				return &P[0];
			};

			/**
			 * get coefficient product
			 */
			const Double* getC2() const {
				return &c2[0];
			};

			/**
			 * get exponent information
			 */
			const Double* getE2() const {
				return &e2[0];
			};

			/**
			 * get factor array
			 */
			const Double* getFac() const {
				return &fac[0];
			};

			/**
			 * get the difference of exp
			 */
			const Double* getExpDiff() const {
				return &e2diff[0];
			};

			/**
			 *  whether the shells are inversed
			 */
			bool inverseShells() const { return shellInversed; };

			/**
			 * whether the shell i is pure?
			 */
			bool isRowShellPure() const { return pureShellI; };

			/**
			 * whether the shell j is pure?
			 */
			bool isColShellPure() const { return pureShellJ; };

			/**
			 * whether the shell pair is significant?
			 */
			bool isSig() const { return (np2 > 0); };
	};

	/**
	 * \class AtomShellPair
	 *
	 * AtomShellPair stores the shell pair data in a vector for two given atom shells
	 *
	 */
	class AtomShellPair {

		private:

			//
			// basis set/shell information
			// We store the local basis set order, which is basis set index inside
			// the atom shell
			// also we store the global basis set order, which is indexed across 
			// the whole molecule shell
			//
			UInt rowAtomShellIndex;     ///< the row atom shell index
			UInt colAtomShellIndex;     ///< the col atom shell index
			UInt iGBasOffset;           ///< atom shell i global basis offset - normal basis set order
			UInt jGBasOffset;           ///< atom shell j global basis offset - normal basis set order
			UInt iNBas;                 ///< number of normal basis set for atom shell i
			UInt jNBas;                 ///< number of normal basis set for atom shell j
			UInt iGCarBasOffset;        ///< atom shell i global basis offset - in Cartesian basis set order
			UInt jGCarBasOffset;        ///< atom shell j global basis offset - in Cartesian basis set order
			UInt iNCarBas;              ///< number of cartesian basis set for atom shell i
			UInt jNCarBas;              ///< number of cartesian basis set for atom shell j

			//
			// shell pair data
			//
			UInt nShellPairs;           ///< number of shell pairs
			ShellPairVec spList;        ///< significant shell pair list 

		public:

			/**
			 * constructor for set up enough memory for building shell pairs
			 */
			AtomShellPair(UInt maxNShellPairs, UInt maxNP2, UInt maxNL):rowAtomShellIndex(-1),
			colAtomShellIndex(-1),iGBasOffset(-1),jGBasOffset(-1),iNBas(-1),jNBas(-1),
			iGCarBasOffset(-1),jGCarBasOffset(-1),iNCarBas(-1),jNCarBas(-1),nShellPairs(0){
				spList.reserve(maxNShellPairs);
				ShellPair sp(maxNP2,maxNL);
				for(UInt iSP=0; iSP<maxNShellPairs; iSP++) spList.push_back(sp);
			};

			/**
			 * initilize the shell pair data for the two given atom shell
			 * we note that the input shell order is fixed
			 * rs is for row atom shell, and cs is for column atom shell
			 * the shell pair is like this :(rs,cs|
			 * however, inside the spList the shell pair data may be inversed
			 * due to the requirement of integral file
			 * \param rs, cs   : input row and col shell data for the whole molecule
			 * \param infor    : significant atom shell pair information
			 * \param thresh   : threshold value to set up shell pair data
			 */
			void init(const MolShell& rs, const MolShell& cs, 
					const SigAtomShellPairInfor& infor, const Double& thresh);

			/**
			 * default destructor
			 */
			~AtomShellPair() { };

			/**
			 * the number of significant shell pairs
			 */
			UInt getNShellPairs() const { return nShellPairs; };

			/**
			 * get the sigfinicant shell pair infor
			 */
			const ShellPair& getShellPair(const UInt& iSP) const {
				return spList[iSP];
			};

			/**
			 * debug printing 
			 */
			void print(UInt iprint) const;

			/**
			 * get global basis set offset - normal order
			 */
			void getGlobalBasOffSet(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iGBasOffset;
				colIndex = jGBasOffset;
			};

			/**
			 * get global basis set offset - Cartesian order
			 */
			void getGlobalCarBasOffSet(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = iGCarBasOffset;
				colIndex = jGCarBasOffset;
			};

			/**
			 * get basis set numbers in Cartesian format
			 */
			void getNCarBas(UInt& rowNBas, UInt& colNBas) const {
				rowNBas = iNCarBas;
				colNBas = jNCarBas;
			};

			/**
			 * get basis set numbers in normal format
			 */
			void getNBas(UInt& rowNBas, UInt& colNBas) const {
				rowNBas = iNBas;
				colNBas = jNBas;
			};

			/**
			 * get the row shell and col shell index 
			 */
			void getAtomShellIndex(UInt& rowIndex, UInt& colIndex) const {
				rowIndex = rowAtomShellIndex;
				colIndex = colAtomShellIndex;
			};

			/**
			 * test that whether the comASP is same with the shell pair data defined here
			 * this is only a debug function
			 */
#ifdef GINTS_PHI_DEBUG
		//	bool same(const UInt& maxNP2, const UInt& maxNL, const Double& thresh, 
		//			const MolShell& rs, const MolShell& cs, const CompactShellPairData_Phi& comASP) const;
#endif

			/**
			 * compare that whether the two atom shell pair are same?
			 * we just use index
			 *
			 * Be careful to use this function. we do not check that
			 * whether the two atom shell pairs comparing inside
			 * are all from same shell data source, you need to check 
			 * it before use.
			 */
			bool operator==(const AtomShellPair& asp) const {
				if (rowAtomShellIndex == asp.rowAtomShellIndex && 
						colAtomShellIndex == asp.colAtomShellIndex) {
					return true;
				}; 
				return false;
			};

	};

}


#endif

