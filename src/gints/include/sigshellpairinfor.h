/**
 * \file    sigshellpairinfor.h
 * \brief   describing the significant shell pairs arragement information
 * \author  Fenglai Liu
 */
#ifndef SIGSHELLPAIRINFOR_H
#define SIGSHELLPAIRINFOR_H
#include "libgen.h"
#include "cptrans.h"
#include<vector>
#include "tbb/scalable_allocator.h"
#include "tbb/tbb.h"
#include "scalablevec.h"
#include "shell.h"
#include "shellprop.h"
#ifdef GINTS_PHI_DEBUG
#include "compactsigmolshellpairinfor.h" 
#endif

namespace shellsize {
	class AtomShellSize;
}

namespace gintsinfor {
	class GIntsInfor;
	class GIntJobInfor;
}

namespace denmtrxinfor {
	class DenMtrxInfor;
}

namespace sigshellpairinfor {

	using namespace shell;
	using namespace shellprop;
	using namespace shellsize;
	using namespace gintsinfor;
	using namespace denmtrxinfor;
	using namespace cptrans;
#ifdef GINTS_PHI_DEBUG
	using namespace compactsigmolshellpairinfor; 
#endif
	using namespace std;
	using namespace tbb;

	//
	// forward declare
	//
	class SigShellPairInfor;
	class SigAtomShellPairInfor;

	// redefine the vector for shell infor array and atom shell infor array
	typedef std::vector<SigShellPairInfor,tbb::scalable_allocator<SigShellPairInfor> > SigShellPairInforVec;
	typedef std::vector<SigAtomShellPairInfor,tbb::scalable_allocator<SigAtomShellPairInfor> > SigAtomShellPairInforVec;

	/**
	 * \class SigShellPairInfor
	 *
	 * SigShellPairInfor set up sigfinicant shell pairs information for two given shells
	 * This class is designed to wrap up the data for significant shell pair
	 *
	 */
	class SigShellPairInfor {

		private:

			UInt rowShellIndex;       ///< the row shell local index(within atom shell)
			UInt colShellIndex;       ///< the col shell local index(within atom shell)
			Double intBoundary;       ///< boundary value for assessing the significant integrals

		public:

			///
			/// constructor
			/// we note that for integral boundary value, it's set to one
			/// in default, which means this shell pair is important in calculation.
			///
			SigShellPairInfor(const Shell& is, const Shell& js):rowShellIndex(is.getLocalShellIndex()),
			colShellIndex(js.getLocalShellIndex()),intBoundary(ONE) { };

			///
			/// constructor
			/// 
			/// this constructor is used to retrieve data from the compact sig shell pair infor 
			///
			SigShellPairInfor(const UInt& rsIndex, const UInt& csIndex, 
					const Double& intBoundVal):rowShellIndex(rsIndex),
			colShellIndex(csIndex),intBoundary(intBoundVal) { };

			/**
			 * assign operator
			 */
			SigShellPairInfor& operator=(const SigShellPairInfor& sp) {
				if (this == &sp) {
					return *this;
				}else{
					rowShellIndex = sp.rowShellIndex;
					colShellIndex = sp.colShellIndex;
					intBoundary   = sp.intBoundary;
					return *this;
				}
			};

			///
			/// destructor
			///
			~SigShellPairInfor() { };

			///
			/// return the row shell index
			///
			UInt getRowShellIndex() const { return rowShellIndex; };

			///
			/// return the col shell index
			///
			UInt getColShellIndex() const { return colShellIndex; };

			///
			/// return the integral boundary value
			///
			Double getIntBoundaryVal() const { return intBoundary; };

			///
			/// this is used to set the integral boundary value
			///
			void setIntBoundaryVal(const Double& val) {
				intBoundary = val;
			};
	};

	/**
	 * \class SigAtomShellPairInfor
	 *
	 * SigAtomShellPairInfor set up sigfinicant shell pairs information for two
	 * given atom shells
	 *
	 * we note that inside the atom shell pair data, the shell pairs are arranged
	 * according to the angular momentum code. For example, firstly we will search
	 * the (S,S) shell pairs in the atom shell pair, then (SP,S), (P,S) etc.
	 *
	 * therefore the same angular code shell pair infor are aggregated together.
	 */
	class SigAtomShellPairInfor {

		private:

			UInt rowAtomShellIndex;                     ///< the row atom shell index
			UInt colAtomShellIndex;                     ///< the col atom shell index
			UInt nSigSP;                                ///< the number of significant shell pairs data in sigSPList
			Double intBoundary;                         ///< boundary value for assessing the significant integrals
			UInt lenC2;                                 ///< if we have composite shell pair, the total length of c2 pair
			UInt lenNP2;                                ///< the length of primitive pair array
			UInt nAngCodeInfor;                         ///< the actual length of nagCodeInfor array

			///
			/// this is the array defines the angular momentum type for this given atom shell pair
			/// data. The angular moentum code is defined in the SHELL_PAIR_ORDER_ARRAY, we will
			/// record all of angular momentum code in this atom shell pair.
			///
			UInt angCodeInfor[MAX_SHELL_PAIR_NUMBER];

			///
			/// for each shell pair angular moemtum code defined in the above array, 
			/// we record the corresponding shell pair offset in the sigSPList and 
			/// the number of shell pairs corresponding to this angular momentum
			///
			UInt angCodeSPOffset[2*MAX_SHELL_PAIR_NUMBER];

			///
			/// this vector stores the whole significant shell pair list
			///
			SigShellPairInforVec sigSPList;       

			/**
			 * form the primitive pairs data for running 
			 * integral calculation later
			 * \param is, js: input shell data, should have L(is) >= L(js)
			 * \param A,  B : centers of is and js
			 * \return iexp2: sum of exponential factors
			 * \return icoe2: multiplication of coef. factors
			 * \return fbra : prefactor for integrals
			 * \return P    : new center by combining two Gaussian primitives
			 */
			void formPrimPairData(const Shell& is, const Shell& js,
					const Double* A, const Double* B, DoubleVec& iexp2,
					DoubleVec& icoe2, DoubleVec& fbra, DoubleVec& P) const;

			/**
			 * normalize the overlap integrals with basis set normalization factor
			 * here the row/col dimension must be same with overlap integrals
			 * \param is, js: input shell data
			 * \return rowScaleNormVec: row shell normalization factors, generated inside
			 * \return colScaleNormVec: col shell normalization factors, generated inside
			 * \return ov             : normalized overlap integrals
			 */
			void normOV(const Shell& is, const Shell& js, 
					DoubleVec& rowScaleNormVec, DoubleVec& colScaleNormVec, 
					DoubleVec& ov) const;

		public:

			/**
			 * constructor for forming the significant shell pair
			 * \param  rowAtomShell: input row atom shell
			 * \param  colAtomShell: input col atom shell
			 * \param  rowAtomShellSize: input row atom shell size data
			 * \param  colAtomShellSize: input col atom shell size data
			 * \param  threshold   : threshold value to determinging the significance of shell pairs 
			 * \param  sameShell   : whether the two input shell data are from same MolShell?
			 */
			SigAtomShellPairInfor(const AtomShell& rowAtomShell, const AtomShell& colAtomShell, 
					const AtomShellSize& rowAtomShellSize, const AtomShellSize& colAtomShellSize,
					const Double& threshold, bool sameShell);

			/**
			 * default constructor
			 */
			SigAtomShellPairInfor():rowAtomShellIndex(-1),colAtomShellIndex(-1),nSigSP(0),
			intBoundary(ONE),lenC2(0),lenNP2(0),nAngCodeInfor(0) { };

			/**
			 * constructor used to build an empty SigAtomShellPairInfor and hold the memory
			 * of sigSPList
			 */
			SigAtomShellPairInfor(UInt nMaxSP):rowAtomShellIndex(-1),colAtomShellIndex(-1),nSigSP(0),
			intBoundary(ONE),lenC2(0),lenNP2(0),nAngCodeInfor(0) { 
				sigSPList.reserve(nMaxSP);
				SigShellPairInfor sp(-1,-1,ONE);
				for(UInt i=0; i<nMaxSP; i++) sigSPList.push_back(sp);
			};

			/**
			 * copy all of significant shell pair data and omitting the insig ones according to the 
			 * input threshold value
			 */
			void copySigData(const SigAtomShellPairInfor& infor, const Double& thresh);

			/**
			 * build the SigAtomShellPairInfor from the compact SigAtomShellPairInfor
			 *
			 * we note, that in this function we will form the data according to the compact ones
			 */
#ifdef GINTS_PHI_DEBUG
			void init(const CompactSigAtomShellPairInfor& aspInfor, const CompactSigShellPairInfor* spInfor);
#endif

			/**
			 * define the operator <
			 * the larger that the integral boundary is, the smaller it's
			 * through this way we will keep all of significant shell pair
			 * information on the head part, and the insig ones at tail
			 */
			bool operator<(const SigAtomShellPairInfor& infor) const {
				return (intBoundary > infor.intBoundary);
			};

			/**
			 * default destructor
			 */
			~SigAtomShellPairInfor() { };

			/**
			 * the number of significant shell pairs
			 */
			UInt getNSigShellPairs() const { return nSigSP; };

			/**
			 * get the length of coefficient pair array
			 * if no composite shell, it will be 0
			 */
			UInt getLenC2() const { return lenC2; };

			/**
			 * get the length of primitive pair array
			 */
			UInt getLenNP2() const { return lenNP2; };

			/**
			 * get the number of angular momentum code in the array of angCodeInfor
			 */
			UInt getNAngCode() const { return nAngCodeInfor; }; 

			/**
			 * get the ang code for the index
			 */
			UInt getAngCode(const UInt& index) const { return angCodeInfor[index]; }; 

			/**
			 * get the shell pair offset data sorted by the angular momentum code
			 */
			void getSPOffsetData(const UInt& index, UInt& offset, UInt& nSP) const { 
				offset = angCodeSPOffset[2*index  ]; 
				nSP    = angCodeSPOffset[2*index+1]; 
			}; 

			/**
			 * whether the atom shell pair is insignificant?
			 */
			bool isInsig() const { return (nSigSP == 0); };

			/**
			 * get the sigfinicant shell pair infor
			 */
			const SigShellPairInfor& getSigShellPairInfor(const UInt& iSP) const {
				return sigSPList[iSP];
			};

			/**
			 * get row atom shell index
			 */
			UInt getRowAtomShellIndex() const { return rowAtomShellIndex; };

			/**
			 * get col atom shell index
			 */
			UInt getColAtomShellIndex() const { return colAtomShellIndex; };

			///
			/// return the integral boundary value
			///
			Double getIntBoundaryVal() const { return intBoundary; };

			///
			/// this is used to set the integral boundary value for shell pair data
			///
			void setIntBoundaryVal(const UInt& iSP, const Double& val) { 
				sigSPList[iSP].setIntBoundaryVal(val);
			};

			///
			/// this is used to set the integral boundary value for whole atom shell pair
			///
			void setIntBoundaryVal(const Double& val) { 
				intBoundary = val;
			};

			///
			/// judge the negligibility of atom shell quartet(<12|34>)
			/// by using the integral boundary value
			/// 
			///
			bool isNegligible(const SigAtomShellPairInfor& spInfor, const Double& thresh) const {
				Double boundary = intBoundary*spInfor.getIntBoundaryVal();
				if (boundary<thresh) return true;
				return false;
			};

			///
			/// judge the negligibility of atom shell quartet(<12|34>)
			/// by using the integral boundary value together with pMax
			/// 
			///
			bool isNegligible(const SigAtomShellPairInfor& spInfor, const Double& pMax,
					const Double& thresh) const {
				Double boundary = intBoundary*spInfor.getIntBoundaryVal()*pMax;
				if (boundary<thresh) return true;
				return false;
			};

			///
			/// judge the negligibility of shell quartet(<12|34>)
			/// by using the integral boundary value
			/// \param iSP     the shell pair index for this infor
			/// \param jSP     the shell pair index for infor passing in
			/// \param spInfor the other part of shell pair information for forming shell quartet
			/// \param thresh  threshold value
			///
			bool isNegligible(const UInt& iSP, const UInt& jSP, 
					const SigAtomShellPairInfor& spInfor, const Double& thresh) const {
				Double boundary = sigSPList[iSP].getIntBoundaryVal();
				const SigShellPairInfor& sp = spInfor.getSigShellPairInfor(jSP);
				boundary *= sp.getIntBoundaryVal();
				if (boundary<thresh) return true;
				return false;
			};

			///
			/// judge the negligibility of shell quartet(<12|34>)
			/// by using the integral boundary value digest with the 
			/// PMax value
			/// \param iSP     the shell pair index for this infor
			/// \param jSP     the shell pair index for infor passing in
			/// \param spInfor the other part of shell pair information for forming shell quartet
			/// \param pMax    the absoluate maximum density matrix value
			/// \param thresh  threshold value
			///
			bool isNegligible(const UInt& iSP, const UInt& jSP, 
					const SigAtomShellPairInfor& spInfor, const Double& pMax,
					const Double& thresh) const {
				Double boundary = sigSPList[iSP].getIntBoundaryVal();
				const SigShellPairInfor& sp = spInfor.getSigShellPairInfor(jSP);
				boundary *= sp.getIntBoundaryVal()*pMax;
				if (boundary<thresh) return true;
				return false;
			};

			/**
			 * bool operator == 
			 * used for debugging purpose
			 */
			bool operator==(const SigAtomShellPairInfor& asp) const; 

			/**
			 * debug printing 
			 */
			void print() const;
	};

	/**
	 * \class SigMolShellPairInfor
	 *
	 * SigMolShellPairInfor set up sigfinicant shell pairs information for two
	 * given molecular shells
	 *
	 */
	class SigMolShellPairInfor {

		private:

			//
			// store the code of the shells as ID
			// so that to test whether shell pair is identical with
			// another 
			//
			UInt rowShellCode;     ///< row shell code
			UInt colShellCode;     ///< row shell code

			//
			// firstly we need to set up some information
			// so that we can build raw shell pair 
			//
			UInt maxL;             ///< maxL for the shells in the shell pair
			UInt maxNShellPairs;   ///< maximum number of shell pairs per atom pair
			UInt maxNP2;           ///< maximum number of primitive pairs for any shell pair
			UInt maxLContraction;  ///< maximum L contraction degree
			UInt maxASPNP2;        ///< maximum length of primitie pairs for atom shell pair data
			UInt maxASPLenC2;      ///< maximum length of c2 array for an atom shell pair

			//
			// secondly, we may also need some information which is needed in
			// digesting the raw integrals etc.
			// we note, that the following maximum dimension data is set up 
			// based on the Cartesian shell information
			//
			UInt maxNBasShell;     ///< maximum number of basis set for shell
			UInt maxNBasAtom;      ///< maximum number of basis set for atom shell

			//
			// this is the batch information 
			// if we set up the integral boundary value and do sorting,
			// we will set up this information 
			//
			// we note that the batches is set up in this way:
			// [0.1,     +infinity)
			// [0.01,    0.1)
			// [0.001,   0.01)
			// [0.0001,  0.001)
			// [0.00001, 0.0001)
			// etc.
			//
			UIntVec batchSigAtomShellPairInfor;      ///< divide the sigAtomSPList into batches
			DoubleVec  batchThreshLimits;            ///< for each batch, it's lower limit value stored

			//
			// now it's the data
			//
			SigAtomShellPairInforVec sigAtomSPList;  ///< significant atom shell pair list 

			///
			/// perform cauchy-schwarz(CS) inequality to get integral boundary for the given
			/// list. We note that the input list is not the sigAtomSPList, therefore
			/// this function is constant
			///
			void setIntBoundary(const GIntJobInfor& infor, const MolShell& bra1, 
					const MolShell& bra2, const DenMtrxInfor& denInfor,
					SigAtomShellPairInforVec& spList) const;

			///
			/// form the significant atom shell pair list to the list
			/// we note that the list just contains the sig atom shell pairs
			/// which characterized by overlap integral
			///
			/// for some jobs, like ERI calculation we may need further 
			/// squeezing the atom shell pairs data by CS integral boundary check
			///
			void formSigAtomSPList(const MolShell& rs, const MolShell& cs, const GIntJobInfor& infor,
					SigAtomShellPairInforVec& list);

			///
			/// this function is used to set up the batch sig atom shell pair information
			/// array according to the integral boundary value
			///
			void setupBatchSigAtomShellPairInfor(const GIntJobInfor& ginfor);

		public:

			/**
			 * constructor for forming the significant shell pair information
			 * this one goes without density matrix information
			 * \param  rs, cs : input row/column molecular shell 
			 * \param  ginfor : gints job information
			 */
			SigMolShellPairInfor(const MolShell& rs, const MolShell& cs, const GIntJobInfor& ginfor);

			/**
			 * constructor for forming the significant shell pair information
			 * this one goes with density matrix information
			 * \param  rs, cs : input row/column molecular shell 
			 * \param  ginfor : gints job information
			 */
			SigMolShellPairInfor(const MolShell& rs, const MolShell& cs, const DenMtrxInfor& denMtrxInfor,
					const GIntJobInfor& ginfor);

			/**
			 * empty constructor - just build an empty object
			 *
			 * later we will use init function to construct the data if we need it
			 */
			SigMolShellPairInfor():rowShellCode(-1),colShellCode(-1),maxL(-1),maxNShellPairs(0),  
			maxNP2(0),maxLContraction(-1),maxASPNP2(0),maxASPLenC2(0),maxNBasShell(0),maxNBasAtom(0) { };    

			/**
			 * forming the significant shell pair information, this is the same function of the constructor
			 * this one goes without density matrix information
			 * \param  rs, cs : input row/column molecular shell 
			 * \param  ginfor : gints job information
			 */
			void init(const MolShell& rs, const MolShell& cs, const GIntJobInfor& ginfor);

			/**
			 * forming the significant shell pair information
			 * this one goes with density matrix information
			 * \param  rs, cs : input row/column molecular shell 
			 * \param  ginfor : gints job information
			 */
			void init(const MolShell& rs, const MolShell& cs, const DenMtrxInfor& denMtrxInfor,
					const GIntJobInfor& ginfor);

			/**
			 * default destructor
			 */
			~SigMolShellPairInfor() { };

			/**
			 * the number of significant shell pairs
			 */
			UInt getNSigAtomShellPairs() const { return sigAtomSPList.size(); };

			/**
			 * get the vector length for the vector inputed for function
			 * convertMatrixIntoSigSPVec
			 * \param rs     row shell corresponding to the shell pair information
			 * \param cs     column shell corresponding to the shell pair information
			 * \param type   the matrix is in Cartesian, normal or Pure basis function type?
			 */
			UInt getSigSPVecLen(const MolShell& rs, const MolShell& cs, const UInt& type) const;

			/**
			 * get the s2 form vector offset for the given row atom shell index and col atom shell index
			 * \param type   the matrix is in Cartesian, normal or Pure basis function type?
			 */
			UInt getS2VecOffset(const MolShell& rs, const MolShell& cs, 
					const UInt& rowShellIndex, const UInt& colShellIndex, const UInt& type) const;

			/**
			 * the function is the driver function to generate the input s2 form vector(v) from matrix M, 
			 * or to transform the data from v to input matrix M. Here inside we can also do the C2P/P2C
			 * and basis set scaling etc. work
			 *
			 * the working functor is TBB_MtrxToS2Vec. 
			 *
			 * \param toVec  whether to generate vector(true), or to move the data from vector to matrix(false)?
			 * \param rs     row shell corresponding to the shell pair information
			 * \param cs     column shell corresponding to the shell pair information
			 * \param M      data matrix
			 * \param v      the s2 form data vector
			 */
			void convertSigSPVec(bool toVec, const GlobalInfor& globInfor, 
					const MolShell& rs, const MolShell& cs, const UInt& transWork, const UInt& scaleWork, 
					const UInt& matrixStatus, Mtrx& M, DoubleVec& v) const;

			/**
			 * get the sigfinicant shell pair infor
			 */
			const SigAtomShellPairInfor& getSigAtomShellPairInfor(const UInt& iSP) const {
				return sigAtomSPList[iSP];
			};

			/**
			 * sometimes we also need the non-const version of fetch
			 */
			SigAtomShellPairInfor& getSigAtomShellPairInfor(const UInt& iSP) {
				return sigAtomSPList[iSP];
			};

			///
			/// squeeze the raw shell pair list and form the final one
			/// this is used, for example; in the cauchy-schwarz(CS) 
			/// boundary check
			///
			void squeezeSPList(const SigAtomShellPairInforVec& list, const GIntJobInfor& ginfor);

			///
			/// get row shell code
			///
			UInt getRowShellCode() const { return rowShellCode; };

			///
			/// get col shell code
			///
			UInt getColShellCode() const { return colShellCode; };

			///
			/// return the maxNP2
			///
			UInt getMaxNP2() const { return maxNP2; };

			///
			/// return the maxNP2 for an atom shell pair
			///
			UInt getMaxAtomSPNP2() const { return maxASPNP2; };

			///
			/// return the maxLenC2 for an atom shell pair
			///
			UInt getMaxAtomSPLenC2() const { return maxASPLenC2; };

			///
			/// return the maxLContraction
			///
			UInt getMaxNL() const { return maxLContraction; };

			///
			/// return the max number of shell pairs
			///
			UInt getMaxNSP() const { return maxNShellPairs; };

			///
			/// return the max L
			///
			UInt getMaxL() const { return maxL; };

			///
			/// return the maximum basis set number for shell
			///
			UInt getMaxNBasisForShell() const { return maxNBasShell; };

			///
			/// return the maximum basis set number for atom shell
			///
			UInt getMaxNBasisForAtomShell() const { return maxNBasAtom; };

			///
			/// whether two shell pairs are same with each other?
			///
			bool operator==(const SigMolShellPairInfor& infor) const{
				return ((rowShellCode == infor.rowShellCode) && 
						(colShellCode == infor.colShellCode));
			}; 

			///
			/// get the number of batches 
			///
			UInt getNBatch() const { return batchSigAtomShellPairInfor.size()/2; };

			///
			/// for the given batch, how many significant atom shell pairs we have
			///
			/// we note, that the input batch index iBatch, should be in the range 
			/// of [0, nBatch)
			///
			UInt getNSigAtomSPInBatch(UInt iBatch) const {
				return (batchSigAtomShellPairInfor[2*iBatch+1] - batchSigAtomShellPairInfor[2*iBatch]);
			};

			///
			/// get the starting sig atom shell pair index for the given batch
			///
			UInt getStartingSigAtomSPIndex(UInt iBatch) const {
				return batchSigAtomShellPairInfor[2*iBatch];
			};

			///
			/// get the end sig atom shell pair index for the given batch
			///
			UInt getEndingSigAtomSPIndex(UInt iBatch) const {
				return batchSigAtomShellPairInfor[2*iBatch+1];
			};

			///
			/// get the lower limit integral boundary value for the current batch
			///
			Double getBatchLowerLimitBoundary(UInt iBatch) const {
				return batchThreshLimits[iBatch];
			};

			/**
			 * debug printing 
			 */
			void print(UInt level = 0) const;

			/**
			 * print out statistic information based on the integral boundary
			 */
			void statPrint(const Double& thresh, const MolShell& rs, const MolShell& cs) const;

			///
			/// this is a debugging function for the compact sig molecule shell pair infor
			///
#ifdef GINTS_PHI_DEBUG
			bool testCompactSigMolSPInfor(const CompactSigMolShellPairInfor& spInfor) const; 
#endif

	};

	///
	/// this is the functor to perform cauchy-schwarz(CS) inequality
	/// evaluation for ERI
	///
	class TBB_CSIntBoundary {

		private:

			//
			// information to do the CS integral boundary
			//
			Double thresh;                      ///< threshold value for generating shell pair
			Double intThresh;                   ///< CS boundary threshold value
			const GIntJobInfor& jobInfor;       ///< job information

			//
			// input shell data and shell pair information
			// we are going to compute the integral boundary value over the sp list
			//
			const DenMtrxInfor& denInfor;       ///< general density matrix information used for boundary set up 
			const MolShell&   bra1;             ///< shell data for bra1
			const MolShell&   bra2;             ///< shell data for bra2
			const SigMolShellPairInfor& infor;  ///< the significant shell pair information
			SigAtomShellPairInforVec& spList;   ///< significant atom shell pair list(output) 

			//
			// these data is used for perform integral scaling
			//
			CPTransData scale;                  ///< this is used to perform basis set scaling
			Double scale_S[1];                  ///< constant for S shell
			Double scale_P[3];                  ///< constant for P shell
			Double scale_SP[4];                 ///< constant for SP shell

			//
			// a small function to get convert vector
			//
			const Double* getScaleVec(const UInt& lmin, const UInt& lmax) const {
				if (lmax == 0) return scale_S;
				if (lmax == 1 && lmin == 1) return scale_P;
				if (lmax == 1 && lmin == 0) return scale_SP;
				const DoubleVec& vec = scale.getConvertVec(lmax);
				return &vec[0];
			};

		public:

			///
			/// constructor for building cs integral boundary
			///
			TBB_CSIntBoundary(const GIntJobInfor& ginfor, const DenMtrxInfor& denInfor0,
					const MolShell& b1, const MolShell& b2,
					const SigMolShellPairInfor& infor0, SigAtomShellPairInforVec& list);

			///
			/// destructor
			///
			~TBB_CSIntBoundary() { };

			///
			/// operator to do the parallel work
			///
			void operator()(const blocked_range<UInt>& r) const;
	};

	///
	/// this is the functor to perform matrix to "s2 form" vector transformation or verse.
	///
	/// "s2 form" vector is composed by atom shell block data, the order to organize it is 
	/// same to the shell pair order in SigMolShellPairInfor. Therefore by using the shell
	/// data and the shell pair information we can convert the matrix into a "S2 form"
	/// vector.
	///
	/// we note that the s2 form vector only contains the data for significant atom shell pairs.
	/// this is the feature that makes s2 form vector much more compact for the data. 
	///
	/// pay attention to the input basType. If no C2P/P2C work is done, then we assume that the 
	/// input vector and the matrix are belonging to the same basis set type defined here.
	///
	class TBB_MtrxToS2Vec {

		private:

			// setting of the work
			// please go to cptrans.h to see the definitions of the macro
			// we just follow the setting over there
			bool doMtrxToVec;                   ///< whether the work is to matrix to vector?
			UInt basType;                       ///< if no C2P or P2C work, we need to know the matrix in Cartesian or pure?
			UInt transStatus;                   ///< whether do C2P or P2C or none?
			UInt scaleStatus;                   ///< whether you scale the Cartesian D, F etc. basis set data?
			UInt matrixStatus;                  ///< with matrix itself or it's transpose?

			// reference data of input
			const MolShell&   rs;               ///< shell data corresponding to the row
			const MolShell&   cs;               ///< shell data corresponding to the column
			const SigMolShellPairInfor& infor;  ///< the significant shell pair information

			// now it's output data
			Mtrx& M;                            ///< the matrix data, could be input or output
			DoubleVec& v;                       ///< the "S2 form" vector, could be input or output

		public:

			///
			/// constructor 
			///
			TBB_MtrxToS2Vec(bool toVec, UInt transWork, UInt scaleWork, UInt mtrxStatus, UInt basType0,
					const MolShell& rs0, const MolShell& cs0, const SigMolShellPairInfor& infor0, 
					Mtrx& M0, DoubleVec& v0):doMtrxToVec(toVec),basType(basType0),transStatus(transWork),
			scaleStatus(scaleWork),matrixStatus(mtrxStatus),rs(rs0),cs(cs0),infor(infor0),M(M0),v(v0) { };

			///
			/// destructor
			///
			~TBB_MtrxToS2Vec() { };

			///
			/// operator to do the parallel work
			///
			void operator()(const blocked_range<UInt>& r) const;
	};
}

#endif

