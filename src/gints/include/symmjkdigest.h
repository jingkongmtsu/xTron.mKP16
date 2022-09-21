/**
 * \file   symmjkdigest.h
 * \brief  doing the digestion for symmetrical J/K(coulomb/exchange) integrals
 * \author Fenglai Liu
 */
#ifndef SYMMJKDIGEST_H
#define SYMMJKDIGEST_H
#include "libgen.h"
#include "matrix.h"
#include "blockmatrix.h"
using namespace matrix;
using namespace blockmatrix;

namespace shellpair {
	class AtomShellPair;
	class ShellPair;
}

namespace denmtrxinfor {
	class DenMtrxInfor;
}

namespace shell {
	class MolShell;
}

namespace integralinfor {
	class SingleIntegralInfor;
}

namespace gintsinfor {
	class GIntJobInfor;
}

namespace symmjkdigest {

	using namespace shell;
	using namespace shellpair;
	using namespace gintsinfor;
	using namespace denmtrxinfor;
	using namespace integralinfor;

	/////////////////////////////////////////////////////////////////////////////////
	//   digestion for exchange type of integrals                                  //
	//   P_{nu,eta}   (mu,nu|lambda,eta) = K(mu,lambda) P24*(12|34)=K13            //
	//   P_{nu,lambda}(mu,nu|lambda,eta) = K(mu,eta)    P23*(12|34)=K14            //
	//   P_{mu,eta}   (mu,nu|lambda,eta) = K(nu,lambda) P14*(12|34)=K23            //
	//   P_{mu,lambda}(mu,nu|lambda,eta) = K(nu,eta)    P13*(12|34)=K14            //
	//                                                                             //
	//   Here below, we divide all of digestions into the following groups,        //
	//   and each group may have one or multiple cases inside.                     //
	//                                                                             //
	//   this is for symmetrical K digestion, which  means; mu, nu, lambda,        //
	//   eta are from same molecule shell data. the group set up and following     //
	//   digestion situation are based on the integral symmetry, that we have:     //
	//   (mu,nu|lambda,eta) = (nu,mu|lambda,eta)                                   //
	//   (mu,nu|lambda,eta) = (mu,nu|eta,lambda)                                   //
	//   (mu,nu|lambda,eta) = (lambda,eta|mu,nu)                                   //
	//   etc.                                                                      //
	//                                                                             //
	//   for the digestions, there are three situations:                           //
	//                                                                             //
	//   1 DUPLICATE_K_DIGESTION:                                                  //
	//   this means that this result block is duplicate. For example, in the       //
	//   case1111, K13,K14,K23,K24 are all same. We only need to use one of        //
	//   them to update the result matrix.  Therefore, K14,K23,K24 are all         //
	//   "duplicate".                                                              //
	//                                                                             //
	//   2  NORMAL_K_DIGESTION:                                                    //
	//   this means that this result block is unique in the digestion, and         //
	//   normal updating will be performed for this result block.                  //
	//                                                                             // 
	//   3  ADDITIONAL_DIGESTION:                                                  //
	//   for case case1213 K11 needs additional digestion. This is because that    //
	//   integral (12|13) and (13|12) is same, therefore if (12|13) is in          //
	//   calculation then (13|12) is not calculated anymore. However, in terms of  //
	//   integral (12|13) for K(1,1) the shell 2 and shell 3 will go over all of   //
	//   basis sets therefore both (12|13) and (13|12) needed to be digested.      //
	//   Thus, for K11 we need to digest with both P23 (for integral (12|13))      //
	//   and P32 (for integral (13|12)) so that to include all of results          //
	//                                                                             //
	/////////////////////////////////////////////////////////////////////////////////

	///
	///group 0
	///the case that should not exist for normal K digestion
	///
	const UInt illegalKDigestCase = 0;

	///
	/// group 1
	/// group that four shells different with each other
	/// (12|34)
	///
	const UInt case1234 = 1;

	///
	/// group 2
	/// group that three shells are different from each other,
	/// and the final shell is same with one of the three shells.
	/// the same shells are not paired in bra/ket
	///
	/// there are four cases in this group, this is case (12|13)
	/// (12|31) etc.
	///
	/// we have additional digestion for K11 
	/// both P23 and P32 needed to be counted in 
	///
	const UInt case1213 = 2;
	const UInt case1231 = 3;
	const UInt case1223 = 4;
	const UInt case1232 = 5;

	///
	/// group 3
	///
	/// group that three shells are different from each other,
	/// and the final shell is same with one of the three shells.
	/// the same shells are paired in bra/ket cases
	///
	/// this is the case of (11|23) and (12|33)
	/// we only take the K12 and K13 one time
	/// the other pair is redundant
	///
	const UInt case1123 = 6;
	const UInt case1233 = 7;

	///
	/// group 4
	///
	/// group that there are three shells same with each other
	/// and the other one is different from the three shells
	///
	/// this is the case (11|12) etc.
	/// 
	/// we only take K11 and K12/K21, however, there's additional
	/// digestion for K11, on both P12 and P21
	///
	const UInt case1112 = 8;
	const UInt case1121 = 9;
	const UInt case1211 = 10;
	const UInt case2111 = 11;

	///
	/// group 5
	///
	/// in this group bra and ket is composed by same shells
	/// but the shell pairs are different from each other 
	///
	/// the shell quartet is like (11|22)
	///
	/// we only take one K12, all of others are redundant
	///
	const UInt case1122 = 12;

	///
	/// group 6
	///
	/// group that bra and ket are same
	/// but the shell pairs is composed by different shells
	///
	/// shell quartet is as (12|12)
	///
	/// we note, that the case (21|12), (12|21) do not exist
	/// for a given shell pair, if shells pair of |12) exists 
	/// then the pair of |21) does not exsit. Therefore (21|12)
	/// and (12|21) will not be in the calculation because
	/// of the symmetry character of shell quartets
	///
	/// we only take K11, K22 and one K12 (the other one is redundant) 
	///
	const UInt case1212 = 13;

	///
	/// group 7
	/// one shell cases
	/// we only take one K11, others are redundant
	///
	const UInt case1111 = 14;

	///
	/// this is the driver class to perform symmetrical JK digestion
	/// for energy/fock matrix calculation. 
	///
	/// This class is used to digest the input raw integrals on shell quartets,
	/// into the global result matrix(input), or into the result atom shell block 
	///
	class SymmJKDigest {

		protected:

			//
			// general information
			//
			bool withLocResult;///< whether we do digest to local result and update AtomFockMtrx result?
			bool isCloseShell; ///< whether this is close shell system (beta density == alpha density)?
			UInt nSpin;        ///< number of spin states
			Double JThresh;    ///< threshold value for J digestion work
			Double KThresh;    ///< threshold value for K digestion work
			Double kCoef;      ///< the linear coefficients before the K matrix/exchange energy, for hybrid functional
			UInt job;          ///< the job name
			UInt order;        ///< the job order

			//
			// basic information for the given shell quartet (mu,nu|lambda,eta)
			// this is also the shell quartet (12|34)
			//
			// we note that the "local offset" is based on atom shell. This is 
			// the basis set offset for this shell in the corresponding atom shell
			//
			UInt mu;           ///< global Cartesian basis set offset for shell mu
			UInt nu;           ///< global Cartesian basis set offset for shell nu
			UInt lam;          ///< global Cartesian basis set offset for shell lambda
			UInt eta;          ///< global Cartesian basis set offset for shell eta
			UInt lmu;          ///< lobal  Cartesian basis set offset for shell mu
			UInt lnu;          ///< local  Cartesian basis set offset for shell nu
			UInt llam;         ///< local  Cartesian basis set offset for shell lambda
			UInt leta;         ///< local Cartesian basis set offset for shell eta
			UInt nmu;          ///< number of Cartesian basis sets for shell mu
			UInt nnu;          ///< number of Cartesian basis sets for shell nu
			UInt nlam;         ///< number of Cartesian basis sets for shell lambda
			UInt neta;         ///< number of Cartesian basis sets for shell eta

			//
			// K digestion, the case number as well as block form result and P
			// P are block form of density matrix, and K are block form of results
			//
			UInt caseNum;      ///< given by mu,nu,lam,eta; which case is it (see above macros)
			Mtrx aP13;         ///< P13*(12|34) = K24 for alpha state
			Mtrx aP14;         ///< P14*(12|34) = K23 for alpha state
			Mtrx aP23;         ///< P23*(12|34) = K14 for alpha state
			Mtrx aP24;         ///< P24*(12|34) = K13 for alpha state
			BlockMtrx aK13;    ///< P24*(12|34) = K13 for alpha state
			BlockMtrx aK14;    ///< P23*(12|34) = K14 for alpha state
			BlockMtrx aK23;    ///< P14*(12|34) = K23 for alpha state
			BlockMtrx aK24;    ///< P13*(12|34) = K24 for alpha state
			Mtrx bP13;         ///< P13*(12|34) = K24 for beta  state
			Mtrx bP14;         ///< P14*(12|34) = K23 for beta  state
			Mtrx bP23;         ///< P23*(12|34) = K14 for beta  state
			Mtrx bP24;         ///< P24*(12|34) = K13 for beta  state
			BlockMtrx bK13;    ///< P24*(12|34) = K13 for beta  state
			BlockMtrx bK14;    ///< P23*(12|34) = K14 for beta  state
			BlockMtrx bK23;    ///< P14*(12|34) = K23 for beta  state
			BlockMtrx bK24;    ///< P13*(12|34) = K24 for beta  state

			//
			// scratch data for J/K digestion, form the block form of result matrix for K13 etc.
			// this is used when the integral calculation will be dispersed into atom shell dimension blocks
			//
			Mtrx P12;          ///< P12*(12|34) = J34
			Mtrx P34;          ///< P34*(12|34) = J12
			BlockMtrx J12;     ///< P34*(12|34) = J12
			BlockMtrx J34;     ///< P12*(12|34) = J34

			///
			/// J digestion factor determiniation for Fock/Energy calculation
			/// 
			/// -  shell pairs diagonally on the bra/ket side, contribution is  1 
			///
			/// -  offdiagonal shell pairs on the bra/ket side, contribution is 2 
			///
			/// -  if bra is same with ket, then one side digestion is duplicate
			/// in default, we alwasy take the bra side digestion, and get rid of 
			/// ket side if it's redundant
			///
			/// braFac is related to bra side digestion, that is to say; P12 -> J34
			/// ketFac is related to ket side digestion, that is to say, P34 -> J12
			///
			void setJDigestFac(Double& braFac, Double& ketFac) const {
				braFac = TWO;  
				ketFac = TWO;  
				if (mu == nu)   braFac = ONE;
				if (lam == eta) ketFac = ONE;
				if (mu == lam && nu == eta) {
					ketFac = ZERO;
				}
			};

			///
			/// this function is to determine the K digestion factor for digestion
			/// on k13, k14 etc. for Fock matrix building and energy calculation
			///
			/// therefore, we should have totally 4 factors
			///
			void setKDigestFac(Double& k13Fac, Double& k14Fac, Double& k23Fac, Double& k24Fac) const {

				//
				// K13: P24*(12|34) = K13
				// we do all of cases, for case1213, case1112, case1211 we need to double
				// the digestion since ADDITIONAL_DIGESTION needed for K13
				//
				k13Fac = ONE;
				if (caseNum == case1213 || caseNum == case1112 || caseNum == case1211) k13Fac = TWO;

				//
				// K14: P23*(12|34) = K14
				// case1233, case1211, case2111, case1122, case1111 are redundant here
				// so their factor is ZERO. 
				// for case1231 and case1121 digestion needs ADDITIONAL_DIGESTION for K14
				// all of others are normal ones, namely
				// case1234 case1213 case1223 case1232 case1123 case1112 case1212
				//
				k14Fac = ONE;
				if (caseNum == case1231 || caseNum == case1121) {
					k14Fac = TWO;
				}else if (caseNum == case1233 || caseNum == case1211 || caseNum == case2111 
						|| caseNum == case1122 || caseNum == case1111) {
					k14Fac = ZERO;
				}

				//
				// K23: P14*(12|34) = K23
				// case1123, case1121, case1112, case1122, case1212, case1111 are redundant here
				// so their factor is ZERO. That's the default value
				// for case1223 and case2111 digestion needs ADDITIONAL_DIGESTION for K23
				// all of others are normal ones, namely
				// case1234 case1213 case1231 case1232 case1233 case1211
				//
				k23Fac = ZERO;
				if (caseNum == case1234 || caseNum == case1213 || 
						caseNum == case1231 || caseNum == case1232 ||
						caseNum == case1233 || caseNum == case1211) {
					k23Fac = ONE;
				}else if (caseNum == case1223 || caseNum == case2111) {
					k23Fac = TWO;
				}

				//
				// K24: P13*(12|34) = K24
				// case1123, case1233, case1112, case1121, case1211, case2111, case1122 
				// and case1111 are redundant here
				// so their factor is ZERO. That's the default value
				//
				// for case1232 digestion needs ADDITIONAL_DIGESTION for K24
				// all of others are normal ones, namely
				// case1234, case1213, case1231, case1223, case1212 
				//
				k24Fac = ZERO;
				if (caseNum == case1234 || caseNum == case1213 || caseNum == case1231 || caseNum == case1223 ||
						caseNum == case1212) {
					k24Fac = ONE;
				}else if (caseNum == case1232) {
					k24Fac = TWO;
				}
			};

		public:

			///
			/// constructor, doing allocation for memory for the matrix inside
			///
			SymmJKDigest(const GIntJobInfor& infor, UInt maxNBasis); 

			///
			/// destructor
			///
			~SymmJKDigest() { };

			///
			/// initialize the information based on the input shell pair data
			/// and initialize the block form of density matrix data
			///
			void initShellInfor(bool switchSP, const ShellPair& braSP, const ShellPair& ketSP, 
					const Mtrx& PA, const Mtrx& PB);

			///
			///  this is the driver function to perform symmetrical K digestion(only)
			///  over shell quartet. The result after digestion only stores locally. 
			///  \param intArray    raw integral array
			///
			void symmKIntsDigest(const Double* intArray);

			///
			///  this is the working function to perform symmetrical J digestion
			///  over shell quartet. The result after digestion only stores locally.
			///  \param intArray    raw integral array
			void symmJIntsDigest(const Double* intArray);

			///
			///  this is the working function to perform symmetrical JK digestion
			///  together over shell quartet. The result after digestion only 
			///  stores locally.
			///  \param intArray    raw integral array
			void symmJKIntsDigest(const Double* intArray);

			///
			/// this function is used to update the local data to the global Fock matrix FA and FB
			///
			void updateFockMatrix(Mtrx& FA, Mtrx& FB) const;

			///
			/// if the result is digested locally, we will update the results
			/// into the fock matrix pieces in atom shell dimension. So here
			/// we will give the reference of result
			///
			const BlockMtrx& getK13(UInt iSpin) const {
				if (iSpin == 1) return bK13;
				return aK13;
			};

			///
			/// if the result is digested locally, we will update the results
			/// into the fock matrix pieces in atom shell dimension. So here
			/// we will give the reference of result
			///
			const BlockMtrx& getK14(UInt iSpin) const {
				if (iSpin == 1) return bK14;
				return aK14;
			};

			///
			/// if the result is digested locally, we will update the results
			/// into the fock matrix pieces in atom shell dimension. So here
			/// we will give the reference of result
			///
			const BlockMtrx& getK23(UInt iSpin) const {
				if (iSpin == 1) return bK23;
				return aK23;
			};

			///
			/// if the result is digested locally, we will update the results
			/// into the fock matrix pieces in atom shell dimension. So here
			/// we will give the reference of result
			///
			const BlockMtrx& getK24(UInt iSpin) const {
				if (iSpin == 1) return bK24;
				return aK24;
			};

			///
			/// if the result is digested locally, we will update the results
			/// into the fock matrix pieces in atom shell dimension. So here
			/// we will give the reference of result
			///
			const BlockMtrx& getJ12() const { return J12; };

			///
			/// if the result is digested locally, we will update the results
			/// into the fock matrix pieces in atom shell dimension. So here
			/// we will give the reference of result
			///
			const BlockMtrx& getJ34() const { return J34; };
	};

	///
	/// here let's define a simple constant, because the J/K digestion
	/// are designed for 4 body integrals, therefore the deriv body is 4
	///
	const UInt SYMM_JK_DERIV_NBODY = 4;

	///
	/// this is the driver class to perform symmetrical JK digestion
	/// for integral derivatives
	///
	class SymmJKDerivDigest : public SymmJKDigest {

		private:

			//
			// these are the information of atom centers for each position
			//
			UInt bra1Pos;         ///< the atom center index for bra1 position
			UInt bra2Pos;         ///< the atom center index for bra2 position
			UInt ket1Pos;         ///< the atom center index for ket1 position
			UInt ket2Pos;         ///< the atom center index for ket2 position

			//
			// raw integral information
			//
			DoubleVec jVec;       ///< scratch array for doing J digestion

			//
			// the result of deriv digestion for current shell quartet
			//
			Mtrx derivResult;     ///< the result for deriv. (4,3) for order = 1, (12,12) for order=2
			
			///
			///  this is the working function to perform symmetrical K digestion
			///  over integral derivatives. It performs for each given derivative array,
			///  therefore at the end we return the derivative value for the given derivatives
			///
			///  \param intArray     raw integral array
			///  \return the derivatives value for K digestion
			///
			Double symmKIntsDerivDigest(const Double* intArray) const;

			///
			/// because of the symmetry of the integrals, that is to consider:
			///
			/// (12|34) <=> (12|43) <=> (21|34) <=> (21|43) <=> 
			/// (34|12) <=> (43|12) <=> (34|21) <=> (43|21)
			///
			/// all of above eight integrals are same(repeat integrals). However, in the practical
			/// calculation to avoid the repeat of the integral calculation we only
			/// pick up one of those. This is done in forming the significant shell pair
			/// information and forming the loop of the shell quartet calculation. 
			///
			/// For the derivatives digestion, the J/K derivatives is like this.
			/// For example, the \f$ deriv_{J}\f$ is like:
			///
			/// \f$ \sum_{\mu\nu}\sum_{\lambda\eta}P_{\mu\nu}P_{\lambda\eta}(\mu\nu|\lambda\eta)^{X}\f$
			///
			/// and \f$ deriv_{K}\f$ is like:
			///
			/// \f$ \sum_{\mu\nu}\sum_{\lambda\eta}P_{\mu\lambda}P_{\nu\eta}(\mu\nu|\lambda\eta)^{X}\f$
			///
			/// that's for the first order derivatives, but second derivative is same. 
			///
			/// for derivative result the summation goes over all of index for each loop, and we have
			/// four loops. Therefore the derivatives requires that we need to count in
			/// all of other identical integrals when digesting. Therefore that's the reason
			/// we have the scaling coefficients.
			///
			/// However, the "identical" integrals in terms of digestion are different between J/K. Let's
			/// take case1234 as example. For J digestion, all of above 8 integrals are identical,
			/// therefore in terms of the 1/2 in front of energy, the scaling factor is four.
			///
			/// However for exchange things is different. For each integral \f$(\mu\nu|\lambda\eta)^{X}\f$,
			/// suggesting this is case1234; therefore all of four basis sets are different from each other.
			/// we have two digestion, P13*P24*(12|34) and P24*P13*(12|34). From the above \f$ deriv_{K}\f$
			/// definition, it's clear that the index 1, 3 and index 2, 4 are symmetrical and exchangable,
			/// and the index 2, 3 and index 1, 4 are exchangable. Therefore, we can re-arrange the 
			/// above repeating integrals of (12|34) as:
			///
			/// (12|34) <=> (34|12) <=> (21|43) <=> (43|21) (identical integrals on index 1, 3 and 2, 4)
			///
			/// and (12|43) <=> (34|21) <=> (21|34) <=> (43|12) (identical integrals on index 1, 4 and 2, 3)  
			///
			/// therefore, for case1234 P13*P24 has 4 identical integrals, so as P14*P23. So the scaling 
			/// factor is 4. considering the 1/2 in front of the enegy, the final coefficient is 2.
			///
			/// this is the example for case1234, all of other cases are similar.
			///
			/// we note, that here what we considered is the symmetry of integrals. 
			///
			void findJKDerivCoef(Double& jCoefs, Double& kCoefs) const {

				// let's determine the scaling coefficients for derivatives
				// we take the similar way as shown above.
				//
				// for case1213, case1231, case1223 and case1232; 
				// they are same with case1234. 
				//
				// for case1123, case1233, case1112, case1121,
				// case1211 and case2111, it's clear that there are
				// only "four" identical integrals for J, so after 
				// multiply with 1/2 the J the coefficients is TWO.
				// 
				// But because that P13*P24 is repeating P23*P14 
				// for all of above cases, so for K exchange there
				// are also 4 integrals, too. it's easy to see that the 
				// K coefficient here is ONE. 
				//
				// for case1122, there's only two identical integrals
				// for J, so the J coefficients is 1. For K because 
				// of the P13*P24 = P23*P14, there are also 2 identical
				// integrals therefore coefficients is 0.5.
				//
				// for case1212, there's only four identical integrals
				// for J so the J coefficients is 2. For K digestion,
				// (12|12) <=> (21|21) on index 1,3 and 2, 4;
				// (21|12) <=> (12|21) on index 1,4 and 2, 3.
				// therefore also considering the 1/2 in front of energy,
				// the coefficient should be 1.
				//
				// for case1111, only one identical integral exists.
				// so J coefficient is 0.5. For K because there's repeat
				// on P13*P24 = P23*P14, the coefficient is 0.25.
				//
				jCoefs = TWO;    
				kCoefs = ONE;    
				if (caseNum == case1234 || caseNum == case1213 || 
						caseNum == case1231 || caseNum == case1223 || caseNum == case1232) {
					jCoefs = FOUR;    
					kCoefs = TWO;    
				}else if (caseNum == case1122) {
					jCoefs = ONE;    
					kCoefs = HALF;    
				}else if (caseNum == case1212) {
					jCoefs = TWO;    
					kCoefs = ONE;    
				}else if (caseNum == case1111) {
					jCoefs = HALF;    
					kCoefs = HALF/TWO;    
				}
			};

			///
			/// for the given elementary derivative information, we return the row and column
			/// index for the both local result and global result
			///
			void findRowColIndex(const ElemDerivInfor& elemDeriv, 
					UInt& localRowIndex, UInt& localColIndex, 
					UInt& globalRowIndex, UInt& globalColIndex) const;

		public:

			///
			/// constructor, doing allocation for memory for the matrix inside
			///
			SymmJKDerivDigest(const GIntJobInfor& infor, UInt maxNBasis):SymmJKDigest(infor,maxNBasis),
			bra1Pos(-1),bra2Pos(-1),ket1Pos(-1),ket2Pos(-1),jVec(maxNBasis*maxNBasis,ZERO),
			derivResult(SYMM_JK_DERIV_NBODY*3, SYMM_JK_DERIV_NBODY*3) { 
				if (order == 1) derivResult.reset(SYMM_JK_DERIV_NBODY,3,false);
			};

			///
			/// destructor
			///
			~SymmJKDerivDigest() { };

			///
			/// initialize the derivatives information based on the input shell pair data
			/// we will consider the shell pair switch inside
			///
			/// also we form the block form of density matrices
			///
			void initDerivInfor(bool switchSP, 
					const ShellPair& braSP, const ShellPair& ketSP, const AtomShellPair& braAtomSP,
					const AtomShellPair& ketAtomSP, const Mtrx& PA, const Mtrx& PB);

			///
			/// driver function to digest the integrals derivatives
			///
			/// it combines the raw integral array value with the density matrix (from alphaDen and 
			/// betaDen), finally form the derivatives result.
			///
			/// based on the integral infor, we can figure out the redundant integrals simply,
			/// see the function inside
			///
			/// for order = 1, result is (nAtoms,3)
			///
			/// for order = 2, result dimension is (3*nAtoms,3*nAtoms), the hessian matrix.
			/// Because it's symmetrical, we only need to lower triangular part
			///
			void digestJKInts(const SingleIntegralInfor& infor, const Double* rawIntArray, Mtrx& result); 
	};

};

#endif
