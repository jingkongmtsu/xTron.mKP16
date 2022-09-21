/**
 * \file   integralinfor.h
 * \brief  collecting information for integral engine
 * \author Fenglai Liu 
 */
#ifndef INTEGRALINFOR_H
#define INTEGRALINFOR_H
#include <string>
#include <vector>
#include "libgen.h"
#include "tbb/scalable_allocator.h"
using namespace std;

namespace integralinfor {

	///
	/// this is to define the derivatives information
	/// the derivatives on X, Y or Z?
	///
	/// additionally, for an arbitrary shell quartet
	/// we also define the position code, like BRA1 etc.
	///
	/// NULL_POS is used to characterize the null position
	/// information
	///
	const UInt NO_DERIV       =   0;
	const UInt GINTS_DERIV_X  =   1;
	const UInt GINTS_DERIV_Y  =   2;
	const UInt GINTS_DERIV_Z  =   3;
	const UInt NULL_POS       =   0;
	const UInt BRA1           =   10;
	const UInt BRA2           =   20;
	const UInt KET1           =   30;
	const UInt KET2           =   40;
	const UInt OPERATOR       =   50;

	///
	/// this value stores the maximum number of unique derivatives
	/// when composing a redundant derivatives through translational
	/// invariance on the RHS.
	///
	/// because the integral may have 4 bodies at most, therefore
	/// the first order redundant derivatives can be written as 
	/// sum of three unique derivatives. For second derivative,
	/// each of the term may also be expanded into another three
	/// terms so the total number of derivatives is 9
	///
	const UInt MAX_NUM_RHS_DERIV = 9;

	///
	/// declare the class
	///
	class ElemDerivInfor;

	///
	/// declare it's vector
	///
	typedef std::vector<ElemDerivInfor,tbb::scalable_allocator<ElemDerivInfor> > ElemDerivInforVec;

	///
	/// declare of the redundant derivatives information class
	///
	class RedundantDerivInfor;

	///
	/// this is for the vector of RedundantDerivInfor
	///
	typedef std::vector<RedundantDerivInfor,tbb::scalable_allocator<RedundantDerivInfor> > RedDerivInforVec;

	///
	/// declare of the single integral infor class
	///
	class SingleIntegralInfor;

	///
	/// now this is the vector of single integral infor
	///
	typedef std::vector<SingleIntegralInfor,tbb::scalable_allocator<SingleIntegralInfor> > IntegralInforVec;

	/**
	 * this class stores the elementary derivative information for each derivative component, e.g. xy.
	 * So that's the class name meaning: elementary derivatives information
	 *
	 * this class supports derivatives up to the 2ed order. If you have 
	 * the derivatives larger than 2ed, you need to expand this class 
	 * to have more information
	 *
	 * if this is 1st order derivative, all of information regarding the 
	 * second position is null. Similarly, if this is for energy calculation
	 * all of default information is null
	 *
	 */
	class ElemDerivInfor {

		public:

			///Possible values for the variables here are constants defined at the beginning of this file.
			UInt firstDerivPos;     ///< bra1, bra2, ket1 or ket2 on the first derivatives position? 
			UInt secondDerivPos;    ///< bra1, bra2, ket1 or ket2 on the second derivatives position?
			UInt firstDerivDir;     ///< x, y or z for first derivatives position?
			UInt secondDerivDir;    ///< x, y or z for second derivatives position?

			///
			/// constructor - all of data will be updated later
			///
			ElemDerivInfor():firstDerivPos(NULL_POS),secondDerivPos(NULL_POS),
			firstDerivDir(NO_DERIV),secondDerivDir(NO_DERIV) { };

			///
			/// add in the first pos information
			///
			void update1stPos(UInt pos) { firstDerivPos  = pos; };

			///
			/// add in the second pos information
			///
			void update2edPos(UInt pos) { secondDerivPos = pos; };

			///
			/// add in the first dir information
			///
			void update1stDir(UInt dir) { firstDerivDir  = dir; };

			///
			/// add in the second dir information
			///
			void update2edDir(UInt dir) { secondDerivDir = dir; };

			///
			/// destructor
			///
			~ElemDerivInfor() { };

			///
			/// get the deriv pos for the second derivatives
			///
			void getDerivPos(UInt& pos1, UInt& pos2) const { 
				pos1 = firstDerivPos;  
				pos2 = secondDerivPos;     
			};

			///
			/// get the deriv pos for the first derivatives
			///
			UInt getDerivPos() const { 
				return firstDerivPos;
			};

			///
			/// get the deriv direction for given index
			///
			void getDerivDir(UInt& dir1, UInt& dir2) const { 
					dir1 = firstDerivDir;
					dir2 = secondDerivDir;
			};

			///
			/// get the deriv direction for given index
			///
			UInt getDerivDir() const { 
					return firstDerivDir;
			};

			///
			/// get the job order
			///
			UInt getOrder() const {
				UInt order = 0;
				if (firstDerivPos != NULL_POS) order++;
				if (secondDerivPos != NULL_POS) order++;
				return order;
			};

			///
			/// operator ==.  Used to determine if a derivative is redundant, i.e. can be calculated using
			/// translational invariance.
			///
			bool operator==(const ElemDerivInfor& infor0) const {

				// if all content are same, 
				// both of them are obvviously same
				if (firstDerivPos    == infor0.firstDerivPos &&
						firstDerivDir  == infor0.firstDerivDir &&
						secondDerivPos == infor0.secondDerivPos &&
						secondDerivDir == infor0.secondDerivDir) return true;

				// now consider interchangable
				if (firstDerivPos    == infor0.secondDerivPos &&
						firstDerivDir  == infor0.secondDerivDir &&
						secondDerivPos == infor0.firstDerivPos &&
						secondDerivDir == infor0.firstDerivDir) return true;

				// now return
				return false;
			};

			///
			/// whether the first deriv position is same with the 
			/// second one?
			///
			bool onSamePos() const {
				if (firstDerivPos == secondDerivPos && firstDerivPos != NULL_POS) return true;
				return false;
			};

			///
			/// whether the first deriv direction is same with the 
			/// second one?
			///
			bool onSameDir() const {
				if (firstDerivDir == secondDerivDir && firstDerivDir != NO_DERIV) return true;
				return false;
			};

			///
			/// print out for debug purpose
			///
			void debugPrint() const;
	};

	/**
	 * infor class complement to the elementary derivatives information. used 
	 * to form the redundant derivatives based on translational invariance.  The redundant
	 * quantities to be calculated on the LHS.  The quantities calculated directly are
	 * on the RHS.
	 *
	 * The elementary derivatives is the integral derivatives directly
	 * calculated from integrals, however; the redundant derivatives
	 * is able to be computed from elementary derivatives by the principle
	 * of translational invariance.
	 */
	class RedundantDerivInfor : public ElemDerivInfor {

		public:

			UInt redundantPos;     ///< the redundant position, bra1, bra2, ket1 or ket2?
			UInt rhsDerivLength;   ///< the length of rhs derivatives
			ElemDerivInfor rhsDeriv[MAX_NUM_RHS_DERIV]; ///< the RHS derivatives information
			UInt rhsDerivPos[MAX_NUM_RHS_DERIV];  ///< for each deriv on the rhs, it's position in array
			Double rhsCoefs[MAX_NUM_RHS_DERIV];   ///< for each deriv on the rhs, it's coefficients

			///
			/// constructor 
			///
			RedundantDerivInfor(const UInt& oper, const UInt& pos, 
					const ElemDerivInfor& deriv);

			///
			/// destructor
			///
			~RedundantDerivInfor() { };

			///
			/// update the information rhs derivative positions
			/// with the input vector
			///
			void updatePos(const ElemDerivInforVec& derivVec);

			///
			/// return the given term of rhs derivatives
			///
			const ElemDerivInfor& getRHSDeriv(UInt i) const {
				return rhsDeriv[i];
			};

			///
			/// return the given term of rhs derivatives coefficient
			///
			Double getRHSDerivCoef(UInt i) const { return rhsCoefs[i]; };

			///
			/// return the given term of rhs derivatives position
			///
			UInt getRHSDerivPos(UInt i) const { return rhsDerivPos[i]; };

			///
			/// return the length of rhs
			///
			UInt getRHSLen() const { return rhsDerivLength; };

			///
			/// print out for debug purpose
			///
			void print() const;
	};

	/**
	 * this class is used to store a single integral file information
	 * based on the two classes above. It contains all the data for
	 * one type of LCode.
	 * for example, for the hgp_os engine the information regarding 
	 * the hgp_os_d_d_d_d.cpp.
	 */
	class SingleIntegralInfor {

		private:

			LInt LCode;          ///< the shell quartet code
			UInt oper;           ///< operator name
			UInt order;          ///< derivatives order

			///
			/// whether this integral file use the memory allocator,
			/// in this case the integral calculation will use array form
			/// variables.
			///
			/// if so, how much memory (in sizeof(Double)) it needs
			///
			UInt memAllocLen;    

			///
			/// unique derivatives information, in array form
			///
			ElemDerivInforVec derivInforArray;

			///
			/// redudant derivatives information
			///
			RedDerivInforVec redDerivInforArray;

			///
			/// based on the read in unique derivatives information,
			/// we are able to form the redundant derivatives information
			///
			void formRedundantDerivInfor();

			///
			/// check the 1st derivatives information
			///
			void check1stDeriv() const;

			///
			/// check the 2ed derivatives information
			///
			void check2edDeriv() const;

		public:

			///
			/// constructor
			///
			SingleIntegralInfor(const LInt& code0, const UInt& oper0, const UInt& order0):LCode(code0),
			oper(oper0),order(order0),memAllocLen(0) { };

			///
			/// destructor
			///
			~SingleIntegralInfor() { };

			///
			/// from the given setting file, we initialize memory length
			///
			void initMemInfor(const string& fileName);

			///
			/// from the given setting file, we initialize the deriv information
			/// array
			///
			void initDerivInfor(const string& fileName);

			///
			/// get the LCode
			///
			LInt getLCode() const { return LCode; };

			///
			/// return the memory usage
			///
			UInt getMemInfor() const { return memAllocLen; };

			///
			/// return the length of elemDerivInfor
			///
			UInt getElemDerivInforLen() const { return derivInforArray.size(); };

			///
			/// return the given elemDerivInfor term
			///
			const ElemDerivInfor& getElemDerivInfor(UInt i) const { return derivInforArray[i]; };

			///
			/// return the length of redundant derivatives infor
			///
			UInt getRedDerivInforLen() const { return redDerivInforArray.size(); };

			///
			/// return the term of redundant derivatives infor
			///
			const RedundantDerivInfor& getRedDerivInfor(UInt i) const { return redDerivInforArray[i]; };

			///
			/// print function
			///
			void print() const;
	};

	/**
	 * Finally this class stores all of information related to the integral engine
	 * 
	 * right now it's binded to the hgp_os engine in the gints_engine folder
	 */
	class IntegralInfor {

		private:

			UInt oper;                       ///< operator name
			UInt order;                      ///< derivatives order
			UInt maxL;                       ///< maximum angular momentum supported
			UInt maxMemLenForSP;             ///< maximum memory length for max L = S/P
			UInt maxMemLenForD;              ///< maximum memory length for max L = D
			UInt maxMemLenForF;              ///< maximum memory length for max L = F
			UInt maxMemLenForG;              ///< maximum memory length for max L = G
			UInt maxMemLenForH;              ///< maximum memory length for max L = H
			UInt maxMemLenForI;              ///< maximum memory length for max L = I
			UInt nElemDerivs;                ///< number of elementary derivatives
			UInt nRedDerivs;                 ///< number of redundant derivatives
			IntegralInforVec integralsInfor; ///< it stores every integral's information

		public:

			///
			/// constructor
			///
			IntegralInfor(const UInt& oper0, const UInt& order0);

			///
			/// destructor
			///
			~IntegralInfor() { };

			///
			/// get the specific integral infor by given the LCode
			///
			const SingleIntegralInfor& getIntegralInfor(const LInt& code) const;

			///
			/// return the maximum L supported here
			///
			UInt getMaxL() const { return maxL; };

			///
			/// return the number of redundant derivatives
			///
			UInt getNRedDerivs() const { return nRedDerivs; };

			///
			/// return the number of elementary derivatives
			///
			UInt getNElemDerivs() const { return nElemDerivs; };

			///
			/// return the memory usage for the given shell
			/// 
			/// the input L is the maximum L for the given integrals calculation 
			///
			UInt getMaxMemLen(UInt L) const; 

			///
			/// return the memory usage for the given shell quartet L Code
			///
			UInt getMemUsage(const LInt& code) const;

			///
			/// for debug purpose, print
			///
			void print() const;
	};

}

#endif

