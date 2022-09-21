/**
 * \file    shellprop.h
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef SHELLPROP_H
#define SHELLPROP_H

#include "libgen.h"
#include <string>

namespace shellprop {

	using namespace std;

	/************************************************************************
	 *                          Data for Shell Structure                    *
	 ************************************************************************/
	/**
	 * shell types 
	 * 
	 * By given this shell types, we are able to calculate the number of 
	 * basis set functions for each shell(n is the given angular momentum
	 * number):
	 * 
	 * Cartesian:
	 * \f$\frac{(n+1)(n+2)(n+3)-n(n+1)(n+2)}{6}\f$
	 *
	 * Spherical:
	 * \f$(n+1)^{2}-n^{2}\f$
	 *
	 */

	/**
	 * maximum high angular momentum allowed in the program
	 * G: 4, H: 5, I: 6
	 *
	 * we will impose the check of angular momentum in the 
	 * shell data forming section(see shell.cpp).
	 */
	const UInt MAX_L = 4;

	/**
	 * this shell pair order array defines the 
	 * shell pair angular momentum order
	 *
	 * for two arbitrary shell pairs, it answers 
	 * the question that which one should be ahead,
	 * and which one should be after.
	 *
	 * the data here is same with the shell pair order
	 * defined in cppints program. However, over there
	 * we have angular momentum up to 20 but here we 
	 * only have up to the limit of MAX_L defined above.
	 *
	 * the array is generated from the folder 
	 * shell/util/shell_pair_order. Please see the code
	 * over there.
	 *
	 * below the comment part is the Lmin and Lmax
	 * corresponding to the LCode. For example,
	 * 100001  , //  1   1   0   1
	 * the LCode is 100001, the bra1 shell is P 
	 * (lmin=1 lmax=1); and the bra2 shell is 
	 * SP (0 1).
	 *
	 * important note: in this program we do not have G, H, I etc.
	 * high angular momentum paired with SP shell. This is for 
	 * the reason to reduce the number of integral files.
	 */
	const Int SHELL_PAIR_ORDER_ARRAY[ ]  =  { 
		0,      //  0   0   0   0 
		100,    //  0   1   0   0 
		1,      //  1   1   0   0 
		2,      //  2   2   0   0 
		100100, //  0   1   0   1 
		100001, //  1   1   0   1 
		1001,   //  1   1   1   1 
		3,      //  3   3   0   0 
		100002, //  2   2   0   1 
		1002,   //  2   2   1   1 
		4,      //  4   4   0   0 
		100003, //  3   3   0   1 
		1003,   //  3   3   1   1 
		2002,   //  2   2   2   2 
		1004,   //  4   4   1   1 
		2003,   //  3   3   2   2 
		2004,   //  4   4   2   2 
		3003,   //  3   3   3   3 
		3004,   //  4   4   3   3 
		4004,   //  4   4   4   4 
	};

	///
	/// this is the maximum shell pair number in the 
	/// array of SHELL_PAIR_ORDER_ARRAY
	///
	const UInt MAX_SHELL_PAIR_NUMBER = 20;

	///
	/// according to the above shell pair order array, here
	/// let's define the composite shell pair information
	/// 
	/// firstly, this is how many composite shell pairs
	///
	const UInt MAX_COM_SHELL_PAIR_NUMBER = 5;

	///
	/// this is composite shell pair list
	///
	/// let's note again, we do not have G, H, I etc. high angular momentums
	/// paired with SP shell
	///
	const UInt COMPOSITE_SHELL_PAIR_ORDER_ARRAY[ ]  =  { 
		100,    //  0   1   0   0   (SP|S)
		100100, //  0   1   0   1   (SP|SP)
		100001, //  1   1   0   1   (P|SP)
		100002, //  2   2   0   1   (D|SP)
		100003  //  3   3   0   1   (F|SP)
	};

	/**
	 * this function return the position of the given angular code in 
	 * the array of COMPOSITE_SHELL_PAIR_ORDER_ARRAY
	 */
	inline UInt getCompositeLCodeIndex(const UInt& angCode) {
		for(UInt i=0; i<MAX_COM_SHELL_PAIR_NUMBER; i++) {
			if (COMPOSITE_SHELL_PAIR_ORDER_ARRAY[i] == angCode) return i;
		}
		return -1;
	};

	/**
	 * given the shell name of s, return the lmin and lmax
	 */
	void getL(const string& s, UInt& lmin, UInt& lmax);

	/**
	 * giving the lmin and lmax, we get the shell name
	 */
	string getShellName(const UInt& lmin, const UInt& lmax);

	/**
	 * number of Cartesian type of basis set functions
	 */
	inline UInt getCartBas(const UInt& lmin, const UInt& lmax) {
		return ((lmax+1)*(lmax+2)*(lmax+3)-lmin*(lmin+1)*(lmin+2))/6;
	};

	/**
	 * number of spherical type of basis set functions
	 */
	inline UInt getPureBas(const UInt& lmin, const UInt& lmax) {
		return (lmax+1)*(lmax+1)-lmin*lmin;
	};

	/**
	 * it returns the angular momentum of l,m,n by given the index of basis set 
	 * in the shell. 
	 *
	 * For the given basis set index and we want to know the total angular momentum,
	 * we can just call this getlmn and add all of l,m and n to get L.
	 *
	 * Input:
	 * \param lmin:   minimum angular momentum in the shell
	 * \param lmax:   maximum angular momentum in the shell
	 * \param index:  the basis set index in this shell
	 * 
	 * output:
	 * \param l:      X's angular momentum
	 * \param m:      Y's angular momentum
	 * \param n:      Z's angular momentum
	 *
	 */
	void getlmn(const UInt& lmin, const UInt& lmax, const UInt& index, 
			UInt& l, UInt& m, UInt& n); 

	/************************************************************************
	 *  #####            functions for coding shell                         *
	 *  here the functions below is going to code input shells provided     *
	 *  by the user. Therefore, in a single job batch (defined by the user  *
	 *  input file), one is able to create unique shell code for all of     *
	 *  basis sets defined in the user input file.                          *
	 ************************************************************************/

	/**
	 * by given the section number(from molecule) and it's part(aux shell index),
	 * we are able to form the code for shell
	 *
	 * We note, that the maximum allowed shell number is same with maximum
	 * number of geometries in input file. see molecule.h
	 */
	UInt codeShell(const UInt& section, const UInt& part); 

	/**
	 * by given the shell code, we can retrieve the corresponding section
	 * and part information
	 */
	void decodeShell(const UInt& scode, UInt& section, UInt& part); 

	/************************************************************************
	 *  #####    functions for coding the angular momentum for shell        *
	 *  Since in this program we support composite shells, then for each    *
	 *  shell we have Lmax and Lmin like SP shell. The code for the shell   *
	 *  is an "ID" for the shell, it's unique constant in this program.     *
	 ************************************************************************/

	///
	/// this is the step length between lmin and lmax
	/// the value indicates that we can support L up to 99
	///
	const UInt ANG_UNIT = 100;

	/**
	 * code angular momentum for a shell 
	 */
	inline UInt codeL(UInt lmin, UInt lmax) {
		if (lmin == lmax) {
			return lmax;
		}else{
			return (lmin+lmax*ANG_UNIT);
		}
	};

	/**
	 * decode L from the angular momentum code above
	 */
	inline void decodeL(const UInt& code, UInt& lmin, UInt& lmax) {
		if (code < ANG_UNIT) {
			lmin = code;
			lmax = code;
		}else{
			lmax = code/ANG_UNIT;
			lmin = code - lmax*ANG_UNIT;
		}
	};


	/************************************************************************
	 *  #####  functions for coding the angular momentums for shell pair    *
	 *  Based on the unique shell code above, we can further form the ID    *
	 *  for shell pairs. We note that for the shell pair, the angular       *
	 *  momentum should satisfy the condition that: L(bra1) >= L(bra2)      *
	 *  if both bra1 and bra2's Lmax equal to each other, then we compare   *
	 *  lmin. For more details, please refer to the shell class operator <. *
	 ************************************************************************/

	///
	/// unit for the bra2 position in terms of bra1
	/// this is the basic unit for forming shell pair code
	///
	const UInt LCODE_UNIT_BRA2 = 1000;

	/**
	 * for the input angular momentums, according to the requirement of 
	 * shell i >= shell j, do we need to switch it?
	 * \param  iLmin, iLmax: angular momentum for shell i
	 * \param  jLmin, jLmax: angular momentum for shell j
	 */
	inline bool switchShell(UInt& iLmin, UInt& iLmax, UInt& jLmin, UInt& jLmax) {
		if (jLmax > iLmax || (jLmax == iLmax && jLmin > iLmin)) return true;
		return false;
	};

	/**
	 * code the angular momentum for the shell pairs/two body integrals
	 *
	 *	We note, that this is not return in LInt type. Since it will never
	 *	exceed the size_t range even if size_t is 32 bit
	 *
	 *	We also use it for generating shell pair code too
	 */
	inline UInt codeL(UInt iLmin, UInt iLmax, UInt jLmin, UInt jLmax) {
		if (switchShell(iLmin,iLmax,jLmin,jLmax)) {
			UInt bra1Code = codeL(jLmin,jLmax);
			UInt bra2Code = codeL(iLmin,iLmax);
			UInt code     = bra1Code+bra2Code*LCODE_UNIT_BRA2;
			return code;
		}
		UInt bra1Code = codeL(iLmin,iLmax);
		UInt bra2Code = codeL(jLmin,jLmax);
		UInt code     = bra1Code+bra2Code*LCODE_UNIT_BRA2;
		return code;
	};

	/**
	 * decode the angular momentum for the shell pairs
	 *	i corresponds to row shell, j corresponds to col shell
	 */
	inline void decodeL(const UInt& code, UInt& iLmin, UInt& iLmax, 
			UInt& jLmin, UInt& jLmax) {
		UInt bra2Code = code/LCODE_UNIT_BRA2;
		UInt bra1Code = code-bra2Code*LCODE_UNIT_BRA2;
		decodeL(bra1Code,iLmin,iLmax);
		decodeL(bra2Code,jLmin,jLmax);
	};

	/************************************************************************
	 *  #####    functions for coding the angular momentums for integrals   *
	 *  for a general integral (bra1,bra2|O|ket1,ket2), if bra1 and bra2    * 
	 *  (so as the ket) part switched, the integral will be still same.     *
	 *  Furthermore, if "bra" and "ket" side switched, the integral is      *
	 *  same, too. Such switch operations is symmetrical operations:        *
	 *  (bra1,bra2|O|ket1,ket2) = (bra2,bra1|O|ket1,ket2)                   *
	 *  (bra1,bra2|O|ket1,ket2) = (bra2,bra1|O|ket2,ket1)                   *
	 *  (bra1,bra2|O|ket1,ket2) = (bra1,bra2|O|ket2,ket1)                   *
	 *  (bra1,bra2|O|ket1,ket2) = (bra2,bra1|O|ket2,ket1)                   *
	 *  (bra1,bra2|O|ket1,ket2) = (ket1,ket2|O|bra2,bra1)                   *
	 *  (bra1,bra2|O|ket1,ket2) = (ket2,ket1|O|bra2,bra1)                   *
	 *  (bra1,bra2|O|ket1,ket2) = (ket2,ket1|O|bra1,bra2)                   *
	 *  (bra1,bra2|O|ket1,ket2) = (ket2,ket1|O|bra2,bra1)                   *
	 *                                                                      *
	 *  As a result, we do not need to provide all of angular momentum      *
	 *  combinations for the integral calculation. We only need a unique    *
	 *  selective subset of them where all of other integrals could be      *
	 *  derived from the above symmetrical operations.                      *
	 *                                                                      *
	 *  Here we define the rule to generate the unique subset for integral  *
	 *  files:                                                              *
	 *  rule 1: angular momentum: bra1 >= bra2(if both bra1 and bra2 exists)*  
	 *  rule 2: angular momentum: ket1 >= ket2(if both ket1 and ket2 exists)*
	 *  rule 3: angular momentum sum: bra1 + bra2 >= ket1 + ket2            *
	 *  (if all of four shells exist,only for Lmax)                         *
	 *  rule 4: if bra1 + bra2 == ket1 + ket2(Lmax compare), then           *
	 *  L(bra2) >= L(ket2) if all of four shells exist                      *
	 *  rule 5: if bra1 + bra2 == ket1 + ket2(Lmax compare) and             *
	 *  L(bra2) == L(ket2), then L(bra1) >= L(ket1)                         *
	 *  only the angular momentum codees which satisfy the above rules      *
	 *  has corresponding integral files.                                   *
	 ************************************************************************/

	///
	/// here below we will define the unit to form
	/// the LCode (this is set so that we can support
	/// very large L, e.g. L = 20 etc.)
	///
	/// LCode is representing the code of angular
	/// momentum for a composite shell property
	/// e.g. shell pair etc.
	///
	/// see the function of codeL and decodL for 
	/// more information
	///
	const UInt LCODE_UNIT_KET1 = 1000000;
	const UInt LCODE_UNIT_KET2 = 1000000000;

	/**
	 * This is the inplementation of above rules so that to determine
	 * whether the input two shell pair should be switched or not (bra 
	 * and ket side is switched)
	 *
	 * since in forming shell pairs, rule 1 and rule 2 is already
	 * satisfied, then we only concentrate on rule 3-5.
	 */
	inline bool switchShellPair(UInt& braLmin1, UInt& braLmax1, 
			UInt& braLmin2, UInt& braLmax2, 
			UInt& ketLmin1, UInt& ketLmax1,
			UInt& ketLmin2, UInt& ketLmax2) {

		// consider rule 3
		if (ketLmax1+ketLmax2>braLmax1+braLmax2) return true;

		// consider rule 4 and 5
		if (ketLmax1+ketLmax2==braLmax1+braLmax2) {

			// rule 4
			// compare bra2 and ket2
			// for bra2/ket2, we will take the bigger one as in the bra side
			if (braLmax2 < ketLmax2) return true;
			if (braLmax2 > ketLmax2) return false;

			// now we have braLmax2 == ketLmax2
			// compare the lmin
			if (braLmin2 < ketLmin2) return true;
			if (braLmin2 > ketLmin2) return false;

			// rule 5
			// now in this case, we must have:
			// braLmax1 == ketLmax1
			// braLmax2 == ketLmax2
			// braLmin2 == ketLmin2
			if (braLmin1 < ketLmin1) return true;
		}

		//
		// all of other cases
		//
		return false;
	};

	/**
	 * code the angular momentum for three body integral
	 * we note that shell i should >= shell j here
	 *	the corresponding integral will be (ij|k)
	 */
	inline LInt codeL(UInt iLmin, UInt iLmax, UInt jLmin, UInt jLmax, UInt kLmin, UInt kLmax) {
		LInt bra1Code = codeL(iLmin,iLmax);
		LInt bra2Code = codeL(jLmin,jLmax);
		LInt ket1Code = codeL(kLmin,kLmax);
		LInt code    = bra1Code+bra2Code*LCODE_UNIT_BRA2+ket1Code*LCODE_UNIT_KET1;
		return code;
	};

	/**
	 * code the angular momentum for four body integral
	 * we note that shell i should >= shell j here
	 * we note that shell k should >= shell l here
	 *	the corresponding integral will be (ij|kl)
	 *
	 *	this function is comment out since in practise we 
	 *	should use the next function of codeSPCodes
	 *
	 inline LInt codeL(UInt iLmin, UInt iLmax, UInt jLmin, UInt jLmax, 
	 UInt kLmin, UInt kLmax, UInt lLmin, UInt lLmax) {
	 LInt bra1Code = codeL(iLmin,iLmax);
	 LInt bra2Code = codeL(jLmin,jLmax);
	 LInt ket1Code = codeL(kLmin,kLmax);
	 LInt ket2Code = codeL(lLmin,lLmax);
	 LInt c1       = bra2Code*LCODE_UNIT_BRA2; 
	 LInt c2       = ket1Code*LCODE_UNIT_KET1;
	 LInt c3       = ket2Code*LCODE_UNIT_KET2;
	 return bra1Code + c1 + c2 + c3;
	 };
	 */

	/**
	 * coded two shell pairs code together 
	 * here we do not check codes
	 * be sure that both of the codes should be 
	 * in "unswitched" mode
	 */
	inline LInt codeSPCodes(UInt braShellPairCode, UInt ketShellPairCode) {
		LInt c = ketShellPairCode*LCODE_UNIT_KET1;
		c = c + braShellPairCode;
		return c;
	};

	/************************************************************************
	 *  #####    functions related to using different basis set order       * 
	 ************************************************************************/

	/**
	 * get the basis set index for the given basis set order
	 *
	 * for the basis set order information, please see the file angmomlist.h
	 *
	 * \param L    : the shell angular momentum
	 * \param index: the index of the basis set for the given shell (the order must be in main basis set order)
	 * \param order: which basis set order are you going to?
	 * \return       the basis set index for the shell in terms of the input basis set order
	 */
	UInt getBasIndexForTheBasisSetOrder(const UInt& L, const UInt& index, const UInt& order);

	/************************************************************************
	 *  #####    functions related to normalization factors etc.            * 
	 ************************************************************************/

	/**
	 * the function normalize the Cartesian type of basis sets
	 * for the given l,m,n basis set
	 * \param  l,m,n  : the given basis set angular momentum
	 * \param  nPrim  : number of primitives for c and e
	 * \param  c:       coefficient array for this shell
	 * \param  e:       exponent array for this shell
	 * \return the normalized factor 
	 *
	 */
	Double normCartBasisSets(const UInt& l, const UInt& m, const UInt& n,
			const UInt& nPrim, const Double* c, const Double* e);

	/**
	 * the function normalize the spherical type of basis sets
	 * \param  iShell: shell type
	 * \param  nPrim : number of primitives for c and e
	 * \param  c:      coefficient array for this shell
	 * \param  e:      exponent array for this shell
	 * \return the normalized factor for the given shell
	 *
	 */
	Double normPureBasisSets(const UInt& iShell,const UInt& nPrim, 
			const Double* c, const Double* e);

	/**
	 * this function is to get transformation matrix to transform
	 * the normalized Cartesian basis sets into normalized spherial type 
	 * basis sets
	 * cartesian dimension is normalized with (lx=L,ly=0,lz=0) 
	 * \param L angular momentum, up to 10
	 * \return M: the transformation matrix, in dimension (nCart,nPure)
	 */
	void getC2PFromL00(const UInt& ang, Double* M); 

	/**
	 * this function is to get transformation matrix to transform
	 * the normalized Cartesian basis sets into normalized spherial type 
	 * basis sets
	 * cartesian dimension is normalized with (lx,ly,lz)
	 * \param L angular momentum, up to 10
	 * \return M: the transformation matrix, in dimension (nCart,nPure)
	 */
	void getC2PFromLxLyLz(const UInt& ang, Double* M); 

	/**
	 * this function is to get transformation matrix to transform
	 * the arbitrary Cartesian basis sets into arbitrary spherial type 
	 * basis sets, no normalization factors appended
	 * \param L angular momentum, up to 10
	 * \return M: the transformation matrix, in dimension (nCart,nPure)
	 */
	void getC2PFromUnnormalized(const UInt& ang, Double* M); 

	/**
	 * this function is to get transformation matrix to transform
	 * the normalized spherical type basis sets into normalized 
	 * Cartesian type basis sets data
	 * cartesian dimension is normalized with (lx=L,ly=0,lz=0) 
	 * \param L angular momentum, up to 10
	 * \return M: the transformation matrix, in dimension (nPure,nCart)
	 */
	void getP2CFromL00(const UInt& ang, Double* M); 

	/**
	 * this function is to get transformation matrix to transform
	 * the normalized spherical type basis sets into normalized 
	 * Cartesian type basis sets data
	 * cartesian dimension is normalized with (lx,ly,lz) 
	 * \param L angular momentum, up to 10
	 * \return M: the transformation matrix, in dimension (nPure,nCart)
	 */
	void getP2CFromLxLyLz(const UInt& ang, Double* M); 

	/**
	 * this function is used to get the ratio between the normalization
	 * factor of N_{L00}(lx=L,ly=0,lz=0) and N_{lmn}(lx=l,ly=m,lz=n)
	 * basically, here what we compute is the ratio of N_{lmn}/N_{L00}
	 * N is a vector
	 */
	void getNormScaleVector(const UInt& ang, Double* M); 
}

#endif

