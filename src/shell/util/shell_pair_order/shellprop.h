/**
 * \file    shellprop.h
 * \author  Fenglai Liu and Jing Kong
 *
 * The data below should be in consistent with the data defineed in 
 * shellprop.h in the shell include folder
 */
#ifndef SHELLPROP_H
#define SHELLPROP_H
#include "libgen.h"

namespace shellprop {

	/************************************************************************
	 *                          Data for Shell Structure                    *
	 ************************************************************************/

	/**
	 * the code for angular momentum numbers
	 * associated with the shell names
	 * We give it up to L = 20
	 */
	const Int SHELL_ANG_MOM_CODE[] = {
		0, 100, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
	};

	/**
	 * this value indicates the maximum number of shell names
	 */
	const Int MAX_SHELL_TYPES = 22;

	///
	/// here below we will define the unit to form
	/// the LCode
	///
	/// LCode is representing the code of angular
	/// momentum for a composite shell property
	/// e.g. shell pair etc.
	///
	/// see the function of codeL and decodL for 
	/// more information
	///
	const Int LCODE_UNIT_BRA2 = 1000;
	const Int ANG_UNIT = 100;

	/**
	 * code angular momentum for a shell 
	 * this should match the WORKING_SHELL_CODE array above
	 */
	inline Int codeL(Int lmin, Int lmax) {
		if (lmin == lmax) {
			return lmax;
		}else{
			return (lmin+lmax*ANG_UNIT);
		}
	};

	/**
	 * code the angular momentum for the shell pairs/two body integrals
	 * we note that shell i should >= shell j here
	 *	i corresponds to row shell, j corresponds to col shell, integral will be (i|j)
	 */
	inline Int codeL(Int iLmin, Int iLmax, Int jLmin, Int jLmax) {
		Int bra1Code = codeL(iLmin,iLmax);
		Int bra2Code = codeL(jLmin,jLmax);
		Int code     = bra1Code+bra2Code*LCODE_UNIT_BRA2;
		return code;
	};

	/**
	 * decode L from the angular momentum code above
	 */
	inline void decodeL(const Int& code, Int& lmin, Int& lmax) {
		if (code < ANG_UNIT) {
			lmin = code;
			lmax = code;
		}else{
			lmax = code/ANG_UNIT;
			lmin = code - lmax*ANG_UNIT;
		}
	};

	/**
	 * decode the angular momentum for the shell pairs
	 *	i corresponds to row shell, j corresponds to col shell
	 */
	inline void decodeL(const Int& code, Int& iLmin, Int& iLmax, Int& jLmin, Int& jLmax) {
		Int bra2Code = code/LCODE_UNIT_BRA2;
		Int bra1Code = code-bra2Code*LCODE_UNIT_BRA2;
		decodeL(bra1Code,iLmin,iLmax);
		decodeL(bra2Code,jLmin,jLmax);
	};

}

#endif

