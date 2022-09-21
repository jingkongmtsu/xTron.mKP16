/**
 * \file    angmomlist.h
 * \brief   this file describe the angular momentum data 
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef ANGMOMLIST_H
#define ANGMOMLIST_H
#include "libgen.h"

namespace shellprop {

	//
	// define the name of basis set orders
	//
	const UInt MAIN_BASIS_SET_ORDER   = 1;  ///< this is the global basis set indices array
	const UInt AIMPAC_BASIS_SET_ORDER = 2;  ///< corresponding to AIMPAC basis set indices array


	/**
	 * The GLOBAL_BASIS_SET_INDICES actually describe the basis set index in terms of
	 * the angular momentum used globally in the whole program. Since in the real calculation, 
	 * the spherical functions will be finally transformed into the linear combination of 
	 * Cartesian functions, so here we only need to list the Cartesian function.
	 * 
	 * This array is used to index the angular momentum directly, no need to calculate it
	 * anymore. so far we supports the L=9 so for I type orbital we can even do the second
	 * derivatives.
	 *
	 * This array returns the angular momentum for x, y and z (lmn) by given the global
	 * basis set index of k. For example, the 2,1,0 returns x^2y^1z^0
	 *
	 * the following data is in libint order
	 */
	const UInt MAX_ANG_GLOBAL_BASIS_SET_INDICES = 9;

	const UInt GLOBAL_BASIS_SET_INDICES[ ] = {
		0,  0,  0,  /*L = 0*/
		1,  0,  0,  /*L = 1*/
		0,  1,  0,
		0,  0,  1,
		2,  0,  0,  /*L = 2*/
		1,  1,  0,
		1,  0,  1,
		0,  2,  0,
		0,  1,  1,
		0,  0,  2,
		3,  0,  0,  /*L = 3*/
		2,  1,  0,
		2,  0,  1,
		1,  2,  0,
		1,  1,  1,
		1,  0,  2,
		0,  3,  0,
		0,  2,  1,
		0,  1,  2,
		0,  0,  3,
		4,  0,  0,  /*L = 4*/
		3,  1,  0,
		3,  0,  1,
		2,  2,  0,
		2,  1,  1,
		2,  0,  2,
		1,  3,  0,
		1,  2,  1,
		1,  1,  2,
		1,  0,  3,
		0,  4,  0,
		0,  3,  1,
		0,  2,  2,
		0,  1,  3,
		0,  0,  4,
		5,  0,  0,  /*L = 5*/
		4,  1,  0,
		4,  0,  1,
		3,  2,  0,
		3,  1,  1,
		3,  0,  2,
		2,  3,  0,
		2,  2,  1,
		2,  1,  2,
		2,  0,  3,
		1,  4,  0,
		1,  3,  1,
		1,  2,  2,
		1,  1,  3,
		1,  0,  4,
		0,  5,  0,
		0,  4,  1,
		0,  3,  2,
		0,  2,  3,
		0,  1,  4,
		0,  0,  5,
		6,  0,  0,  /*L = 6*/
		5,  1,  0,
		5,  0,  1,
		4,  2,  0,
		4,  1,  1,
		4,  0,  2,
		3,  3,  0,
		3,  2,  1,
		3,  1,  2,
		3,  0,  3,
		2,  4,  0,
		2,  3,  1,
		2,  2,  2,
		2,  1,  3,
		2,  0,  4,
		1,  5,  0,
		1,  4,  1,
		1,  3,  2,
		1,  2,  3,
		1,  1,  4,
		1,  0,  5,
		0,  6,  0,
		0,  5,  1,
		0,  4,  2,
		0,  3,  3,
		0,  2,  4,
		0,  1,  5,
		0,  0,  6,
		7,  0,  0,  /*L = 7*/
		6,  1,  0,
		6,  0,  1,
		5,  2,  0,
		5,  1,  1,
		5,  0,  2,
		4,  3,  0,
		4,  2,  1,
		4,  1,  2,
		4,  0,  3,
		3,  4,  0,
		3,  3,  1,
		3,  2,  2,
		3,  1,  3,
		3,  0,  4,
		2,  5,  0,
		2,  4,  1,
		2,  3,  2,
		2,  2,  3,
		2,  1,  4,
		2,  0,  5,
		1,  6,  0,
		1,  5,  1,
		1,  4,  2,
		1,  3,  3,
		1,  2,  4,
		1,  1,  5,
		1,  0,  6,
		0,  7,  0,
		0,  6,  1,
		0,  5,  2,
		0,  4,  3,
		0,  3,  4,
		0,  2,  5,
		0,  1,  6,
		0,  0,  7,
		8,  0,  0,  /*L = 8*/
		7,  1,  0,
		7,  0,  1,
		6,  2,  0,
		6,  1,  1,
		6,  0,  2,
		5,  3,  0,
		5,  2,  1,
		5,  1,  2,
		5,  0,  3,
		4,  4,  0,
		4,  3,  1,
		4,  2,  2,
		4,  1,  3,
		4,  0,  4,
		3,  5,  0,
		3,  4,  1,
		3,  3,  2,
		3,  2,  3,
		3,  1,  4,
		3,  0,  5,
		2,  6,  0,
		2,  5,  1,
		2,  4,  2,
		2,  3,  3,
		2,  2,  4,
		2,  1,  5,
		2,  0,  6,
		1,  7,  0,
		1,  6,  1,
		1,  5,  2,
		1,  4,  3,
		1,  3,  4,
		1,  2,  5,
		1,  1,  6,
		1,  0,  7,
		0,  8,  0,
		0,  7,  1,
		0,  6,  2,
		0,  5,  3,
		0,  4,  4,
		0,  3,  5,
		0,  2,  6,
		0,  1,  7,
		0,  0,  8,
		9,  0,  0,  /*L = 9*/
		8,  1,  0,
		8,  0,  1,
		7,  2,  0,
		7,  1,  1,
		7,  0,  2,
		6,  3,  0,
		6,  2,  1,
		6,  1,  2,
		6,  0,  3,
		5,  4,  0,
		5,  3,  1,
		5,  2,  2,
		5,  1,  3,
		5,  0,  4,
		4,  5,  0,
		4,  4,  1,
		4,  3,  2,
		4,  2,  3,
		4,  1,  4,
		4,  0,  5,
		3,  6,  0,
		3,  5,  1,
		3,  4,  2,
		3,  3,  3,
		3,  2,  4,
		3,  1,  5,
		3,  0,  6,
		2,  7,  0,
		2,  6,  1,
		2,  5,  2,
		2,  4,  3,
		2,  3,  4,
		2,  2,  5,
		2,  1,  6,
		2,  0,  7,
		1,  8,  0,
		1,  7,  1,
		1,  6,  2,
		1,  5,  3,
		1,  4,  4,
		1,  3,  5,
		1,  2,  6,
		1,  1,  7,
		1,  0,  8,
		0,  9,  0,
		0,  8,  1,
		0,  7,  2,
		0,  6,  3,
		0,  5,  4,
		0,  4,  5,
		0,  3,  6,
		0,  2,  7,
		0,  1,  8,
		0,  0,  9
	};

	/**
	 * The AIMPAC_BASIS_SET_INDICES describe the basis set index for the AIMPAC
	 * program. It's used to generate the wfn type of file.
	 */
	const UInt MAX_ANG_AIMPAC_BASIS_SET_INDICES = 4;

	const UInt AIMPAC_BASIS_SET_INDICES[ ] = {
		0,  0,  0,  /*L = 0*/
		1,  0,  0,  /*L = 1*/
		0,  1,  0,
		0,  0,  1,
		2,  0,  0,  /* DXX */
		0,  2,  0,  /* DYY */
		0,  0,  2,  /* DZZ */
		1,  1,  0,  /* DXY */
		1,  0,  1,  /* DXZ */
		0,  1,  1,  /* DYZ */
		3,  0,  0,  /* FXXX*/
		0,  3,  0,  /* FYYY*/
		0,  0,  3,  /* FZZZ*/
		2,  1,  0,  /* FXXY*/
		2,  0,  1,  /* FXXZ*/
		0,  2,  1,  /* FYYZ*/
		1,  2,  0,  /* FXYY*/
		1,  0,  2,  /* FXZZ*/
		0,  1,  2,  /* FYZZ*/
		1,  1,  1,  /* FXYZ*/
		4,	 0,  0,  /* GXXXX*/  
		0,  4,  0,	/* GYYYY*/ 
		0,  0,  4,	/* GZZZZ*/ 
		3,	 1,  0,  /* GXXXY*/ 
		3,  0,  1,  /* GXXXZ*/ 
		1,  3,  0,  /* GXYYY*/ 
		0,  3,  1,  /* GYYYZ*/ 
		1,  0,  3,  /* GXZZZ*/ 
		0,  1,  3,  /* GYZZZ*/ 
		2,  2,  0,  /* GXXYY*/ 
		2,  0,  2,  /* GXXZZ*/ 
		0,  2,  2,  /* GYYZZ*/ 
		2,  1,  1,  /* GXXYZ*/ 
		1,  2,  1,  /* GXYYZ*/ 
		1,  1,  2   /* GXYZZ*/ 
	};
}


#endif

