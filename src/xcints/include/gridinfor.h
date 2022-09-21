/**
 * \file    gridinfor.h
 * \brief   describing the constant grid information 
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef GRIDINFOR_H
#define GRIDINFOR_H
#include "libgen.h"

namespace gridinfor {

	/**
	 * valid number of angular points in Lebedev method
	 */
	const UInt VALID_LEBEDEV_GRID_NUMBER [ ] = {
		6,   14,  26,  38,  50,  74,   86,   110,  146,  170,  194,  230,  266,  302,
		350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 
		4334,4802,5294,5810
	};	

	/**
	 * totally there are 32 points allowed for Lebedev grids
	 */
	const UInt NUMBER_VALID_LEBEDEV = 32;

	// prune grid method list
	const UInt NON_PRUNE_GRID       = 0;
	const UInt SG1_GRID             = 1;
	const UInt BAKER_GRID           = 2;

	// choice of grid
	const UInt STANDARD_GRID        = 1;
	const UInt COARSE_GRID          = 2;
	const UInt FINE_GRID            = 3;
}

#endif
