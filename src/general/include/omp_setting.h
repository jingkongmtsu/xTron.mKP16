/**
 * \file   omp_setting.h
 * \author Fenglai Liu and Jing Kong
 * setting function for OpenMP
 */
#ifndef OMP_SETTING_H
#define OMP_SETTING_H
#include<cstdio>
#include<boost/thread.hpp>
#include "libgen.h"
#include "omp.h"

namespace omp_setting {

	///
	/// testing that wether omp is doable 
	/// for the following section?
	/// sections that open mp could be appiled:
	/// mkl math library;
	///
	inline bool omp_doable() {
#ifdef USE_MKL
		return true;
#else
		return false;
#endif
	};

	///
	/// initilize the omp env
	///
	inline void omp_init() {

		// test that if omp is doable?
		if (! omp_doable()) return;

		// get the number of possible threads
		// we assume that each processor take on thread
		Int nProcs = boost::thread::hardware_concurrency();
#ifdef OMP_DEBUG
		printf("in the omp_init, number of procs: %d\n", nProcs);
#endif

		// now set it to parallel
		omp_set_num_threads(nProcs);
	};

	///
	/// turn off the omp and set the number of threads 
	/// back to 1
	///
	inline void omp_turnoff() {

		// test that if omp is doable?
		if (! omp_doable()) return;

		// 
		// we disable the omp by setting the nProcs
		// to be 1 
		//
		omp_set_num_threads(1);
	};

}

#endif

