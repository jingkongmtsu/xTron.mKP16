/**
 * \file    scfmacro.h
 * \brief   job definition for scf convergence procedure
 * \author  Fenglai Liu 
 */
#ifndef SCFMACRO_H
#define SCFMACRO_H
#include<string>
#include "libgen.h"

namespace scfmacro {

	///////////////////////////////////////////////
	// scf convergence algorithm defined here    //
	// range from 1 to 99                        //
	///////////////////////////////////////////////
	const UInt SCF_NULL_JOB  =  0;     ///< job is "null"
	const UInt SCF_EDIIS     =  1;     ///< EDIIS job
	const UInt SCF_ADIIS     =  2;     ///< ADIIS job
	const UInt SCF_DIIS      =  3;     ///< normal Pulay's DIIS job
	const UInt SCF_GDM       =  4;     ///< GDM job

	///////////////////////////////////////////////////
	// constant used by scf convergence algorithm    //
	///////////////////////////////////////////////////
	
	///
	/// for the algorithm to globally search the optimum coefficients used in
	/// EDIIS/ADIIS, we need to have a space limit to use the algorithm. Such
	/// limit is set to 6, because the exhaustive search used in the findCoefs()
	/// function takes step about 1.0^((n-1)*2); so in this case n = 6.
	/// however, we found that for n = 6 the running usually very slow; therefore
	/// we restrict the space limit to be 5
	///
	/// for the exhaustive search, we note that it's precision up to 0.01.
	///
	const UInt SCF_EADIIS_SPACE_LIMIT = 5;

	///
	/// for the algorithm to globally search the optimum coefficients used in
	/// EDIIS/ADIIS, we have an algorithm that is able to search the 
	/// the opt coefficients in precision of 0.001 in findCoefs() of scfeadiis.cpp. 
	/// this requires that the space is as large as 4, the steps took is 
	/// about 1.0^((n-1)*3). so n = 4 gives 1.0^9 steps.
	///
	const UInt SCF_EADIIS_SPACE_LIMIT_EXACT_SOLUTION = 4;

   ///////////////////////////////////////////////
	// define SCF convergence procedure problem  //
	// possibly encounted during SCF process     //
	// also define the corresponding solution    //
	// ranging from 100 to 199                   //
   ///////////////////////////////////////////////
	
	///
	/// in default, scf convergence everything good
	///
	const UInt SCFCONV_NORMAL_STATUS = 100;

	///
	/// in the DIIS controller class, we keep a history record
	/// of all previous index array data for the given job. Sometimes
	/// you may found that the derived new index array just repeats 
	/// the old one in the record. It indicates that the SCF just 
	/// stucks, so we need a solution on that 
	///
	const UInt DIIS_SPACE_REPEAT = 101; 

	///
	/// we may have the problem that the current SCF index is not 
	/// included in the result SCF index space. this can be a problem
	/// needed to be addressed 
	///
	const UInt CURRENT_SCF_INDEX_NOT_IN_SPACE = 102; 

	///
	/// when the SCF status is still far away from the convergence,
	/// it may be trapped into the oscillating status. That is,
	/// the SCF energy changes very little between the previous 
	/// SCF iterations and current SCF iteration, and also the DIIS error
	/// is still large
	///
	const UInt SCF_IN_OSCILLATION_STATUS = 103;

	//////////////////////////////////////////
	// solutions for above problems defined //
	// ranging from 201 to 300              //
	//////////////////////////////////////////
	
	///
	/// for each problem, let's solve the possible solution space
	/// the solution space means, if we meet the problem for the 
	/// 1st, 2ed, 3rd... times, what is the corresponding solution
	/// we have for them?
	///
	const UInt MAXSOL = 3;
	
	///
	/// let's define the NULL type solution
	/// 
	/// it means that there's no solution for the current problem,
	/// so we just stick to what we have.
	/// 
	const UInt NULL_SCF_SOLUTION  = 0;

	///
	/// to enlarge the SCF index space for DIIS/EDIIS/ADIIS etc.
	///
	const UInt ENLARGE_SCF_INDEX_SPACE  = 201; 

	///
	/// restart the scf index space for the given DIIS/EDIIS/ADIIS etc. method
	/// we will start from the scf index which gives the lowest energy in the 
	/// history 
	///
	const UInt RESTART_SCF_INDEX_SPACE  = 202; 

	///
	/// switch to a new SCF convergence method
	/// we will start from the scf index which gives the lowest energy in the 
	/// history 
	///
	const UInt SWITCH_TO_NEW_SCF_CONV_METHOD = 203; 

	///////////////////////////////////////////////
	// common SCF name macro                     //
	// range from 1000 to 2000                   //
	///////////////////////////////////////////////

	///
	/// define MACRO of DIIS_ERROR
	///
	const UInt SCF_DIIS_ERROR = 1000;
	
	///
	/// define MACRO of SCF energy
	///
	const UInt SCF_ENERGY = 1001;


	///////////////////////////////////////////////////
	// functions for job status check etc.           //
	///////////////////////////////////////////////////
	
	/**
	 * is the current job meaninful?
	 */
	inline bool isNullSCFConvJob(const UInt& jobType) {
		if (jobType == SCF_NULL_JOB) return true;
		return false;
	};

	/**
	 * do we use EDIIS?
	 */
	inline bool useEDIIS(const UInt& jobType) {
		if (jobType == SCF_EDIIS) return true;
		return false;
	};

	/**
	 * do we use ADIIS?
	 */
	inline bool useADIIS(const UInt& jobType) {
		if (jobType == SCF_ADIIS) return true;
		return false;
	};

	/**
	 * do we use DIIS?
	 */
	inline bool useDIIS(const UInt& jobType) {
		if (jobType == SCF_DIIS) return true;
		return false;
	};

	/**
	 * is this job the DIIS like procedure?
	 * including DIIS/EDIIS/ADIIS
	 */
	inline bool isDIISLikeJob(const UInt& jobType) {
		if (jobType == SCF_DIIS || jobType == SCF_EDIIS || jobType == SCF_ADIIS) return true;
		return false;
	};

	/**
	 * do we use GDM?
	 */
	inline bool useGDM(const UInt& jobType) {
		if (jobType == SCF_GDM) return true;
		return false;
	};

	/**
	 * get the corresponding scf algorithm by given the string name
	 */
	UInt getSCFAlgorithm(const std::string& algorithm);

	/**
	 * get the name of scf algorithm
	 */
	std::string getSCFAlgName(UInt job);

	/**
	 * for a given string of possible solution,
	 * let's return it's value 
	 */
	UInt getSolutionVal(const std::string& key);

	/**
	 * from the given key we will derive the corresponding problem string name
	 */
	void getProblemName(UInt key, std::string& name);

	/**
	 * from the given key we will derive the corresponding solution string name
	 */
	void getSolutionName(UInt key, std::string& name);
}

#endif

