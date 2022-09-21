/**
 * \file   integraljobs.h
 * \brief  define the analytical/numerical integral job name for the whole program
 * \author Fenglai Liu 
 */
#ifndef INTEGRALJOBS_H
#define INTEGRALJOBS_H

// begin the offload attribute region
//#pragma offload_attribute (push, target (mic))

#include <string>
#include "libgen.h"
using namespace std;

namespace integraljobs {

	//////////////////////////////////////////////////////////////////////
	//                  integral job name definition                    //
	//////////////////////////////////////////////////////////////////////

	//
	// here it defines all of possible analytical/numerical integral jobs
	// starting from 1
	//
	const UInt  NULL_INT_JOB        =  0;     ///< the null integral job

	//
	// here the first batch of jobs are defined for one/two body integrals
	// it ranges from 1 to 200
	//
	const UInt  TWO_BODY_OVERLAP    =  1;     ///< two body overlap integrals
	const UInt  KINETIC             =  2;     ///< kinetic integrals
	const UInt  NUCLEAR_ATTRACTION  =  3;     ///< nulcear attraction integrals
	const UInt  EX_DENSITY          =  4;     ///< exchange energy density 
	const UInt  ORTHOGONAL_MATRIX   =  5;     ///< orthogonal matrix (S^{-1/2}, S is overlap matrix)
	const UInt  MOM_P               =  6;     ///< (psi_{1}|x,y,z|psi_{2}) type of integrals, psi is a basis set function

	//
	// three body integral jobs name, it ranges from 201 to 400
	//
	const UInt  THREE_BODY_OVERLAP  =  201;   ///< three body overlap integrals
	const UInt  THREE_BODY_KINETIC  =  202;   ///< three body kinetic integrals 

	//
	// four body integrals job name, it ranges from 401 to 1999
	//
	const UInt  COULOMB             =  401;   ///< normal Coulomb integrals
	const UInt  EXCHANGE            =  402;   ///< normal Exchange integrals
	const UInt  COULOMB_EXCHANGE    =  403;   ///< normal Coulomb and Exchange integrals together

	//
	// general numerical integral job definition
	// 
	// here we do not have explicit functional information. The functional
	// information is defined in xcfunc class. Here for numerical integral
	// here we just give the very general category.
	//
	// the job ranges from 2000 to 3000
	//
	const UInt GROUND_STATE_DFT     = 2001;   ///< ground state DFT calculation
	const UInt TDDFT                = 2002;   ///< time dependent DFT calculation


	//////////////////////////////////////////////////////////////////////
	//            functions related to the integral job                 //
	//////////////////////////////////////////////////////////////////////

	/**
	 * whether the integral job is well defined?
	 */
	inline bool isNULLIntJob(UInt job) {
		if (job == NULL_INT_JOB) return true;
		return false;
	};

	/**
	 * whether it's DFT job?
	 */
	inline bool isDFTJob(UInt job) {
		if (job == GROUND_STATE_DFT || job == TDDFT) return true;
		return false;
	};

	/**
	 * for the given job, check it's integral order
	 */
	UInt getOperOrder(UInt job);

	/**
	 * checking that whether the job ID is valid
	 * and whether it's consistent with the general information 
	 * already inside
	 */
	bool checkIntJob(UInt intJob);

	/**
	 * this is mainly for debugging use
	 * it's to get the job string name with the input job code
	 */
	string getJobName(UInt intJob);

	/**
	 * this is to get the operator name for the setting file
	 * in the setting folder, basically it's related to the 
	 * memory/derivatives information of analytical integrals
	 */
	string getOperNameInSettingFile(UInt intJob);

	/**
	 * given by the integral job, do we need the molecule geometrical 
	 * data for calculation?
	 */
	inline bool needMolData(UInt intJob) {
		if (intJob == NUCLEAR_ATTRACTION) return true;
		return false;
	};  

	/**
	 * do we do coulomb integrals?
	 */
	inline bool doJ(UInt intJob)  { 
		if (intJob == COULOMB || intJob == COULOMB_EXCHANGE) return true;
		return false;
	};

	/**
	 * do we do exchange integrals?  
	 **/
	inline bool doK(UInt intJob)  { 
		if (intJob == EXCHANGE || intJob == COULOMB_EXCHANGE || intJob == EX_DENSITY) return true;
		return false;
	};

	/**
	 * is it the job only with exchange?
	 */
	inline bool onlyDoK(UInt intJob)  { 
		if (intJob == EXCHANGE || intJob == EX_DENSITY) return true;
		return false;
	};

	/**
	 * is it the job only with coulomb?
	 */
	inline bool onlyDoJ(UInt intJob)  { 
		if (intJob == COULOMB) return true;
		return false;
	};

	/**
	 * do we do coulomb and exchange integrals together?  
	 **/
	inline bool doJKTogether(UInt intJob)  { 
		if (intJob == COULOMB_EXCHANGE) return true;
		return false;
	};

	/**
	 * do we perform ERI calculation?
	 */
	inline bool doERI(UInt intJob)  {
		if (intJob == EXCHANGE || intJob == COULOMB || 
				intJob == COULOMB_EXCHANGE) return true;
		return false;
	};

	/**
	 * whether this is two body integral job: (a|O|b)
	 */
	inline bool twoBodyInts(UInt intJob)  {
		if (intJob == TWO_BODY_OVERLAP  || intJob == KINETIC || intJob == NUCLEAR_ATTRACTION ||
				intJob == EX_DENSITY || intJob == MOM_P) return true;
		return false;
	};

	/**
	 * whether this is three body integral job: (ab|O|c) 
	 */
	inline bool threeBodyInts(UInt intJob)  {
		if (intJob == THREE_BODY_OVERLAP || intJob == THREE_BODY_KINETIC) return true;
		return false;
	};

	/**
	 * whether this is four body integral job: (ab|O|cd)
	 */
	inline bool fourBodyInts(UInt intJob) {
		if (intJob == COULOMB || intJob == EXCHANGE ||
				intJob == COULOMB_EXCHANGE ) return true;
		return false;
	};

	/**
	 * does our integral work need to digest with density matrix?
	 *
	 * usually when the final result needs the digestion with 
	 * density matrix, and the same time the raw integrals itself
	 * are too big to keep on disk/memory, the raw integrals will 
	 * be digested directly with density matrix when it's formed.
	 * For example, the J/K jobs are just like this.
	 *
	 * Therefore, for such integral job to test the significance 
	 * of the shell quartet/integral, we need to consider the  
	 * density matrix information. This is what this function means.
	 *
	 * currently for all of four body integral jobs we need to do this
	 * you need to keep an eye on that
	 */
	inline bool digestWithDensityMatrix(UInt intJob) {
		if (fourBodyInts(intJob)) return true;
		return false;
	};

	///
	/// for numerical JK calculation, for debugging purpose we may need
	/// to calculate the energy density on the batch grids.
	///
	/// this function is to calculate the alpha exchange energy. We note,
	/// the exchange energy density (K variable) is always put in the head
	/// of the array
	///
	inline UInt getKVarPos(UInt intJob, UInt iSpin) {
		if (doK(intJob)) return iSpin;
		return -1;
	};

	///
	/// for numerical JK calculation, for debugging purpose we may need
	/// to calculate the energy density on the batch grids.
	///
	/// this function is to calculate the Coulomb energy(J variable). We note,
	/// it's always after the K variable
	///
	inline UInt getJVarPos(UInt intJob, UInt nSpin) {
		if (doJ(intJob)) {
			UInt pos = 0;
			if (doK(intJob)) {
				pos += 1;
				if (nSpin == 2) pos +=1;
			}
			return pos;
		}
		return -1;
	};
}

// end the offload attribute region
//#pragma offload_attribute (pop)

#endif

