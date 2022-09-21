/**
 * cpp file associated with integraljobs.h
 * \author Fenglai Liu 
 */
#include<string>
#include<cstdio>
#include "excep.h"
#include "integraljobs.h"
using namespace std;
using namespace excep;
using namespace integraljobs;

// begin the offload attribute region
//#pragma offload_attribute (push, target (mic))

UInt integraljobs::getOperOrder(UInt job)
{
	// only the non-DFT jobs has operator order
	if (isDFTJob(job)) {
		string info = "DFT job does not have operator order";
		Excep excep("integraljobs","getOperOrder",EXCEPTION_INPUT_PARAMETER_INVALID,info);
		handleExcep(excep);
	}

	// now let's check the order
	UInt operOrder = -1;
	if (twoBodyInts(job)) {
		operOrder = 2;
	}else if (threeBodyInts(job)) {
		operOrder = 3;
	}else if (fourBodyInts(job)) {
		operOrder = 4;
	}
	return operOrder;
}

bool integraljobs::checkIntJob(UInt intJob) 
{
	// now check the job names
	switch(intJob) {

		// null integral job, it's still a
		// valid integral job
		case  NULL_INT_JOB:       
			return true;
			break;

		// two body integrals
		case  TWO_BODY_OVERLAP:   
			return true;
			break;
		case  KINETIC:            
			return true;
			break;
		case  NUCLEAR_ATTRACTION: 
			return true;
			break;
		case  EX_DENSITY:         
			return true;
			break;
		case  ORTHOGONAL_MATRIX:  
			return true;
			break;
		case  MOM_P:
			return true;
			break;


		// three body integrals 
		case  THREE_BODY_OVERLAP: 
			return true;
			break;
		case  THREE_BODY_KINETIC: 
			return true;
			break;

		// four body integrals 
		case  COULOMB:            
			return true;
			break;
		case  EXCHANGE:           
			return true;
			break;
		case  COULOMB_EXCHANGE:   
			return true;
			break;

		// DFT jobs
		case GROUND_STATE_DFT:
			return true;
			break;
		case TDDFT:
			return true;
			break;

		// default value
		default:
			break;
	}

	// in default we return false
	return false;
}

string integraljobs::getJobName(UInt intJob) 
{
	// now check the job names
	switch(intJob) {

		// null integral job, it's still a
		// valid integral job
		case  NULL_INT_JOB:       
			return "NULL_INTEGRAL_JOB";
			break;

		// two body integrals
		case  TWO_BODY_OVERLAP:   
			return "TWO_BODY_OVERLAP";
			break;
		case  KINETIC:            
			return "TWO_BODY_KINETIC";
			break;
		case  NUCLEAR_ATTRACTION: 
			return "NUCLEAR_ATTRACTION";
			break;
		case  EX_DENSITY:         
			return "EXCHANGE_ENERGY_DENSITY";
			break;
		case  ORTHOGONAL_MATRIX:  
			return "ORTHOGONAL_MATRIX";
			break;
		case  MOM_P:  
			return "MOMENTUM INTEGRAL MATRIX(DIPOLE)";
			break;

		// three body integrals 
		case  THREE_BODY_OVERLAP: 
			return "THREE_BODY_OVERLAP";
			break;
		case  THREE_BODY_KINETIC: 
			return "THREE_BODY_KINETIC";
			break;

		// four body integrals 
		case  COULOMB:            
			return "COULOMB_INTEGRALS";
			break;
		case  EXCHANGE:           
			return "EXCHANGE_INTEGRALS";
			break;
		case  COULOMB_EXCHANGE:   
			return "COULOMB_EXCHANGE_INTEGRALS";
			break;

		// DFT jobs
		case GROUND_STATE_DFT:
			return "GROUND_STATE_DFT";
			break;
		case TDDFT:
			return "TDDFT";
			break;

		// default value
		default:
			break;
	}

	// if we go here it means something wrong
	printf("the given integral job code is %d\n", (Int)intJob);
	string info = "invalid integral job name is given, failed to get it's string name";
	Excep excep("integraljobs","getJobName",EXCEPTION_INPUT_PARAMETER_INVALID,info);
	handleExcep(excep);

	// in default we return false
	return "INVALID_INTEGRAL_JOB_NAME";
}

string integraljobs::getOperNameInSettingFile(UInt oper) 
{
	// get the name used for setting file
	// basically, it's memory setting file 
	// and the derivative setting file
	//
	// currently the setting file is generated from cppints
	// who is generator for the integral codes
	//
	// the return name should be same in accordance with 
	// the setting file in setting folder. If error happens,
	// please check the files in setting folder
	string operName = "OPER_NAME_IS_INVALID";
	if (oper == TWO_BODY_OVERLAP){
		operName = "twobodyoverlap";
	}else if (oper == KINETIC) {
		operName = "kinetic";
	}else if (oper == NUCLEAR_ATTRACTION) {
		operName = "nai";
	}else if (oper == EX_DENSITY) {
		operName = "esp";
	}else if (oper == MOM_P) {
		operName = "mom_p";
	}else if (oper == COULOMB || oper == EXCHANGE || oper == COULOMB_EXCHANGE) {
		operName = "eri";
	}else{
		printf("passed in operator is(check the integraljobs.h for the macro): %d\n", (Int)oper);
		string infor = "failed to identify the operator name for the corresponding setting file";
		Excep excep("integraljobs","getOperNameInSettingFile",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}
	return operName;
}

// end the offload attribute region
//#pragma offload_attribute (pop)

