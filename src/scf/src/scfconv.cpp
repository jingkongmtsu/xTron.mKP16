/**
 * cpp file related to scfconv.h
 * \author fenglai liu 
 */
#include<iostream>
#include<algorithm> 
#include "globalinfor.h"
#include "xcfunc.h"
#include "textread.h"
#include "molecule.h"
#include "parameterparsing.h"
#include "denmtrx.h"
#include "spinmatrix.h"
#include "scfparam.h"
#include "fock.h"
#include "scfconv.h"
using namespace globalinfor;
using namespace xcfunc;
using namespace textread;
using namespace molecule;
using namespace parameterparsing;
using namespace denmtrx;
using namespace spinmatrix;
using namespace scfparam;
using namespace fock;
using namespace scfconv;
using namespace std;

SCFConv::SCFConv(const GlobalInfor& globInfor, const SCFParam& par, const string& algorithm):infor(globInfor),
	convCriteria(par.getConvCriteria()),maxSCFCycles(par.getMaxSCFCycles()),nSpin(par.getNSpin()),
	scfIndex(-1),jobType(SCF_DIIS),firstJob(SCF_DIIS),secondJob(SCF_DIIS),diisError(ZERO),
	secondJobSwitchOnLimit(0.01E0),algorithmReplaceDIIS(SCF_EDIIS),algorithmReplaceEDIIS(SCF_DIIS),
	algorithmReplaceADIIS(SCF_DIIS),scfLinearThresh(1.0E-10),
	energyList(par.getMaxSCFCycles()),diisErrors(par.getMaxSCFCycles()),
	scfStatus(par.getMaxSCFCycles()),scfSolution(par.getMaxSCFCycles()),jobInfor(par.getMaxSCFCycles()),
	FockMtrxList(infor,"scffocklist",par.getSec(),par.getNSpin(),par.useFile()),  
	denMtrxList(infor,"scfdenmtrxlist",par.getSec(),par.getNSpin(),par.useFile()),
	diisController(globInfor,par),scfConvSol(globInfor,par)
{
	// read in parameters
	string input = globInfor.getInputFile();
	UInt section = par.getSec();
	ParameterParsing pp(input,"scfconv",section);
	if (pp.hasAnyParameters()) {

		//
		// below it's job setting section
		//
		string key = algorithm;
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			LineParse l(value);
			if (l.getNPieces() == 1) {
				string tmpValue = l.findValue(0);
				w.capitalize(tmpValue);
				if (tmpValue == "EDIIS") {
					firstJob  = SCF_EDIIS;
					secondJob = SCF_EDIIS;
				}else if (tmpValue == "ADIIS") {
					firstJob  = SCF_ADIIS;
					secondJob = SCF_ADIIS;
				}else if (tmpValue == "DIIS") {
					firstJob  = SCF_DIIS;
					secondJob = SCF_DIIS;
				}else if (tmpValue == "GDM") {
					firstJob  = SCF_GDM;
					secondJob = SCF_GDM;
				}else{
					cout << "the available scf algorithm is available in scfmacro.h" << endl;
					cout << "current invalid scf algorithm is " << tmpValue << endl;
					string infor = "the input scf algorithm method is invalid";
					Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
			}else if (l.getNPieces() == 2) {
				for(UInt i=0; i<2; i++) {
					string tmp = l.findValue(i);
					w.capitalize(tmp);
					UInt job = -1;
					if (tmp == "EDIIS") {
						job = SCF_EDIIS;
					}else if (tmp == "ADIIS") {
						job = SCF_ADIIS;
					}else if (tmp == "DIIS") {
						job = SCF_DIIS;
					}else{
						cout << "the available scf algorithm is available in scfmacro.h" << endl;
						cout << "current invalid scf algorithm is " << tmp << endl;
						string infor = "the input scf algorithm method is invalid";
						Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
						handleExcep(excep);
					}
					if (i == 0) {
						firstJob  = job;
					}else{
						secondJob = job;
					}
				}

			}else{
				string infor = "failed to resolve the input keywords for scf algorithms";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		//
		// for DIIS, if you want to switch to a new job, what job is that?
		//
		key = "new_algorithm_replace_diis";
		if (pp.hasKeyDefined(key)) {
			string value    = pp.getValue(key);
			UInt algorithm  = getSCFAlgorithm(value);
			if (algorithm == SCF_NULL_JOB) {
				string infor = "failed to resolve the method given in new_algorithm_replace_diis";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			algorithmReplaceDIIS = algorithm;
		}

		//
		// for EDIIS, if you want to switch to a new job, what job is that?
		//
		key = "new_algorithm_replace_ediis";
		if (pp.hasKeyDefined(key)) {
			string value    = pp.getValue(key);
			UInt algorithm  = getSCFAlgorithm(value);
			if (algorithm == SCF_NULL_JOB) {
				string infor = "failed to resolve the method given in new_algorithm_replace_ediis";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			algorithmReplaceEDIIS = algorithm;
		}

		//
		// for ADIIS, if you want to switch to a new job, what job is that?
		//
		key = "new_algorithm_replace_adiis";
		if (pp.hasKeyDefined(key)) {
			string value    = pp.getValue(key);
			UInt algorithm  = getSCFAlgorithm(value);
			if (algorithm == SCF_NULL_JOB) {
				string infor = "failed to resolve the method given in new_algorithm_replace_adiis";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			algorithmReplaceADIIS = algorithm;
		}

		// when we switch on the second job?
		key = "threshold_second_job_switch_on";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "invalid value of threshold_second_job_switch_on. not double.";
				Excep excep("SCFConv","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp<ZERO || tmp>ONE) {
				string infor = "invalid threshold_second_job_switch_on value, which is < 0 or > 1";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			secondJobSwitchOnLimit = tmp; 
		}
	}

	// reset job information
	// and all of data archive
	jobType = firstJob;
	energyList.clear();
	diisErrors.clear();
	scfStatus.clear();
	scfSolution.clear();
	jobInfor.clear();

	// for current job, if this is DIIS like job
	// we need to set the job for DIIS controller
	if (isDIISLikeJob(jobType)) {
		diisController.setJob(jobType);
	}
}

SCFConv::SCFConv(const GlobalInfor& globInfor, const SCFParam& par):infor(globInfor),
	convCriteria(par.getConvCriteria()),maxSCFCycles(par.getMaxSCFCycles()),nSpin(par.getNSpin()),
	scfIndex(-1),jobType(SCF_DIIS),firstJob(SCF_DIIS),secondJob(SCF_DIIS),diisError(ZERO),
	secondJobSwitchOnLimit(0.01E0),algorithmReplaceDIIS(SCF_EDIIS),algorithmReplaceEDIIS(SCF_DIIS),
	algorithmReplaceADIIS(SCF_DIIS),scfLinearThresh(1.0E-10),
	energyList(par.getMaxSCFCycles()),diisErrors(par.getMaxSCFCycles()),
	scfStatus(par.getMaxSCFCycles()),scfSolution(par.getMaxSCFCycles()),jobInfor(par.getMaxSCFCycles()),
	FockMtrxList(infor,"scffocklist",par.getSec(),par.getNSpin(),par.useFile()),  
	denMtrxList(infor,"scfdenmtrxlist",par.getSec(),par.getNSpin(),par.useFile()),
	diisController(globInfor,par),scfConvSol(globInfor,par)
{
	// read in parameters
	string input = globInfor.getInputFile();
	UInt section = par.getSec();
	ParameterParsing pp(input,"scfconv",section);
	if (pp.hasAnyParameters()) {

		//
		// below it's job setting section
		//
		string key = "scf_algorithm";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			LineParse l(value);
			if (l.getNPieces() == 1) {
				string tmpValue = l.findValue(0);
				w.capitalize(tmpValue);
				if (tmpValue == "EDIIS") {
					firstJob  = SCF_EDIIS;
					secondJob = SCF_EDIIS;
				}else if (tmpValue == "ADIIS") {
					firstJob  = SCF_ADIIS;
					secondJob = SCF_ADIIS;
				}else if (tmpValue == "DIIS") {
					firstJob  = SCF_DIIS;
					secondJob = SCF_DIIS;
				}else if (tmpValue == "GDM") {
					firstJob  = SCF_GDM;
					secondJob = SCF_GDM;
				}else{
					cout << "the available scf algorithm is available in scfmacro.h" << endl;
					cout << "current invalid scf algorithm is " << tmpValue << endl;
					string infor = "the input scf algorithm method is invalid";
					Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
			}else if (l.getNPieces() == 2) {
				for(UInt i=0; i<2; i++) {
					string tmp = l.findValue(i);
					w.capitalize(tmp);
					UInt job = -1;
					if (tmp == "EDIIS") {
						job = SCF_EDIIS;
					}else if (tmp == "ADIIS") {
						job = SCF_ADIIS;
					}else if (tmp == "DIIS") {
						job = SCF_DIIS;
					}else{
						cout << "the available scf algorithm is available in scfmacro.h" << endl;
						cout << "current invalid scf algorithm is " << tmp << endl;
						string infor = "the input scf algorithm method is invalid";
						Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
						handleExcep(excep);
					}
					if (i == 0) {
						firstJob  = job;
					}else{
						secondJob = job;
					}
				}

			}else{
				string infor = "failed to resolve the input keywords for scf algorithms";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		//
		// for DIIS, if you want to switch to a new job, what job is that?
		//
		key = "new_algorithm_replace_diis";
		if (pp.hasKeyDefined(key)) {
			string value    = pp.getValue(key);
			UInt algorithm  = getSCFAlgorithm(value);
			if (algorithm == SCF_NULL_JOB) {
				string infor = "failed to resolve the method given in new_algorithm_replace_diis";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			algorithmReplaceDIIS = algorithm;
		}

		//
		// for EDIIS, if you want to switch to a new job, what job is that?
		//
		key = "new_algorithm_replace_ediis";
		if (pp.hasKeyDefined(key)) {
			string value    = pp.getValue(key);
			UInt algorithm  = getSCFAlgorithm(value);
			if (algorithm == SCF_NULL_JOB) {
				string infor = "failed to resolve the method given in new_algorithm_replace_ediis";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			algorithmReplaceEDIIS = algorithm;
		}

		//
		// for ADIIS, if you want to switch to a new job, what job is that?
		//
		key = "new_algorithm_replace_adiis";
		if (pp.hasKeyDefined(key)) {
			string value    = pp.getValue(key);
			UInt algorithm  = getSCFAlgorithm(value);
			if (algorithm == SCF_NULL_JOB) {
				string infor = "failed to resolve the method given in new_algorithm_replace_adiis";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			algorithmReplaceADIIS = algorithm;
		}

		// when we switch on the second job?
		key = "threshold_second_job_switch_on";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "invalid value of threshold_second_job_switch_on. not double.";
				Excep excep("SCFConv","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp<ZERO || tmp>ONE) {
				string infor = "invalid threshold_second_job_switch_on value, which is < 0 or > 1";
				Excep excep("SCFConv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			secondJobSwitchOnLimit = tmp; 
		}
	}

	// reset job information
	// and all of data archive
	jobType = firstJob;
	energyList.clear();
	diisErrors.clear();
	scfStatus.clear();
	scfSolution.clear();
	jobInfor.clear();

	// for current job, if this is DIIS like job
	// we need to set the job for DIIS controller
	if (isDIISLikeJob(jobType)) {
		diisController.setJob(jobType);
	}
}

void SCFConv::getLowestEnergyInfor(Double& E, UInt& index) const
{
	E     = ZERO;
	index = 0;
	for(UInt i=0; i<energyList.size(); i++) {
		if (energyList[i]<E) {
			E = energyList[i];
			index = i;
		}
	}
}

void SCFConv::resetJob()
{
	UInt newJob = -1;
	if (useDIIS(jobType)) {
		newJob = algorithmReplaceDIIS;
	}else if (useEDIIS(jobType)) {
		newJob = algorithmReplaceEDIIS;
	}else if (useADIIS(jobType)) {
		newJob = algorithmReplaceADIIS;
	}else{
		string infor = "invalid scf convergence algorithm";
		Excep excep("SCFConv","resetJob",EXCEPTION_SCFCONV_ERROR,infor);
		handleExcep(excep);
	}

	// now set the job
	jobType = newJob;
}

bool SCFConv::checkUpSCFSolutions(UInt status, UInt solution) const
{
	// if the solution is restarting the index space,
	// we need to check whether it has been used before
	// since restarting index space always lead a "good" 
	// status after restarting, so just in case the 
	// solutions are not trapped in this point 
	if (solution == RESTART_SCF_INDEX_SPACE) {

		// check the status and solution array
		if (scfStatus.size() != scfSolution.size() || scfStatus.size() != jobInfor.size()) {
			string infor = "arrays scfStatus, scfSolution and jobInfor does not have same length";
			Excep excep("SCFConv","checkUpSCFSolutions",EXCEPTION_SCFCONV_ERROR,infor);
			handleExcep(excep);
		}

		// now let's check the history
		// do we have the same status in the history?
		UInt scfIndex = maxSCFCycles;
		for(UInt i=0; i<scfStatus.size(); i++) {
			UInt oldStatus   = scfStatus[i];
			UInt oldSolution = scfSolution[i];
			if (oldStatus == status && oldSolution == solution) {
				scfIndex = i;
				break;
			}
		}

		// whether they are corresponding to same algorithm?
		if (scfIndex != maxSCFCycles) {
			UInt algo = jobInfor[scfIndex];
			if (algo == jobType) {
				return false;
			}
		}

		// in default we return true
		return true;
	}

	// in default we return normal good checking
	return true;
}

bool SCFConv::deriveIndexSpace() 
{
	// let's try to derive the scf index array
	// also update the scf status
	UInt status   = diisController.derivIndexSpace(*this,false); 

	// now let's list possible problems
	UInt solution = NULL_SCF_SOLUTION; 
	if (status != SCFCONV_NORMAL_STATUS) {

		// we will try to solve the scf problem in a loop
		// until the status becomes normal
		// or we running out of the solutions
		while(true) {

			// get the current solution
			// if no solution found, we just break
			solution = scfConvSol.getCurrentSolution(status);
			if (solution == NULL_SCF_SOLUTION) {

				// print out the wanring sign
				cout << "current SCF problem is " << status << endl;
				string infor = "SCF convergence is not running in good status, but we ran out of solutions";
				Excep excep("SCFConv","deriveIndexSpace",EXCEPTION_SCFCONV_WARNING,infor);
				handleExcep(excep);

				// then let's finally derive the index array and save it
				UInt newStatus = diisController.derivIndexSpace(*this,true);
				status = newStatus;

				// now let's break
				break;
			}

			// check up the possible solution and status
			// if this is bad, we will try to bypass the 
			// following code and go to next round
			if (! checkUpSCFSolutions(status,solution)) {
				continue;
			}

			// now let's see how we can do
			if (solution == ENLARGE_SCF_INDEX_SPACE) {
				diisController.enlargeSCFIndexSpace();
			}else if (solution == RESTART_SCF_INDEX_SPACE) {
				diisController.resetIndexSpace(*this,solution);
			}else if (solution == SWITCH_TO_NEW_SCF_CONV_METHOD) {
				diisController.resetIndexSpace(*this,solution);
				resetJob();

				// if the new method is still DIIS like method,
				// then the diis controller should know it
				if (isDIISLikeJob(jobType)) {
					diisController.setJob(jobType);
				}
			}else {
				string infor = "the given solution for scf convergence problem is not known";
				Excep excep("SCFConv","deriveIndexSpace",EXCEPTION_SCFCONV_ERROR,infor);
				handleExcep(excep);
			}

			// now let's try to get the new index space
			// if everything good, we will break the loop
			UInt newStatus = diisController.derivIndexSpace(*this,false); 
			if (newStatus == SCFCONV_NORMAL_STATUS) {
				break;
			}else{
				status = newStatus;
			}
		}
	}

	// now let's push the status and solution to record
	scfStatus.push_back(status); 
	scfSolution.push_back(solution); 
}

void SCFConv::convSCF(const UInt& scfIndex0, const Double& E, const DenMtrx& denMtrx, 
		const OneEMtrx& oneEMtrx, SpinMatrix& fock)
{
	// update scf index 
	scfIndex = scfIndex0;

	// if no need to do scf convergence, we just return at the 
	// beginning
	if (noNeedSCFConv()) return;

	// updating to the data archive first
	energyList.push_back(E);
	bool isSymm = true;
	for(UInt iSpin=0; iSpin<fock.getNSpin(); iSpin++) {
		FockMtrxList.storeData(fock.getMtrx(iSpin),isSymm,iSpin);
		denMtrxList.storeData(denMtrx.getMtrx(iSpin),isSymm,iSpin);      
	}

	// update DIIS information
	// also the diis error
	diisController.updateDIISInfor(denMtrx,oneEMtrx,fock);
	diisError = diisController.DIISError();
	diisErrors.push_back(diisError);

	// does the SCF end, so that we do not 
	// need to do anything here?
	// but we need to update the scf status
	// so that to make the SCF print working correct
	if (firstSCFCycle() || scfEnds()) {
		scfStatus.push_back(SCFCONV_NORMAL_STATUS); 
		scfSolution.push_back(NULL_SCF_SOLUTION);
		jobInfor.push_back(jobType);
		return;
	}

	// update the scf index space to include current scf index
	diisController.incrementSCFIndexSpace(*this);

	// if we already have GDM runs
	// return
#ifdef WITH_GDM
	if (firstJob == SCF_GDM || secondJob == SCF_GDM) {
		return;
	}
#endif

	// now let's derive the scf index array for DIIS like procedure
	if (isDIISLikeJob(jobType)) {
		deriveIndexSpace(); 
	}
	
	// now let's derive the new fock matrix 
	if (isDIISLikeJob(jobType)) {
		diisController.convSCF(denMtrx,*this,fock);
	}else{
		string infor = "currently we do not have other non-DIIS like scf algorithms available";
		Excep excep("SCFConv","convSCF",EXCEPTION_SCFCONV_ERROR,infor);
		handleExcep(excep);
	}

	// push the current job into the jobInfor as a record
	jobInfor.push_back(jobType);

	// now let's see whether it's ready for second job switch on?
	if (diisError<secondJobSwitchOnLimit && firstJob != secondJob && jobType == firstJob) {
		jobType = secondJob;
		if (isDIISLikeJob(jobType)) {
			diisController.setJob(jobType);
		}
	}

	// reset the solution index
	scfConvSol.resetSolArrayIndex(); 
	if (scfConvSol.doPrint()) {
		scfConvSol.print();
	}
}

void SCFConv::print() const
{
	cout << "***************************************" << endl;
	cout << "debug information for SCF CONV        "  << endl;
	cout << "***************************************" << endl;
	cout << "standard convergence criteria: " << convCriteria << endl;
	cout << "maximum SCF cycles           : " << maxSCFCycles << endl;
	cout << "spin states                  : " << nSpin << endl;
	cout << "first  SCF algorithm         : " << getSCFAlgName(firstJob) << endl;
	cout << "second SCF algorithm         : " << getSCFAlgName(secondJob) << endl;
	cout << "current SCF algorithm        : " << getSCFAlgName(jobType) << endl;
	cout << "limit val. 2ed job switch on : " << secondJobSwitchOnLimit << endl;
	cout << "current SCF iteration index  : " << scfIndex << endl;
	cout << "DIIS error for current scf   : " << diisError << endl;

	// now let's print out the scf history that scfconv stores
	// we only do it when we need the SCF process
	if ((! noNeedSCFConv()) && (! firstSCFCycle())) {

		// check the dimension
		// all of them should equal with each other
		if (energyList.size() != diisErrors.size()) {
			string infor = "something wrong here, energyList size does not equal to diisErrors array";
			Excep excep("SCFConv","print",EXCEPTION_SCFCONV_ERROR,infor);
			handleExcep(excep);
		}
		if (energyList.size() != scfStatus.size()) {
			cout << "energy list size " << energyList.size() << endl;
			cout << "scf status  size " << scfStatus.size() << endl;
			string infor = "something wrong here, energyList size does not equal to scfStatus array";
			Excep excep("SCFConv","print",EXCEPTION_SCFCONV_ERROR,infor);
			handleExcep(excep);
		}
		if (energyList.size() != scfSolution.size()) {
			cout << "energy list size " << energyList.size() << endl;
			cout << "scfSolution size " << scfSolution.size() << endl;
			string infor = "something wrong here, energyList size does not equal to scfSolution array";
			Excep excep("SCFConv","print",EXCEPTION_SCFCONV_ERROR,infor);
			handleExcep(excep);
		}
		if (energyList.size() != jobInfor.size()) {
			cout << "energy list size " << energyList.size() << endl;
			cout << "jobInfor    size " << jobInfor.size() << endl;
			string infor = "something wrong here, energyList size does not equal to jobInfor array";
			Excep excep("SCFConv","print",EXCEPTION_SCFCONV_ERROR,infor);
			handleExcep(excep);
		}

		// now let's print out
		cout << "SCF history iteration status data " << endl;
		printf("%-4s %-16s %-16s %-10s %-35s %-35s\n", "iter", "energy", "DIIS error", "algorithm", 
				"status", "solution");
		string status;
		string solution;
		for(UInt i=0; i<energyList.size(); i++) {
			getProblemName(scfStatus[i],status);
			getSolutionName(scfSolution[i],solution);
			string algo = getSCFAlgName(jobInfor[i]);
			printf("%-4d %-16.10f %-16.10f %-10s %-35s %-35s\n", (Int)i, energyList[i], 
					diisErrors[i], algo.c_str(),status.c_str(), solution.c_str());
		}
	}
}
