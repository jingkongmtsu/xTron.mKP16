/**
 * cpp file related to diiscontroller.h
 * \author fenglai liu 
 */
#include<iostream>
#include<algorithm> 
#include<cmath> 
#include "globalinfor.h"
#include "textread.h"
#include "parameterparsing.h"
#include "spinmatrix.h"
#include "scfeadiis.h"
#include "scfconv.h"
#include "scfparam.h"
#include "diiscontroller.h"
using namespace globalinfor;
using namespace textread;
using namespace parameterparsing;
using namespace spinmatrix;
using namespace scfeadiis;
using namespace scfdiis;
using namespace scfconv;
using namespace scfparam;
using namespace std;

DIISController::DIISController(const GlobalInfor& globInfor, 
		const SCFParam& par):printSCFEADIIS(false),printSCFDIIS(false), 
	printController(0),nOscillationCycles(5),oscillationThresh(1.0E-4),
	oscillationErrorThresh(1.0E-3),maxSCFCycles(par.getMaxSCFCycles()),
	currentJob(SCF_NULL_JOB),maxPrecEADIIS(3),useExpensiveEADIISSolver(false),
	maxEADIISSpace(4),maxDIISSpace(10),DIISSpaceIncNumber(5),EADIISSpaceIncNumber(1),
	scfVecSpaceChoice(SCF_ENERGY),diisError(ZERO),
	indexArray(maxSCFCycles),historyIndexArray(maxSCFCycles*(maxSCFCycles+1)/2),
	historyIndexRecord(2*maxSCFCycles),coefs(maxSCFCycles),scfIndexSpace(maxSCFCycles),
	scfDIIS(par)
{
	// now let's clear the content of the array
	indexArray.clear();
	historyIndexArray.clear();
	historyIndexRecord.clear();
	scfIndexSpace.clear();
	coefs.clear();

	// now let's read in the information
	string input = globInfor.getInputFile();
	UInt section = par.getSec();
	ParameterParsing pp(input,"scfconv",section);
	if (pp.hasAnyParameters()) {

		// define the number of SCF cycles for checking oscillation status
		string key = "number_scf_for_checking_oscillation";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "number_scf_for_checking_oscillation can not be processed, not an integer";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp<=2 || tmp>maxSCFCycles) {
				string infor = "invalid number_scf_for_checking_oscillation value, " 
					"which is <=2 or > max scf cycles";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			nOscillationCycles = tmp;
		}

		// define the energy difference for checking wether it's in oscillation status
		key = "energy_threshold_for_checking_oscillation";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "energy_threshold_for_checking_oscillation can not be processed, not double";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp>ONE) {
				string infor = "invalid energy_threshold_for_checking_oscillation value, " 
					"can not be larger than one";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			oscillationThresh = tmp;
		}

		// define the error creteria value for checking wether it's in oscillation status, not converged
		key = "error_threshold_for_checking_oscillation";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "error_threshold_for_checking_oscillation can not be processed, not double";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp>ONE) {
				string infor = "invalid error_threshold_for_checking_oscillation value, " 
					"can not be larger than one";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			oscillationErrorThresh = tmp;
		}

		// maximum dimension space for DIIS
		key = "max_space_diis";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "max_space_diis can not be processed, not an integer";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp<=1 || tmp>maxSCFCycles) {
				string infor = "invalid max_space_diis value, which is <=1 or > max scf cycles";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			maxDIISSpace = tmp;
		}

		// incremental space value if you want to enlarge DIIS space
		key = "diis_space_inc_number";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "diis_space_inc_number can not be processed, not an integer";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (maxDIISSpace+tmp>maxSCFCycles) {
				string infor = "invalid diis_space_inc_number value, the result DIIS space > max scf cycles";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			DIISSpaceIncNumber = tmp;
		}

		// do we go with expensive EADIIS solver? It uses the same algorithm, but it garantees
		// the precision to third digit
		key = "use_expensive_eadiis_solver";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				useExpensiveEADIISSolver = true;
			}else if (value == "FALSE" || value == "F") {
				useExpensiveEADIISSolver = false;
			}else{
				string infor = "only true or false allowed to process use_expensive_eadiis_solver";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// return the resilt precision for EDIIS/ADIIS procedure
		key = "result_precision_eadiis";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "result_precision_eadiis can not be processed, not an integer";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			maxPrecEADIIS = tmp;
		}

		// maximum dimension space for EDIIS/ADIIS
		key = "max_space_eadiis";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "max_space_eadiis can not be processed, not an integer";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp<=1) {
				string infor = "invalid max_space_eadiis value, which is <=1";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}else if (tmp > SCF_EADIIS_SPACE_LIMIT) {
				cout << "SCF_EADIIS_SPACE_LIMIT is " << SCF_EADIIS_SPACE_LIMIT << endl;
				string infor = "too large max_space_eadiis value, the space dimension"
					"should be < SCF_EADIIS_SPACE_LIMIT; else EADIIS calculation will be very expensive";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);

			}
			maxEADIISSpace = tmp;
		}

		// incremental space value if you want to enlarge EDIIS/ADIIS space
		key = "eadiis_space_inc_number";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "eadiis_space_inc_number can not be processed, not an integer";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (maxEADIISSpace+tmp>SCF_EADIIS_SPACE_LIMIT) {
				cout << "SCF_EADIIS_SPACE_LIMIT is " << SCF_EADIIS_SPACE_LIMIT << endl;
				string infor = "too large eadiis_space_inc_number value, the result enlarged space dimension"
					"should be < SCF_EADIIS_SPACE_LIMIT; else EADIIS calculation will be very expensive";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			EADIISSpaceIncNumber = tmp;
		}

		// for forming the scf vector space, do we use diis error or the energy?
		key = "scf_vector_space_choice";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			if (value == "DIIS_ERROR") {
				scfVecSpaceChoice = SCF_DIIS_ERROR;
			}else if (value == "SCF_ENERGY") {
				scfVecSpaceChoice = SCF_ENERGY;
			}else{
				string infor = "invalid scf_vector_space_choice, should be either DIIS_ERROR or SCF_ENERGY";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		//
		// debug choice section
		//

		// whether to do debug printing for eadiis?
		key = "scfeadiis";
		if (pp.hasKeyDefined(key)) {
			printSCFEADIIS = true; 
		}

		// whether to do debug printing for diis?
		key = "scfdiis";
		if (pp.hasKeyDefined(key)) {
			printSCFDIIS = true; 
		}

		// whether to do debug printing for this class?
		key = "diiscontroller";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "can not transfer the debug value for diis_controller to integer";
				Excep excep("DIISController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp == 1 || tmp == 2) {
				printController = tmp; 
			}else{
				cout << "diis_controller debug options" << endl;
				cout << "0:  no print out for diis_controller class itself, default option" << endl;
				cout << "1:  call print() function to output normal printing" << endl;
				cout << "2:  debug printing to output more information" << endl;
				string infor = "to debug diis_controller the given value is invalid, must be 1 or 2";
				Excep excep("DIISController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}
	}
}

void DIISController::resetSCFIndexSpace(const SCFConv& scfconv)
{
	//
	// the general scheme is to take the lowest energy/diis error
	// from the history as the first element after resetting 
	// the whole scf index space.
	//
	// here we need to note that we do not take the first
	// SCF iteration result. This is because the first iteration
	// may have very low energy (for example, when using atomic guess);
	// which is far away from the real one
	//
	// we need to check the length of original data (diis error should
	// has same length with energy). We need to make sure that there's 
	// enough data to pick up
	//
	const DoubleVec& energyList = scfconv.getHistEnergyData();
	const DoubleVec& diisErrors = scfconv.getHistDIISErrorData();
	if (energyList.size() != diisErrors.size()) {
		string infor = "energy result length does not equal to diis error result length";
		Excep excep("DIISController","resetSCFIndexSpace",EXCEPTION_DIIS_CONTROLLER_ERROR,infor);
		handleExcep(excep);
	}
	if (energyList.size() <= 1) {
		string infor = "there's only 1 or no SCF energy and we can not reset scfIndexSpace vector";
		Excep excep("DIISController","resetSCFIndexSpace",EXCEPTION_DIIS_CONTROLLER_ERROR,infor);
		handleExcep(excep);
	}

	// now clear the whole content of scf index history
	scfIndexSpace.clear();
	
	// let's find out the lowest one in the history and
	// this will be the new starting point for 
	// restarting/switching algorithm
	UInt index  = -1;
	if (scfVecSpaceChoice == SCF_ENERGY) {
		Double minE = energyList[1];
		for(UInt i=1; i<energyList.size(); i++) {
			if (minE > energyList[i]) {
				minE  = energyList[i];
				index = i;
			}
		}
	}else{
		Double minError = diisErrors[1];
		for(UInt i=1; i<diisErrors.size(); i++) {
			if (minError > diisErrors[i]) {
				minError = diisErrors[i];
				index = i;
			}
		}
	}

	// now let's push back the index
	scfIndexSpace.push_back(index);
}

void DIISController::updateIndexSpace(const SCFConv& scfconv, UIntVec& scfIndexVec) const 
{
	//
	// here this function picks up the nth lowest energy/diis_error
	// so that to form the space to interpolate/extropolate the SCF
	// result Fock matrix.
	//
	// the pick up process is like this. All of results are chosen
	// from the scf index space. We will pick up the scf index
	// that has lowest energy or diis error.
	//
	
	// let's see how much current vectors we should have
	// usually it should increase 1 based on the last SCF iteration
	// that means, the current scf should have +1 vector space
	// than the last scf iteration
	UInt nSCFVec = indexArray.size() + 1;
	if (useADIIS(currentJob) || useEDIIS(currentJob)) {
		if (nSCFVec>maxEADIISSpace) nSCFVec = maxEADIISSpace;

		// check the expensive solver dimension
		// if we use expensive solver for EDIIS, make sure that it's less than 
		// the limit
		if(useExpensiveEADIISSolver && nSCFVec>SCF_EADIIS_SPACE_LIMIT_EXACT_SOLUTION) {
			nSCFVec = SCF_EADIIS_SPACE_LIMIT_EXACT_SOLUTION;
			cout << "limit value is " << SCF_EADIIS_SPACE_LIMIT_EXACT_SOLUTION << endl;
			string infor = "because you required to use the expensive solver for EDIIS/ADIIS, "
				"the scf vector space is reset to limit value";
			Excep excep("DIISController","checkSCFIndexSpace",EXCEPTION_SCFCONV_WARNING,infor);
		}else if (nSCFVec>SCF_EADIIS_SPACE_LIMIT){
			nSCFVec = SCF_EADIIS_SPACE_LIMIT;
			cout << "limit value is " << SCF_EADIIS_SPACE_LIMIT << endl;
			string infor = "for EDIIS/ADIIS there's a limit value for "
				"the number of scf index space involved, so we reset to use the limit value";
			Excep excep("DIISController","checkSCFIndexSpace",EXCEPTION_SCFCONV_WARNING,infor);
		}	

	}else if (useDIIS(currentJob)) {
		if (nSCFVec>maxDIISSpace) nSCFVec = maxDIISSpace;
	}else{
		string infor = "the current job is invalid";
		Excep excep("DIISController","updateIndexSpace",EXCEPTION_DIIS_CONTROLLER_ERROR,infor);
		handleExcep(excep);
	}

	// also check with the input list size
	// since valList records all of available SCF iteration
	// data, we should not have calculated nSCFVec
	// larger than that
	if (nSCFVec>scfIndexSpace.size()) {
		string infor = "the calculated nSCFVec > available scf index number, invalid situation";
		Excep excep("DIISController","updateIndexSpace",EXCEPTION_DIIS_CONTROLLER_ERROR,infor);
		handleExcep(excep);
	}

	// let's get the a copy of scf index space
	IntVec scfIndexList(scfIndexSpace.size());
	for(UInt i=0; i<scfIndexSpace.size(); i++) scfIndexList[i] = (Int)scfIndexSpace[i];

	// initialize the scf index vector
	// let's pick up the scf index corresponding to the 
	// lowest one in energy or diis error
	scfIndexVec.assign(nSCFVec,0);
	for(UInt i=0; i<nSCFVec; i++) {
		Double minVal    = LARGE_FLOAT_NUMBER;
		UInt minIndex    = 0;  // this is the index for scfIndexList corresponding to minVal
		UInt scfIndex    = 0;  // this is the result SCF index corresponding to minVal
		for(UInt iSCF=0; iSCF<scfIndexList.size(); iSCF++) {
			Int index = scfIndexList[iSCF];
			if (index < 0) continue;
			Double val = scfconv.getHistEnergy(index);
			if (scfVecSpaceChoice == SCF_DIIS_ERROR) val = scfconv.getHistDIISError(index);
			if (val<minVal) {
				minVal    = val;
				minIndex  = iSCF;
				scfIndex  = (UInt)index;
			}
		}

		// now let's push the value into the result
		scfIndexVec[i]         = scfIndex;
		scfIndexList[minIndex] = -1;
	}

	// finally, it's better that we sort the index order
	// we still follow the original order of SCF index
	std::stable_sort(scfIndexVec.begin(), scfIndexVec.end());

	// debug print out
	if (printController == 2) {
		cout << "the generated scf space vector" << endl;
		printf("%-5s  %-16s\n", "index", "scf index in vec");
		for(UInt i=0; i<nSCFVec; i++) {
			printf("%-5d  %-16d\n", (Int)i, (Int)scfIndexVec[i]);
		}
	}
}

UInt DIISController::checkSCFIndexSpace(const SCFConv& conv, const UIntVec& newSCFVecSpace) const
{

	// let's see what kind of problem we have for the new
	// scf index space
	
	//
	// possible problem: the new scf index array repeating
	// the old one in the archive
	//
	// because all of scf index are formed with the original order
	// thefore we are able to compare side by side
	//
	UInt nIters = historyIndexRecord.size()/2;
	for(UInt iIter=0; iIter<nIters; iIter++) {
		UInt iBegin = historyIndexRecord[2*iIter+0];
		UInt iLen   = historyIndexRecord[2*iIter+1];

		// whether the current index space has different dimension
		// with the given history?
		if (iLen != newSCFVecSpace.size()) continue;

		// now let's check every element
		bool same = true;
		for(UInt i=0; i<iLen; i++) {
			if (newSCFVecSpace[i] != historyIndexArray[iBegin+i]) {
				same = false;
				break;
			}
		}

		// if we have same result, we need to check it
		// raise up a warning sign in SCF
		if (same) {
			cout << "the generated space vector is same with " << iIter << " in history data" << endl;
			string infor = "the current SCF vector space for DIIS/EDIIS etc. "
				"is same with an old one in archive";
			Excep excep("DIISController","checkSCFIndexSpace",EXCEPTION_SCFCONV_WARNING,infor);
			handleExcep(excep);
			return DIIS_SPACE_REPEAT;
		}
	}

	//
	// possible problem: to check that whether SCF is in oscillation status
	//
	// we note that this problem maybe related to the problem that current
	// SCF index is not in the result space solution. It's better to check
	// this problem first, then to check that problem
	//
	bool checkOsc = true;
	if (checkOsc) {
		UInt currentSCFIndex = conv.getCurrentSCFIndex();
		bool hasOscillation = true;

		// we need to have last n consecutive SCF cycles for checking
		// also we only check the situation that current diis error is significantly large
		// we compare the history energy with the current energy to see the difference
		if (diisError>oscillationErrorThresh && currentSCFIndex>nOscillationCycles) {
			Double E = conv.getCurrentEnergy();
			const DoubleVec& energyList = conv.getHistEnergyData();
			for(UInt i=currentSCFIndex-1; i>currentSCFIndex-nOscillationCycles; i--) {
				Double diff = fabs(E-energyList[i]);
				if (diff>oscillationThresh) {
					hasOscillation = false;
					break;
				}
			}
		}else{
			hasOscillation = false;
		}

		// if we really in the Oscillation status, then we need to report
		// this problem
		if (hasOscillation) {
			string infor = "SCF seems in the oscillation status while SCF is still far away from convergence";
			Excep excep("DIISController","checkSCFIndexSpace",EXCEPTION_SCFCONV_WARNING,infor);
			handleExcep(excep);
			return SCF_IN_OSCILLATION_STATUS;
		}
	}
	
	//
	// possible problem: the current SCF Index is not appearing in
	// the result space solution. Because the energy is picking up
	// from lowest to second-lowest and so on, so this implies that
	// the current SCF index has high energy/diis error than all 
	// of scf cycles in the space solution
	//
	// for the situation that we may restart the current index space,
	// we do not do the check for the first run. Since the first run
	// will pick up the scf index which gives the lowest energy/diis error
	//
	if (newSCFVecSpace.size()>1) {
		UInt currentSCFIndex = conv.getCurrentSCFIndex();
		bool hasCurrentSCFIndex = false;
		for(UInt i=0; i<newSCFVecSpace.size(); i++) {
			if (currentSCFIndex == newSCFVecSpace[i]) {
				hasCurrentSCFIndex = true;
				break;
			}
		}
		if (! hasCurrentSCFIndex) {
			string infor = "the current SCF index is not in the derived scf vector space";
			Excep excep("DIISController","checkSCFIndexSpace",EXCEPTION_SCFCONV_WARNING,infor);
			handleExcep(excep);
			return CURRENT_SCF_INDEX_NOT_IN_SPACE;
		}
	}
	
	// finally let's return
	return SCFCONV_NORMAL_STATUS;
}

void DIISController::updateSpaceArchive(const UIntVec& newSCFVecSpace)
{
	// firstly update the record
	UInt beginPos = historyIndexArray.size();
	UInt len = newSCFVecSpace.size();
	historyIndexRecord.push_back(beginPos);
	historyIndexRecord.push_back(len);

	// now let's add in data
	for(UInt i=0; i<newSCFVecSpace.size(); i++) {
		historyIndexArray.push_back(newSCFVecSpace[i]);
	}

	// debug print out
	if (printController == 2) {
		UInt nIters = historyIndexRecord.size()/2;
		cout << "current history archive of scf vector space(including current one)" << endl;
		cout << "there are totally " << nIters << " vector space data" << endl;

		// print out the data
		for(UInt iIter=0; iIter<nIters; iIter++) {
			cout << endl;
			cout << "*************" << endl;
			cout << "the " << iIter << " vector space data" << endl;
			cout << "*************" << endl;
			UInt iBegin = historyIndexRecord[2*iIter+0];
			UInt iLen   = historyIndexRecord[2*iIter+1];
			printf("%-5s  %-16s\n", "index", "scf index in vec");
			for(UInt i=0; i<iLen; i++) {
				printf("%-5d  %-16d\n", (Int)i, (Int)historyIndexArray[iBegin+i]);
			}
		}
	}
}

void DIISController::updateDIISInfor(const DenMtrx& denMtrx, 
		const OneEMtrx& oneEMtrx, const SpinMatrix& fock)
{
	scfDIIS.updateDIISInfor(oneEMtrx,denMtrx,fock);
	diisError = scfDIIS.diisError();
}

UInt DIISController::derivIndexSpace(const SCFConv& conv, bool saveResult) 
{
	// deriv the scf index space 
	UIntVec scfIndexVec(maxSCFCycles);
	updateIndexSpace(conv,scfIndexVec);
	UInt status = checkSCFIndexSpace(conv,scfIndexVec);

	// if everything is good, we will write it into the history
	if (status == SCFCONV_NORMAL_STATUS || saveResult) {
		updateSpaceArchive(scfIndexVec);
		indexArray = scfIndexVec;
	}

	// now return
	return status;
}

void DIISController::incrementSCFIndexSpace(const SCFConv& conv)
{
	// now let's take care of the scf index space
	// we should push the current scf index into the space now
	UInt currentSCFIndex = conv.getCurrentSCFIndex();
	bool hasIt = false;
	for(UInt i=0; i<scfIndexSpace.size(); i++) {
		if (scfIndexSpace[i] == currentSCFIndex) {
			hasIt = true;
			break;
		}
	}
	if(hasIt) {
		string info = "why the scf index space already has current SCF index?";
		Excep excep("DIISController","incrementSCFIndexSpace",EXCEPTION_DIIS_CONTROLLER_ERROR,info);
		handleExcep(excep);
	}
	scfIndexSpace.push_back(currentSCFIndex);
}

void DIISController::convSCF(const DenMtrx& denMtrx, const SCFConv& conv, SpinMatrix& fock)
{
	// now let's derive the coefficients based on the 
	// index array
	coefs.assign(getNTerms(),ZERO);
	if (useEDIIS(currentJob) || useADIIS(currentJob)) {
		SCFEADIIS eadiis(*this);
		eadiis.minimization(*this,conv,denMtrx,fock);
		eadiis.updateCoefs(coefs);
		if (printSCFEADIIS) {
			eadiis.print();
		}
	}else if (useDIIS(currentJob)) {
		scfDIIS.updateCoefs(indexArray,coefs);
		if (printSCFDIIS) {
			scfDIIS.print();
		}
	}else{
		string infor = "the current job is not known to do SCF convergence job";
		Excep excep("DIISController","convSCF",EXCEPTION_DIIS_CONTROLLER_ERROR,infor);
		handleExcep(excep);
	}

	// double check the length of coefficients and index array
	if (indexArray.size() != coefs.size()) {
		string info = "the input scf index array size does not equal to coefficients array size "; 
		Excep excep("DIISController","convSCF",EXCEPTION_DIIS_CONTROLLER_ERROR,info);
		handleExcep(excep);
	}

	// now retrieve data from history to interpolate/extropolate
	// Fock matrix
	const HistDataMan& FockMtrxList = conv.getHistFockData();
	for(UInt iSpin=0; iSpin<fock.getNSpin(); iSpin++) {

		// get the Fock matrix, and clean it
		Mtrx& F  = fock.getMtrx(iSpin);
		F.set(ZERO);

		// set up an empty matrix for tmp Fock
		UInt nRow = F.getRow();
		UInt nCol = F.getCol();
		Mtrx f(nRow,nCol);

		// now get the history Fock matrix
		// the real history data, fockMtrxList is same or
		// longer than the length of data we want to do
		// extrapolation/interpolation
		for(UInt iFock=0; iFock<indexArray.size(); iFock++) {

			// get the index corresponding to iFock
			UInt index = indexArray[iFock];

			// dimension information
			// we also need to check it
			UInt rowDim = FockMtrxList.rowDim(index,iSpin);
			UInt colDim = FockMtrxList.colDim(index,iSpin);
			if (rowDim != nRow || colDim != nCol) {
				string info1 = "the history fock matrix index is: " + boost::lexical_cast<string>(index) + " ";
				string info2 = " result fock dimension is: " + boost::lexical_cast<string>(nRow) + " " + 
					boost::lexical_cast<string>(nCol);
				string info3 = " the given history fock in dimension of: " + 
					boost::lexical_cast<string>(rowDim) + " " + boost::lexical_cast<string>(colDim);
				string infor = info1 + info2 + info3;
				Excep excep("DIISController","convSCF",EXCEPTION_HISTORY_FOCK_DIMENSION_ERROR,infor);
				handleExcep(excep);
			}

			// load in data
			f.set(ZERO);
			FockMtrxList.retrieveData(f,index,iSpin);

			// scale Fock and add it into the result
			Double c = coefs[iFock];
			f.scale(c);
			F.add(f);
		}
	}

	// finally, if it's required to print out we just print out
	if (printController>0) {
		print();
	}
}

void DIISController::print() const
{
	cout << "***************************************" << endl;
	cout << "debug information for DIISController   " << endl;
	cout << "***************************************" << endl;
	cout << "current DIIS like algorithm    : " << currentJob << endl;
	cout << "result precision digit EADIIS  : " << maxPrecEADIIS << endl;
	cout << "current DIIS error             : " << diisError << endl;
	cout << "current DIIS space limit       : " << maxDIISSpace << endl;
	cout << "current EADIIS space limit     : " << maxEADIISSpace << endl;
	cout << "increment for DIIS space       : " << DIISSpaceIncNumber << endl;
	cout << "increment for EADIIS space     : " << EADIISSpaceIncNumber << endl;
	if (scfVecSpaceChoice == SCF_ENERGY) {
		cout << "SCF vector space choice        : scf energy" << endl;
	}else{
		cout << "SCF vector space choice        : scf diis error " << endl;
	}

	// now print out scf index space
	cout << "print out the scf index space: the result index array will be picked up from this space " << endl;
	printf("%-12s  %-12s\n", "index", "scf index");
	for(UInt i=0; i<scfIndexSpace.size(); i++) {
		printf("%-12d  %-12d\n", (Int)i, (Int)scfIndexSpace[i]);
	}
	printf("\n");

	// the result index array and corresponding coefficients
	if (useADIIS(currentJob) || useEDIIS(currentJob)) {
		if (useExpensiveEADIISSolver) {
			cout << "the base EADIIS solution(first round) coef. is accurate as 0.001" << endl;
		}else{
			cout << "the base EADIIS solution(first round) coef. is accurate as 0.01" << endl;
		}
		cout << "for more information, please check the findCoefs() function in scfeadiis.cpp" << endl;
		cout << "EDIIS/ADIIS scf index vector space, and it's coeff.: " << endl;
	}else if (useDIIS(currentJob)) {
		cout << "DIIS scf index vector space, and it's coeff.: " << endl;
	}
	printf("%-12s  %-12s  %-12s\n", "index", "scf index", "coefs");
	for(UInt i=0; i<getNTerms(); i++) {
		printf("%-12d  %-12d  %-12.6f\n", (Int)i, (Int)indexArray[i], coefs[i]);
	}
	printf("\n");
}
