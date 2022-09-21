/**
 * cpp file related to scfintscontroller.h
 * \author fenglai liu 
 */
#include <cstdio>
#include "blas.h"
#include "parameterparsing.h"
#include "textread.h"
#include "globalinfor.h"
#include "molecule.h"
#include "spinmatrix.h"
#include "scfconv.h"
#include "scfparam.h"
#include "scfintscontroller.h"
using namespace blas;
using namespace parameterparsing;
using namespace textread;
using namespace globalinfor;
using namespace molecule;
using namespace spinmatrix;
using namespace scfconv;
using namespace scfparam;
using namespace scfintscontroller;

SCFIntsController::SCFIntsController(const SCFParam& param):scfOption(NO_CONTROLLER),
	initialCriteria(1.0E-2),csThreshRatio(1.0E-6),gints4DSPThreshRatio(1.0E-6),
	gints4DThreshRatio(1.0E-5),xcintsThreshRatio(1.0E-6),criteria(1.0E-2),
	inputGIntsInfor(param.getGIntsInfor()),workGIntsInfor(inputGIntsInfor),
	inputXCIntsInfor(param.getXCIntsInfor()),workXCIntsInfor(inputXCIntsInfor)
{
	const GlobalInfor& infor = param.getGlobalInfor();
	string input = infor.getInputFile();
	UInt section = param.getSec();
	ParameterParsing pp(input,"scfintscontroller",section);
	if (pp.hasAnyParameters()) {

		// initial criteria 
		string key = "initial_criteria";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "invalid initial_criteria value set";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			if (tmp<ZERO || tmp>ONE) {
				string infor = "invalid initial_criteria value, which is < 0 or > 1";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			initialCriteria = tmp;
		}

		// integral controller option
		key = "integral_controller_option";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "USER") {
				scfOption = AFTER_CONVERGE_WITH_USER_CHOICE; 
			}else if (value == "STEPWIZE") {
				scfOption = AFTER_CONVERGE_WITH_STEPWIZE_PROGRESS; 
			}else if (value == "NO_CONTROLLER") {
				scfOption = NO_CONTROLLER; 
			}else{
				string infor = "invalid option given in integral_controller_option";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// ratios set by the user
		key = "gints4d_thresh_ratio";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "invalid gints4d_thresh_ratio value set";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			if (tmp<ZERO || tmp>ONE) {
				string infor = "invalid gints4d_thresh_ratio value, which is < 0 or > 1";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			gints4DThreshRatio = tmp;
		}

		// ratios set by the user
		key = "gints4d_shell_pair_thresh_ratio";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "invalid shell_pair_thresh_ratio value set";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			if (tmp<ZERO || tmp>ONE) {
				string infor = "invalid gints4d_shell_pair_thresh_ratio value, which is < 0 or > 1";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			gints4DSPThreshRatio = tmp;
		}

		// ratios set by the user
		key = "xcints_thresh_ratio";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "invalid xcints_thresh_ratio value set";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			if (tmp<ZERO || tmp>ONE) {
				string infor = "invalid xcints_thresh_ratio value, which is < 0 or > 1";
				Excep excep("SCFIntsController","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			xcintsThreshRatio = tmp;
		}
	}

	// if the number of SCF = 1
	// then we do not do anything
	// and disable the future reset of 
	// integral threshold, too
	UInt nSCFCycles = param.getMaxSCFCycles();
	if (nSCFCycles == 1) {
		scfOption = NO_CONTROLLER; 
		return;
	}

	// finally, initilize the integral work center
	// according to the initial criteria
	bool withInit = true;
	setIntsInfor(withInit); 
}

void SCFIntsController::setIntsInfor(bool withInitCriteria) 
{
	//
	// now let's set the integral information for gints and xcints
	// 
	// our strategy is like this, the reset threshold value should 
	// not exceed the user set/default value(which should be accurate
	// enough).
	//
	Double c = criteria;
	if (withInitCriteria) c= initialCriteria;

	// initilize the workGIntsInfor
	if (withSCFIntsController()) {

		// setting initial value of shell pair
		Double thresh = c*gints4DSPThreshRatio;
		if (thresh>inputGIntsInfor.spGInts4DThreshVal()) {
			workGIntsInfor.changeThresh(GINTS4D_SHELL_PAIR_THRESH,thresh);
		}

		// setting initial value of gints4D
		thresh = c*gints4DThreshRatio;
		if (thresh>inputGIntsInfor.int4DThreshVal()) {
			workGIntsInfor.changeThresh(GINTS4D_THRESH,thresh);
		}
	}
}

