/**
 * cpp file related to scfmacro.h
 * \author fenglai liu 
 */
#include<iostream>
#include "excep.h"
#include "scfmacro.h"
using namespace excep;
using namespace scfmacro;
using namespace std;

UInt scfmacro::getSCFAlgorithm(const std::string& algorithm)
{
	if (algorithm == "EDIIS") {
		return SCF_EDIIS;
	}else if (algorithm == "ADIIS") {
		return SCF_ADIIS;
	}else if (algorithm == "DIIS") {
		return SCF_DIIS;
	}else if (algorithm == "GDM") {
		return SCF_GDM;
	}else{
		cout << "the input algorithm is not matching with what we currently have" << endl;
		cout << "the available scf algorithms are" << endl;
		cout << "DIIS  EDIIS  ADIIS" << endl;
	}

	// now let's return
	return SCF_NULL_JOB;
}

std::string scfmacro::getSCFAlgName(UInt algorithm)
{
	if (algorithm == SCF_EDIIS) {
		return "EDIIS";
	}else if (algorithm == SCF_ADIIS) {
		return "ADIIS";
	}else if (algorithm == SCF_DIIS) {
		return "DIIS";
	}else if (algorithm == SCF_GDM) {
		return "GDM";
	}else{
		cout << "the input algorithm is not matching with what we currently have" << endl;
		cout << "the available scf algorithms are" << endl;
		cout << "DIIS  EDIIS  ADIIS" << endl;
	}

	// now let's return
	return "SCF_JOB_INVALID";
}

UInt scfmacro::getSolutionVal(const std::string& key)
{
	// let's list all possible solutions
	if (key == "ENLARGE_SCF_INDEX_SPACE") {
		return ENLARGE_SCF_INDEX_SPACE;
	}else if (key == "RESTART_SCF_INDEX_SPACE") {
		return RESTART_SCF_INDEX_SPACE;
	}else if (key == "SWITCH_TO_NEW_SCF_CONV_METHOD") {
		return SWITCH_TO_NEW_SCF_CONV_METHOD;
	}else{

		// print out error information
		cout << "your input solution keyword for solving scf convergence problem is invalid" << endl;
		cout << "here are what's available" << endl;
		cout << "ENLARGE_SCF_INDEX_SPACE: enlarge the scf index space for DIIS/EDIIS etc. method" << endl;
		cout << "RESTART_SCF_INDEX_SPACE: restart the scf index space for DIIS/EDIIS etc. method" << endl;
		cout << "SWITCH_TO_NEW_SCF_CONV_METHOD: switch to a different scf congergence algorithm"  << endl;
	}

	// now let's return
	return NULL_SCF_SOLUTION;
}

void scfmacro::getProblemName(UInt key, std::string& name)
{
	// let's list all possible problem
	name = "NORMAL";
	if (key != SCFCONV_NORMAL_STATUS) {
		if (key == DIIS_SPACE_REPEAT) {
			name = "DIIS_SPACE_REPEAT";
		}else if (key == CURRENT_SCF_INDEX_NOT_IN_SPACE) {
			name = "CURRENT_SCF_INDEX_NOT_IN_SPACE";
		}else if (key == SCF_IN_OSCILLATION_STATUS) {
			name = "SCF_IN_OSCILLATION_STATUS";
		}else{
			std::string infor = "the input integer value is invalid";
			Excep excep("scfmacro","getProblemName",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}
	}
}

void scfmacro::getSolutionName(UInt key, std::string& name)
{
	// let's list all possible solutions
	name = "NO_SOLUTION_NEEDED";
	if (key != NULL_SCF_SOLUTION) {
		if (key == ENLARGE_SCF_INDEX_SPACE) {
			name = "ENLARGE_SCF_INDEX_SPACE";
		}else if (key == RESTART_SCF_INDEX_SPACE) {
			name = "RESTART_SCF_INDEX_SPACE";
		}else if (key == SWITCH_TO_NEW_SCF_CONV_METHOD) {
			name = "SWITCH_TO_NEW_SCF_CONV_METHOD";
		}else{
			std::string infor = "the input integer value is invalid";
			Excep excep("scfmacro","getSolutionName",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}
	}
}

