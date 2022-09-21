/**
 * cpp file related to scfconvsol.h
 * \author fenglai liu 
 */
#include<cstdio>
#include<iostream>
#include<algorithm> 
#include "excep.h"
#include "globalinfor.h"
#include "textread.h"
#include "parameterparsing.h"
#include "scfparam.h"
#include "scfconvsol.h"
using namespace excep;
using namespace globalinfor;
using namespace textread;
using namespace parameterparsing;
using namespace scfparam;
using namespace scfconvsol;
using namespace std;

SCFConvSol::SCFConvSol(const GlobalInfor& globInfor, const SCFParam& par):printSCFConvSol(false),
	spaceVecRepeatIndex(0),currentScfIndexNotInsideIndex(0),scfOscillationIndex(0)
{
	//
	// set up default solution for each scf possible problem
	//
	
	// problems for space vector repeating
	for(UInt i=0; i<MAXSOL; i++) spaceVecRepeat[i] = NULL_SCF_SOLUTION;
	for(UInt i=0; i<MAXSOL; i++) {
		if (i == 0) spaceVecRepeat[i] = ENLARGE_SCF_INDEX_SPACE;
		if (i == 1) spaceVecRepeat[i] = RESTART_SCF_INDEX_SPACE;
		if (i == 2) spaceVecRepeat[i] = SWITCH_TO_NEW_SCF_CONV_METHOD;
	}

	// problem for current scf index is not in the result scf vector space
	for(UInt i=0; i<MAXSOL; i++) currentScfIndexNotInside[i] = NULL_SCF_SOLUTION;
	for(UInt i=0; i<MAXSOL; i++) {
		if (i == 0) currentScfIndexNotInside[i] = ENLARGE_SCF_INDEX_SPACE;
		if (i == 1) currentScfIndexNotInside[i] = RESTART_SCF_INDEX_SPACE;
		if (i == 2) currentScfIndexNotInside[i] = SWITCH_TO_NEW_SCF_CONV_METHOD;
	}
	
	// problem for solving scf oscillation problem
	for(UInt i=0; i<MAXSOL; i++) scfOscillation[i] = NULL_SCF_SOLUTION;
	for(UInt i=0; i<MAXSOL; i++) {
		if (i == 0) scfOscillation[i] = ENLARGE_SCF_INDEX_SPACE;
		if (i == 1) scfOscillation[i] = SWITCH_TO_NEW_SCF_CONV_METHOD;
	}
	
	// read in parameters
	string input = globInfor.getInputFile();
	UInt section = par.getSec();
	ParameterParsing pp(input,"scfconv",section);
	if (pp.hasAnyParameters()) {

		// if you have trouble for diis space repeating, how you handle it?
		string key = "diis_space_repeat_trouble";
		if (pp.hasKeyDefined(key)) {

			// now read in
			string value = pp.getValue(key);
			LineParse l(value);
			for(UInt i=0; i<MAXSOL; i++) spaceVecRepeat[i] = NULL_SCF_SOLUTION;
			for(UInt i=0; i<l.getNPieces(); i++) {
				string tmpValue = l.findValue(i);
				UInt solution = scfmacro::getSolutionVal(tmpValue);
				if (solution == NULL_SCF_SOLUTION) {
					string infor = "invalid solution for diis_space_repeat_trouble keyword";
					Excep excep("SCFConvSol","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
				if (i<MAXSOL) {
					spaceVecRepeat[i] = solution;
				}else{
					cout << "only maximum " << MAXSOL << "solutions allowed to solving "
						"the scf convergence problem" << endl;
					string infor = "too many solutions provided for diis_space_repeat_trouble problem";
					Excep excep("SCFConvSol","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
			}
		}

		// if you have trouble for energy going up, how you handle it?
		key = "current_scf_index_not_in_space";
		if (pp.hasKeyDefined(key)) {

			// now read in
			string value = pp.getValue(key);
			LineParse l(value);
			for(UInt i=0; i<MAXSOL; i++) currentScfIndexNotInside[i] = NULL_SCF_SOLUTION;
			for(UInt i=0; i<l.getNPieces(); i++) {
				string tmpValue = l.findValue(i);
				UInt solution = scfmacro::getSolutionVal(tmpValue);
				if (solution == NULL_SCF_SOLUTION) {
					string infor = "invalid solution for current_scf_index_not_in_space keyword";
					Excep excep("SCFConvSol","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
				if (i<MAXSOL) {
					currentScfIndexNotInside[i] = solution;
				}else{
					cout << "only maximum " << MAXSOL << "solutions allowed to solving "
						"the scf convergence problem" << endl;
					string infor = "too many solutions provided for current_scf_index_not_in_space problem";
					Excep excep("SCFConvSol","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
			}
		}

		// if you have trouble for SCF oscillation status, how you handle it?
		key = "scf_in_oscillation_status";
		if (pp.hasKeyDefined(key)) {

			// now read in
			string value = pp.getValue(key);
			LineParse l(value);
			for(UInt i=0; i<MAXSOL; i++) scfOscillation[i] = NULL_SCF_SOLUTION;
			for(UInt i=0; i<l.getNPieces(); i++) {
				string tmpValue = l.findValue(i);
				UInt solution = scfmacro::getSolutionVal(tmpValue);
				if (solution == NULL_SCF_SOLUTION) {
					string infor = "invalid solution for scf_in_oscillation_status keyword";
					Excep excep("SCFConvSol","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
				if (i<MAXSOL) {
					scfOscillation[i] = solution;
				}else{
					cout << "only maximum " << MAXSOL << "solutions allowed to solving "
						"the scf convergence problem" << endl;
					string infor = "too many solutions provided for scf_in_oscillation_status problem";
					Excep excep("SCFConvSol","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
			}
		}

		// whether we print out the scf convergence solution class?
		key = "scfconvsol";
		if (pp.hasKeyDefined(key)) {
			printSCFConvSol = true; 
		}
	}
}

UInt SCFConvSol::getCurrentSolution(UInt problem) 
{
	// let's see the problem
	// get the corresponding solution 
	// and increment the solution index
	UInt solution = -1;
	if (problem == DIIS_SPACE_REPEAT) {
		if (spaceVecRepeatIndex>=MAXSOL) {
			solution = NULL_SCF_SOLUTION;
		}else{
			solution = spaceVecRepeat[spaceVecRepeatIndex];
			spaceVecRepeatIndex++;  
		}
	}else if (problem == CURRENT_SCF_INDEX_NOT_IN_SPACE) {
		if (currentScfIndexNotInsideIndex>=MAXSOL) {
			solution = NULL_SCF_SOLUTION;
		}else{
			solution = currentScfIndexNotInside[currentScfIndexNotInsideIndex];
			currentScfIndexNotInsideIndex++;  
		}
	}else if (problem == SCF_IN_OSCILLATION_STATUS) {
		if (scfOscillationIndex>=MAXSOL) {
			solution = NULL_SCF_SOLUTION;
		}else{
			solution = scfOscillation[scfOscillationIndex];
			scfOscillationIndex++;  
		}
	}else {
		cout << "the input scf convergence status is " << problem << endl;
		string infor = "invalid scf convergence status provided here, check scfmacro.h for status definition";
		Excep excep("SCFConvSol","getCurrentSolution",EXCEPTION_SCFCONV_ERROR,infor);
		handleExcep(excep);
	}

	// now return the solution
	return solution;
}

void SCFConvSol::print() const
{
	cout << "**************************************************" << endl;
	cout << "debug information for SCF CONV SOLUTIONS class " << endl;
	cout << "**************************************************" << endl;

	cout << "print out the solution array for DIIS space repeating problem" << endl;
	string sol;
	printf("%-5s %-35s\n", "index", "solution");
	for(UInt i=0; i<MAXSOL; i++) {
		getSolutionName(spaceVecRepeat[i],sol);
		printf("%-5d %-35s\n", (Int)i, sol.c_str());
	}
	printf("\n");

	cout << "print out the solution array for \"current SCF index not in "
		"result scf vector\" problem" << endl;
	printf("%-5s %-35s\n", "index", "solution");
	for(UInt i=0; i<MAXSOL; i++) {
		getSolutionName(currentScfIndexNotInside[i],sol);
		printf("%-5d %-35s\n", (Int)i, sol.c_str());
	}
	printf("\n");

	cout << "print out the solution array for \" SCF in scfOscillation status\" problem" << endl;
	printf("%-5s %-35s\n", "index", "solution");
	for(UInt i=0; i<MAXSOL; i++) {
		getSolutionName(scfOscillation[i],sol);
		printf("%-5d %-35s\n", (Int)i, sol.c_str());
	}
	printf("\n");
}
