/**
 * \file    scfparam.cpp
 * \brief   collecting SCF parameters 
 * \author  Fenglai Liu 
 */
#include <cstdio>
#include <boost/filesystem.hpp>
#include "parameterparsing.h"
#include "textread.h"
#include "molecule.h"
#include "shell.h"
#include "excep.h"
#include "scfparam.h"
using namespace boost::filesystem;
using namespace parameterparsing;
using namespace textread;
using namespace shell;
using namespace excep;
using namespace molecule;
using namespace scfparam;

SCFParam::SCFParam(const GlobalInfor& globInfor, 
		const Molecule& mol):gInfor(globInfor),gIntInfor(globInfor,mol),xcIntInfor(globInfor,mol),
	xcFunc(globInfor.getInputFile(),mol.getSec()),scfGuessMethod(SCF_GUESS_CORE),
	moPath("NONE"),denMtrxPath("NONE"),initDenMtrxPath("NONE"),guessDenMtrxPath("NONE"),
	withFile(true),withPostSCF(true),isCloseShell(mol.isCloseShell()),
	nSpin(mol.getNSpin()),convCriteria(1.0E-5),maxSCFCycles(50),section(mol.getSec()),
	doDipoleMOMPostSCF(true),doVDWPostSCF(false),
	keepSCFSilent(false),doJKSeparately(false),printSCFTiming(false),saveJKMtrx(false),saveXCMtrx(false),
	printInitDenMtrx(false),printDenMtrx(false),printResultDenMtrx(false),
	printCoreMtrx(false),printJKMtrx(false),printXCMtrx(false),
	printCalFockMtrx(false),printSCFFockMtrx(false),printSCFConv(false),printMO(0)         
{
	string input = globInfor.getInputFile();
	UInt section = mol.getSec();
	ParameterParsing pp(input,"scf",section);
	if (pp.hasAnyParameters()) {

		// threshold value for SCF convergence
		string key = "scf_convergence";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "invalid scf_convergence value set";
				Excep excep("SCFParam","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			convCriteria = tmp;
		}

		// SCF cycles
		key = "max_scf_cycles";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "max_scf_cycles can not be processed, not an integer";
				Excep excep("SCFParam","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp > 0) {
				maxSCFCycles = tmp;
			}else{
				string infor = "invalid number of maximum SCF cycles set";
				Excep excep("SCFParam","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// do we use file to keep result in SCF?
		key = "keep_result_in_memory";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				withFile = false;
			}
		}

		// scf guess method
		key = "scf_guess";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "READ") {
				scfGuessMethod = SCF_GUESS_READ;
			}else if(value == "CORE") {
				scfGuessMethod = SCF_GUESS_CORE;
			}else if(value == "ATOMIC") {
				scfGuessMethod = SCF_GUESS_ATOMIC;
			}else{
				string infor = "invalid SCF guess method";
				Excep excep("SCFParam","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// do you want to do dipole moment calculation in the post SCF calculation?
		key = "do_dipole_moment";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				doDipoleMOMPostSCF = true;
			}else if (value == "FALSE" || value == "F") {
				doDipoleMOMPostSCF = false;
			}else{
				string infor = "invalid choice for do_dipole_moment, only true or false allowed";
				Excep excep("SCFParam","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// do you want to do vdw calculation in post SCF?
		key = "do_vdw";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				doVDWPostSCF = true;
			}else if (value == "FALSE" || value == "F") {
				doVDWPostSCF = false;
			}else{
				string infor = "invalid choice for do_vdw, only true or false allowed";
				Excep excep("SCFParam","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// read in mo path for reading guess
		key = "guess_mo_path";
		if (pp.hasKeyDefined(key)) {
			moPath = pp.getValue(key);
		}

		// print final density matrix into the given path
		key = "result_density_matrix_path";
		if (pp.hasKeyDefined(key)) {
			denMtrxPath = pp.getValue(key);
		}

		// the guess density matrix path
		// this is for reading the initial density matrix
		key = "guess_density_matrix_path";
		if (pp.hasKeyDefined(key)) {
			guessDenMtrxPath = pp.getValue(key);
		}

		// the initial density matrix path
		// this is for writing the initial density matrix
		key = "initial_density_matrix_path";
		if (pp.hasKeyDefined(key)) {
			initDenMtrxPath = pp.getValue(key);
		}

		// do we need to print out mo when SCF ends?
		key = "print_mo";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "print_mo can not be processed, not an integer";
				Excep excep("SCFParam","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp >=0 && tmp<= 3) {
				printMO = tmp;
			}else{
				printf("print_mo == 0: no mo data print out\n");
				printf("print_mo == 1: only mo energy data print out\n");
				printf("print_mo == 2: only mo matrix print out\n");
				printf("print_mo == 3: both mo energy and matrix data print out\n");
				string infor = "invalid value of print_mo, only 0, 1, 2 and 3 allowed";
				Excep excep("SCFParam","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// do you want to do coulomb/exchange apart?
		key = "do_jk_apart";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				doJKSeparately = true;
			}
		}

		// do you want to save coulomb/exchange matrix on disk?
		key = "save_jk_matrix";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				saveJKMtrx = true;
			}
		}

		// do you want to save xc matrix on disk?
		key = "save_xc_matrix";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				saveXCMtrx = true;
			}
		}

		// do you want to print out the timing data for SCF?
		key = "print_scf_timing";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				printSCFTiming = true;
			}
		}

		//
		// now below is for debugging print out
		//

		// whether to print out xc matrix during the SCF?
		key = "xcmtrx";
		if (pp.hasKeyDefined(key)) {
			printXCMtrx = true; 
		}

		// whether to print out jk matrix during the SCF?
		key = "jkmtrx";
		if (pp.hasKeyDefined(key)) {
			printJKMtrx = true; 
		}

		// whether to print out core matrix?
		key = "coremtrx";
		if (pp.hasKeyDefined(key)) {
			printCoreMtrx = true; 
		}

		// whether to print out initial density matrix?
		key = "initial_denmtrx";
		if (pp.hasKeyDefined(key)) {
			printInitDenMtrx = true; 
		}

		// whether to print out result density matrix?
		key = "result_denmtrx";
		if (pp.hasKeyDefined(key)) {
			printResultDenMtrx = true; 
		}

		// whether to print out density matrix?
		key = "denmtrx";
		if (pp.hasKeyDefined(key)) {
			printDenMtrx = true; 
		}

		// whether to print out calculated fock matrix?
		key = "fock";
		if (pp.hasKeyDefined(key)) {
			printCalFockMtrx = true; 
		}

		// whether to print out fock matrix for scf convergence?
		key = "scf_fock";
		if (pp.hasKeyDefined(key)) {
			printSCFFockMtrx = true; 
		}

		// whether to print out scfconv?
		key = "scfconv";
		if (pp.hasKeyDefined(key)) {
			printSCFConv = true; 
		}
	}

	// post processing for reading guess
	if (scfGuessMethod == SCF_GUESS_READ) {

		// you need to provide some file
		if(moPath == "NONE") {
			string infor = "no mo files information is given for reading the scf guess";
			Excep excep("SCFParam","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}

		// for the given file, it should exist
		path p(moPath.c_str());
		if (! exists(p)) {
			string info = "the corresponding mo file folder does not exist: " + p.string();
			Excep excep("SCFParam","constructor",EXCEPTION_FILE_MISSING,info);
			handleExcep(excep);
		}
	}

	// post processing section for 
	// gintsInfor and xcFunc
	gIntInfor.updateFromXCFunc(xcFunc);
}

