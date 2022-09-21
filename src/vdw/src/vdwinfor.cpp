/**
 * \file    vdwinfor.cpp
 * \brief   top class to mornitor vdw calculation
 * \author  Fenglai Liu and Jing Kong
 */
#include "textread.h"
#include "parameterparsing.h"
#include "globalinfor.h"
#include "excep.h"
#include "shell.h"
#include "molecule.h"
#include "matrix.h"
#include "denmtrx.h"
#include "xcintsinfor.h"
#include "xdm.h"
#include "vdwinfor.h"
using namespace textread;
using namespace parameterparsing;
using namespace globalinfor;
using namespace excep;
using namespace shell;
using namespace molecule;
using namespace matrix;
using namespace denmtrx;
using namespace xcintsinfor;
using namespace xdm;
using namespace vdwinfor;

VDWInfor::VDWInfor(const GlobalInfor& infor, const Molecule& mol):order(3),
	withNonAdditiveB3(false),energy(ZERO),a1(0.83E0), a2(1.55E0*ANGSTROM_TO_BOHR),
	funcName("VDWBR89"),printXDMTiming(false)
{ 
	// deal with input key words
	UInt section = mol.getSec();
	string input = infor.getInputFile();
	ParameterParsing pp(input,"vdw",section);
	if (! pp.hasAnyParameters()) return;

	// parameter of a1 and a2 in xdm expression
	// in default all given in angstrom
	WordConvert w;
	string key = "xdm_parameters";
	if (pp.hasKeyDefined(key)) {
		string value = pp.getValue(key);
		LineParse l(value);
		if (l.getNPieces() == 2) {

			// a1
			Double tmp = ZERO;
			string tmpValue = l.findValue(0);
			if (!w.toDouble(tmpValue,tmp)) {
				string infor = "we can not process xdm a1 value. not double.";
				Excep excep("VDWInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			a1 = tmp;

			// a2
			tmp = ZERO;
			tmpValue = l.findValue(1);
			if (!w.toDouble(tmpValue,tmp)) {
				string infor = "we can not process xdm a2 value. not double.";
				Excep excep("VDWInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			a2 = tmp*ANGSTROM_TO_BOHR;
		}else{
			string infor = "the xdm parameters should have and only have two data fields, you failed on that";
			Excep excep("VDWInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}
	}

	// the functional used for calculating exchange hole
	key = "xdm_exchange_hole";
	if (pp.hasKeyDefined(key)) {
		string value = pp.getValue(key);
		funcName = value;
	}

	// the functional used for calculating exchange hole
	key = "xdm_non_addative_b3";
	if (pp.hasKeyDefined(key)) {
		string value = pp.getValue(key);
		w.capitalize(value);
		if (value == "TRUE" || value == "T") {
			withNonAdditiveB3 = true;
		}else if (value == "FALSE" || value == "F") {
			withNonAdditiveB3 = false;
		}else{
			string infor = "only true or false allowed to process xdm_non_addative_b3";
			Excep excep("VDWInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}
	}

	// do you want to print out timiing data for VDW part?
	key = "print_xdm_timing";
	if (pp.hasKeyDefined(key)) {
		string value = pp.getValue(key);
		w.capitalize(value);
		if (value == "TRUE" || value == "T") {
			printXDMTiming = true;
		}else if (value == "FALSE" || value == "F") {
			printXDMTiming = false;
		}else{
			string infor = "only true or false allowed to process print_xdm_timing";
			Excep excep("VDWInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}
	}
}

void VDWInfor::doXDM(const Molecule& mol,const MolShell& ms,const XCIntJobInfor& xcinfor,const DenMtrx& den)
{
	XDM xdm(*this,ms,xcinfor);
	xdm.energy(*this,ms,mol,den,printXDMTiming);
	energy = xdm.getEnergy();
}



