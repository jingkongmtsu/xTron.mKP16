/**
 * cpp file associated with gintsinfor.h
 * \author Fenglai Liu 
 */
#include <cmath>
#include <boost/lexical_cast.hpp>
#include "parameterparsing.h"
#include "textread.h"
#include "excep.h"
#include "globalinfor.h"
#include "molecule.h"
#include "xcfunc.h"
#include "integraljobs.h"
#include "gintsinfor.h"
using namespace parameterparsing;
using namespace textread;
using namespace excep;
using namespace globalinfor;
using namespace molecule;
using namespace xcfunc;
using namespace integraljobs;
using namespace gintsinfor;

////////////////////////////////////////////////////////////////////
//                       #### GIntsInfor ####                     // 
////////////////////////////////////////////////////////////////////

GIntsInfor::GIntsInfor(const GlobalInfor& infor, const Molecule& mol):ginfor(infor),
	closeShell(mol.isCloseShell()),nSpin(mol.getNSpin()),doScaleforK(false),kCoeff(ONE),
	spGInts2DThresh(1.0E-12),spGInts4DThresh(1.0E-12),gints2DThresh(1.0E-12),
	gints2DThreshDeriv(1.0E-12),gints4DThresh(1.0E-10),gints4DThreshDeriv(1.0E-10),
	nAtomBlocksInList(50),linearDependenceThresh(1.0E-6)
{
	// now deal with input key words
	string input = infor.getInputFile();
	UInt section = mol.getSec();
	ParameterParsing pp(input,"gints",section);
	if (pp.hasAnyParameters()) {

		// threshold value in generating shell pair
		// for gints2d
		string key = "shellpair_threshold_gints2d";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for shell pair. not double.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			spGInts2DThresh = tmp;
		}

		// threshold value in generating shell pair
		// for gints4d
		key = "shellpair_threshold_gints4d";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for shell pair. not double.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			spGInts4DThresh = tmp;
		}

		key = "gints2d_threshold";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for gints2d. not double.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			gints2DThresh = tmp;
		}

		key = "gints4d_threshold";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for gints4d. not double.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			gints4DThresh = tmp;
		}

		key = "gints2deriv_threshold";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for gints2d derivatives. not double.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			gints2DThreshDeriv = tmp;
		}

		key = "gints4deriv_threshold";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for gints4d derivatives. not double.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			gints4DThreshDeriv = tmp;
		}

		key = "linear_dependence_threshold";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for basis set linear dependence. not double.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			linearDependenceThresh = tmp;
		}

		key = "number_atom_block_results_in_list";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "we can not process value for number_atom_block_results_in_list. not UInt.";
				Excep excep("GIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			nAtomBlocksInList = tmp;
		}
	}
}

void GIntsInfor::updateFromXCFunc(const XCFunc& xcfunc) 
{
	// do you have hybrid functional in xcfunc?
	if (xcfunc.isHybrid()) {
		kCoeff = xcfunc.HFCoeff();
		if (fabs(kCoeff-ONE)<THRESHOLD_MATH) {
			doScaleforK = false;
		}else{
			doScaleforK = true;
		}
	}
}

void GIntsInfor::changeThresh(UInt i, Double thresh) {
	if (i == GINTS4D_SHELL_PAIR_THRESH) {
		spGInts4DThresh = thresh;
	}else if (i == GINTS2D_SHELL_PAIR_THRESH) {
		spGInts2DThresh = thresh;
	}else if (i == GINTS2D_THRESH) {
		gints2DThresh = thresh;
	}else if (i == GINTS4D_THRESH) {
		gints4DThresh = thresh;
	}else{
		string infor = "the given thresh mode can not be recognized: " + boost::lexical_cast<string>(i);
		Excep excep("GIntJobInfor","changeThresh",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}
}

////////////////////////////////////////////////////////////////////
//                      #### GIntJobInfor ####                    // 
////////////////////////////////////////////////////////////////////
GIntJobInfor::GIntJobInfor(const GIntsInfor& infor, 
		const UInt& intJob0, UInt jobOrder0):GIntsInfor(infor),intJob(intJob0),
	jobOrder(jobOrder0),engineInfor(intJob,jobOrder)
{
	bool success = checkIntJob(intJob);
	if (! success) {
		string info1 = "the job name is " + getJobName(intJob);
		string info2 = " the job order is " + boost::lexical_cast<string>(jobOrder);
		string infor = info1 + info2;
		Excep excep("GIntJobInfor","constructor",EXCEPTION_GINTS_INVALID_INTEGRAL_JOB,infor);
		handleExcep(excep);
	}
}

