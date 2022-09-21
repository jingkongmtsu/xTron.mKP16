/**
 * \file    fractional spin information file
 * \brief   collecting the parameters for performing the fractional spin calculation
 * \author  Fenglai Liu and Jing Kong
 */
#include <boost/filesystem.hpp>
#include <cmath>
#include <cstdio>
#include "parameterparsing.h"
#include "textread.h"
#include "excep.h"
#include "molecule.h"
#include "globalinfor.h"
#include "fracspininfor.h"
using namespace boost::filesystem;
using namespace parameterparsing;
using namespace globalinfor;
using namespace textread;
using namespace excep;
using namespace molecule;
using namespace fracspininfor;

FracSpinInfor::FracSpinInfor(const GlobalInfor& globInfor, const Molecule& mol):nAlphaScaledMO(0),
	alphaMOBeginIndex(-1),nBetaScaledMO(0),betaMOBeginIndex(-1),newMol(mol)
{
	// get the input keyword
	string input = globInfor.getInputFile();
	UInt section = mol.getSec();
	ParameterParsing pp(input,"frac_spin",section);
	if (pp.hasAnyParameters()) {

		// alpha fractional spin information : mo begin index
		string key = "alpha_frac_infor_mo_begin_index";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "the alpha mo beginning index can not be processed";
				Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			alphaMOBeginIndex = tmp;
		}

		// alpha fractional spin information : number of mo
		key = "alpha_frac_infor_nmo";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "the number of alpha mo for scaling can not be processed";
				Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			nAlphaScaledMO = tmp; 
		}

		// alpha fractional spin information : the scaled value 
		key = "alpha_frac_infor_scale_value";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			LineParse l(value);
			for(UInt i=0; i<l.getNPieces(); i++) {
				string val = l.findValue(i);
				Double tmp = ZERO;
				if (!w.toDouble(val,tmp)) {
					string infor = "the scale value for alpha mo can not be processed";
					Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
				if (tmp<ZERO) {
					string infor = "the scale value for alpha mo must be > 0";
					Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
				alphaMOScalVals.push_back(sqrt(tmp));
			}
		}

		// beta fractional spin information : mo begin index
		key = "beta_frac_infor_mo_begin_index";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "the beta mo beginning index can not be processed";
				Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			betaMOBeginIndex = tmp;
		}

		// beta fractional spin information : number of mo
		key = "beta_frac_infor_nmo";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "the number of beta mo for scaling can not be processed";
				Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			nBetaScaledMO = tmp; 
		}

		// beta fractional spin information : scaled value
		key = "beta_frac_infor_scale_value";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			LineParse l(value);
			for(UInt i=0; i<l.getNPieces(); i++) {
				string val = l.findValue(i);
				Double tmp = ZERO;
				if (!w.toDouble(val,tmp)) {
					string infor = "the scale value for beta mo can not be processed";
					Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
				if (tmp<ZERO) {
					string infor = "the scale value for beta mo must be > 0";
					Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
				betaMOScalVals.push_back(sqrt(tmp));
			}
		}
	}else{
		return;
	}

	//
	// check the number of scaled value and nmo
	// they should equal with each other
	//
	if (betaMOScalVals.size() != nBetaScaledMO) {
		string infor = "the number of scale values for beta does not equal to the number of scaled mo";
		Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}
	if (alphaMOScalVals.size() != nAlphaScaledMO) {
		string infor = "the number of scale values for alpha does not equal to the number of scaled mo";
		Excep excep("FracSpinInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	//
	// now let's reset the number of alpha and beta electrons
	// since the scaling is done from top to down order for mo,
	// then the given MO starting index denotes the number of electrons
	// for the alpha/beta spin
	//
	// this is because the scaled MO always on the top of all MO
	//
	if (nAlphaScaledMO > 0) {
		UInt newAlpha = alphaMOBeginIndex+1;
		newMol.changeEleNum(0,newAlpha);
	}
	if (nBetaScaledMO > 0) {
		UInt newBeta  = betaMOBeginIndex+1;
		newMol.changeEleNum(1,newBeta);
	}

	// finally set the whole molecule into unrestricted case
	if (nAlphaScaledMO > 0 && nBetaScaledMO > 0) {
		newMol.forceUNRestricted();
	}
}

