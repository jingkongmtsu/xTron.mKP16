/**
 * CPP files corresponding to the xcvarinfor.h
 * \author  Fenglai liu and Jing kong
 */
#include <iostream>
#include <string>
#include "excep.h"
#include "xcvarinfor.h"
using namespace excep;
using namespace xcvarinfor;
using namespace std;

string xcvarinfor::getVarName(const UInt& v) {
	if (v == RA) return "RA";
	if (v == RB) return "RB";
	if (v == GAA) return "GAA";
	if (v == GAB) return "GAB";
	if (v == GBB) return "GBB";
	if (v == DAX) return "DAX";
	if (v == DAY) return "DAY";
	if (v == DAZ) return "DAZ";
	if (v == DBX) return "DBX";
	if (v == DBY) return "DBY";
	if (v == DBZ) return "DBZ";
	if (v == TA) return "TA";
	if (v == TB) return "TB";
	if (v == LA) return "LA";
	if (v == LB) return "LB";
	if (v == EXA) return "EXA";
	if (v == EXB) return "EXB";
	return "NONE";
}

const UInt* xcvarinfor::getDFTVarArray (const UInt& var) {
	if (var == RHO) {
		return RHO_ARRAY;
	}else if (var == GRHO) {
		return GRHO_ARRAY;
	}else if (var == GAMMA) {
		return GAMMA_ARRAY;
	}else if (var == TAU) {
		return TAU_ARRAY;
	}else if (var == LAP) {
		return LAP_ARRAY;
	}else if (var == EXRHO) {
		return EXRHO_ARRAY;
	}

	// something error here
	// we just return null pointer here
	string infor = "the input var name is invalid: " + getVarName(var);
	Excep excep("xcvarinfor","getDFTVarArray",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
	handleExcep(excep);
	return NULL;
}

UInt xcvarinfor::getDFTVarArrayLength (const UInt& var) {
	if (var == GAMMA) return 3; // GAA GAB GBB
	if (var == GRHO)  return 6; // DA and DB on X Y Z
	return 2; // one alpha var and one beta var
}

bool xcvarinfor::isBetaDFTVar(const UInt& var) {
	if (var == RA || var == DAX || var == DAY || var == DAZ ||
			var == GAA || var == TA || var == LA || var == EXA) {
		return false;
	}else{
		return true;
	}
}

UInt xcvarinfor::getCounterpartVar(const UInt& v) {
	if (v == RA) return RB;
	if (v == RB) return RA;
	if (v == GAA) return GBB;
	if (v == GAB) return GAB;
	if (v == GBB) return GAA;
	if (v == DAX) return DBX;
	if (v == DAY) return DBY;
	if (v == DAZ) return DBZ;
	if (v == DBX) return DAX;
	if (v == DBY) return DAY;
	if (v == DBZ) return DAZ;
	if (v == TA) return TB;
	if (v == TB) return TA;
	if (v == LA) return LB;
	if (v == LB) return LA;
	if (v == EXA) return EXB;
	if (v == EXB) return EXA;

	// now it's something error
	string infor = "the input var name is invalid: " + getVarName(v);
	Excep excep("xcvarinfor","getCounterpartVar",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
	handleExcep(excep);
	return 0;
}

bool xcvarinfor::isGammaVar(const UInt& var) {
	if ( var == GAA || var == GAB || var == GBB ) {
		return true;
	}else{
		return false;
	}
}

bool xcvarinfor::isDRhoVar(const UInt& var) {
	if ( var == DAX || var == DAY || var == DAZ ||
			var == DBX || var == DBY || var == DBZ ) {
		return true;
	}else {
		return false;
	}
}

bool xcvarinfor::doesDRhoVarMatchGammaVar(const UInt& gvar, const UInt& dvar) {
	if (dvar == DAX || dvar == DAY || dvar == DAZ) {
		if (gvar == GAA || gvar == GAB) return true;
		return false;
	}else if(dvar == DBX || dvar == DBY || dvar == DBZ) {
		if (gvar == GBB || gvar == GAB) return true;
		return false;
	}else {
		return false;
	}
}

Double xcvarinfor::getFacGammaIntoDRho(const UInt& gammaVar) {
	if (gammaVar == GAB) return 1.0E0;
	return 2.0E0;
}

UInt xcvarinfor::getVarGammaIntoDRho(const UInt& gammaVar, const UInt& dVar) {
	if (gammaVar == GAB) {
		if (dVar == DAX) return DBX;
		if (dVar == DAY) return DBY;
		if (dVar == DAZ) return DBZ;
		if (dVar == DBX) return DAX;
		if (dVar == DBY) return DAY;
		if (dVar == DBZ) return DAZ;

		// something wrong here
		string infor = "the input gamma var is GAB, dvar name is invalid: " + getVarName(dVar);
		Excep excep("xcvarinfor","getVarGammaIntoDRho",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
		return 0;
	}else{
		return dVar;
	}
}

