/**
 * CPP files corresponding to the xcvar.h
 * \author  Fenglai liu and Jing kong
 */
#include <iostream>
#include "xcfunc.h"
#include "excep.h"
#include "integraljobs.h"
#include "xcvarinfor.h"
#include "xcvar.h"
using namespace xcfunc;
using namespace excep;
using namespace integraljobs;
using namespace xcvarinfor;
using namespace xcvar;

void XCVar::setDerivOrder(const XCFunc& xcfunc, const UInt& job, const UInt& jobOrder)
{
	// set the variable, functional derivatives order
	// according to job and job order
	if (job == GROUND_STATE_DFT) {
		if (jobOrder == 0) {
			maxVarDerivOrder  = 0;   
			maxFunDerivOrder  = 1;   
		}else if (jobOrder == 1) {
			maxVarDerivOrder  = 1;   
			maxFunDerivOrder  = 1;   
		}else if (jobOrder == 2) {
			maxVarDerivOrder  = 2;   
			maxFunDerivOrder  = 2;   
		}
	}else if (job == TDDFT) {
		if (jobOrder == 0) {
			maxVarDerivOrder  = 0;   
			maxFunDerivOrder  = 2;   
		}else if (jobOrder == 1) {
			maxVarDerivOrder  = 1;   
			maxFunDerivOrder  = 3;   
		}else if (jobOrder == 2) {
			maxVarDerivOrder  = 2;   
			maxFunDerivOrder  = 4;   
		}
	}

	// this is for the energy
	maxBasisDerivOrder = 0;
	if (xcfunc.hasVarTau() || xcfunc.hasVarGRho()) {
		maxBasisDerivOrder = 1;
	}
	if (xcfunc.hasVarLap()) {
		maxBasisDerivOrder = 2;
	}

	// for derivative calculation, 
	// we also need to add in derivative order
	maxBasisDerivOrder += jobOrder;
}

XCVar::XCVar(const XCFunc& xcfunc, const UInt& nSpin, const UInt& job, 
		const UInt& jobOrder):maxVarDerivOrder(-1), 
	maxFunDerivOrder(-1),maxBasisDerivOrder(-1),
	totalAlphaVar(0),totalVarNum(0),totalNumVarWithGamma(0)
{
	// set the variable and functional derivatives order
	// we need the maxBasisDerivOrder for all of Exc/non-Exc calculations
	// therefore in xcvar we still call this function then
	// check whether this is normal DFT job
	setDerivOrder(xcfunc,job,jobOrder);

	// for non-DFT job, we just return
	if (! isDFTJob(job)) return;

	// whether we have beta variables?
	bool hasBetaVar = true;
	if (nSpin == 1) hasBetaVar = false;

	// get the variable types information
	varTypes.reserve(MAX_DFTVAR_TYPES);
	if (xcfunc.hasVarRho()) {
		varTypes.push_back(RHO);
	}
	if (xcfunc.hasVarGRho()) {
		varTypes.push_back(GRHO);
	}
	if (xcfunc.hasVarTau()) {
		varTypes.push_back(TAU);
	}
	if (xcfunc.hasVarLap()) {
		varTypes.push_back(LAP);
	}
	if (xcfunc.hasVarExRho()) {
		varTypes.push_back(EXRHO);
	}

	// now let's count the variables
	// we note that for variables with gamma,
	// we do not distinguish the close/open 
	// shell case -since in the fortran code,
	// we will calculate them together
	totalAlphaVar = 0;
	totalNumVarWithGamma = 0;
	for(UInt i=0; i<varTypes.size(); i++) {
		if (varTypes[i] == GRHO) {
			totalAlphaVar += 3;
			totalNumVarWithGamma += 3;
		}else{
			totalAlphaVar += 1;
			totalNumVarWithGamma += 2;
		}
	}
	totalVarNum = totalAlphaVar;
	if (hasBetaVar) totalVarNum = 2*totalAlphaVar;

	// finally, let's set the functional derivatives information
	// first for non-gamma variable list
	// second for gamma variable list
	for(UInt order=1; order<=maxFunDerivOrder; order++) {
		setVarDeriv(order, false);
		setVarDeriv(order, true);
	}
}

void XCVar::setVarDeriv(const UInt& order, bool withGamma) 
{
	// for the d1Vars etc. if we do not have beta var
	// then in the list the beta var should not exist
	// however, since functional derivatives calculation
	// always has both alpha var and beta var, therefore
	// for creating variable list for gamma we always
	// including beta var
	bool noBetaVar = false;
	if (!withGamma && !hasBetaVar()) noBetaVar = true;

	// now doing the job for each order
	UInt varTypesNum = varTypes.size();
	if (order == 1) {

		// now let's set the variable information
		// we will put it into the d1Vars
		// so the information is used in both 
		// variable generation as well as 1st 
		// functional derivatives
		UInt n = totalVarNum;
		if (withGamma) {
			n = getNumVarWithGamma();
			d1GVars.reserve(n);      
		}else{
			d1Vars.reserve(n);      
		}
		for(UInt ivar=0; ivar<varTypesNum; ivar++) {
			UInt vType = varTypes[ivar];
			if (vType == GRHO && withGamma) vType = GAMMA;
			const UInt* vars = getDFTVarArray(vType);
			for(UInt i=0; i<getDFTVarArrayLength(vType); i++) {
				UInt v = vars[i];
				if (noBetaVar && isBetaDFTVar(v)) continue;
				if(withGamma) {
					d1GVars.push_back(v);
				}else{
					d1Vars.push_back(v);
				}
			}
		}
	}else if (order == 2) {

		UInt n = getNumFuncDeriv(order,withGamma);
		if (withGamma) {
			d2GVars.reserve(2*n);      
		}else{
			d2Vars.reserve(2*n);      
		}

		// now we step into var types counting
		for(UInt iNVarType=0; iNVarType<varTypesNum; iNVarType++) {
			for(UInt jNVarType=0; jNVarType<=iNVarType; jNVarType++) {
				UInt iVType = varTypes[iNVarType];
				UInt jVType = varTypes[jNVarType];
				if (iVType == GRHO && withGamma) iVType = GAMMA;
				if (jVType == GRHO && withGamma) jVType = GAMMA;

				if (iVType == jVType) {
					const UInt* vars = getDFTVarArray(iVType);
					for(UInt i=0; i<getDFTVarArrayLength(iVType); i++) {
						for(UInt j=0; j<=i; j++) {
							UInt iv = vars[i];
							UInt jv = vars[j];
							if (noBetaVar && isBetaDFTVar(iv)) continue;
							if (noBetaVar && isBetaDFTVar(jv)) continue;
							if(withGamma) {
								d2GVars.push_back(jv);
								d2GVars.push_back(iv);
							}else{
								d2Vars.push_back(jv);
								d2Vars.push_back(iv);
							}
						}
					}
				}else{
					const UInt* iVars = getDFTVarArray(iVType);
					const UInt* jVars = getDFTVarArray(jVType);
					for(UInt i=0; i<getDFTVarArrayLength(iVType); i++) {
						for(UInt j=0; j<getDFTVarArrayLength(jVType); j++) {
							UInt iv = iVars[i];
							UInt jv = jVars[j];
							if (noBetaVar && isBetaDFTVar(iv)) continue;
							if (noBetaVar && isBetaDFTVar(jv)) continue;
							if(withGamma) {
								d2GVars.push_back(jv);
								d2GVars.push_back(iv);
							}else{
								d2Vars.push_back(jv);
								d2Vars.push_back(iv);
							}
						}
					}
				}
			}
		}
	}else if (order == 3) {

		UInt n = getNumFuncDeriv(order,withGamma);
		if (withGamma) {
			d3GVars.reserve(3*n);      
		}else{
			d3Vars.reserve(3*n);      
		}

		// now we step into var types counting
		for(UInt iNVarType=0; iNVarType<varTypesNum; iNVarType++) {
			for(UInt jNVarType=0; jNVarType<=iNVarType; jNVarType++) {
				for(UInt kNVarType=0; kNVarType<=jNVarType; kNVarType++) {
					UInt iVType = varTypes[iNVarType];
					UInt jVType = varTypes[jNVarType];
					UInt kVType = varTypes[kNVarType];
					if (iVType == GRHO && withGamma) iVType = GAMMA;
					if (jVType == GRHO && withGamma) jVType = GAMMA;
					if (kVType == GRHO && withGamma) kVType = GAMMA;

					// here we note that we do not need to judge
					// whether iVType == kVType
					// in this case, it must have jVType == kVType
					if (iVType == jVType || jVType == kVType) {
						const UInt* iVars = getDFTVarArray(iVType);
						const UInt* jVars = getDFTVarArray(jVType);
						const UInt* kVars = getDFTVarArray(kVType);
						for(UInt i=0; i<getDFTVarArrayLength(iVType); i++) {
							UInt jNVar = getDFTVarArrayLength(jVType) - 1;
							if (iVType == jVType) jNVar = i;
							for(UInt j=0; j<=jNVar; j++) {
								UInt kNVar = getDFTVarArrayLength(kVType) - 1;
								if (jVType == kVType) kNVar = j;
								for(UInt k=0; k<=kNVar; k++) {
									UInt iv = iVars[i];
									UInt jv = jVars[j];
									UInt kv = kVars[k];
									if (noBetaVar && isBetaDFTVar(iv)) continue;
									if (noBetaVar && isBetaDFTVar(jv)) continue;
									if (noBetaVar && isBetaDFTVar(kv)) continue;
									if(withGamma) {
										d3GVars.push_back(kv);
										d3GVars.push_back(jv);
										d3GVars.push_back(iv);
									}else{
										d3Vars.push_back(kv);
										d3Vars.push_back(jv);
										d3Vars.push_back(iv);
									}
								}
							}
						}
					}else{
						const UInt* iVars = getDFTVarArray(iVType);
						const UInt* jVars = getDFTVarArray(jVType);
						const UInt* kVars = getDFTVarArray(kVType);
						for(UInt i=0; i<getDFTVarArrayLength(iVType); i++) {
							for(UInt j=0; j<getDFTVarArrayLength(jVType); j++) {
								for(UInt k=0; k<getDFTVarArrayLength(kVType); k++) {
									UInt iv = iVars[i];
									UInt jv = jVars[j];
									UInt kv = kVars[k];
									if (noBetaVar && isBetaDFTVar(iv)) continue;
									if (noBetaVar && isBetaDFTVar(jv)) continue;
									if (noBetaVar && isBetaDFTVar(kv)) continue;
									if(withGamma) {
										d3GVars.push_back(kv);
										d3GVars.push_back(jv);
										d3GVars.push_back(iv);
									}else{
										d3Vars.push_back(kv);
										d3Vars.push_back(jv);
										d3Vars.push_back(iv);
									}
								}
							}
						}
					}
				}
			}
		}
	}else{
		string infor = "currently we only support order <=3 for doing functional derivatives";
		Excep excep("XCVar","setVarDeriv", EXCEPTION_XCINTS_INVALID_FUNC_DERIV_ORDER, infor);
		handleExcep(excep);
	}
}

bool XCVar::hasGRho() const {
	for(UInt pos=0; pos<totalVarNum; pos++) {
		if (d1Vars[pos] == DAX || d1Vars[pos] == DAY ||
				d1Vars[pos] == DAZ || d1Vars[pos] == DBX ||
				d1Vars[pos] == DBY || d1Vars[pos] == DBZ
			) return true;
	}
	return false;
}

bool XCVar::hasTau() const {
	for(UInt pos=0; pos<totalVarNum; pos++) {
		if (d1Vars[pos] == TA || d1Vars[pos] == TB ) return true;
	}
	return false;
}

bool XCVar::hasLap() const {
	for(UInt pos=0; pos<totalVarNum; pos++) {
		if (d1Vars[pos] == LA || d1Vars[pos] == LB ) return true;
	}
	return false;
}

bool XCVar::hasExRho() const {
	for(UInt pos=0; pos<totalVarNum; pos++) {
		if (d1Vars[pos] == EXA || d1Vars[pos] == EXB ) return true;
	}
	return false;
}

bool XCVar::hasRho() const {
	for(UInt pos=0; pos<totalVarNum; pos++) {
		if (d1Vars[pos] == RA || d1Vars[pos] == RB ) return true;
	}
	return false;
}

UInt XCVar::getNumFuncDeriv(const UInt& order, bool withGamma) const {
	UInt n = totalVarNum;
	if (withGamma) n = getNumVarWithGamma();
	if (order == 1) {
		return n;
	}else if (order == 2) {
		return (n*(n+1))/2;
	}else if (order == 3) {
		return (n*(n+1)*(n+2))/6;
	}else{
		string infor = "currently we only support order <=3 for doing functional derivatives";
		Excep excep("XCVar","getNumFuncDeriv", EXCEPTION_XCINTS_INVALID_FUNC_DERIV_ORDER, infor);
		handleExcep(excep);
		return -1;
	}
}

UInt XCVar::getVar(const UInt& i, bool withGamma) const { 

	// since XCVar may be empty, so there's possible no variable
	// information available. We check the situation here
	if (d1GVars.size() == 0) {
		string infor = "there's no any variable defined in the XCVar class";
		Excep excep("XCVar","getVar",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
		return -1;
	}

	// now return the stuff
	if (withGamma) return d1GVars[i];
	return d1Vars[i];
}

UInt XCVar::getVarPos(const UInt& var, bool withGamma) const {
	// if the given variable is in beta, and we do not
	// have beta; we will return it's alpha counterpart
	if (withGamma) {
		UInt v = var;
		for(UInt pos=0; pos<totalNumVarWithGamma; pos++) {
			if (d1GVars[pos] == v) return pos;
		}
		string infor = "fail to find the position for given var: " + getVarName(var);
		Excep excep("XCVar","getVarPos",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
		return -1;
	}else{
		UInt v = var;
		if(!hasBetaVar() && isBetaDFTVar(var)) v = getCounterpartVar(var);
		for(UInt pos=0; pos<totalVarNum; pos++) {
			if (d1Vars[pos] == v) return pos;
		}
		string infor = "fail to find the position for given var: " + getVarName(var);
		Excep excep("XCVar","getVarPos",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
		return -1;
	}
} 

UInt XCVar::getVarPos(const UInt& var1, const UInt& var2, bool withGamma) const {
	// if the given variable is in beta, and we do not
	// have beta; we will return it's alpha counterpart
	// additionally, we consider the situation that:
	// F''_(x,y) = F''_{y,x}
	if (withGamma) {
		UInt v1 = var1;
		UInt v2 = var2;
		for(UInt pos=0; pos<getNumFuncDeriv(2,withGamma); pos++) {
			if (d2GVars[2*pos] == v1 && d2GVars[2*pos+1] == v2) return pos;
			if (d2GVars[2*pos] == v2 && d2GVars[2*pos+1] == v1) return pos;
		}
		string infor = "fail to find the position for given vars: " + getVarName(v1) + getVarName(v2);
		Excep excep("XCVar","getVarPos",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
		return -1;
	}else{
		UInt v1 = var1;
		UInt v2 = var2;
		if(!hasBetaVar() && isBetaDFTVar(v1)) v1 = getCounterpartVar(var1);
		if(!hasBetaVar() && isBetaDFTVar(v2)) v2 = getCounterpartVar(var2);
		for(UInt pos=0; pos<getNumFuncDeriv(2,withGamma); pos++) {
			if (d2Vars[2*pos] == v1 && d2Vars[2*pos+1] == v2) return pos;
			if (d2Vars[2*pos] == v2 && d2Vars[2*pos+1] == v1) return pos;
		}
		string infor = "fail to find the position for given vars: " + getVarName(v1) + getVarName(v2);
		Excep excep("XCVar","getVarPos",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
		return -1;
	}
} 

bool XCVar::do1stGradHalfVar() const {

	// if the maxVarDerivOrder is not defined, we need to 
	// issue an error
	if (maxVarDerivOrder == static_cast<UInt>(-1)) {
		string infor = "maxVarDerivOrder is not defined, so that this function is meaningless to call";
		Excep excep("XCVar","do1stGradHalfVar",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
	}


	// for gradient variable, if we only have rho; then we do not
	// need to do it 
	if (maxVarDerivOrder == 0) {

		// for energy calculation
		// if we have tau. lap etc. we must do it
		// else we do not need it
		if (hasTau() || hasLap()) {
			return true;
		}else{
			return false;
		}
	} else {
		string infor = "for the case maxVarDerivOrder > 0 it's not implemented yet";
		Excep excep("XCVar","do1stGradHalfVar",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
		handleExcep(excep);
	}

	return false;

}

bool XCVar::do2edGradHalfVar() const {

	// if the maxVarDerivOrder is not defined, we need to 
	// issue an error
	if (maxVarDerivOrder == static_cast<UInt>(-1)) {
		string infor = "maxVarDerivOrder is not defined, so that this function is meaningless to call";
		Excep excep("XCVar","do2edGradHalfVar",EXCEPTION_XCINTS_INVALID_XCVAR,infor);
		handleExcep(excep);
	}

	// for the order 0 DFT variable, in default we do not need to
	// do the 2ed gradient half variable

	// do we do second derivatives?
	if (maxVarDerivOrder >= 2) {
		return true;
	}

	// for 1st order gradient variable, we should consider
	// the situation
	if (maxVarDerivOrder == 1) {
		if (hasTau() || hasLap()) {
			return true;
		}else{
			return false;
		}
	}

	return false;

}

void XCVar::setupVarInfor(Int* infor) const
{
	//
	// we note, that here the code is strictly following
	// the variable sequence:
	// rho->grho->tau->lap->exrho
	// each of them has a position in the infor array
	// we judge wether it exists in the given array
	//
	for(UInt i=0; i<MAX_DFTVAR_TYPES; i++) infor[i] = -1;
	if (hasRho())  infor[0] = 1;
	if (hasGRho()) infor[1] = 1;
	if (hasTau())  infor[2] = 1;
	if (hasLap())  infor[3] = 1;
	if (hasExRho())infor[4] = 1;
}

bool XCVar::neglectVar(const UInt& var) const
{
	// note that we only do neglecting work
	// if we do not have beta variable
	if (hasBetaVar()) return false;

	// now real job
	if (isGammaVar(var)) {
		// for close shell, we do not count in the GBB, but GAB should be in count
		if (var == GBB) return true;
	}else{
		// for close shell, we do not count in the beta part
		if (isBetaDFTVar(var)) return true;
	}

	// default return status
	return false;
}

void XCVar::print() const {

	cout << "*******************************************" << endl;
	cout << "*            XCVar Results                *" << endl;
	cout << "*******************************************" << endl;

	// if no variable defined, just return
	if (d1Vars.size() == 0) return;

	// d1vars
	cout << "d1Vars data " << endl; 
	for(UInt i=0; i<d1Vars.size(); i++) {
		UInt v = d1Vars[i];
		string vname = getVarName(v);
		cout << vname << endl;
	}
	cout << endl;

	// d1gvars
	cout << "d1GVars data " << endl; 
	for(UInt i=0; i<d1GVars.size(); i++) {
		UInt v = d1GVars[i];
		string vname = getVarName(v);
		cout << vname << endl;
	}
	cout << endl;

	if (d2Vars.size() > 0) {

		// d2vars
		cout << "d2Vars data " << endl; 
		UInt n = getNumFuncDeriv(2,false);
		for(UInt i=0; i<n; i++) {
			UInt v1 = d2Vars[2*i+0];
			UInt v2 = d2Vars[2*i+1];
			string v1name = getVarName(v1);
			string v2name = getVarName(v2);
			cout << v1name  << "_" << v2name << endl;
		}
		cout << endl;

		// d2GVars
		cout << "d2GVars data " << endl; 
		n = getNumFuncDeriv(2,true);
		for(UInt i=0; i<n; i++) {
			UInt v1 = d2GVars[2*i+0];
			UInt v2 = d2GVars[2*i+1];
			string v1name = getVarName(v1);
			string v2name = getVarName(v2);
			cout << v1name  << "_" << v2name << endl;
		}
		cout << endl;
	}

	if (d3Vars.size() > 0) {

		// d3vars
		cout << "d3Vars data " << endl; 
		UInt n = getNumFuncDeriv(3,false);
		for(UInt i=0; i<n; i++) {
			UInt v1 = d3Vars[3*i+0];
			UInt v2 = d3Vars[3*i+1];
			UInt v3 = d3Vars[3*i+2];
			string v1name = getVarName(v1);
			string v2name = getVarName(v2);
			string v3name = getVarName(v3);
			cout << v1name  << "_" << v2name  << "_" << v3name << endl;
		}
		cout << endl;

		// d3gvars
		cout << "d3GVars data " << endl; 
		n = getNumFuncDeriv(3,true);
		for(UInt i=0; i<n; i++) {
			UInt v1 = d3GVars[3*i+0];
			UInt v2 = d3GVars[3*i+1];
			UInt v3 = d3GVars[3*i+2];
			string v1name = getVarName(v1);
			string v2name = getVarName(v2);
			string v3name = getVarName(v3);
			cout << v1name  << "_" << v2name  << "_" << v3name << endl;
		}
		cout << endl;
	}
}

