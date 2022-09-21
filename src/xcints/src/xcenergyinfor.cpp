/**
 * CPP files corresponding to the xcenergyinfor.h
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include<cstdio>
#include<cmath>
#include "blas1.h"
#include "excep.h"
#include "xcfunc.h"
#include "xcintsinfor.h"
#include "batchgrid.h"
#include "xcenergyinfor.h"
using namespace excep;
using namespace blas;
using namespace xcfunc;
using namespace xcintsinfor;
using namespace batchgrid;
using namespace xcenergyinfor;
using namespace std;

XCEnergyInfor::XCEnergyInfor(const XCIntJobInfor& infor, 
		const XCFunc& xcfunc):nGrids(0),nEComponents(0) 
{
	// if it's none, which is the default value; then we do not do it
	// and return
	if (! infor.doXCEnergyProfile()) return;

	// now let's see the functional
	// it's not all of XC functional can do the energy profiling,
	// we can perform the job for the following functionals
	// if the given exchange correlation name is not in the below
	// list, we will issue a warning message and move on
	string exFuncName = xcfunc.getEXName();

	// let's see that whether it's required for xc energy profiling
	// because we do not know the batch grid number, therefore
	// the batchEList will be initialized later
	if (exFuncName == "KP14"  || exFuncName == "KP14_MNDPAR1" || exFuncName == "KP14_MNDPAR2") {

		// let's initialize the data
		nEComponents = 4; // B13 correlation OPP and PAR, KP14C two components
		xcEList.assign(nEComponents,ZERO);
		xcENameList.reserve(nEComponents);
		xcENameList.push_back(XC_E_PROFILE_KP14C_U_STATIC_OPP);
		xcENameList.push_back(XC_E_PROFILE_KP14C_U_STATIC_PAR);
		xcENameList.push_back(XC_E_PROFILE_B13COOR_OPP);
		xcENameList.push_back(XC_E_PROFILE_B13COOR_PAR);
	}
}

UInt XCEnergyInfor::getTheNameCode(const string& exname, UInt iFunc) const
{
	if (exname == "B13COOR_OPP") return XC_E_PROFILE_B13COOR_OPP;
	if (exname == "B13COOR_PAR") return XC_E_PROFILE_B13COOR_PAR;
	if (exname == "KP14C") {
		if (iFunc == 1) return XC_E_PROFILE_KP14C_U_STATIC_OPP;
		if (iFunc == 2) return XC_E_PROFILE_KP14C_U_STATIC_PAR;
	}
	return -1;
}

void XCEnergyInfor::getName(const UInt& code, string& name) const
{
	name = "INVALID_XC_NAME_TO_RETURN";
	if (code == XC_E_PROFILE_B13COOR_OPP) name = "B13COOR_OPP";
	if (code == XC_E_PROFILE_B13COOR_PAR) name = "B13COOR_PAR";
	if (code == XC_E_PROFILE_KP14C_U_STATIC_OPP) name = "KP14C_U_STATIC_OPP";
	if (code == XC_E_PROFILE_KP14C_U_STATIC_PAR) name = "KP14C_U_STATIC_PAR";
	return;
}

void XCEnergyInfor::init(const BatchGrid& batchGrid)
{
	nGrids = batchGrid.getNGrids();
	batchEList.assign(nGrids*nEComponents,ZERO);
}

void XCEnergyInfor::print() const 
{
	if (nEComponents == 0) return;
	string params = "";
	cout << "$ParamOptLSF ";
	for(UInt iFunc=0; iFunc<xcENameList.size(); iFunc++) {

		// get the functional component name
		string name;
		getName(xcENameList[iFunc],name);

		// double check, it should not be zero, or close to zero
		Double exc = xcEList[iFunc];
		params = params + name + " ";
		printf("%-16.9f ", exc);
	}
	cout << "$end" << endl;
	cout << "#" << params << endl;

}

void XCEnergyInfor::formBatchXCEnergy(const BatchGrid& batchGrid)
{
	for(UInt iFunc=0; iFunc<xcENameList.size(); iFunc++) {
		const Double* wts = batchGrid.getGridWts();
		const Double* funval = &batchEList[iFunc*nGrids];
		Double exc = vdot(funval,wts,nGrids);
		xcEList[iFunc] += exc;
	}
}

void XCEnergyInfor::updateEComponent(const Double* F, const string& funcName, UInt iFunc)
{
	// firstly let's search the component name
	// if we can not get it, we just return
	UInt pos = xcENameList.size();
	for(UInt iName=0; iName<xcENameList.size(); iName++) {
		const UInt& name = xcENameList[iName];
		if (getTheNameCode(funcName,iFunc) == name) {
			pos = iName;
			break;
		}
	}

	// do we get the name?
	if (pos == xcENameList.size()) {
		return;
	}

	// now copy the functional value data to the given vector
	Double* E = &batchEList[nGrids*pos];
	for(UInt i=0; i<nGrids; i++) {
		E[i] = F[i];
	}
}

