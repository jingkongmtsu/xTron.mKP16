/**
 * CPP files corresponding to the oddelec.h
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include<cstdio>
#include "excep.h"
#include "shell.h"
#include "xcintsinfor.h"
#include "batchvar.h"
#include "xcvar.h"
#include "denmtrx.h"
#include "batchgrid.h"
#include "functionallist.h"
#include "oddelec.h"
using namespace excep;
using namespace shell;
using namespace xcintsinfor;
using namespace batchvar;
using namespace xcvar;
using namespace batchgrid;
using namespace oddelec;
using namespace std;

ODDElec::ODDElec(const MolShell& ms, const XCIntJobInfor& infor) 
{
	if (infor.doOddElec()) {
		oddElecSum.assign(ms.getNAtomShells(),ZERO);
	}
}

void ODDElec::doOddElecPopulation(const XCIntJobInfor& xcinfor, const BatchGrid& bg, 
		const BatchVar& bvar, const XCVar& xcvar,const DenMtrx& den) 
{
	// get the number of electrons
	Int NA = den.getNEle(0);
	Int NB = den.getNEle(1);
	bool isSingleEleSystem = den.isSingleEleSystem();

	// get the density variables for alpha state
	UInt varPos = -1;
	if (xcvar.hasRho()) varPos = xcvar.getVarPos(RA);
	const Double* oriRhoA  = bvar.getVar(varPos);
	varPos = -1;
	if (xcvar.hasGRho()) varPos = xcvar.getVarPos(DAX);
	const Double* oriDRA   = bvar.getVar(varPos);
	varPos = -1;
	if (xcvar.hasTau()) varPos = xcvar.getVarPos(TA);
	const Double* oriTA    = bvar.getVar(varPos);
	varPos = -1;
	if (xcvar.hasLap()) varPos = xcvar.getVarPos(LA);
	const Double* oriLA    = bvar.getVar(varPos);
	varPos = -1;
	if (xcvar.hasExRho()) varPos = xcvar.getVarPos(EXA);
	const Double* oriEXA   = bvar.getVar(varPos);

	// before get beta variables, we need to prepare 
	// empty density variable for single electron case
	// 
	// this is because the functional valculation will
	// use both alpha and beta variables, however; single
	// electron system like H atom they do not have beta
	// part - or in the other word, beta part is all zero
	// therefore here we need to make a special treatment
	DoubleVec betaVar;
	DoubleVec betaGVar;
	UInt ng = bg.getNGrids();
	if (isSingleEleSystem) {
		betaVar.assign(ng,ZERO);
		betaGVar.assign(3*ng,ZERO);
	}
	
	// now get the beta variables
	const Double* oriRhoB  = NULL;
	const Double* oriDRB   = NULL;
	const Double* oriTB    = NULL;
	const Double* oriLB    = NULL;
	const Double* oriEXB   = NULL;
	if (isSingleEleSystem) {
		if (xcvar.hasRho()) oriRhoB  = &betaVar[0];
		if (xcvar.hasGRho()) oriDRB  = &betaGVar[0];
		if (xcvar.hasTau()) oriTB    = &betaVar[0];
		if (xcvar.hasLap()) oriLB    = &betaVar[0];
		if (xcvar.hasExRho()) oriEXB = &betaVar[0];
	}else{
		varPos = -1;
		if (xcvar.hasRho()) varPos = xcvar.getVarPos(RB);
		oriRhoB  = bvar.getVar(varPos);
		varPos = -1;
		if (xcvar.hasGRho()) varPos = xcvar.getVarPos(DBX);
		oriDRB   = bvar.getVar(varPos);
		varPos = -1;
		if (xcvar.hasTau()) varPos = xcvar.getVarPos(TB);
		oriTB    = bvar.getVar(varPos);
		varPos = -1;
		if (xcvar.hasLap()) varPos = xcvar.getVarPos(LB);
		oriLB    = bvar.getVar(varPos);
		varPos = -1;
		if (xcvar.hasExRho()) varPos = xcvar.getVarPos(EXB);
		oriEXB   = bvar.getVar(varPos);
	}

	// we need to perform some treatment to the inputs 
	// variables first
#ifdef WITH_SINGLE_PRECISION

	// we need to copy the floating point variables into double precision ones
	// by getting the bvar's first posision, we will copy all of data from
	// single precision to double precision
	UInt len = xcvar.getVarNum()*bg.getNGrids();
	DoublePrecVec temp_var(len);
	const Double* var0 = bvar.getVarPos(0);
	for(UInt i=0; i<len; i++) {
		temp_var[i] = static_cast<double>(var0[i]);
	}
	
	// set up single electron case as shown above
	DoublePrecVec temp_betaVar;
	DoublePrecVec temp_betaGVar;
	if (isSingleEleSystem) {
		temp_betaVar.assign(ng,ZERO);
		temp_betaGVar.assign(3*ng,ZERO);
	}

	// now fetch the variable data
	// alpha rho 
	UInt newVarPos = -1;
	if (xcvar.hasRho()) newVarPos = xcvar.getVarPos(RA);
	const double* rhoA  = NULL;
	if (newVarPos >= 0) rhoA  = &temp_var[newVarPos*ng];

	// beta rho
	const double* rhoB  = NULL;
	if (isSingleEleSystem) {
		if (xcvar.hasRho()) rhoB  = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasRho()) newVarPos = xcvar.getVarPos(RB);
		if (newVarPos >= 0) rhoB  = &temp_var[newVarPos*ng];
	}

	// DRA
	newVarPos = -1;
	if (xcvar.hasGRho()) newVarPos = xcvar.getVarPos(DAX);
	const double* DRA   = NULL;
	if (newVarPos >= 0) DRA   = &temp_var[newVarPos*ng];

	// DRB  
	const double* DRB   = NULL;
	if (isSingleEleSystem) {
		if (xcvar.hasGRho()) DRB  = &temp_betaGVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasGRho()) newVarPos = xcvar.getVarPos(DBX);
		if (newVarPos >= 0) DRB   = &temp_var[newVarPos*ng];
	}

	// TA
	newVarPos = -1;
	if (xcvar.hasTau()) newVarPos = xcvar.getVarPos(TA);
	const double* TA    = NULL;
	if (newVarPos >= 0) TA    = &temp_var[newVarPos*ng];

	// TB
	const double* TB    = NULL;
	if (isSingleEleSystem) {
		if (xcvar.hasTau()) TB    = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasTau()) newVarPos = xcvar.getVarPos(TB);
		if (newVarPos >= 0) TB    = &temp_var[newVarPos*ng];
	}

	// LA
	newVarPos = -1;
	if (xcvar.hasLap()) newVarPos = xcvar.getVarPos(LA);
	const double* LA    = NULL;
	if (newVarPos >= 0) LA    = &temp_var[newVarPos*ng];

	// LB
	const double* LB    = NULL;
	if (isSingleEleSystem) {
		if (xcvar.hasLap()) LB    = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasLap()) newVarPos = xcvar.getVarPos(LB);
		if (newVarPos >= 0) LB    = &temp_var[newVarPos*ng];
	}

	// EXA
	newVarPos = -1;
	if (xcvar.hasExRho()) newVarPos = xcvar.getVarPos(EXA);
	const double* EXRA  = NULL;
	if (newVarPos >= 0) EXRA  = &temp_var[newVarPos*ng];

	// EXB
	const double* EXRB  = NULL;
	if (isSingleEleSystem) {
		if (xcvar.hasExRho()) EXRB = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasExRho()) newVarPos = xcvar.getVarPos(EXB);
		if (newVarPos >= 0) EXRB  = &temp_var[newVarPos*ng];
	}

#else

	const double* rhoA  = oriRhoA;
	const double* rhoB  = oriRhoB;
	const double* DRA   = oriDRA;
	const double* DRB   = oriDRB;
	const double* TA    = oriTA;
	const double* TB    = oriTB;
	const double* LA    = oriLA;
	const double* LB    = oriLB;
	const double* EXRA  = oriEXA;
	const double* EXRB  = oriEXB;

#endif

	// set up possible hirshfeld weights for odd electron population
	// right now it's all zero
	DoublePrecVec hirWeights(ng,ZERO);
	const double* hir  = &hirWeights.front();

	// now let's call fortran function to do the calculation
	Double tol  = xcinfor.getFuncTol();
	Int ng_     = static_cast<Int>(ng);
	double tol_ = static_cast<double>(tol);
	double a1   = xcinfor.getOddElecPar(0);
	double a2   = xcinfor.getOddElecPar(1);
	DoublePrecVec TF(ng);
	becke05_odd_electron(&ng_,&NA,&NB,&tol_,hir,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,EXRA,EXRB,&a1,&a2,&TF[0]);

	// now we do the summation
	Double pop = ZERO;
	const Double* wts = bg.getGridWts();
	for(UInt i=0; i<ng; i++) {
		pop += TF[i]*wts[i];
	}

	// finally let's assign the result to the given center
	UInt parentAtom = bg.AtomInCurrentBatch();
	if (oddElecSum.size()<=parentAtom) {
		string infor = "oddElecSum length is less than the given atom index: " + boost::lexical_cast<string>(parentAtom); 
		Excep excep("ODDElec","doOddElecPopulation",EXCEPTION_VECTOR_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
	oddElecSum[parentAtom] += pop;
}

void ODDElec::add(const XCIntJobInfor& infor, const ODDElec& oddElec) 
{
	// do we just return?
	if (! infor.doOddElec()) {
		return;
	}

	// check the size
	if (oddElecSum.size() != oddElec.oddElecSum.size()) {
		string infor = "oddElecSum length does not equal to the input one"; 
		Excep excep("ODDElec","add",EXCEPTION_VECTOR_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	// now adding
	for(UInt i=0; i<oddElecSum.size(); i++) {
		oddElecSum[i] += oddElec.oddElecSum[i];
	}
}

void ODDElec::printResults(const MolShell& ms) const
{
	cout << "*******************************************" << endl;
	cout << "*    Results of ODD Elec. Calculation     *" << endl;
	cout << "*******************************************" << endl;
	cout << "Odd electron population summary for each atom " << endl;
	printf ("%-12s  %-16s\n", "atom index", "population value");
	UInt nAtoms = ms.getNAtomShells();
	for(UInt oriIndex=0; oriIndex<nAtoms; oriIndex++) {
		printf ("%-12d  %-8.8f\n", (Int)(oriIndex+1), (double)oddElecSum[oriIndex]);
	}
	Double sum = ZERO;
	for(UInt i=0; i<oddElecSum.size(); i++) {
		sum += oddElecSum[i];
	}
	printf ("Odd electron population summary for all of atoms is: %-12.8f\n", sum);
}
