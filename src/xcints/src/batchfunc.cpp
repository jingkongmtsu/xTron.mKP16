/**
 * CPP files corresponding to the batchfunc.h
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include "excep.h"
#include "shell.h"
#include "xcintsinfor.h"
#include "xcfunc.h"
#include "batchvar.h"
#include "xcvar.h"
#include "batchgrid.h"
#include "functionallist.h"
#include "sigatombasis.h"
#include "matrix.h"
#include "blas.h"
#include "blas1.h"
#include "denmtrx.h"
#include "hirshfeld.h"
#include "xcenergyinfor.h"
#include "batchfunc.h"
using namespace excep;
using namespace shell;
using namespace xcintsinfor;
using namespace sigatombasis;
using namespace batchvar;
using namespace xcfunc;
using namespace xcvar;
using namespace batchgrid;
using namespace matrix;
using namespace blas;
using namespace denmtrx;
using namespace hirshfeld; 
using namespace xcenergyinfor; 
using namespace batchfunc;
using namespace std;

BatchFunc::BatchFunc(const BatchVar& bvar, const DenMtrx& den, 
		const XCFunc& xcfunc, const XCVar& xcvar, 
		const XCIntJobInfor& infor):ng(bvar.getNGrids()),nDen(1),nAlpha(den.getNEle(0)),
	nBeta(den.getNEle(1)),nEXFunc(0),nECFunc(0),
	nFunc(xcfunc.getNFuncValues()),nD1FWithGamma(xcvar.getNumVarWithGamma()),
	nD1F(xcvar.getVarNum()),tol(infor.getFuncTol()),hirWeights(ng,ZERO),funvals(ng*nFunc,ZERO),d1F(ng,nD1F)
{
	// currently we do not suppor non-linear functional
	if (! xcfunc.isLinear()) {
		string infor = "so far the support for non-linear functional is not available yet";
		Excep excep("BatchFunc","constructor",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
		handleExcep(excep);
	}

	// now let's see how many exchange/correlation functional
	if (xcfunc.hasCorrelation()) nECFunc = 1;
	nEXFunc = nFunc - nECFunc;

	//
	// determine the number of density for open shell case
	//
	// how many density we should actually use in evaluating the functional 
	// derivatives? From one hand, we could refer to the situation that whether
	// we have the "beta" variable. However, even if this is open shell
	// calculation, it's possible that both of two alpha and beta densities
	// are still same. therefore, the actual running result for functional
	// derivatives is equivalent to the close shell case.
	// 
	// Therefore, we do not use beta variable idea. Instead, we go to compare 
	// the difference between alpha rho and beta rho - since for all of 
	// functionals we always has rho as variable, and the difference for the 
	// variable, is absolutely from the input density matrix.
	// 
	// if RA = RB, no matter open/close shell cases, we will make nDen = 1
	// in the functional calculation. Otherwise nDen = 2.
	//
	if (xcvar.hasBetaVar()) {

		// get the density
		UInt varPos = xcvar.getVarPos(RA);
		const Double* rhoA  = bvar.getVar(varPos);
		varPos = xcvar.getVarPos(RB);
		const Double* rhoB  = bvar.getVar(varPos);

		// define the absolute difference between the rhoA and rhoB
		Double MAE = ZERO;
		for(UInt i=0; i<ng; i++) {
			MAE += fabs(rhoA[i]-rhoB[i]);
		}
		MAE /= ng;

		// now let's see the number of density
		nDen = 2;
		if (MAE < tol) nDen = 1;
	}

	// an exceptional case: for single electron case,
	// the nDen should be two, where one state DFT variable
	// are all zero
	if (den.isSingleEleSystem()) {
		nDen = 2;
	}
}

void BatchFunc::getBatchExc(const BatchGrid& grid, Double& exc) const
{
	exc = ZERO;
	const Double* wts = grid.getGridWts();
	for(UInt iFunc=0; iFunc<nFunc; iFunc++) {
		const Double* funval = &funvals[iFunc*ng];
		exc += vdot(funval,wts,ng);
	}
}

void BatchFunc::doFuncDeriv(const BatchVar& bvar, const XCFunc& xcfunc, 
		const XCVar& xcvar, XCEnergyInfor& xcEList) 
{
	// set up tmp array - since in most cases we need to scale
	// the functional derivatives result, therefore we use 
	// tmp array to hold each calculation result
	// we note that this tmp array must be in double precision
	UInt nTempFunc = (nEXFunc>=nECFunc ? nEXFunc : nECFunc);
	UInt tmpDataLength = (nD1FWithGamma+nTempFunc)*ng;
	DoublePrecVec TF_array(tmpDataLength);
	double* TF   = &TF_array.front();
	double* TD1F = TF + nTempFunc*ng;

	// functional derivatives array in terms of gamma
	UInt gFlen = nD1FWithGamma*ng;
	DoubleVec gF(gFlen,ZERO);

	// set up the infor array used in functional calculation
	Int infor[MAX_DFTVAR_TYPES];
	xcvar.setupVarInfor(infor);

	// Grab all the variables to be passed to fortran routines.
	const Double *rhoA(0), *rhoB(0), *DRA(0), *DRB(0), *TA(0), *TB(0), *LA(0), *LB(0), 
	              *EXRA(0), *EXRB(0); 
	bvar.varForFort(rhoA, rhoB, DRA, DRB, TA, TB, LA, LB, EXRA, EXRB, xcvar);

	//Hirshfeld weights.
	const double* hirw  = &hirWeights.front();

	// next some treatments to the input integers etc.
	// the things input into fortran must be signed integer
	Int ng_     = static_cast<Int>(ng);
	Int nDen_   = static_cast<Int>(nDen);
	double tol_ = static_cast<double>(tol);

	// setting for the TPSS correlation functional
	// for ordinary case set use_lambda form functional as 0
	Int use_lambda = 0;

	// now let's go to the real calculation part
	// firstly it's exchange
	if (nEXFunc > 0) {

		for(UInt i=0; i<xcfunc.getLenEx(); i++) {

			// prepare
			const string& exname = xcfunc.getEXName(i);
			if (exname == "HF") continue;
			TF_array.assign(tmpDataLength,ZERO);

			// go to each pieces
			if (exname == "SLATER") {
				functional_s(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,TF,TD1F);
			}else if (exname == "PW86X") {
				functional_pw86x(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TF,TD1F);
			}else if (exname == "BECKE88") {
				functional_becke88(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TF,TD1F);
			}else if (exname == "PBEX") {
				functional_pbex(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TF,TD1F);
			}else if (exname == "PW91X") {
				functional_pw91x(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TF,TD1F);
			}else if (exname == "TPSSX") {
				functional_tpssx(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TA,TB,TF,TD1F);
			}else if (exname == "SCANX") {
				functional_scanx(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TA,TB,TF,TD1F);
			}else if (exname == "VDWBR89") {
				br89b2_vdwx(TF,TD1F,infor,rhoA,rhoB,DRA,DRB,LA,LB,TA,TB,&ng_,&nDen_);
			}else if (exname == "BR89X") {
				br89xx(TF,TD1F,infor,rhoA,rhoB,DRA,DRB,LA,LB,TA,TB,&ng_,&nDen_);
			}else if (exname == "B05_NDOP") {

				// let's set the values P 
				double p = static_cast<double>(xcfunc.getBecke05PVal());

				// now call the functional
				becke05_ndop(infor,&nDen_,&ng_,&tol_,&p,hirw,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,EXRA,EXRB,TF,TD1F);
			}else if (exname == "B05_NDPAR") {

				// let's set the values P and Q
				double p = static_cast<double>(xcfunc.getBecke05PVal());
				double q = static_cast<double>(xcfunc.getBecke05QVal());

				// NA and NB
				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
				becke05_ndpar(infor,&nDen_,&ng_,&NA,&NB,&tol_,&p,&q,hirw,
						rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,EXRA,EXRB,TF,TD1F);
			}else{
				string infor = "this exchange functional is not available yet: " + exname;
				Excep excep("BatchFunc","doFuncDeriv",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
				handleExcep(excep);
			}

			// do we do energy profiling?
			if (xcEList.doXCEnergyProfiling()) {
				xcEList.updateEComponent(TF,exname);
			}

			// scale the functional derivatives in case that it's go with tau/exrho with 1/2 
			// factor
			scalFuncDeriv(exname,xcvar,TD1F);

			// scale the functional and its derivatives
			Double excoe  = xcfunc.getEXFactor(i);
			Double* F = &funvals.front();
			for(UInt ip=0; ip<nEXFunc*ng; ip++) {
				F[ip] += excoe*TF[ip];
			}
			Double* D1F = &gF.front();
			for(UInt ip=0; ip<nD1FWithGamma*ng; ip++) {
				D1F[ip] += excoe*TD1F[ip];
			}
		}
	}

	// second it's correlation
	if (nECFunc > 0) {

		for(UInt i=0; i<xcfunc.getLenEc(); i++) {

			// prepare
			const string& ecname = xcfunc.getECName(i);
			TF_array.assign(tmpDataLength,ZERO);

			// real calculation
			if (ecname == "PW92C") {
				functional_pw92c(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,TF,TD1F);
			}else if (ecname == "VWN5") {
				functional_vwn5(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,TF,TD1F);
			}else if (ecname == "VWN1RPA") {
				functional_vwn1rpa(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,TF,TD1F);
			}else if (ecname == "PW86C") {

				// for pw86x, things is different
				// we use GAA etc. array directly
				DoublePrecVec GAA(ng);
				DoublePrecVec GAB(ng);
				DoublePrecVec GBB(ng);
				const double* DRAX = DRA;
				const double* DRAY = DRA+ng;
				const double* DRAZ = DRA+2*ng;
				const double* DRBX = DRB;
				const double* DRBY = DRB+ng;
				const double* DRBZ = DRB+2*ng;
				for(UInt ip=0; ip<ng; ip++) {
					GAA[ip] = DRAX[ip]*DRAX[ip]+DRAY[ip]*DRAY[ip]+DRAZ[ip]*DRAZ[ip];
					GAB[ip] = DRAX[ip]*DRBX[ip]+DRAY[ip]*DRBY[ip]+DRAZ[ip]*DRBZ[ip];
					GBB[ip] = DRBX[ip]*DRBX[ip]+DRBY[ip]*DRBY[ip]+DRBZ[ip]*DRBZ[ip];
				}
				functional_pw86c(infor,&ng_,&nDen_,rhoA,rhoB,&GAA.front(),
						&GAB.front(),&GBB.front(),TF,TD1F);

			}else if (ecname == "PW91C") {
				functional_pw91c(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TF,TD1F);
			}else if (ecname == "PBEC") {
				functional_pbec(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TF,TD1F);
			}else if (ecname == "LYP") {
				functional_lyp(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TF,TD1F);
			}else if (ecname == "TPSSC") {
				functional_tpssc(infor,&ng_,&nDen_,&tol_,&use_lambda,rhoA,rhoB,DRA,DRB,TA,TB,TF,TD1F);
			}else if (ecname == "SCANC") {
				functional_scanc(infor,&ng_,&nDen_,&tol_,&use_lambda,rhoA,rhoB,DRA,DRB,TA,TB,TF,TD1F);
			}else if (ecname == "PSTS_ND") {
				functional_psts(infor,&ng_,&nDen_,&tol_,rhoA,rhoB,DRA,DRB,TA,TB,EXRA,EXRB,TF,TD1F);
			}else if (ecname == "VDWBR89") {
				br89b2_vdwx(TF,TD1F,infor,rhoA,rhoB,DRA,DRB,LA,LB,TA,TB,&ng_,&nDen_);
			}else if (ecname == "BR94COOR_PAR") {
                                double p = static_cast<double>(xcfunc.getBecke05PVal());
                                double q = static_cast<double>(xcfunc.getBecke05QVal());

				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
			br94coor_par(infor,&nDen_,&NA,&NB,&ng_,&tol_,&p,&q,hirw,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,EXRA,EXRB,TF,TD1F);
			}else if (ecname == "BR94COOR_OP") {
                                double p = static_cast<double>(xcfunc.getBecke05PVal());

				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
			br94coor_op(infor,&nDen_,&NA,&NB,&ng_,&tol_,&p,hirw,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,EXRA,EXRB,TF,TD1F);
			}else if (ecname == "B13COOR_PAR") {

				// let's set the values P and Q
				double p = static_cast<double>(xcfunc.getBecke05PVal());
                                double q = static_cast<double>(xcfunc.getBecke05QVal());
				Int method = static_cast<Int>(xcfunc.getB13CorrMethod());

				// NA and NB
				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
				b13coor_par(&method,infor,&nDen_,&NA,&NB,&ng_,&tol_,&p,&q,hirw,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,
						EXRA,EXRB,TF,TD1F);
			}else if (ecname == "B13COOR_OPP") {

				// let's set the values P and Q
				double p = static_cast<double>(xcfunc.getBecke05PVal());
				Int method = static_cast<Int>(xcfunc.getB13CorrMethod());

				// NA and NB
				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
				b13coor_opp(&method,infor,&nDen_,&NA,&NB,&ng_,&tol_,&p,hirw,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,
						EXRA,EXRB,TF,TD1F);
			}else if (ecname == "B13STRONG_AC2") {
				// let's set the values P and Q
				double p = static_cast<double>(xcfunc.getBecke05PVal());
				double q = static_cast<double>(xcfunc.getBecke05QVal());
				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
				b13strong_ac2(infor,&nDen_,&ng_,&NA,&NB,&tol_,&p,&q,hirw,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,
						EXRA,EXRB,TF,TD1F);

				// see the find difference for strong component
				//double delta_   = 0.01E0;
				//double criteria = 0.00001E0;
				//find_diff_check_deriv1_array(infor,&ng_,&tol_,&delta_,
				//		&criteria,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,EXRA,EXRB);
			}else if (ecname == "B13STRONG_AC3") {
				// let's set the values P and Q
				double p = static_cast<double>(xcfunc.getBecke05PVal());
				double q = static_cast<double>(xcfunc.getBecke05QVal());
				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
				b13strong_ac3(infor,&nDen_,&ng_,&NA,&NB,&tol_,&p,&q,hirw,rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,
						EXRA,EXRB,TF,TD1F);
			}else if (ecname == "KP14C") {


				// KP14 depends on three parameters
				Int NA = static_cast<Int>(nAlpha);
				Int NB = static_cast<Int>(nBeta);
				double b = static_cast<double>(xcfunc.getKP14BVal());
				double alpha = static_cast<double>(xcfunc.getKP14AlphaVal());
				double ndparCoef = static_cast<double>(xcfunc.getKP14CNDPARVal());
				double ndparCap = static_cast<double>(xcfunc.getKP14CNDPARCapVal());
				double p = static_cast<double>(xcfunc.getBecke05PVal());
				double q = static_cast<double>(xcfunc.getBecke05QVal());
				Int ndparMethod = static_cast<Int>(xcfunc.getB05NDParMethod());

				// set up an array for collecting the energy decomposition results
				DoublePrecVec TF1_array(2*ng);
				double* TF1   = &TF1_array.front();

				// because it uses the b13 correlation functional
				// so it also needs the NA and NB
				kp14ec(&ndparMethod,infor,&nDen_,&ng_,&NA,&NB,&b,&alpha,&ndparCoef,&ndparCap,&tol_,&p,&q,hirw,
				       rhoA,rhoB,DRA,DRB,TA,TB,LA,LB,EXRA,EXRB,TF,TF1,TD1F);

				// possible energy profiling code for KP14C functional
				if (xcEList.doXCEnergyProfiling()) {
					xcEList.updateEComponent(TF1,ecname,1);
					xcEList.updateEComponent(&TF1[ng],ecname,2);
				}
			}else{
				string infor = "this correlation functional is not available yet: " + ecname;
				Excep excep("BatchFunc","doFuncDeriv",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
				handleExcep(excep);
			}

			// do we do energy profiling?
			if (xcEList.doXCEnergyProfiling() && ecname != "KP14C") {
				xcEList.updateEComponent(TF,ecname);
			}

			// scale the functional derivatives in case that it's go with tau/exrho with 1/2 
			// factor
			scalFuncDeriv(ecname,xcvar,TD1F);

			// scale the functional and its derivatives
			Double eccoe  = xcfunc.getECFactor(i);
			Double* F = &funvals.front() + nEXFunc*ng;
			for(UInt ip=0; ip<nECFunc*ng; ip++) {
				F[ip] += eccoe*TF[ip];
			}
			Double* D1F = &gF.front();
			for(UInt ip=0; ip<nD1FWithGamma*ng; ip++) {
				D1F[ip] += eccoe*TD1F[ip];
			}
		}
	}

	// finally, we do post work
	const Double* D1F = &gF.front();
	transformGammaDerivIntoDRhoDeriv(D1F,xcfunc,xcvar,bvar);
}

void BatchFunc::scalFuncDeriv(const string& name, const XCVar& xcvar, Double* TD1F) const
{
	// let's see whether the functional needs to be scaled
	// because the tau/exrho uses the factor of 1/2
	bool forTau   = false;
	bool forExRho = false;
	if (name == "PSTS_ND") {
		forExRho = true;
	}

	// if we do not do any scaling work, just return
	if (!forExRho && !forTau) return;

	// let's scale the tau variable
	// just in case the tau variable is going with 1/2
	if (forTau) {
		for(UInt i=0; i<xcvar.getNumVarWithGamma(); i++) {
			bool gammaVar = true;
			UInt var = xcvar.getVar(i,gammaVar);
			if (var == TA || var == TB) {
				UInt varPos = xcvar.getVarPos(var,true);
				vscal(&TD1F[varPos*ng],HALF,ng);
			}
		}
	}

	// similarly, let's scale the exchange energy density
	// just in case the variable is going with 1/2
	if (forExRho) {
		for(UInt i=0; i<xcvar.getNumVarWithGamma(); i++) {
			bool gammaVar = true;
			UInt var = xcvar.getVar(i,gammaVar);
			if (var == EXA || var == EXB) {
				UInt varPos = xcvar.getVarPos(var,true);
				vscal(&TD1F[varPos*ng],HALF,ng);
			}
		}
	}
}

Double BatchFunc::b13StrongFreeAtomGammaBRVal(const UInt& Z) const 
{
	// choose the value
	// this is derived from the table IV for the given B13Strong paper
	switch(Z) {
		case 1:
			return ZERO;
		case 2:
			return -0.627E0;
		case 3:
			return -0.157E0;
		case 4:
			return -0.225E0;
		case 5:
			return -0.12E0;
		case 6:
			return -0.081E0;
		case 7:
			return -0.066E0;
		case 8:
			return -0.161E0;
		case 9: 
			return -0.207E0;
		case 10:
			return -0.249E0;
		case 11:
			return -0.152E0;
		case 12:
			return -0.123E0;
		case 13:
			return -0.067E0;
		case 14:
			return -0.022E0;
		case 15:
			return  0.016E0;
		case 16:
			return  0.034E0;
		case 17:
			return  0.051E0;
		case 18:
			return  0.066E0;
		default:
			break;
	}

	// we should not be here
	// so just raise up error infor
	// and return
	cout << "the atomic number passed in is: " << Z << endl;
	string infor = "only atomic number 1-18 has the gamma BR adjustment parameter so far";
	Excep excep("BatchFunc","b13StrongFreeAtomGammaBRVal",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
	handleExcep(excep);
	return ZERO;
}

void BatchFunc::formBR89EXHoleWeights(const MolShell& ms, const SigAtomBasis& sigAtomBasis, 
		const Hirshfeld& hirshfeld) 
{
	// this is to form the weights modified Br89 exchange hole 
	// here the implementation is just following the paper below:
	// J. Chem. Phys. 138, 074109 (2013)
	// Density functionals for static, dynamical, and strong correlation
	// Axel D. Becke

	// firstly let's see whether this is free atom?
	// for free atom, according to the equation C1 in appendix C
	if (ms.getNAtomShells() == 1) {

		// get the atomic number for this atom
		// because it's only one atom, therefore
		const AtomShell& atom = ms.getAtomShell(0);
		UInt Z = atom.getAtomic();

		// now get the weights
		// and assign
		Double g = b13StrongFreeAtomGammaBRVal(Z);
		hirWeights.assign(hirWeights.size(),g);
	}else{

		// initialize the weights
		hirWeights.assign(hirWeights.size(),ZERO);

		// now let's see sig atom list
		const UIntVec& sigAtomList = sigAtomBasis.getSigAtoms();
		for(UInt iAtom=0; iAtom<sigAtomBasis.getNSigAtoms(); iAtom++) {

			// get the gamma BR value for this sig atom
			UInt atomIndex = sigAtomList[iAtom];
			const AtomShell& as = ms.getAtomShell(atomIndex);
			UInt Z = as.getAtomic();
			Double g = b13StrongFreeAtomGammaBRVal(Z);

			// now let's combine with hirshfeld weights
			const Double* hirshfeldWeights = hirshfeld.getWts(iAtom);
			for(UInt iG=0; iG<hirWeights.size(); iG++) {
				hirWeights[iG] += g*hirshfeldWeights[iG];
			}
		}

		// debug
		// let's see the weights value
		//printf("debugging for forming BR89 hole weights after combined with Hirshfeld weights\n");
		//for(UInt i=0; i<hirWeights.size(); i++) {
		//	printf("for grid %4d weights is %-12.7f\n", (Int)i, hirWeights[i]);
		//}
	}
}

void BatchFunc::transformGammaDerivIntoDRhoDeriv(const Double* TD1F,
		const XCFunc& xcfunc, const XCVar& xcvar, const BatchVar& bvar)
{
	for(UInt i=0; i<xcvar.getNumVarWithGamma(); i++) {
		bool gammaVar = true;
		UInt var = xcvar.getVar(i,gammaVar);
		if (xcvar.neglectVar(var)) continue;
		if (!isGammaVar(var)) {
			UInt oldVarPos = xcvar.getVarPos(var,true);
			UInt newVarPos = xcvar.getVarPos(var);
			vcopy(&TD1F[oldVarPos*ng],d1F.getPtr(0,newVarPos),ng);
		}else{
			for(UInt j=0; j<xcvar.getVarNum(); j++) {
				UInt newVar = xcvar.getVar(j);
				if(isDRhoVar(newVar) && doesDRhoVarMatchGammaVar(var,newVar)) {
					Double fac = getFacGammaIntoDRho(var);
					UInt oldVarPos = xcvar.getVarPos(var,true);
					UInt newVarPos = xcvar.getVarPos(newVar);
					UInt newVar2   = getVarGammaIntoDRho(var, newVar);
					UInt newVarPos2= xcvar.getVarPos(newVar2);
					const Double* v =  bvar.getVar(newVarPos2);
					vaxyaddz(&TD1F[oldVarPos*ng],v,fac,d1F.getPtr(0,newVarPos),ng);
				}
			}
		}
	}
}

void BatchFunc::print() {

	cout << "*******************************************" << endl;
	cout << "*                BatchFunc                *" << endl;
	cout << "*******************************************" << endl;
	Mtrx excMtrx(ng,nEXFunc+nECFunc);
	vcopy(&funvals.front(), excMtrx.getPtr(), (nEXFunc+nECFunc)*ng);
	excMtrx.print("functional values for exchange and correlation");
	d1F.print("1st functional derivatives");
}

