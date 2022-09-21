/**
 * CPP files corresponding to the batchvar.h
 * \author  Fenglai Liu and Jing Kong
 */
#include<iostream>
#include "shell.h"
#include "xcintsinfor.h"
#include "batchbasis.h"
#include "batchgrid.h"
#include "xcvar.h"
#include "integraljobs.h"
#include "halfjkrho.h"
#include "dftderivinfor.h"
#include "blas.h"
#include "blas1.h"
#include "denmtrx.h"
#include "dftmatrix.h"
#include "batchvar.h"
using namespace shell;
using namespace xcintsinfor;
using namespace batchbasis;
using namespace batchgrid;
using namespace xcvar;
using namespace integraljobs;
using namespace halfjkrho;
using namespace blas;
using namespace denmtrx;
using namespace dftmatrix;
using namespace batchvar;
using namespace std;

void BatchVar::buildHalfVar(const BatchBasis& basis, const Mtrx& denM, 
		const XCVar& xcvar, MtrxVec& halfVar) const 
{
	// P_munu*phi_mu
	// this is default and we always do that
	const Mtrx& phi0 = basis.getPhi(0);
	Mtrx& halfVar0   = halfVar[0];
	bool transposePhi = false;
	bool transposeDen = false;
	halfVar0.mult(phi0,denM,transposePhi,transposeDen,ONE,ZERO);

	// P_munu*phi_mu_x' etc.
	bool doVar1 = xcvar.do1stGradHalfVar();
	if (doVar1) {
		for(UInt order=0; order<N_DERIV_1; order++) {
			UInt deriv = XC_DERIV_ORDER_1[order];
			const Mtrx& phi1 = basis.getPhi(deriv);
			Mtrx& halfVar1   = halfVar[order+1];
			halfVar1.mult(phi1,denM,transposePhi,transposeDen,ONE,ZERO);
		}
	}

	// P_munu*phi_mu_x'x' etc.
	bool doVar2 = xcvar.do2edGradHalfVar();
	if (doVar2) {

		// compute the offset for halfvar
		UInt offset = 1;
		if (doVar1) {
			offset = N_DERIV_1+1;
		}

		// combination
		for(UInt order=0; order<N_DERIV_2; order++) {
			UInt deriv = XC_DERIV_ORDER_2[order];
			const Mtrx& phi2 = basis.getPhi(deriv);
			Mtrx& halfVar2   = halfVar[order+offset];
			halfVar2.mult(phi2,denM,transposePhi,transposeDen,ONE,ZERO);
		}
	}
}

void BatchVar::buildVar(const BatchBasis& basis, const MtrxVec& halfVar, 
					const UInt& iSpin, const XCVar& xcvar) 
{
	// variable with order 0
	for(UInt ivar=0; ivar<xcvar.getVarNum(); ivar++) {

		// variable information
		UInt var    = xcvar.getVar(ivar);
		UInt varPos = xcvar.getVarPos(var);

		// determine whether we do alpha var or beta var?
		if (iSpin == 0 && isBetaDFTVar(var)) continue;
		if (iSpin == 1 && ! isBetaDFTVar(var)) continue;
		
		//rho
		if (var == RA || var == RB) {
			const Mtrx& phi   = basis.getPhi(0);
			const Mtrx& halfV = halfVar[0];
			Double* Rho       = &vars[varPos*nGrids];
			for(UInt i=0; i<nSigBasis; i++) {
				vmuladd(phi.getPtr(0,i),halfV.getPtr(0,i),Rho,nGrids);
			}
		}

		//gradient rho
		//order for half var is in X, Y, Z
		// Drho_x = D(rho)_x = 2*sum_munu P_munu (D(phi_mu)_x*phi_nu)
		if (isDRhoVar(var)) {

			// determine the var
			UInt deriv = XC_DERIV_X;
			if (var == DAY || var == DBY) {
				deriv = XC_DERIV_Y;
			} else if (var == DAZ || var == DBZ) {
				deriv = XC_DERIV_Z;
			}

			// now let's do it
			const Mtrx& phi   = basis.getPhi(deriv);
			const Mtrx& halfV = halfVar[0];
			Double* GRho      = &vars[varPos*nGrids];
			for(UInt i=0; i<nSigBasis; i++) {
				vmuladd(phi.getPtr(0,i),halfV.getPtr(0,i),GRho,nGrids);
			}
			vscal(GRho,TWO,nGrids);
		}

		//tau 
		if (var == TA || var == TB) {
			for(UInt dorder=0; dorder<N_DERIV_1; dorder++) {
				UInt deriv = XC_DERIV_ORDER_1[dorder];
				const Mtrx& phi   = basis.getPhi(deriv);
				const Mtrx& halfV = halfVar[dorder+1];
				Double* Tau       = &vars[varPos*nGrids];
				for(UInt i=0; i<nSigBasis; i++) {
					vmuladd(phi.getPtr(0,i),halfV.getPtr(0,i),Tau,nGrids);
				}
			}
		}

		//lap
		if (var == LA || var == LB) {

			// copy the tau variable if we have it
			if (xcvar.hasTau()) {
				UInt v = TA;
				if (var == LB) v = TB;
				UInt newVarPos = xcvar.getVarPos(v);
				Double* Tau = &vars[newVarPos*nGrids];
				Double* Lap = &vars[varPos*nGrids];
				vcopy(Tau,Lap,nGrids);
			}else{
				// copy the tau code above to pretend that we are 
				// calculating tau variable - but actually this is 
				// for Lap
				for(UInt dorder=0; dorder<N_DERIV_1; dorder++) {
					UInt deriv = XC_DERIV_ORDER_1[dorder];
					const Mtrx& phi   = basis.getPhi(deriv);
					const Mtrx& halfV = halfVar[dorder+1];
					Double* Tau       = &vars[varPos*nGrids];
					for(UInt i=0; i<nSigBasis; i++) {
						vmuladd(phi.getPtr(0,i),halfV.getPtr(0,i),Tau,nGrids);
					}
				}
			}

			// now we add another component - D2(phi)*phi
			// they are from XX, YY and ZZ
			Double* Lap = &vars[varPos*nGrids];
			for(UInt i=0; i<3; i++) {

				// set the deriv
				UInt deriv = XC_DERIV_XX;
				if (i == 1) deriv = XC_DERIV_YY;
				if (i == 2) deriv = XC_DERIV_ZZ;

				// now let's do the combination
				const Mtrx& phi2  = basis.getPhi(deriv);
				const Mtrx& halfV = halfVar[0];
				for(UInt i=0; i<nSigBasis; i++) {
					vmuladd(phi2.getPtr(0,i),halfV.getPtr(0,i),Lap,nGrids);
				}
			}

			// finally, scale the lap
			vscal(Lap,TWO,nGrids);
		}
	}	
}

void BatchVar::buildExRho(const MolShell& ms, const SigAtomBasis& sigList, 
		const DenMtrx& denxxx, const BatchBasis& basis, const HalfJKRho& halfJKRho,
		const XCVar& xcvar) 
{
	// prepare the density matrix in sig order
	// only one dimension is changed into sig order
	DFTMatrix den(ms,sigList,denMtrx.getNSpin());
	for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
		den.intoSigOrder(iSpin,denMtrx.getMtrx(iSpin),sigList);
	}

	// prepare the half var
	Mtrx halfVar(nGrids,nSigBasis);

	// now let's begin our work
	for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {

		// combine density matrix with halfExRho
		const SpinMatrix& hExRho = halfJKRho.getHalfExRho();
		bool transpose = false;
		halfVar.mult(hExRho.getMtrx(iSpin),den.getMtrx(iSpin),transpose,transpose,ONE,ZERO);

		// now combine with phi
		for(UInt ivar=0; ivar<xcvar.getVarNum(); ivar++) {

			// variable information
			// we find ex var
			UInt var    = xcvar.getVar(ivar);
			if (var != EXA && var != EXB) continue;
			UInt varPos = xcvar.getVarPos(var);

			// determine whether we do alpha var or beta var?
			if (iSpin == 0 && isBetaDFTVar(var)) continue;
			if (iSpin == 1 && ! isBetaDFTVar(var)) continue;
		
			// digestion
			const Mtrx& phi = basis.getPhi(0);
			Double* EXRho   = &vars[varPos*nGrids];
			for(UInt i=0; i<nSigBasis; i++) {
				vmuladd(phi.getPtr(0,i),halfVar.getPtr(0,i),EXRho,nGrids);
			}
		}
	}
}

BatchVar::BatchVar(const SigAtomBasis& sigList, const BatchBasis& basis, const XCVar& xcvar, 
		const DenMtrx& denM):nGrids(basis.getNGrids()),
	nSigBasis(basis.getNSigBasis()),vars(xcvar.getVarNum()*nGrids,ZERO), denMtrx(denM)
{
	// pre-build the halfVar array
	// in default we always has rho combined
	UInt halfVarSize = 1;
	if (xcvar.do1stGradHalfVar()) halfVarSize += N_DERIV_1;
	if (xcvar.do2edGradHalfVar()) halfVarSize += N_DERIV_2;
	MtrxVec halfVar(halfVarSize);
	for(UInt i=0; i<halfVarSize; i++) {
		halfVar[i].init(nGrids,nSigBasis);
	}

	// form significant order density matrix
	// this is used for all of variables calculation
	// 
	// for the matrix whose two dimensions are both 
	// in sig order, it's used everythere excep
	// exchange rho 
	//
	// therefore, exchange rho will be calculated 
	// later and indepently
	DFTMatrix den(sigList,sigList,denMtrx.getNSpin());
	for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
		den.intoSigOrder(iSpin,denMtrx.getMtrx(iSpin),sigList,sigList);
	}

	// build variable
	UInt nSpin = den.getNSpin();
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// build half var
		buildHalfVar(basis,den.getMtrx(iSpin),xcvar,halfVar);

		// build the whole variable
		buildVar(basis,halfVar,iSpin,xcvar);
	}
}

BatchVar::BatchVar(const BatchBasis& basis, const XCVar& xcvar, 
		const DFTMatrix& den, const DenMtrx& denM):nGrids(basis.getNGrids()),
	nSigBasis(basis.getNSigBasis()),vars(xcvar.getVarNum()*nGrids,ZERO),
  denMtrx(denM)
{
	// pre-build the halfVar array
	// in default we always has rho combined
	UInt halfVarSize = 1;
	if (xcvar.do1stGradHalfVar()) halfVarSize += N_DERIV_1;
	if (xcvar.do2edGradHalfVar()) halfVarSize += N_DERIV_2;
	MtrxVec halfVar(halfVarSize);
	for(UInt i=0; i<halfVarSize; i++) {
		halfVar[i].init(nGrids,nSigBasis);
	}

	// build variable
	UInt nSpin = den.getNSpin();
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// build half var
		buildHalfVar(basis,den.getMtrx(iSpin),xcvar,halfVar);

		// build the whole variable
		buildVar(basis,halfVar,iSpin,xcvar);
	}
}

void BatchVar::buildRIVar(const BatchBasis& basis, const XCVar& xcvar, 
		const Double* den, const UInt& iSpin) 
{
	for(UInt ivar=0; ivar<xcvar.getVarNum(); ivar++) {

		// variable information
		UInt var    = xcvar.getVar(ivar);
		UInt varPos = xcvar.getVarPos(var);

		// determine whether we do alpha var or beta var?
		if (iSpin == 0 && isBetaDFTVar(var)) continue;
		if (iSpin == 1 && ! isBetaDFTVar(var)) continue;

		// rho
		if (var == RA || var == RB) {
			const Mtrx& phi = basis.getPhi(0);
			Double* Rho = &vars[varPos*nGrids];
			mmul(phi.getPtr(),den,Rho,nGrids,1,nSigBasis,ONE,ZERO);
		}

		// gradient rho
		if (isDRhoVar(var)) {

			// determine the var
			UInt deriv = XC_DERIV_X;
			if (var == DAY || var == DBY) {
				deriv = XC_DERIV_Y;
			} else if (var == DAZ || var == DBZ) {
				deriv = XC_DERIV_Z;
			}
			const Mtrx& phi1  = basis.getPhi(deriv);
			Double* GRho = &vars[varPos*nGrids];
			mmul(phi1.getPtr(),den,GRho,nGrids,1,nSigBasis,ONE,ZERO);
		}
	}
}

void BatchVar::smoothDFTVar(const XCVar& xcvar) 
{

	// currently we just set all of negative density to ZERO
	for(UInt ivar=0; ivar<xcvar.getVarNum(); ivar++) {
		UInt var    = xcvar.getVar(ivar);
		UInt varPos = xcvar.getVarPos(var);
		if (var == RA || var == RB) {
			Double* Rho = &vars[varPos*nGrids];
			for(UInt i=0; i<nGrids; i++) {
				if(Rho[i] < 0) Rho[i] = ZERO; 
			}
		}
	}
}

BatchVar::BatchVar(const BatchBasis& basis, const XCVar& xcvar, 
		const DFTVect& vec, const DenMtrx& denM):nGrids(basis.getNGrids()),
	nSigBasis(basis.getNSigBasis()),vars(xcvar.getVarNum()*nGrids,ZERO),
  denMtrx(denM)
{
	// build variable
	UInt nSpin = vec.getNSpin();;
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		const Double* denVec = vec.getConstVect(iSpin);
		buildRIVar(basis,xcvar,denVec,iSpin);
	}

	// we need to smooth it after produce it 
	// since RI density may be less than 0
	smoothDFTVar(xcvar);
}

BatchVar::BatchVar(const MolShell& ms, const SigAtomBasis& sigList, 
		const BatchBasis& basis, const HalfJKRho& halfJKRho,
		const XCIntJobInfor& infor, const DenMtrx& denM):nGrids(basis.getNGrids()),
	nSigBasis(basis.getNSigBasis()),vars(infor.getJKVarNum()*nGrids,ZERO),
  denMtrx(denM)
{
	// firstly let's see whether we have exchange energy density
	UInt intJob = infor.getIntJob();
	if (infor.doNumericalK()) {

		// prepare the density matrix in sig order
		// only one dimension is changed into sig order
		DFTMatrix den(ms,sigList,denMtrx.getNSpin());
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			den.intoSigOrder(iSpin,denMtrx.getMtrx(iSpin),sigList);
		}

		// prepare the half var
		Mtrx halfVar(nGrids,nSigBasis);

		// now let's begin our work
		for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {

			// combine density matrix with halfExRho
			const SpinMatrix& hExRho = halfJKRho.getHalfExRho();
			bool transpose = false;
			halfVar.mult(hExRho.getMtrx(iSpin),den.getMtrx(iSpin),transpose,transpose,ONE,ZERO);

			// digestion
			UInt varPos     = getKVarPos(intJob,iSpin);
			const Mtrx& phi = basis.getPhi(0);
			Double* EXRho   = &vars[varPos*nGrids];
			for(UInt i=0; i<nSigBasis; i++) {
				vmuladd(phi.getPtr(0,i),halfVar.getPtr(0,i),EXRho,nGrids);
			}
		}
	}

	// now let's see Coulomb energy variable
	if (infor.doNumericalJ()) {

		// firstly convert the density matrix
		// both alpha and beta are together
		DFTMatrix den(sigList,sigList,1);
		for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {
			den.intoSigOrder(0,denMtrx.getMtrx(iSpin),sigList,sigList);
		}

		// set up halfRho
		Mtrx halfRho(nGrids,nSigBasis);

		// P_munu*phi_mu
		const Mtrx& phi0  = basis.getPhi(0);
		const Mtrx& den0  = den.getMtrx(0);
		bool transposePhi = false;
		bool transposeDen = false;
		halfRho.mult(phi0,den0,transposePhi,transposeDen,ONE,ZERO);

		// now form rho
		DoubleVec rho(nGrids,ZERO);
		for(UInt i=0; i<nSigBasis; i++) {
			vmuladd(phi0.getPtr(0,i),halfRho.getPtr(0,i),&rho[0],nGrids);
		}

		// now let's combine the result
		UInt pos   = getJVarPos(intJob,infor.getNSpin());
		Double* eJ = &vars[pos*nGrids];
		const DoubleVec& halfEJ = halfJKRho.getHalfCouRhoVec();
		for(UInt i=0; i<nGrids; i++) {
			eJ[i] += halfEJ[i]*rho[i];
		}
	}
}

Double BatchVar::getNEle(const UInt& iSpin, const XCVar& xcvar,
		const BatchGrid& bg) const
{
	// get the variable position
	UInt v = RA;
	if (iSpin == 1) v = RB;
	UInt varPos = xcvar.getVarPos(v);

	// estimate the number of electrons
	const Double* var = getVar(varPos);
	const Double* wts = bg.getGridWts();
	Double nEle = vdot(var,wts,nGrids);
	return nEle;
}

Double BatchVar::getEX(const UInt& iSpin, const XCVar& xcvar,
		const BatchGrid& bg) const
{
	// get the variable position
	UInt v = EXA;
	if (iSpin == 1) v = EXB;
	UInt varPos = xcvar.getVarPos(v);

	// estimate the energy
	// we note, energy is with 1/2
	const Double* var = getVar(varPos);
	const Double* wts = bg.getGridWts();
	Double ex = HALF*vdot(var,wts,nGrids);
	return ex;
}

Double BatchVar::getEX(const UInt& iSpin, const XCIntJobInfor& infor, const BatchGrid& bg) const
{
	// estimate the energy
	// we note, energy is with 1/2
	UInt intJob       = infor.getIntJob();
	UInt varPos       = getKVarPos(intJob,iSpin);
	const Double* var = getVar(varPos);
	const Double* wts = bg.getGridWts();
	Double ex = HALF*vdot(var,wts,nGrids);
	return ex;
}

Double BatchVar::getEJ(const XCIntJobInfor& infor, const BatchGrid& bg) const
{
	// estimate the energy
	// we note, energy is with 1/2
	UInt intJob       = infor.getIntJob();
	UInt varPos       = getJVarPos(intJob,infor.getNSpin());
	const Double* var = getVar(varPos);
	const Double* wts = bg.getGridWts();
	Double ec = HALF*vdot(var,wts,nGrids);
	return ec;
}

void BatchVar::print(const XCVar& xcvar) const
{

	cout << "*******************************************" << endl;
	cout << "*                 BatchVar                *" << endl;
	cout << "*******************************************" << endl;

	// print each variable we support
	if (xcvar.hasRho()) {
		UInt pos = xcvar.getVarPos(RA);
		Mtrx alphaRho(nGrids,1);
		vcopy(&vars[pos*nGrids],alphaRho.getPtr(),nGrids);
		alphaRho.print("Alpha Rho");
		UInt pos2 = xcvar.getVarPos(RB);
		if (pos2 != pos) {
			Mtrx betaRho(nGrids,1);
			vcopy(&vars[pos2*nGrids],betaRho.getPtr(),nGrids);
			betaRho.print("Beta Rho");
		}
	}

	if (xcvar.hasGRho()) {
		UInt pos = xcvar.getVarPos(DAX);
		Mtrx alphaRho(nGrids,3);
		vcopy(&vars[pos*nGrids],alphaRho.getPtr(),3*nGrids);
		alphaRho.print("DAX, DAY and DAZ");
		UInt pos2 = xcvar.getVarPos(DBX);
		if (pos2 != pos) {
			Mtrx betaRho(nGrids,3);
			vcopy(&vars[pos2*nGrids],betaRho.getPtr(),3*nGrids);
			betaRho.print("DBX, DBY and DBZ");
		}
	}

	if (xcvar.hasTau()) {
		UInt pos = xcvar.getVarPos(TA);
		Mtrx alphaRho(nGrids,1);
		vcopy(&vars[pos*nGrids],alphaRho.getPtr(),nGrids);
		alphaRho.print("Alpha Tau");
		UInt pos2 = xcvar.getVarPos(TB);
		if (pos2 != pos) {
			Mtrx betaRho(nGrids,1);
			vcopy(&vars[pos2*nGrids],betaRho.getPtr(),nGrids);
			betaRho.print("Beta Tau");
		}
	}

	if (xcvar.hasLap()) {
		UInt pos = xcvar.getVarPos(LA);
		Mtrx alphaRho(nGrids,1);
		vcopy(&vars[pos*nGrids],alphaRho.getPtr(),nGrids);
		alphaRho.print("Alpha Laplacian");
		UInt pos2 = xcvar.getVarPos(LB);
		if (pos2 != pos) {
			Mtrx betaRho(nGrids,1);
			vcopy(&vars[pos2*nGrids],betaRho.getPtr(),nGrids);
			betaRho.print("Beta Laplacian");
		}
	}

	if (xcvar.hasExRho()) {
		UInt pos = xcvar.getVarPos(EXA);
		Mtrx alphaRho(nGrids,1);
		vcopy(&vars[pos*nGrids],alphaRho.getPtr(),nGrids);
		alphaRho.print("Alpha Exchange Energy Density");
		UInt pos2 = xcvar.getVarPos(EXB);
		if (pos2 != pos) {
			Mtrx betaRho(nGrids,1);
			vcopy(&vars[pos2*nGrids],betaRho.getPtr(),nGrids);
			betaRho.print("Beta Exchange Energy Density");
		}
	}
}

void BatchVar::print(const XCIntJobInfor& infor) const
{
	cout << "*******************************************" << endl;
	cout << "*                 BatchVar                *" << endl;
	cout << "*******************************************" << endl;

	UInt intJob = infor.getIntJob();

	// whether we do the Coulomb
	if (doJ(intJob)) {
		UInt pos = getJVarPos(intJob,infor.getNSpin());
		Mtrx JRho(nGrids,1);
		vcopy(&vars[pos*nGrids],JRho.getPtr(),nGrids);
		JRho.print("Coulomb energy density");
	}

	// whether we do the exchange
	if (doK(intJob)) {

		// alpha 
		UInt pos = getKVarPos(intJob,0);
		Mtrx KRho(nGrids,1);
		vcopy(&vars[pos*nGrids],KRho.getPtr(),nGrids);
		KRho.print("Alpha Exchange energy density");
		
		// beta
		if (infor.getNSpin() == 2) {
			pos = getKVarPos(intJob,1);
			vcopy(&vars[pos*nGrids],KRho.getPtr(),nGrids);
			KRho.print("Beta Exchange energy density");
		}
	}
}

void BatchVar::xhole(const XCVar& xcvar, vector<DoubleVec>& xholeVal) const
{
	UInt nSpin = xholeVal.size();
	vector<UInt> posEx(nSpin);
	posEx[0] = xcvar.getVarPos(EXA);
	posEx[nSpin-1] = xcvar.getVarPos(EXB);
	vector<UInt> posDen(nSpin);
	posDen[0] = xcvar.getVarPos(RA);
	posDen[nSpin-1] = xcvar.getVarPos(RB);
	for (UInt is = 0; is < nSpin; is++)
	{
		//rho is exchange energy density!
    Mtrx rho(nGrids,1);
    vcopy(&vars[posEx[is]*nGrids],rho.getPtr(),nGrids);
    //rho.print("Exchange Energy Density");
    Mtrx den(nGrids,1);
    vcopy(&vars[posDen[is]*nGrids],den.getPtr(),nGrids);
    //den.print("density");
		for ( UInt ig = 0; ig < nGrids; ig++ )
		{
			//cout<<"THE NGRID IS "<<nGrids << endl;
			if ( den(ig,0) > 1E-10 )
			{
				//xholeVal[is][ig] = rho(ig,0); //For comparing with x density.
				xholeVal[is][ig] = rho(ig,0)/den(ig,0);  //For hole calculation.
			}
			else
			{
				xholeVal[is][ig] =0;
			}
		  //cout <<"alphaRho = "<<alphaRho(ig,0)<<" alphaRhoDen = "<<alphaRhoDen(ig,0)<<endl;
		  //cout <<"the ["<<ig<<"] xhole in batchvar is"<<xholeVal[ig]<<endl;
  	}
	}
}


void BatchVar::varForFort(const Double*& rhoA, const Double*& rhoB, 
     const Double*& DRA, const Double*& DRB, const Double*& TA, 
     const Double*& TB, const Double*& LA, const Double*& LB, 
     const Double*& EXRA, const Double*& EXRB, const XCVar& xcvar) const
{
	UInt ng = nGrids;
	// get the density variables for alpha state
	UInt varPos = -1;
	if (xcvar.hasRho()) varPos = xcvar.getVarPos(xcvarinfor::RA);
	const Double* oriRhoA  = getVar(varPos);
	varPos = -1;
	if (xcvar.hasGRho()) varPos = xcvar.getVarPos(xcvarinfor::DAX);
	const Double* oriDRA   = getVar(varPos);
	varPos = -1;
	if (xcvar.hasTau()) varPos = xcvar.getVarPos(xcvarinfor::TA);
	const Double* oriTA    = getVar(varPos);
	varPos = -1;
	if (xcvar.hasLap()) varPos = xcvar.getVarPos(xcvarinfor::LA);
	const Double* oriLA    = getVar(varPos);
	varPos = -1;
	if (xcvar.hasExRho()) varPos = xcvar.getVarPos(xcvarinfor::EXA);
	const Double* oriEXA   = getVar(varPos);

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
	if (denMtrx.isSingleEleSystem()) {
		betaVar.assign(ng,ZERO);
		betaGVar.assign(3*ng,ZERO);
	}
	
	// now get the beta variables
	const Double* oriRhoB  = NULL;
	const Double* oriDRB   = NULL;
	const Double* oriTB    = NULL;
	const Double* oriLB    = NULL;
	const Double* oriEXB   = NULL;
	if (denMtrx.isSingleEleSystem()) {
		if (xcvar.hasRho()) oriRhoB  = &betaVar[0];
		if (xcvar.hasGRho()) oriDRB  = &betaGVar[0];
		if (xcvar.hasTau()) oriTB    = &betaVar[0];
		if (xcvar.hasLap()) oriLB    = &betaVar[0];
		if (xcvar.hasExRho()) oriEXB = &betaVar[0];
	}else{
		varPos = -1;
		if (xcvar.hasRho()) varPos = xcvar.getVarPos(xcvarinfor::RB);
		oriRhoB  = getVar(varPos);
		varPos = -1;
		if (xcvar.hasGRho()) varPos = xcvar.getVarPos(xcvarinfor::DBX);
		oriDRB   = getVar(varPos);
		varPos = -1;
		if (xcvar.hasTau()) varPos = xcvar.getVarPos(xcvarinfor::TB);
		oriTB    = getVar(varPos);
		varPos = -1;
		if (xcvar.hasLap()) varPos = xcvar.getVarPos(xcvarinfor::LB);
		oriLB    = getVar(varPos);
		varPos = -1;
		if (xcvar.hasExRho()) varPos = xcvar.getVarPos(xcvarinfor::EXB);
		oriEXB   = getVar(varPos);
	}

	// functional derivatives array in terms of gamma
	// UInt gFlen = nD1FWithGamma*ng;
	// DoubleVec gF(gFlen,ZERO);

	// set up tmp array - since in most cases we need to scale
	// the functional derivatives result, therefore we use 
	// tmp array to hold each calculation result
	// we note that this tmp array must be in double precision
	// UInt nTempFunc = (nEXFunc>=nECFunc ? nEXFunc : nECFunc);
	// UInt tmpDataLength = (nD1FWithGamma+nTempFunc)*ng;
	// DoublePrecVec TF_array(tmpDataLength);
	// double* TF   = &TF_array.front();
	// double* TD1F = TF + nTempFunc*ng;

	// set up the infor array used in functional calculation
	Int infor[MAX_DFTVAR_TYPES];
	xcvar.setupVarInfor(infor);

	// we need to perform some treatment to the inputs 
	// variables first
#ifdef WITH_SINGLE_PRECISION

	// we need to copy the floating point variables into double precision ones
	// by getting the bvar's first posision, we will copy all of data from
	// single precision to double precision
	UInt len = xcvar.getVarNum()*bg.getNGrids();
	DoublePrecVec temp_var(len);
	const Double* var0 = getVarPos(0);
	for(UInt i=0; i<len; i++) {
		temp_var[i] = static_cast<double>(var0[i]);
	}
	
	// set up single electron case as shown above
	DoublePrecVec temp_betaVar;
	DoublePrecVec temp_betaGVar;
	if (denMtrx.isSingleEleSystem()) {
		temp_betaVar.assign(ng,ZERO);
		temp_betaGVar.assign(3*ng,ZERO);
	}

	// now fetch the variable data
	// alpha rho 
	UInt newVarPos = -1;
	if (xcvar.hasRho()) newVarPos = xcvar.getVarPos(RA);
	rhoA  = NULL;
	if (newVarPos >= 0) rhoA  = &temp_var[newVarPos*ng];

	// beta rho
	rhoB  = NULL;
	if (denMtrx.isSingleEleSystem()) {
		if (xcvar.hasRho()) rhoB  = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasRho()) newVarPos = xcvar.getVarPos(RB);
		if (newVarPos >= 0) rhoB  = &temp_var[newVarPos*ng];
	}

	// DRA
	newVarPos = -1;
	if (xcvar.hasGRho()) newVarPos = xcvar.getVarPos(DAX);
	DRA   = NULL;
	if (newVarPos >= 0) DRA   = &temp_var[newVarPos*ng];

	// DRB  
	DRB   = NULL;
	if (denMtrx.isSingleEleSystem()) {
		if (xcvar.hasGRho()) DRB  = &temp_betaGVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasGRho()) newVarPos = xcvar.getVarPos(DBX);
		if (newVarPos >= 0) DRB   = &temp_var[newVarPos*ng];
	}

	// TA
	newVarPos = -1;
	if (xcvar.hasTau()) newVarPos = xcvar.getVarPos(TA);
	TA    = NULL;
	if (newVarPos >= 0) TA    = &temp_var[newVarPos*ng];

	// TB
	TB    = NULL;
	if (denMtrx.isSingleEleSystem()) {
		if (xcvar.hasTau()) TB    = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasTau()) newVarPos = xcvar.getVarPos(TB);
		if (newVarPos >= 0) TB    = &temp_var[newVarPos*ng];
	}

	// LA
	newVarPos = -1;
	if (xcvar.hasLap()) newVarPos = xcvar.getVarPos(LA);
	LA    = NULL;
	if (newVarPos >= 0) LA    = &temp_var[newVarPos*ng];

	// LB
	LB    = NULL;
	if (denMtrx.isSingleEleSystem()) {
		if (xcvar.hasLap()) LB    = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasLap()) newVarPos = xcvar.getVarPos(LB);
		if (newVarPos >= 0) LB    = &temp_var[newVarPos*ng];
	}

	// EXA
	newVarPos = -1;
	if (xcvar.hasExRho()) newVarPos = xcvar.getVarPos(EXA);
	EXRA  = NULL;
	if (newVarPos >= 0) EXRA  = &temp_var[newVarPos*ng];

	// EXB
	EXRB  = NULL;
	if (denMtrx.isSingleEleSystem()) {
		if (xcvar.hasExRho()) EXRB = &temp_betaVar[0];
	}else{
		newVarPos = -1;
		if (xcvar.hasExRho()) newVarPos = xcvar.getVarPos(EXB);
		if (newVarPos >= 0) EXRB  = &temp_var[newVarPos*ng];
	}

	// possible hirshfeld weights
	DoublePrecVec temp_hirWeights(ng);
	for(UInt i=0; i<ng; i++) temp_hirWeights[i] = hirWeights[i];
	const double* hirw  = &temp_hirWeights.front();

#else

	rhoA  = oriRhoA;
	rhoB  = oriRhoB;
	DRA   = oriDRA;
	DRB   = oriDRB;
	TA    = oriTA;
	TB    = oriTB;
	LA    = oriLA;
	LB    = oriLB;
	EXRA  = oriEXA;
	EXRB  = oriEXB;

#endif
}

