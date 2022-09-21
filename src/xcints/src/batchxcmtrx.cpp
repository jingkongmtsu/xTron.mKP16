/**
 * CPP files corresponding to the batchxcmtrx.h
 * \author  Fenglai Liu and Jing Kong
 */
#include "shell.h"
#include "batchbasis.h"
#include "batchgrid.h"
#include "batchfunc.h"
#include "xcfunc.h"
#include "xcvar.h"
#include "blas.h"
#include "dftderivinfor.h"
#include "halfjkrho.h"
#include "xcintsinfor.h"
#include "sigatombasis.h"
#include "batchxcmtrx.h"
using namespace shell;
using namespace batchbasis;
using namespace batchgrid;
using namespace batchfunc;
using namespace halfjkrho;
using namespace xcfunc;
using namespace xcvar;
using namespace blas;
using namespace xcintsinfor;
using namespace sigatombasis;
using namespace batchxcmtrx;

/***************************************************************************
 *    ####                 BatchXCMtrx  Class                              *
 ***************************************************************************/
BatchXCMtrx::BatchXCMtrx(const SigAtomBasis& sigList, const BatchBasis& basis, 
		const BatchFunc& func, const BatchGrid& grid, 
		const XCVar& xcvar, const XCIntsInfor& infor):DFTMatrix(sigList,sigList,infor.getNSpin())
{
	// first, build the halfFxc result
	UInt nBas   = basis.getNSigBasis();
	UInt nGrids = grid.getNGrids();
	UInt nHalfFxc = 1;
	if (xcvar.hasTau() || xcvar.hasLap()) nHalfFxc += 3;
	MtrxVec halfFxc(nHalfFxc);
	for(UInt i=0; i<nHalfFxc; i++) {
		halfFxc[i].init(nGrids,nBas);
	}

	// work for each spin state
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		if (iSpin == 1) {
			for(UInt i=0; i<nHalfFxc; i++) {
				halfFxc[i].set(ZERO);
			}
		}
		buildHalfFxcMtrx(basis,func,grid,xcvar,iSpin,halfFxc);
		buildFxcMtrx(basis,xcvar,iSpin,halfFxc);
	}
}

void BatchXCMtrx::buildHalfFxcMtrx(const BatchBasis& basis, const BatchFunc& func,
		const BatchGrid& grid, const XCVar& xcvar, const UInt& iSpin, MtrxVec& halfFxc) const 
{
	// prepare work
	UInt nGrids = grid.getNGrids();
	UInt nSigBas= basis.getNSigBasis();

	// firstly, we will deal with rho, gradient rho and lap
	// here we do D1F*Phi
	for(UInt iVar=0; iVar<xcvar.getVarNum(); iVar++) {

		// functional derivatives
		UInt var    = xcvar.getVar(iVar);
		UInt varPos = xcvar.getVarPos(var);
		const Double* D1F = func.getD1F(varPos);

		// checking the spin state
		if (iSpin == 0 && isBetaDFTVar(var)) continue;
		if (iSpin == 1 && ! isBetaDFTVar(var)) continue;

		// rho: F'(r)*phi(r)
		// here ultiply 0.5 is because we will do Fock=Fock+Fock^{T} 
		if(var == RA || var == RB) {
			const Mtrx& phi = basis.getPhi(0);
			Mtrx& Fxc = halfFxc[0];
			for(UInt i=0; i<nSigBas; i++) {
				vaxyaddz(D1F,phi.getPtr(0,i),0.5E0,Fxc.getPtr(0,i),nGrids);
			}
		}

		// gradient rho: F'(r)*phi'(r)
		if(isDRhoVar(var)) {
			UInt deriv = XC_DERIV_X;
			if (var == DAY || var == DBY) deriv = XC_DERIV_Y;
			if (var == DAZ || var == DBZ) deriv = XC_DERIV_Z;
			const Mtrx& phi = basis.getPhi(deriv);
			Mtrx& Fxc = halfFxc[0];
			for(UInt i=0; i<nSigBas; i++) {
				vmuladd(D1F,phi.getPtr(0,i),Fxc.getPtr(0,i),nGrids);
			}
		}

		// lap: F'(r)*phi''(r)(laplacian gradient) 
		if(var == LA || var == LB) {
			for(UInt order=0; order< 3; order++) {
				UInt deriv = XC_DERIV_XX;
				if (order == 1) deriv = XC_DERIV_YY;
				if (order == 2) deriv = XC_DERIV_ZZ;
				const Mtrx& phi = basis.getPhi(deriv);
				Mtrx& Fxc = halfFxc[0];
				for(UInt i=0; i<nSigBas; i++) {
					vmuladd(D1F,phi.getPtr(0,i),Fxc.getPtr(0,i),nGrids);
				}
			}
		}
	}

	// now we multiply the above result with the weights
	const Double* wts = grid.getGridWts();
	Mtrx& Fxc = halfFxc[0];
	for(UInt i=0; i<nSigBas; i++) {
		vmul(wts,Fxc.getPtr(0,i),Fxc.getPtr(0,i),nGrids);
	}

	// secondly, we deal with the tau and lap
	// here we do D1F*Phi'
	if (xcvar.hasTau() || xcvar.hasLap()) {
		for(UInt iVar=0; iVar<xcvar.getVarNum(); iVar++) {

			// get the var
			UInt var = xcvar.getVar(iVar);
			if(iSpin == 0) {
				if (isBetaDFTVar(var)) continue;
			}else{
				if(! isBetaDFTVar(var)) continue;
			}

			// let's do thing, avoid all of other variables except Lap and tau
			if(var != TA && var != TB && var != LA && var != LB) continue;

			// let's do tau and lap
			// the difference is only that for Tau we need to have factor 0.5
			// which same reason for rho
			// for lap factor is set to one
			Double fac = ONE;
			if(var == TA || var == TB) {
				fac = HALF;
			}

			// now let's combine the D1F with the phi'
			UInt varPos = xcvar.getVarPos(var);
			const Double* D1F = func.getD1F(varPos);
			for(UInt order=0; order< N_DERIV_1; order++) {
				UInt deriv = XC_DERIV_ORDER_1[order];
				const Mtrx& phi = basis.getPhi(deriv);
				Mtrx& Fxc = halfFxc[order+1];
				for(UInt i=0; i<nSigBas; i++) {
					vaxyaddz(D1F,phi.getPtr(0,i),fac,Fxc.getPtr(0,i),nGrids);
				}
			}
		}

		// now we multiply the above result with the weights
		const Double* wts = grid.getGridWts();
		for(UInt order=0; order< N_DERIV_1; order++) {
			Mtrx& Fxc = halfFxc[order+1];
			for(UInt i=0; i<nSigBas; i++) {
				vmul(wts,Fxc.getPtr(0,i),Fxc.getPtr(0,i),nGrids);
			}
		}
	}
}

void BatchXCMtrx::buildFxcMtrx(const BatchBasis& basis,const XCVar& xcvar, 
	const UInt& iSpin, const MtrxVec& halfFxc)
{
	// set the result 
	Mtrx& Fxc = getMtrx(iSpin);

	// firstly, combine the first section data in halfFxc with phi value
	const Mtrx& phi = basis.getPhi(0);
	const Mtrx& halfFxc0 = halfFxc[0];
	Fxc.mult(phi,halfFxc0,true,false,ONE,ZERO);

	// secondly, check whether we have second section of halfFxc 
	// for this section, we combine with phi derivatives value
	if (xcvar.hasTau() || xcvar.hasLap()) {
		for(UInt order=0; order< N_DERIV_1; order++) {
			UInt deriv = XC_DERIV_ORDER_1[order];
			const Mtrx& phi  = basis.getPhi(deriv);
			const Mtrx& hFxc = halfFxc[order+1];
			Fxc.mult(phi,hFxc,true,false,ONE,ONE);
		}
	}

	// finally, add its transpose so that to keep fock matrix
	// symmetrized
	// this is performed later after xcints derived the final
	// matrix data, see the function of doMtrx, before updating
	// the global results
	//Fxc.addTranspose();
}

/***************************************************************************
 *    ####                 BatchEXRhoMtrx  Class                           *
 ***************************************************************************/
BatchEXRhoMtrx::BatchEXRhoMtrx(const MolShell& ms, const SigAtomBasis& sigList, 
		const BatchBasis& basis, const HalfJKRho& halfJKRho,
		const BatchFunc& func, const BatchGrid& grid, 
		const XCVar& xcvar, const XCIntJobInfor& infor):DFTMatrix(ms,sigList,infor.getNSpin())
{
	// if no exrho exist, just return
	if (! xcvar.hasExRho()) return;

	// we build the FXC in terms of exchange rho here
	UInt nSigBas= basis.getNSigBasis();
	UInt nGrids = grid.getNGrids();
	Mtrx halfFxc(nGrids,nSigBas);
	const SpinMatrix& halfExRho = halfJKRho.getHalfExRho(); 
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {

		// functional derivatives
		UInt var = EXA;
		if (iSpin == 1) var = EXB;
		UInt varPos = xcvar.getVarPos(var);
		const Double* D1F = func.getD1F(varPos);

		//
		// this part we will combine the functional derivatives
		// with the phi value. Because we do add transpose finally
		// before updating the global result (see xcints.cpp, the 
		// function of doMtrx), therefore in general we need to do
		// halfFxc = F'(r)*phi(r)*0.5
		//
		// However, the exchange energy density in our program
		// is defined as:
		// E_{x}(r) = \sum_{mu,nu,lambda,eta}P_{mu,nu}P_{lambda,eta}
		//            \phi_{mu}(r)\phi_{lambda}(r)\phi_{nu}(r')\phi_{eta}(r')
		// Here it's worthy to note that comparing with exchange energy
		// formula, it does not have 1/2. Therefore, if you sum up
		// all of exchange energy density, it's 2 times bigger than
		// the exchange energy:
		// E_{x} = 1/2*(\sum_{r}w(r)E_{x}(r))
		// w(r) is the weights per grid point.
		//
		// When you do partial derivatives for E_{x}(r) in terms of 
		// P_{mu,nu}, because P_{mu,nu} and P_{lambda,eta} are 
		// symmetrical, it will generate a factor of 2 in the result.
		// This is just like when you do partical derivative for E_{x}
		// in terms of P_{mu,nu}, the exchange Fock matrix part does 
		// not have the 1/2 anymore. 
		//
		// As a result, the 0.5 needed by addTranspose is cancelled
		// with the factor of 2 generating in doing partial derivatives.
		// This is what this long notes about.
		//
		const Mtrx& phi = basis.getPhi(0);
		const Double* wts = grid.getGridWts();
		for(UInt i=0; i<nSigBas; i++) {
			vmul(D1F,phi.getPtr(0,i),halfFxc.getPtr(0,i),nGrids);
			vmul(wts,halfFxc.getPtr(0,i),halfFxc.getPtr(0,i),nGrids);
		}

		// now we need to combine the halfFxc and 
		// the halfExRho together
		const Mtrx& hExRho = halfExRho.getMtrx(iSpin);
		Mtrx& Fxc = getMtrx(iSpin);
		Fxc.mult(hExRho,halfFxc,true,false,ONE,ZERO);
	}
}

/***************************************************************************
 *    ####                 BatchXCVect  Class                              *
 ***************************************************************************/
BatchXCVect::BatchXCVect(const SigAtomBasis& sigList, const BatchBasis& basis, 
		const BatchFunc& func, const BatchGrid& grid, 
		const XCVar& xcvar, const XCIntsInfor& infor):DFTVect(sigList,infor.getNSpin())
{
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		buildXCVect(basis,func,grid,xcvar,iSpin);
	}
}

void BatchXCVect::buildXCVect(const BatchBasis& basis, const BatchFunc& func,
		const BatchGrid& grid, const XCVar& xcvar, const UInt& iSpin) 
{
	// prepare work
	UInt nGrids  = grid.getNGrids();
	UInt nBas    = basis.getNSigBasis();
	UInt offset  = iSpin*nData;
	Double* Fxc  = &vectList[offset];

	// real work
	const Double* wts = grid.getGridWts();
	for(UInt iVar=0; iVar<xcvar.getVarNum(); iVar++) {

		// functional derivatives
		UInt var    = xcvar.getVar(iVar);
		UInt varPos = xcvar.getVarPos(var);
		const Double* D1F = func.getD1F(varPos);

		// rho
		if(var == RA || var == RB) {
			const Mtrx& phi = basis.getPhi(0);
			for(UInt i=0; i<nBas; i++) {
				Fxc[i] += vdot3(wts,D1F,phi.getPtr(0,i),nGrids);
			}
		} else if (isDRhoVar(var)) {
			UInt deriv = XC_DERIV_X;
			if (var == DAY || var == DBY) {
				deriv = XC_DERIV_Y;
			}
			if (var == DAZ || var == DBZ) {
				deriv = XC_DERIV_Z;
			}
			const Mtrx& phi1 = basis.getPhi(deriv);
			for(UInt i=0; i<nBas; i++) {
				Fxc[i] += vdot3(wts,D1F,phi1.getPtr(0,i),nGrids);
			}
		} 
	}
}

