/**
 * CPP files corresponding to the halfjkrho.h
 * \author  Fenglai Liu
 */
#include <iostream>
#include <cstdio>
#include "excep.h"
#include "xcvar.h"
#include "shellprop.h"
#include "batchgrid.h"
#include "shell.h"
#include "cptrans.h"
#include "blas.h"
#include "matrix.h"
#include "sigatombasis.h"
#include "xcintsinfor.h"
#include "denmtrx.h"
#include "dftmatrix.h"
#include "espints.h"
#include "integraljobs.h"
#include "batchbasis.h"
#include "sigshellpairinfor.h"
#include "halfjkrho.h"
using namespace excep;
using namespace xcvar;
using namespace batchgrid;
using namespace shell;
using namespace cptrans;
using namespace shellprop;
using namespace blas;
using namespace matrix;
using namespace sigatombasis;
using namespace xcintsinfor;
using namespace denmtrx;
using namespace dftmatrix;
using namespace espints;
using namespace integraljobs;
using namespace batchbasis;
using namespace sigshellpairinfor;
using namespace halfjkrho;
using namespace std;

bool HalfJKRho::hasHalfExRho() const 
{
	for(UInt iSpin=0; iSpin<halfExRho.getNSpin(); iSpin++) {
		const Mtrx& hExRho = halfExRho.getMtrx(iSpin);
		if (hExRho.getRow() == 0 || hExRho.getCol() == 0) {
			return false;
		}
	}
	return true;
}

HalfJKRho::HalfJKRho(const XCIntJobInfor& infor, const XCVar& xcvar,
		const BatchGrid& grid, const MolShell& molShl, 
		const SigAtomBasis& sigList, const DoubleVec& cartDen):nGrids(grid.getNGrids()),
	nSigBasis(sigList.getNSigBasis()),halfExRho(infor.getNSpin()) 
{
	// now let's see whether we have the exchange energy density 
	if (infor.doNumericalK() || xcvar.hasExRho()) {
		halfExRho.init(nGrids,molShl.getNBas());
	}

	// whether we have the Coulomb energy density 
	if (infor.doNumericalJ()) {
		halfCouRhoVec.assign(nGrids,ZERO);
		UInt len = cartDen.size(); 
		halfCouRhoVecMatrix.assign(len,ZERO);
	}
}

void HalfJKRho::doHalfJKRho(const XCIntJobInfor& infor, const XCVar& xcvar, const MolShell& ms, 
		const BatchGrid& grid, const SigAtomBasis& sigList, const BatchBasis& basis,
		const DenMtrx& den, const DoubleVec& cartDen, const SigMolShellPairInfor& spData)
{
	// let's check the half JK rho defined before proceed
	// if nothing defined we just returned
	if (! hasHalfExRho() && ! hasHalfCouRho()) {
		return;
	}

	// do we do exchange matrix or producing the exrho?
	SpinMatrix PmunuPhinu(den.getNSpin());
	if (infor.doNumericalK() || xcvar.hasExRho()) {

		// initialize the result
		PmunuPhinu.init(ms.getNCarBas(),nGrids);

		// we need to convert the density matrix into the sig order
		// the result sigDen is in dimension of (nSigBas,nWholeBas)
		DFTMatrix sigDen(sigList,ms,den.getNSpin());
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			sigDen.intoSigOrder(iSpin,den.getMtrx(iSpin),sigList); 
		}

		// now let's appened the C2P transformation as well as basis set
		// scaling into the sigDen
		// after this step sigDen is in (nSigBas,nWholeCartBas)
		const GlobalInfor& ginfor = infor.getGlobalInfor();
		CPTransBasisNorm trans(ginfor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_COL);
		for(UInt iSpin=0; iSpin<sigDen.getNSpin(); iSpin++) {
			Mtrx& D = sigDen.getMtrx(iSpin);
			trans.transform(ms,ms,D);
		}

		// let's make a local copy of Phi0, and transpose it
		// so that the later matrix multiplication can be easier
		// after this step Phi0 is in dimension of (nSigBas,nG)
		Mtrx Phi0(basis.getPhi(0));
		Phi0.transpose(false);

		// now let's combine sigDen, which is Pmunu, with the Phinu value
		// which is the data in  Phi0
		// PmunuPhinu = sigDen*Phi
		// the result dimension is (nWholeCartBas,nGrids)
		for(UInt iSpin=0; iSpin<sigDen.getNSpin(); iSpin++) {
			const Mtrx& D = sigDen.getMtrx(iSpin);
			Mtrx& T = PmunuPhinu.getMtrx(iSpin);
			T.mult(D,Phi0,true,false,ONE,ZERO);
		}

		/*
		 * we abaondoned the use of pMaxInfor
		 * we tested it for ala10 with 6-31G**; it seems that only %20 percentage
		 * of time is saved. this does not worthy the screening effort
		 *
		// form the pMaxInfor
		// this is (nShell,nGrid) 
		// we also count in the effects of the weights so that to 
		// maximum the screening effects
		/
		bool formPmaxInfor = true;
		Mtrx PmaxInfor(ms.getNShell(),nGrids);
		if (formPmaxInfor) {

		// so let's first combine the weights with the PmunuPhinu
		// we do it within a braket because it's only locally used
		// finally make it positive
		SpinMatrix PmunuPhinu2(PmunuPhinu);
		const Double* w = grid.getGridWts(); 
		for(UInt iSpin=0; iSpin<PmunuPhinu2.getNSpin(); iSpin++) {
		Mtrx& T  = PmunuPhinu2.getMtrx(iSpin);
		UInt len = T.getRow(); 
		for(UInt i=0; i<T.getCol(); i++) {
		Double wts = w[i];
		vscal(T.getPtr(0,i),wts,len); 
		vabs(T.getPtr(0,i),len); 
		}
		}

		// now let's try to get PmaxInfor
		// single spin state is easy
		// multiple spin state we need to compare both
		// alpha and beta
		if (PmunuPhinu2.getNSpin() == 1) {
		const Mtrx& T  = PmunuPhinu2.getMtrx(0);

		// loop over grids
		for(UInt iG=0; iG<nGrids; iG++) {
		const Double* t0 = T.getPtr(0,iG);
		Double*       p0 = PmaxInfor.getPtr(0,iG);

		// now let's loop over shell to derive the result
		UInt shellIndex = 0;
		for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {
		const AtomShell& atomShell = ms.getAtomShell(iAtom);
		for(UInt iShell=0; iShell<atomShell.getNShell(); iShell++) {
		const Shell& shell = atomShell.getShell(iShell);
		UInt nCarBas = shell.getNCarBas();
		UInt offset  = shell.getBasisIndex(0,TYPE_CART);
		Double pMax  = ZERO;
		for(UInt j=0; j<nCarBas; j++) {
		Double v = t0[offset+j];
		if (v>pMax) pMax = v;
		}
		p0[shellIndex] = pMax;
		shellIndex++;
		}
		}
		}
		}else{
		const Mtrx& T0  = PmunuPhinu2.getMtrx(0);
		const Mtrx& T1  = PmunuPhinu2.getMtrx(1);

		// loop over grids
		for(UInt iG=0; iG<nGrids; iG++) {
		const Double* t0 = T0.getPtr(0,iG);
		const Double* t1 = T1.getPtr(0,iG);
		Double*       p0 = PmaxInfor.getPtr(0,iG);

		// now let's loop over shell to derive the result
		UInt shellIndex = 0;
		for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {
			const AtomShell& atomShell = ms.getAtomShell(iAtom);
			for(UInt iShell=0; iShell<atomShell.getNShell(); iShell++) {
				const Shell& shell = atomShell.getShell(iShell);
				UInt nCarBas = shell.getNCarBas();
				UInt offset  = shell.getBasisIndex(0,TYPE_CART);
				Double pMax  = ZERO;
				for(UInt j=0; j<nCarBas; j++) {
					Double v0 = t0[offset+j];
					Double v1 = t1[offset+j];
					Double v  = v0 > v1 ? v0 : v1;
					if (v>pMax) pMax = v;
				}
				p0[shellIndex] = pMax;
				shellIndex++;
			}
		}
}
}
}
*/
	}

	// do we do Coulomb matrix?
	// if we do it we need to get the density scaled with weights
	// so that we can use it to form the halfCouRhoVecMatrix
	DoubleVec wRho;
	if (infor.doNumericalJ()) {

		// now let's form the density
		wRho.assign(nGrids,ZERO);

		// firstly convert the density matrix
		// however, we need the spin resolved one
		// so we add alpha and beta together
		DFTMatrix denMtrx(sigList,sigList,1);
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			denMtrx.intoSigOrder(0,den.getMtrx(iSpin),sigList,sigList);
		}

		// set up halfRho
		// P_munu*phi_mu
		Mtrx halfRho(nGrids,nSigBasis);
		const Mtrx& phi0  = basis.getPhi(0);
		const Mtrx& den0  = denMtrx.getMtrx(0);
		bool transposePhi = false;
		bool transposeDen = false;
		halfRho.mult(phi0,den0,transposePhi,transposeDen,ONE,ZERO);

		// now form rho
		for(UInt i=0; i<nSigBasis; i++) {
			vmuladd(phi0.getPtr(0,i),halfRho.getPtr(0,i),&wRho[0],nGrids);
		}

		// debug
		// this should be the real rho
		/*
		for(UInt i=0; i<nGrids; i++) {
			printf("for %d  rho value %-14.8f\n", (Int)i, wRho[i]);
		}
		*/

		// finally we need to conbine the rho with the weights
		const Double* w = grid.getGridWts();
		vmul(w,&wRho[0],&wRho[0],nGrids); 
		
		/*
		for(UInt i=0; i<nGrids; i++) {
			printf("for %d  weight with rho value %-14.8f\n", (Int)i, wRho[i]);
		}
		*/

		// we need to consider the close shell case
		// for alpha and beta are same, we need to scale it with two
		// however, for the Coulomb energy there's also a 1/2 in front
		// considering the two factors:
		// 1 for the close shell and number of spin is 1, the 1/2 
		// cancels the effect of alpha  + beta;
		// 2 for the open shell case(nSpin == 2), or the case that 
		// nSpin == 1 but not close shell case(such as single electron);
		// we need to scale 1/2
		//
		// all in all, for not close shell (case 2), we scale 1/2
		//
		if (! infor.closeShell()) {
			vscal(&wRho[0],HALF,nGrids); 
		}
	}

	// extract the grid points out
	DoubleVec pts(nGrids*3);
	const Double* xyz = grid.getGridCoord(); 
	for(UInt i=0; i<nGrids; i++) {
		pts[3*i  ] = xyz[3*i  ];
		pts[3*i+1] = xyz[3*i+1];
		pts[3*i+2] = xyz[3*i+2];
	}

	// set the job
	UInt job = infor.getIntJob();
	if (isDFTJob(job) && xcvar.hasExRho()) {
		job = EX_DENSITY;
	}

	// get the gints infor
	const GIntsInfor& gintsInfor = infor.getGIntsInfor();

	// now let's form the half exrho/Coulomb rho
	ESPInts excrho(gintsInfor,job);
	excrho.doMtrx(infor,ms,spData,pts,PmunuPhinu,cartDen,wRho,
			halfExRho,halfCouRhoVec,halfCouRhoVecMatrix);

	// debug halfCouRhoVecMatrix
	//printf("the vector matrix\n");
	//for(UInt i=0; i<halfCouRhoVecMatrix.size(); i++) {
	//	if (fabs(halfCouRhoVecMatrix[i])>1.0E-14) {
	//		printf("%i  %-16.10f\n", (Int)i, halfCouRhoVecMatrix[i]);
	//	}
	//}
}

void HalfJKRho::print() const {

	cout << "*******************************************" << endl;
	cout << "*         HalfJKRho  Results              *" << endl;
	cout << "*******************************************" << endl;

	if (hasHalfExRho()) {
		halfExRho.print("half exchange energy density");
	}
	cout << endl;

	if (hasHalfCouRho()) {
		printf("half Coulomb energy density vector\n");
		for(UInt i=0; i<halfCouRhoVec.size(); i++) {
			printf("%i  %-16.10f\n", (Int)i, halfCouRhoVec[i]);
		}
	}
	cout << endl;
}


void HalfJKRho::doHalfJKRhoXHole(const XCIntJobInfor& infor, const XCVar& xcvar, const MolShell& ms, 
		const BatchGrid& grid, const SigAtomBasis& sigList, const BatchBasis& basis,
			    const DenMtrx& den, const DoubleVec& cartDen, const SigMolShellPairInfor& spData, Double sValue)//yw
{
	// let's check the half JK rho defined before proceed
	// if nothing defined we just returned
	if (! hasHalfExRho() && ! hasHalfCouRho()) {
		return;
	}

	// do we do exchange matrix or producing the exrho?
	SpinMatrix PmunuPhinu(den.getNSpin());
	if (infor.doNumericalK() || xcvar.hasExRho()) {

		// initialize the result
		PmunuPhinu.init(ms.getNCarBas(),nGrids);

		// we need to convert the density matrix into the sig order
		// the result sigDen is in dimension of (nSigBas,nWholeBas)
		DFTMatrix sigDen(sigList,ms,den.getNSpin());
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			sigDen.intoSigOrder(iSpin,den.getMtrx(iSpin),sigList); 
		}

		// now let's appened the C2P transformation as well as basis set
		// scaling into the sigDen
		// after this step sigDen is in (nSigBas,nWholeCartBas)
		const GlobalInfor& ginfor = infor.getGlobalInfor();
		CPTransBasisNorm trans(ginfor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_COL);
		for(UInt iSpin=0; iSpin<sigDen.getNSpin(); iSpin++) {
			Mtrx& D = sigDen.getMtrx(iSpin);
			trans.transform(ms,ms,D);
		}

		// let's make a local copy of Phi0, and transpose it
		// so that the later matrix multiplication can be easier
		// after this step Phi0 is in dimension of (nSigBas,nG)
		Mtrx Phi0(basis.getPhi(0));
		Phi0.transpose(false);

		// now let's combine sigDen, which is Pmunu, with the Phinu value
		// which is the data in  Phi0
		// PmunuPhinu = sigDen*Phi
		// the result dimension is (nWholeCartBas,nGrids)
		for(UInt iSpin=0; iSpin<sigDen.getNSpin(); iSpin++) {
			const Mtrx& D = sigDen.getMtrx(iSpin);
			Mtrx& T = PmunuPhinu.getMtrx(iSpin);
			T.mult(D,Phi0,true,false,ONE,ZERO);
		}

	}

	// do we do Coulomb matrix?
	// if we do it we need to get the density scaled with weights
	// so that we can use it to form the halfCouRhoVecMatrix
	DoubleVec wRho;
	if (infor.doNumericalJ()) {

		// now let's form the density
		wRho.assign(nGrids,ZERO);

		// firstly convert the density matrix
		// however, we need the spin resolved one
		// so we add alpha and beta together
		DFTMatrix denMtrx(sigList,sigList,1);
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			denMtrx.intoSigOrder(0,den.getMtrx(iSpin),sigList,sigList);
		}

		// set up halfRho
		// P_munu*phi_mu
		Mtrx halfRho(nGrids,nSigBasis);
		const Mtrx& phi0  = basis.getPhi(0);
		const Mtrx& den0  = denMtrx.getMtrx(0);
		bool transposePhi = false;
		bool transposeDen = false;
		halfRho.mult(phi0,den0,transposePhi,transposeDen,ONE,ZERO);

		// now form rho
		for(UInt i=0; i<nSigBasis; i++) {
			vmuladd(phi0.getPtr(0,i),halfRho.getPtr(0,i),&wRho[0],nGrids);
		}

		// debug
		// this should be the real rho
		/*
		for(UInt i=0; i<nGrids; i++) {
			printf("for %d  rho value %-14.8f\n", (Int)i, wRho[i]);
		}
		*/

		// finally we need to conbine the rho with the weights
		const Double* w = grid.getGridWts();
		vmul(w,&wRho[0],&wRho[0],nGrids); 
		
		/*
		for(UInt i=0; i<nGrids; i++) {
			printf("for %d  weight with rho value %-14.8f\n", (Int)i, wRho[i]);
		}
		*/

		// we need to consider the close shell case
		// for alpha and beta are same, we need to scale it with two
		// however, for the Coulomb energy there's also a 1/2 in front
		// considering the two factors:
		// 1 for the close shell and number of spin is 1, the 1/2 
		// cancels the effect of alpha  + beta;
		// 2 for the open shell case(nSpin == 2), or the case that 
		// nSpin == 1 but not close shell case(such as single electron);
		// we need to scale 1/2
		//
		// all in all, for not close shell (case 2), we scale 1/2
		//
		if (! infor.closeShell()) {
			vscal(&wRho[0],HALF,nGrids); 
		}
	}

	// extract the grid points out
	DoubleVec pts(nGrids*3);
	const Double* xyz = grid.getGridCoord(); 
	for(UInt i=0; i<nGrids; i++) {
		pts[3*i  ] = xyz[3*i  ];
		pts[3*i+1] = xyz[3*i+1];
		pts[3*i+2] = xyz[3*i+2];
	}

	// set the job
	UInt job = infor.getIntJob();
	if (isDFTJob(job) && xcvar.hasExRho()) {
		job = EX_DENSITY;
	}

	// get the gints infor
	const GIntsInfor& gintsInfor = infor.getGIntsInfor();

	// now let's form the half exrho/Coulomb rho
	ESPInts excrho(gintsInfor,job);
	excrho.doMtrxXHole(infor,ms,spData,pts,PmunuPhinu,cartDen,wRho,
		      halfExRho,halfCouRhoVec,halfCouRhoVecMatrix, sValue);//yw

	// debug halfCouRhoVecMatrix
	//printf("the vector matrix\n");
	//for(UInt i=0; i<halfCouRhoVecMatrix.size(); i++) {
	//	if (fabs(halfCouRhoVecMatrix[i])>1.0E-14) {
	//		printf("%i  %-16.10f\n", (Int)i, halfCouRhoVecMatrix[i]);
	//	}
	//}
}

