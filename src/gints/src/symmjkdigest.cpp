/**
 * functions and classes associated with symmjkdigest.h
 * \author Fenglai Liu
 */
#include<cmath>
#include "gintsinfor.h"
#include "shellprop.h"
#include "matrix.h"
#include "shell.h"
#include "blas.h"
#include "blas1.h"
#include "shellpair.h"
#include "denmtrxinfor.h"
#include "integraljobs.h"
#include "integralinfor.h"
#include "gintsinfor.h"
#include "symmjkdigest.h"
using namespace gintsinfor;
using namespace shellprop;
using namespace matrix;
using namespace shell;
using namespace blas;
using namespace shellpair;
using namespace denmtrxinfor;
using namespace integraljobs;
using namespace integralinfor;
using namespace gintsinfor;
using namespace symmjkdigest;

SymmJKDigest::SymmJKDigest(const GIntJobInfor& infor, 
		UInt maxNBasis):withLocResult(infor.useBlockResult()),
	isCloseShell(infor.isCloseShell()),nSpin(infor.getNSpin()),
	JThresh(infor.pickUpThresh()),KThresh(infor.pickUpThresh()),
	kCoef(infor.getKCoeff()),job(infor.getIntJob()),order(infor.getJobOrder()),
	mu(0),nu(0),lam(0),eta(0),lmu(0),lnu(0),llam(0),leta(0),
	nmu(0),nnu(0),nlam(0),neta(0),caseNum(illegalKDigestCase),
	aP13(maxNBasis,maxNBasis),aP14(maxNBasis,maxNBasis),
	aP23(maxNBasis,maxNBasis),aP24(maxNBasis,maxNBasis),        
	aK13(maxNBasis,maxNBasis),aK14(maxNBasis,maxNBasis),
	aK23(maxNBasis,maxNBasis),aK24(maxNBasis,maxNBasis),        
	bP13(maxNBasis,maxNBasis),bP14(maxNBasis,maxNBasis),
	bP23(maxNBasis,maxNBasis),bP24(maxNBasis,maxNBasis),        
	bK13(maxNBasis,maxNBasis),bK14(maxNBasis,maxNBasis),
	bK23(maxNBasis,maxNBasis),bK24(maxNBasis,maxNBasis),        
	P12(maxNBasis,maxNBasis),P34(maxNBasis,maxNBasis),
	J12(maxNBasis,maxNBasis),J34(maxNBasis,maxNBasis) { } 

void SymmJKDigest::initShellInfor(bool switchSP, const ShellPair& braSP, const ShellPair& ketSP, 
		const Mtrx& PA, const Mtrx& PB)
{
	// set the shell information first
	if (switchSP) {
		ketSP.getGlobalCarBasOffSet(mu,nu);
		braSP.getGlobalCarBasOffSet(lam,eta);
		ketSP.getLocalCarBasOffSet(lmu,lnu);
		braSP.getLocalCarBasOffSet(llam,leta);
		ketSP.getNCarBas(nmu,nnu);
		braSP.getNCarBas(nlam,neta);
	}else{
		braSP.getGlobalCarBasOffSet(mu,nu);
		ketSP.getGlobalCarBasOffSet(lam,eta);
		braSP.getLocalCarBasOffSet(lmu,lnu);
		ketSP.getLocalCarBasOffSet(llam,leta);
		braSP.getNCarBas(nmu,nnu);
		ketSP.getNCarBas(nlam,neta);
	}

	// now let's determine the case number for K digestion
	if (doK(job) || order > 0) {
		caseNum = illegalKDigestCase;
		if (mu != nu && mu != lam && mu != eta && nu != lam && nu != eta && lam != eta) {
			caseNum = case1234;
		}else if (mu != nu && mu == lam && mu != eta && nu != eta) {
			caseNum = case1213;
		}else if (mu != nu && mu != lam && mu == eta && nu != lam) {
			caseNum = case1231;
		}else if (nu != mu && nu == lam && nu != eta && mu != eta) {
			caseNum = case1223;
		}else if (nu != mu && nu != lam && nu == eta && mu != lam) {
			caseNum = case1232;
		}else if (mu == nu && mu != lam && mu != eta && lam != eta) {
			caseNum = case1123;
		}else if (lam == eta && lam != mu && lam != nu && mu != nu) {
			caseNum = case1233;
		}else if (mu == nu && mu == lam && mu != eta) {
			caseNum = case1112;
		}else if (mu == nu && mu == eta && mu != lam) {
			caseNum = case1121;
		}else if (mu != nu && mu == lam && mu == eta) {
			caseNum = case1211;
		}else if (mu != nu && nu == lam && nu == eta) {
			caseNum = case2111;
		}else if (mu == nu && lam == eta && mu != lam) {
			caseNum = case1122;
		}else if (mu == lam && nu == eta && mu != nu) {
			caseNum = case1212;
		}else if (mu == nu && mu == lam && mu == eta) {
			caseNum = case1111;
		}else {
			string info = "symmetry K digestion fail to find a corresponding case";
			Excep excep("SymmJKDigest","initShellInfor",EXCEPTION_GINTS_ILLEGAL_K_DIGESTION,info);
			handleExcep(excep);
		}
	}

	// let's initialize the density matrix blocks
	bool reform = false;
	if (doK(job)) {
		aP13.reset(nmu,nlam,reform); 
		aP14.reset(nmu,neta,reform);      
		aP23.reset(nnu,nlam,reform);     
		aP24.reset(nnu,neta,reform);    
		aP13.copyFromMatrix(mu,lam,0,0,nmu,nlam,PA);
		aP14.copyFromMatrix(mu,eta,0,0,nmu,neta,PA);
		aP23.copyFromMatrix(nu,lam,0,0,nnu,nlam,PA);
		aP24.copyFromMatrix(nu,eta,0,0,nnu,neta,PA);
		if (nSpin == 2) {
			bP13.reset(nmu,nlam,reform); 
			bP14.reset(nmu,neta,reform);      
			bP23.reset(nnu,nlam,reform);     
			bP24.reset(nnu,neta,reform);    
			bP13.copyFromMatrix(mu,lam,0,0,nmu,nlam,PB);
			bP14.copyFromMatrix(mu,eta,0,0,nmu,neta,PB);
			bP23.copyFromMatrix(nu,lam,0,0,nnu,nlam,PB);
			bP24.copyFromMatrix(nu,eta,0,0,nnu,neta,PB);
		}
	}
	if (doJ(job)) {

		// initialize
		P12.reset(nmu,nnu,reform);    
		P34.reset(nlam,neta,reform);    

		// P12
		for(UInt j=0; j<nnu; j++) {
			vcopy(PA.getPtr(mu,nu+j),P12.getPtr(0,0+j),nmu);
		}
		if (nSpin == 2) {
			for(UInt j=0; j<nnu; j++) {
				vaxpy(PB.getPtr(mu,nu+j),P12.getPtr(0,0+j),ONE,nmu);
			}
		}else if (isCloseShell) {
			vscal(P12.getPtr(),TWO,nmu*nnu);
		}

		// P34
		for(UInt j=0; j<neta; j++) {
			vcopy(PA.getPtr(lam,eta+j),P34.getPtr(0,0+j),nlam);
		}
		if (nSpin == 2) {
			for(UInt j=0; j<neta; j++) {
				vaxpy(PB.getPtr(lam,eta+j),P34.getPtr(0,0+j),ONE,nlam);
			}
		}else if (isCloseShell) {
			vscal(P34.getPtr(),TWO,nlam*neta);
		}
	}

	// reset the block result if needed
	if (order == 0) {

		// let's see which offset we should use
		// the local one or global one
		UInt muIndex  = mu;
		UInt nuIndex  = nu;
		UInt lamIndex = lam;
		UInt etaIndex = eta;
		if (withLocResult) {
			muIndex  = lmu;
			nuIndex  = lnu;
			lamIndex = llam;
			etaIndex = leta;
		}

		// now reset the block matrix result
		if (doK(job)) {
			aK13.init(nmu,nlam,muIndex,lamIndex); 
			aK14.init(nmu,neta,muIndex,etaIndex);      
			aK23.init(nnu,nlam,nuIndex,lamIndex);     
			aK24.init(nnu,neta,nuIndex,etaIndex);    
			if (nSpin == 2) {
				bK13.init(nmu,nlam,muIndex,lamIndex); 
				bK14.init(nmu,neta,muIndex,etaIndex);      
				bK23.init(nnu,nlam,nuIndex,lamIndex);     
				bK24.init(nnu,neta,nuIndex,etaIndex);    
			}
		}
		if (doJ(job)) {
			J12.init(nmu,nnu,muIndex,nuIndex);    
			J34.init(nlam,neta,lamIndex,etaIndex);    
		}
	}
}

void SymmJKDigest::symmJKIntsDigest(const Double* intArray) 
{
	// now set the J/K digestion factor
	Double braFac, ketFac;
	setJDigestFac(braFac,ketFac);
	Double k13Fac,k14Fac,k23Fac,k24Fac;
	setKDigestFac(k13Fac,k14Fac,k23Fac,k24Fac);

	// get the scale coefficients for exchange
	Double coefs = MINUS_ONE;
	if (fabs(kCoef)>ZERO) {
		coefs = coefs*kCoef;
	}

	// now we scale all of K digestion factor with coefs
	k13Fac = k13Fac*coefs;
	k14Fac = k14Fac*coefs;
	k23Fac = k23Fac*coefs;
	k24Fac = k24Fac*coefs;

	// single spin state
	if (nSpin == 1) {
		UInt d=0;
		for(UInt l=0; l<neta; l++) {
			for(UInt k=0; k<nlam; k++) {

				// set for J digestion
				Double sum = ZERO;
				Double p34Val = P34(k,l);

				// loop over bra part
				for(UInt j=0; j<nnu; j++) {
					for(UInt i=0; i<nmu; i++) {
						Double v = intArray[d];

						// exchange digestion
						aK13(i,k) += k13Fac*v*aP24(j,l);  
						aK14(i,l) += k14Fac*v*aP23(j,k);
						aK23(j,k) += k23Fac*v*aP14(i,l);
						aK24(j,l) += k24Fac*v*aP13(i,k);

						// coulomb part
						J12(i,j)  += ketFac*v*p34Val;
						sum       += v*P12(i,j);
						d++;
					}
				}
				J34(k,l) = braFac*sum;
			}
		}
	}else{
		UInt d=0;
		for(UInt l=0; l<neta; l++) {
			for(UInt k=0; k<nlam; k++) {

				// set for J digestion
				Double sum = ZERO;
				Double p34Val = P34(k,l);

				// loop over bra part
				for(UInt j=0; j<nnu; j++) {
					for(UInt i=0; i<nmu; i++) {
						Double v = intArray[d];

						// exchange digestion
						// alpha spin components
						aK13(i,k) += k13Fac*v*aP24(j,l);  
						aK14(i,l) += k14Fac*v*aP23(j,k);
						aK23(j,k) += k23Fac*v*aP14(i,l);
						aK24(j,l) += k24Fac*v*aP13(i,k);
						
						// beta spin components
						bK13(i,k) += k13Fac*v*bP24(j,l);  
						bK14(i,l) += k14Fac*v*bP23(j,k);
						bK23(j,k) += k23Fac*v*bP14(i,l);
						bK24(j,l) += k24Fac*v*bP13(i,k);

						// coulomb part
						J12(i,j)  += ketFac*v*p34Val;
						sum       += v*P12(i,j);
						d++;
					}
				}
				J34(k,l) = braFac*sum;
			}
		}
	}
}

void SymmJKDigest::symmKIntsDigest(const Double* intArray) 
{
	// now set the J/K digestion factor
	Double k13Fac,k14Fac,k23Fac,k24Fac;
	setKDigestFac(k13Fac,k14Fac,k23Fac,k24Fac);

	// get the scale coefficients for exchange
	Double coefs = MINUS_ONE;
	if (fabs(kCoef)>ZERO) {
		coefs = coefs*kCoef;
	}

	// now we scale all of K digestion factor with coefs
	k13Fac = k13Fac*coefs;
	k14Fac = k14Fac*coefs;
	k23Fac = k23Fac*coefs;
	k24Fac = k24Fac*coefs;

	// single spin state
	if (nSpin == 1) {
		UInt d=0;
		for(UInt l=0; l<neta; l++) {
			for(UInt k=0; k<nlam; k++) {
				for(UInt j=0; j<nnu; j++) {
					for(UInt i=0; i<nmu; i++) {
						Double v = intArray[d];
						aK13(i,k) += k13Fac*v*aP24(j,l);  
						aK14(i,l) += k14Fac*v*aP23(j,k);
						aK23(j,k) += k23Fac*v*aP14(i,l);
						aK24(j,l) += k24Fac*v*aP13(i,k);
						d++;
					}
				}
			}
		}
	}else{
		UInt d=0;
		for(UInt l=0; l<neta; l++) {
			for(UInt k=0; k<nlam; k++) {
				for(UInt j=0; j<nnu; j++) {
					for(UInt i=0; i<nmu; i++) {
						Double v = intArray[d];
						aK13(i,k) += k13Fac*v*aP24(j,l);  
						aK14(i,l) += k14Fac*v*aP23(j,k);
						aK23(j,k) += k23Fac*v*aP14(i,l);
						aK24(j,l) += k24Fac*v*aP13(i,k);
						bK13(i,k) += k13Fac*v*bP24(j,l);  
						bK14(i,l) += k14Fac*v*bP23(j,k);
						bK23(j,k) += k23Fac*v*bP14(i,l);
						bK24(j,l) += k24Fac*v*bP13(i,k);
						d++;
					}
				}
			}
		}
	}
}

void SymmJKDigest::symmJIntsDigest(const Double* intArray) 
{
	// now set the J/K digestion factor
	Double braFac, ketFac;
	setJDigestFac(braFac,ketFac);

	// now make J12 and J34
	UInt d=0;
	for(UInt l=0; l<neta; l++) {
		for(UInt k=0; k<nlam; k++) {

			// set for J digestion
			Double sum = ZERO;
			Double p34Val = P34(k,l);

			// loop over bra part
			for(UInt j=0; j<nnu; j++) {
				for(UInt i=0; i<nmu; i++) {
					Double v  = intArray[d];
					J12(i,j) += ketFac*v*p34Val;
					sum      += v*P12(i,j);
					d++;
				}
			}
			J34(k,l) = braFac*sum;
		}
	}
}

void SymmJKDigest::updateFockMatrix(Mtrx& FA, Mtrx& FB) const
{
	// finally digest them into the global results of FA, FB
	if (nSpin == 1) {

		// exchange
		if (doK(job)) {
			FA.updateFromMatrix(0,0,mu,lam,nmu,nlam,ONE,aK13);
			FA.updateFromMatrix(0,0,nu,lam,nnu,nlam,ONE,aK23);
			FA.updateFromMatrix(0,0,mu,eta,nmu,neta,ONE,aK14);
			FA.updateFromMatrix(0,0,nu,eta,nnu,neta,ONE,aK24);
		}

		// coulomb
		if (doJ(job)) {
			FA.updateFromMatrix(0,0,mu,nu,nmu,nnu,ONE,J12);
			FA.updateFromMatrix(0,0,lam,eta,nlam,neta,ONE,J34);
		}
	}else{
		// exchange for alpha and beta
		if (doK(job)) {
			FA.updateFromMatrix(0,0,mu,lam,nmu,nlam,ONE,aK13);
			FA.updateFromMatrix(0,0,nu,lam,nnu,nlam,ONE,aK23);
			FA.updateFromMatrix(0,0,mu,eta,nmu,neta,ONE,aK14);
			FA.updateFromMatrix(0,0,nu,eta,nnu,neta,ONE,aK24);
			FB.updateFromMatrix(0,0,mu,lam,nmu,nlam,ONE,bK13);
			FB.updateFromMatrix(0,0,nu,lam,nnu,nlam,ONE,bK23);
			FB.updateFromMatrix(0,0,mu,eta,nmu,neta,ONE,bK14);
			FB.updateFromMatrix(0,0,nu,eta,nnu,neta,ONE,bK24);
		}

		// coulomb for alpha and beta
		if (doJ(job)) {
			FA.updateFromMatrix(0,0,mu,nu,nmu,nnu,ONE,J12);
			FA.updateFromMatrix(0,0,lam,eta,nlam,neta,ONE,J34);
			FB.updateFromMatrix(0,0,mu,nu,nmu,nnu,ONE,J12);
			FB.updateFromMatrix(0,0,lam,eta,nlam,neta,ONE,J34);
		}
	}
}

///////////////////////////////////////////////////////////////////////////
//            @@@@       digestion for JK derivatives                    //
///////////////////////////////////////////////////////////////////////////
Double SymmJKDerivDigest::symmKIntsDerivDigest(const Double* intArray) const
{
	// get the scale coefficients for exchange
	Double coefs = MINUS_ONE;
	if (fabs(kCoef)>ZERO) {
		coefs = coefs*kCoef;
	}

	// initialize the derivatives value
	Double deriv = ZERO;
	UInt d=0;
	Double v = ZERO;
	if (nSpin == 1) {
		if (isCloseShell) coefs *= TWO; // here we consider the alpha+beta together
		for(UInt l=0; l<neta; l++) {
			for(UInt k=0; k<nlam; k++) {
				for(UInt j=0; j<nnu; j++) {
					for(UInt i=0; i<nmu; i++) {
						v = coefs*intArray[d];
						deriv += v*aP24(j,l)*aP13(i,k);  
						deriv += v*aP23(j,k)*aP14(i,l);
						d++;
					}
				}
			}
		}
	}else{
		for(UInt l=0; l<neta; l++) {
			for(UInt k=0; k<nlam; k++) {
				for(UInt j=0; j<nnu; j++) {
					for(UInt i=0; i<nmu; i++) {
						v = coefs*intArray[d];
						deriv += v*(aP24(j,l)*aP13(i,k)+bP24(j,l)*bP13(i,k));  
						deriv += v*(aP23(j,k)*aP14(i,l)+bP23(j,k)*bP14(i,l));
						d++;
					}
				}
			}
		}
	}
	return deriv;
}

void SymmJKDerivDigest::initDerivInfor(bool switchSP, 
		const ShellPair& braSP, const ShellPair& ketSP, const AtomShellPair& braAtomSP,
		const AtomShellPair& ketAtomSP, const Mtrx& PA, const Mtrx& PB) 
{
	// initialize all of information for symmjkdist object
	initShellInfor(switchSP,braSP,ketSP,PA,PB);

	// initialize the atom index from atom shell pair
	UInt bra1AtomIndex = -1;
	UInt bra2AtomIndex = -1;
	braAtomSP.getAtomShellIndex(bra1AtomIndex,bra2AtomIndex);
	UInt ket1AtomIndex = -1;
	UInt ket2AtomIndex = -1;
	ketAtomSP.getAtomShellIndex(ket1AtomIndex,ket2AtomIndex);

	// let's initialize the data
	// this is for the normal (12|34)
	bra1Pos = bra1AtomIndex;
	bra2Pos = bra2AtomIndex;
	ket1Pos = ket1AtomIndex;
	ket2Pos = ket2AtomIndex;
	if (switchSP) {
		if (braSP.inverseShells()) {
			if (ketSP.inverseShells()) {
				// (43|21)
				bra1Pos = ket2AtomIndex;
				bra2Pos = ket1AtomIndex;
				ket1Pos = bra2AtomIndex;
				ket2Pos = bra1AtomIndex;
			}else{
				// (34|21)
				bra1Pos = ket1AtomIndex;
				bra2Pos = ket2AtomIndex;
				ket1Pos = bra2AtomIndex;
				ket2Pos = bra1AtomIndex;
			}
		}else{
			if (ketSP.inverseShells()) {
				// (43|12)
				bra1Pos = ket2AtomIndex;
				bra2Pos = ket1AtomIndex;
				ket1Pos = bra1AtomIndex;
				ket2Pos = bra2AtomIndex;
			}else{
				// (34|12)
				bra1Pos = ket1AtomIndex;
				bra2Pos = ket2AtomIndex;
				ket1Pos = bra1AtomIndex;
				ket2Pos = bra2AtomIndex;
			}
		}	
	}else{
		if (braSP.inverseShells()) {
			if (ketSP.inverseShells()) {
				// (21|43)
				bra1Pos = bra2AtomIndex;
				bra2Pos = bra1AtomIndex;
				ket1Pos = ket2AtomIndex;
				ket2Pos = ket1AtomIndex;
			}else{
				// (21|34)
				bra1Pos = bra2AtomIndex;
				bra2Pos = bra1AtomIndex;
				ket1Pos = ket1AtomIndex;
				ket2Pos = ket2AtomIndex;
			}
		}else{
			if (ketSP.inverseShells()) {
				// (12|43)
				bra1Pos = bra1AtomIndex;
				bra2Pos = bra2AtomIndex;
				ket1Pos = ket2AtomIndex;
				ket2Pos = ket1AtomIndex;
			}
		}	
	}

	// now for the J, initialize scratch array
	if (doJ(job)) {
		jVec.assign(nmu*nnu,ZERO);
	}

	// initialize the local result
	derivResult.set(ZERO);
}

void SymmJKDerivDigest::findRowColIndex(const ElemDerivInfor& elemDeriv, 
		UInt& localRowIndex, UInt& localColIndex, UInt& globalRowIndex, UInt& globalColIndex) const 
{
	// now let's derive the row and col
	if (order == 1) {

		// get the pos and dir
		UInt pos = elemDeriv.getDerivPos();
		UInt dir = elemDeriv.getDerivDir();

		// set the derivPos and derivDir
		// 0 is for BRA1, 1 is for BRA2 etc. until KET2
		UInt derivPos = 0;
		UInt derivDir = 0; 
		UInt globalDerivPos = bra1Pos;
		if (pos == BRA2) {
			derivPos = 1;
			globalDerivPos = bra2Pos;
		}else if (pos == KET1) {
			derivPos = 2;
			globalDerivPos = ket1Pos;
		}else if (pos == KET2) {
			derivPos = 3;
			globalDerivPos = ket2Pos;
		}
		if (dir == GINTS_DERIV_Y) {
			derivDir = 1;
		}else if (dir == GINTS_DERIV_Z) {
			derivDir = 2;
		}

		// now update row and col index
		localRowIndex  = derivPos;
		localColIndex  = derivDir;
		globalRowIndex = globalDerivPos;
		globalColIndex = derivDir;

	}else if (order == 2) {

		// get the pos and dir
		UInt pos1,pos2,dir1,dir2; 
		elemDeriv.getDerivPos(pos1,pos2);
		elemDeriv.getDerivDir(dir1,dir2);

		// set the derivPos and derivDir
		UInt derivPos1 = 0;
		UInt derivDir1 = 0; 
		UInt derivPos2 = 0;
		UInt derivDir2 = 0; 
		UInt globalDerivPos1 = bra1Pos;
		UInt globalDerivPos2 = bra1Pos;

		// update first deriv infor
		if (pos1 == BRA2) {
			derivPos1 = 1;
			globalDerivPos1 = bra2Pos;
		}else if (pos1 == KET1) {
			derivPos1 = 2;
			globalDerivPos1 = ket1Pos;
		}else if (pos1 == KET2) {
			derivPos1 = 3;
			globalDerivPos1 = ket2Pos;
		}
		if (dir1 == GINTS_DERIV_Y) {
			derivDir1 = 1;
		}else if (dir1 == GINTS_DERIV_Z) {
			derivDir1 = 2;
		}

		// update second deriv infor
		if (pos2 == BRA2) {
			derivPos2 = 1;
			globalDerivPos2 = bra2Pos;
		}else if (pos2 == KET1) {
			derivPos2 = 2;
			globalDerivPos2 = ket1Pos;
		}else if (pos2 == KET2) {
			derivPos2 = 3;
			globalDerivPos2 = ket2Pos;
		}
		if (dir2 == GINTS_DERIV_Y) {
			derivDir2 = 1;
		}else if (dir2 == GINTS_DERIV_Z) {
			derivDir2 = 2;
		}

		// now get the index
		// because row and col are identical,
		// therefore we only stuff the lower triangular part
		localRowIndex  = 3*derivPos1+derivDir1;
		localColIndex  = 3*derivPos2+derivDir2;
		globalRowIndex = 3*globalDerivPos1+derivDir1;
		globalColIndex = 3*globalDerivPos2+derivDir2;
		if (localRowIndex<localColIndex) {
			UInt tmp      = localRowIndex;
			localRowIndex = localColIndex;
			localColIndex = tmp;
		}
		if (globalRowIndex<globalColIndex) {
			UInt tmp       = globalRowIndex;
			globalRowIndex = globalColIndex;
			globalColIndex = tmp;
		}
	}
}

void SymmJKDerivDigest::digestJKInts(const SingleIntegralInfor& infor,
		const Double* rawIntArray, Mtrx& result) 
{
	// set up the scaling coefficients
	// for the first order deriv
	Double kScal = ONE;
	Double jScal = ONE;
	findJKDerivCoef(jScal,kScal);

	// perform digestion for the raw integrals
	UInt intLen  = nmu*nnu*nlam*neta;
	UInt rowABCD = nmu*nnu;
	UInt colABCD = nlam*neta;
	for(UInt iDeriv=0; iDeriv<infor.getElemDerivInforLen(); iDeriv++) {

		// find the row col index
		// check whether this is integral deriv on redundant position?
		UInt rowIndex  = -1;
		UInt colIndex  = -1;
		UInt gRowIndex = -1;
		UInt gColIndex = -1;
		const ElemDerivInfor& intInfor = infor.getElemDerivInfor(iDeriv);
		findRowColIndex(intInfor,rowIndex,colIndex,gRowIndex,gColIndex);
		//printf("the raw integrals %f for deriv %d\n", rawIntArray[iDeriv*intLen], iDeriv); 

		// now exchange
		Double derivK  = ZERO;
		if (doK(job)) {
			derivK = kScal*symmKIntsDerivDigest(&rawIntArray[iDeriv*intLen]);
			//printf("the result without scaling factor %f\n", derivK/kScal); 
			//printf("the raw integrals after into K digestion %f\n", derivK); 
		}

		// coulomb
		// coulomb part is simply P12*Int(nmu*nnu,nlam*neta)*P34
		Double derivJ  = ZERO;
		if(doJ(job)) {

			// digestion for bra part first
			const Double* intArray = &rawIntArray[iDeriv*intLen];
			for(UInt j=0; j<colABCD; j++) {
				Double v = vdot(P12.getPtr(),&intArray[j*rowABCD],rowABCD);
				jVec[j]  = v;
			}

			// now combine the jVec with ket part digestion
			derivJ = jScal*vdot(P34.getPtr(),&jVec.front(),colABCD);
			//printf("the scaling factor %f\n", jScal); 
			//printf("the result without scaling factor %f\n", derivJ/jScal); 
			//printf("the raw integrals after into J digestion %f\n", derivJ); 
		}

		// now set the total result
		Double deriv = derivJ + derivK;
		derivResult(rowIndex,colIndex) += deriv; 

		// update global result
      // for order = 2, we also need to consider the scaling for the 
      // derivatives like (da/dx, b|da/dx,d). We should have two in front of it
      Double scal = ONE;
      if (order == 2 && gRowIndex == gColIndex) {
         // whether the deriv is on different position?
         UInt pos1,pos2; 
         intInfor.getDerivPos(pos1,pos2);
         if (pos1 != pos2) scal = TWO;
      }
		result(gRowIndex,gColIndex) += scal*deriv;	
	}

	// perform digestion for the derived integrals - which are the redundant ones
	for(UInt iDeriv=0; iDeriv<infor.getRedDerivInforLen(); iDeriv++) {

		// get the corresponding redundant integral infor
		const RedundantDerivInfor& redIntInfor = infor.getRedDerivInfor(iDeriv); 

		// find the row col index. Row is the index of atoms.  Col is for the components
		// of derivatives, e.g. xz.
		UInt rowIndex = -1;
		UInt colIndex = -1;
		UInt gRowIndex = -1;
		UInt gColIndex = -1;
		findRowColIndex(redIntInfor,rowIndex,colIndex,gRowIndex,gColIndex);

		// now for this derivatives, loop over rhs of elementary derivatives
		UInt nRHS = redIntInfor.getRHSLen();
		Double deriv = ZERO;
		for(UInt i=0; i<nRHS; i++) {
			Double coef = redIntInfor.getRHSDerivCoef(i);
			const ElemDerivInfor& intInfor = redIntInfor.getRHSDeriv(i);
			UInt iRow  = -1;
			UInt iCol  = -1;
			UInt iGRow = -1;
			UInt iGCol = -1;
			findRowColIndex(intInfor,iRow,iCol,iGRow,iGCol);
			deriv += coef*derivResult(iRow,iCol);
		}
		//printf("the red integrals of %d after the red integral digestion %f\n", iDeriv, deriv); 

		/* for debug purpose
		 * compare with direct calculation
		// reset the redundant integral array
		vset(redIntArray,ZERO,intLen);

		// if the nRHS is 3, we have a faster function to do it
      if (nRHS == 3) {

         // coefficients 
         Double c1 = redIntInfor.getRHSDerivCoef(0);
         Double c2 = redIntInfor.getRHSDerivCoef(1);
         Double c3 = redIntInfor.getRHSDerivCoef(2);

         // integral 
         UInt pos1 = redIntInfor.getRHSDerivPos(0);
         UInt pos2 = redIntInfor.getRHSDerivPos(1);
         UInt pos3 = redIntInfor.getRHSDerivPos(2);
         const Double* p1 = &rawIntArray[pos1*intLen];
         const Double* p2 = &rawIntArray[pos2*intLen];
         const Double* p3 = &rawIntArray[pos3*intLen];

         // now let's do it
         v3axpy(c1,c2,c3,p1,p2,p3,&redIntArray[0],intLen); 
      }else{
         for(UInt i=0; i<nRHS; i++) {
            Double coef = redIntInfor.getRHSDerivCoef(i);
            UInt intPos = redIntInfor.getRHSDerivPos(i);
            vaxpy(&rawIntArray[intPos*intLen],&redIntArray[0],coef,intLen,false); 
         }
      }   
      printf("the red integrals %f for deriv %d\n", redIntArray[0], iDeriv); 

      // now exchange
      // we note that we may need to scale the exchange
      Double derivK  = ZERO;
      if (doK(job)) {
         derivK = kScal*symmKIntsDerivDigest(&redIntArray[0]);
         //printf("the result without scaling factor %f\n", derivK/kScal); 
         printf("the red integrals after into K digestion %f\n", derivK); 
      }   

      // coulomb
      // coulomb part is simply P12*Int(nmu*nnu,nlam*neta)*P34
      Double derivJ  = ZERO;
      if(doJ(job)) {

         // digestion for bra part first
         for(UInt j=0; j<colABCD; j++) {
            Double v = vdot(P12.getPtr(),&redIntArray[j*rowABCD],rowABCD);
            jVec[j]  = v;
         }

         // now combine the jVec with ket part digestion
         derivJ = jScal*vdot(P34.getPtr(),&jVec.front(),colABCD);
         printf("the red integrals after into J digestion %f\n", derivJ); 
      }   

		// compare
		if (fabs(deriv - (derivJ+derivK))>0.000001E0) {
			printf("the red integrals are different between two ways calculation %d\n", iDeriv); 
			printf("first way %f, old way %f\n", deriv, derivJ+derivK); 
		}
		*/

		// we do not need to update the local results
		// we do directly the global result for this case
      // for order = 2, we also need to consider the scaling for the 
      // derivatives like (da/dx, b|da/dx,d). We should have two in front of it
      Double scal = ONE;
      if (order == 2 && gRowIndex == gColIndex) {
         // whether the deriv is on different position?
         UInt pos1,pos2; 
         redIntInfor.getDerivPos(pos1,pos2);
         if (pos1 != pos2) scal = TWO;
      }
		result(gRowIndex,gColIndex) += scal*deriv;	
	}
}
