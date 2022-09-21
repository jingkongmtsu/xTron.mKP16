/**
 * functions and classes associated with gints4deriv.h
 * working class to handle the 4D analytical integrals derivatives calculation
 * \author Fenglai Liu 
 */
#include<cstdio>
#include<cmath>
#include "blas.h"
#include "blas1.h"
#include "excep.h"
#include "shell.h"
#include "matrix.h"
#include "denmtrx.h"
#include "molecule.h"
#include "shellpair.h"
#include "shellprop.h"
#include "spinmatrix.h"
#include "gintsinfor.h"
#include "hgp_os_ints.h"
#include "localmemscr.h"
#include "integraljobs.h"
#include "denmtrxinfor.h"
#include "symmjkdigest.h"
#include "sigshellpairinfor.h"
#include "gints4deriv.h"
using namespace blas;
using namespace excep;
using namespace shell;
using namespace matrix;
using namespace denmtrx;
using namespace molecule;
using namespace shellpair;
using namespace shellprop;
using namespace spinmatrix;
using namespace gintsinfor;
using namespace localmemscr;
using namespace integraljobs;
using namespace denmtrxinfor;
using namespace symmjkdigest;
using namespace sigshellpairinfor;
using namespace gints4deriv;

TBB_GInts4Deriv::TBB_GInts4Deriv(bool spStatus, 
		const UInt& braBatchIndex0, const UInt& ketBatchIndex0, 
		const Molecule& mol, const GIntJobInfor& ginfor0, 
		const SigMolShellPairInfor& bra, const SigMolShellPairInfor& ket,
		const MolShell& b1, const MolShell& b2, const MolShell& k1, const MolShell& k2, 
		const Mtrx& alpha, const Mtrx& beta,const DenMtrxInfor& denMtrxInfor0):sameShellPair(spStatus),
	nAtoms(mol.getNAtoms()),nSpin(ginfor0.getNSpin()),thresh(ginfor0.pickUpThresh()),
	spThresh(ginfor0.pickUpSPThresh()),braBatchIndex(braBatchIndex0),ketBatchIndex(ketBatchIndex0),
	bra1(b1),bra2(b2),ket1(k1),ket2(k2),braInfo(bra),ketInfo(ket),ginfor(ginfor0),
	alphaDen(alpha),betaDen(beta),denMtrxInfor(denMtrxInfor0)
{ 
	// initialize the result
	if (ginfor.getJobOrder() == 1) {
		derivData.init(nAtoms,3);
	}else if (ginfor.getJobOrder() == 2) {
		derivData.init(nAtoms*3,nAtoms*3);
	}
}

TBB_GInts4Deriv::TBB_GInts4Deriv(const TBB_GInts4Deriv& tbb_gints4deriv, 
		split):sameShellPair(tbb_gints4deriv.sameShellPair),nAtoms(tbb_gints4deriv.nAtoms),
	nSpin(tbb_gints4deriv.nSpin),thresh(tbb_gints4deriv.thresh),spThresh(tbb_gints4deriv.spThresh),
	braBatchIndex(tbb_gints4deriv.braBatchIndex),ketBatchIndex(tbb_gints4deriv.ketBatchIndex),
	bra1(tbb_gints4deriv.bra1),bra2(tbb_gints4deriv.bra2),ket1(tbb_gints4deriv.ket1),
	ket2(tbb_gints4deriv.ket2),braInfo(tbb_gints4deriv.braInfo),ketInfo(tbb_gints4deriv.ketInfo),
	ginfor(tbb_gints4deriv.ginfor),alphaDen(tbb_gints4deriv.alphaDen),betaDen(tbb_gints4deriv.betaDen), 
	denMtrxInfor(tbb_gints4deriv.denMtrxInfor)
{ 
	// initialize the result
	if (ginfor.getJobOrder() == 1) {
		derivData.init(nAtoms,3);
	}else if (ginfor.getJobOrder() == 2) {
		derivData.init(nAtoms*3,nAtoms*3);
	}
}

void TBB_GInts4Deriv::operator()(const blocked_range<UInt>& r) 
{
	//
	// setup scratch data for ERI calculation
	// 1  atom shell pair (bra and ket)
	// 2  heap memory management if necessary
	// 3  integral array, for elementary derivatives and 
	//    additional one for redundant derivatives
	// 4  jk digestion class
	//

	// shell pair data
	UInt maxNP2 = braInfo.getMaxNP2();
	UInt maxNL  = braInfo.getMaxNL();
	UInt maxNSP = braInfo.getMaxNSP();
	AtomShellPair bra(maxNSP,maxNP2,maxNL);
	maxNP2 = ketInfo.getMaxNP2();
	maxNL  = ketInfo.getMaxNL();
	maxNSP = ketInfo.getMaxNSP();
	AtomShellPair ket(maxNSP,maxNP2,maxNL);

	// heap memory tool
	UInt maxLBra = braInfo.getMaxL();
	UInt maxLKet = ketInfo.getMaxL();
	UInt maxL    = maxLBra>=maxLKet ? maxLBra : maxLKet;
	UInt heapMemLen = 0;
	if (ginfor.useHeapMem(maxL)) {
		heapMemLen = ginfor.setMemSize(maxL);
	}
	LocalMemScr scr(heapMemLen);

	// integral derivatives array for elementary/redundant derivatives
	UInt maxBraNCarBas = braInfo.getMaxNBasisForShell(); 
	UInt maxKetNCarBas = ketInfo.getMaxNBasisForShell(); 
	UInt maxIntLen     = maxBraNCarBas*maxBraNCarBas*maxKetNCarBas*maxKetNCarBas;
	UInt nElemDeriv    = ginfor.getNElemDeriv();
	DoubleVec intArray(maxIntLen*nElemDeriv);

	// symm jk digestion
	UInt maxNCarBas = maxBraNCarBas>maxKetNCarBas ? maxBraNCarBas : maxKetNCarBas;
	SymmJKDerivDigest jkDigest(ginfor,maxNCarBas);

	// now real work begins
	for( UInt n=r.begin(); n!=r.end(); ++n ) {

		//
		//get the atom shell pair data information
		//
		UInt rowAtomSPIndex = -1;
		UInt colAtomSPIndex = -1;
		UInt nBraSigAtomSP  = braInfo.getNSigAtomSPInBatch(braBatchIndex);
		UInt braBeginIndex  = braInfo.getStartingSigAtomSPIndex(braBatchIndex);
		UInt ketBeginIndex  = ketInfo.getStartingSigAtomSPIndex(ketBatchIndex);
		if (sameShellPair && braBatchIndex == ketBatchIndex) {
			UInt rowIndex = -1;
			UInt colIndex = -1;
			getRowColIndexFromPackIndex(nBraSigAtomSP,n,rowIndex,colIndex); 
			rowAtomSPIndex = braBeginIndex + rowIndex;
			colAtomSPIndex = ketBeginIndex + colIndex;
		}else{
			UInt colIndex  = n/nBraSigAtomSP;
			UInt rowIndex  = n-colIndex*nBraSigAtomSP;
			rowAtomSPIndex = braBeginIndex + rowIndex;
			colAtomSPIndex = ketBeginIndex + colIndex;
		}
		const SigAtomShellPairInfor& sigBraAtomSPInfor = braInfo.getSigAtomShellPairInfor(rowAtomSPIndex); 
		const SigAtomShellPairInfor& sigKetAtomSPInfor = ketInfo.getSigAtomShellPairInfor(colAtomSPIndex); 

		//
		// whether the atom shell quartet is neligible?
		// for only J job, things is a bit of different. Because 
		// the integral boundary, in this case will contain the Pmunu value,
		// the other cases that the integral boundary value is only the integral itself
		//
		bool isJImportant = false;
		bool isKImportant = false;
		UInt intJob = ginfor.getIntJob();
		if (doJ(intJob)) {
			if (onlyDoJ(intJob)) {
				if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,thresh)) isJImportant = true;
			}else{
				if(denMtrxInfor.isJDerivSig(thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isJImportant = true;
			}
		}
		if (doK(intJob)) {
			if (denMtrxInfor.isKDerivSig(0,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isKImportant = true;
			if (nSpin == 2) {
				if (! isKImportant && denMtrxInfor.isKDerivSig(1,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) {
					isKImportant = true;
				}
			}
		}
		if (! isJImportant && ! isKImportant) continue;

		//
		//form atom shell pair data
		//
		bra.init(bra1,bra2,sigBraAtomSPInfor,spThresh);
		ket.init(ket1,ket2,sigKetAtomSPInfor,spThresh);

		//
		//generate integral result
		//
		UInt iNSP     = bra.getNShellPairs(); 
		UInt jNSP     = ket.getNShellPairs(); 
		for(UInt jSP=0; jSP<jNSP; jSP++) {
			UInt iSPStart = 0;
			if (sameShellPair && bra == ket) iSPStart = jSP;
			for(UInt iSP=iSPStart; iSP<iNSP; iSP++) {

				// get shell pair data
				// and begin to prepare data
				const ShellPair& braSP = bra.getShellPair(iSP);
				const ShellPair& ketSP = ket.getShellPair(jSP);
				//braSP.print(3);
				//ketSP.print(3);

				// also perform screening test using cauchy-schwarz creteria
				// for the given shell quartet
				bool isJSig = false;
				bool isKSig = false;
				Double pMax = ZERO;
				if (doJ(intJob)) {
					Double jPPairMax = denMtrxInfor.getJPPairMax(braSP,ketSP);
					if (onlyDoJ(intJob)) {
						if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,thresh)) isJSig = true;
					}else{
						if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,jPPairMax,thresh)) isJSig = true;
					}
					if (jPPairMax>pMax) pMax = jPPairMax;
				}
				if (doK(intJob)) {
					Double kPPairMax = denMtrxInfor.getKPPairMax(braSP,ketSP);
					if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,kPPairMax,thresh)) isKSig = true;
					if (kPPairMax>pMax) pMax = kPPairMax;
				}

				// if both J and K insigfinicant, we just continue
				if (! isJSig && ! isKSig) continue;

				// angular momentum
				UInt iLmin = -1;
				UInt iLmax = -1;
				UInt jLmin = -1;
				UInt jLmax = -1;
				braSP.getL(iLmin,iLmax,jLmin,jLmax);
				UInt maxL1 = braSP.getLmax();
				UInt kLmin = -1;
				UInt kLmax = -1;
				UInt lLmin = -1;
				UInt lLmax = -1;
				ketSP.getL(kLmin,kLmax,lLmin,lLmax);
				UInt maxL2 = ketSP.getLmax();
				bool switchSP = switchShellPair(iLmin,iLmax,jLmin,jLmax,kLmin,kLmax,lLmin,lLmax);
				UInt maxL  = maxL1 >= maxL2 ? maxL1 : maxL2;

				// LCode
				LInt LCode = 0;
				UInt braLCode = braSP.getLCode();
				UInt ketLCode = ketSP.getLCode();
				if (switchSP) {
					LCode = codeSPCodes(ketLCode,braLCode);
				}else{
					LCode = codeSPCodes(braLCode,ketLCode);
				}

				// primitive pairs data
				UInt inp2          = braSP.getNP2();
				UInt jnp2          = ketSP.getNP2();
				const Double* ic2  = braSP.getC2();
				const Double* ie2  = braSP.getE2();
				const Double* A    = braSP.getA();
				const Double* B    = braSP.getB();
				const Double* P    = braSP.getP();
				const Double* fbra = braSP.getFac();
				const Double* jc2  = ketSP.getC2();
				const Double* je2  = ketSP.getE2();
				const Double* C    = ketSP.getA();
				const Double* D    = ketSP.getB();
				const Double* Q    = ketSP.getP();
				const Double* fket = ketSP.getFac();
				const Double* ie2diff= braSP.getExpDiff();
				const Double* je2diff= ketSP.getExpDiff();

				// initilize the integral
				UInt iNCarBas = 0; 
				UInt jNCarBas = 0;
				braSP.getNCarBas(iNCarBas,jNCarBas);
				UInt kNCarBas = 0; 
				UInt lNCarBas = 0;
				ketSP.getNCarBas(kNCarBas,lNCarBas);
				UInt intLen = iNCarBas*jNCarBas*kNCarBas*lNCarBas;
				vset(&intArray.front(),ZERO,intLen*nElemDeriv);

				// now ERI calculation
				// right now we just set for normal ERI, so omega = 0
				Double omega = ZERO;
				if (switchSP) {
					if (ginfor.getJobOrder() == 1) {
						hgp_os_eri_d1(LCode,maxL,jnp2,inp2,pMax,omega,jc2,je2,je2diff,fket,Q,C,D,
								ic2,ie2,ie2diff,fbra,P,A,B,&intArray.front(),scr);
					}else if (ginfor.getJobOrder() == 2) {
						hgp_os_eri_d2(LCode,maxL,jnp2,inp2,pMax,omega,jc2,je2,je2diff,fket,Q,C,D,
								ic2,ie2,ie2diff,fbra,P,A,B,&intArray.front(),scr);
					}
				}else{
					if (ginfor.getJobOrder() == 1) {
						hgp_os_eri_d1(LCode,maxL,inp2,jnp2,pMax,omega,ic2,ie2,ie2diff,fbra,P,A,B,
								jc2,je2,je2diff,fket,Q,C,D,&intArray.front(),scr);
					}else if (ginfor.getJobOrder() == 2) {
						hgp_os_eri_d2(LCode,maxL,inp2,jnp2,pMax,omega,ic2,ie2,ie2diff,fbra,P,A,B,
								jc2,je2,je2diff,fket,Q,C,D,&intArray.front(),scr);
					}
				}

				// now get the integral infor and let's prepare for digestion
				const SingleIntegralInfor& intInfor = ginfor.getIntegralInfor(LCode);
				if (intInfor.getMemInfor() > 0) {
					scr.reset();
				}

				// let's judge that whether we should do digestion
				bool doDigestion = false;
				for(UInt i=0; i<intLen*nElemDeriv; i++) {
					Double val = pMax*fabs(intArray[i]);
					if (val>thresh) {
						doDigestion = true;
						break;
					}
				}

				// now let's go to the job
				if (! doDigestion) continue;

				// now initialize the jk digestion information
				// and digest raw integrals into result
				jkDigest.initDerivInfor(switchSP,braSP,ketSP,bra,ket,alphaDen,betaDen);
				jkDigest.digestJKInts(intInfor,&intArray.front(),derivData);
			}
		}
	}
}

void GInts4Deriv::do4DAODeriv(const SigMolShellPairInfor& braInfor, const SigMolShellPairInfor& ketInfor,
	const MolShell& bra1, const MolShell& bra2, 
	const MolShell& ket1, const MolShell& ket2, 
	const Mtrx& alphaDen, const Mtrx& betaDen, const DenMtrxInfor& denMtrxInfor, Mtrx& result) const
{
	//
	// setup scratch data for ERI calculation
	// 1  atom shell pair (bra and ket)
	// 2  heap memory management if necessary
	// 3  integral array
	//

	// shell pair data
	UInt maxNP2 = braInfor.getMaxNP2();
	UInt maxNL  = braInfor.getMaxNL();
	UInt maxNSP = braInfor.getMaxNSP();
	AtomShellPair bra(maxNSP,maxNP2,maxNL);
	maxNP2 = ketInfor.getMaxNP2();
	maxNL  = ketInfor.getMaxNL();
	maxNSP = ketInfor.getMaxNSP();
	AtomShellPair ket(maxNSP,maxNP2,maxNL);

	// heap memory tool
	UInt maxLBra = braInfor.getMaxL();
	UInt maxLKet = ketInfor.getMaxL();
	UInt maxL    = maxLBra>=maxLKet ? maxLBra : maxLKet;
	UInt heapMemLen = 0;
	if (ginfor.useHeapMem(maxL)) {
		heapMemLen = ginfor.setMemSize(maxL);
	}
	LocalMemScr scr(heapMemLen);

	// integral derivatives array for elementary/redundant derivatives
	UInt maxBraNCarBas = braInfor.getMaxNBasisForShell(); 
	UInt maxKetNCarBas = ketInfor.getMaxNBasisForShell(); 
	UInt maxIntLen     = maxBraNCarBas*maxBraNCarBas*maxKetNCarBas*maxKetNCarBas;
	UInt nElemDeriv    = ginfor.getNElemDeriv();
	DoubleVec intArray(maxIntLen*nElemDeriv);

	// symm jk digestion
	UInt maxNCarBas = maxBraNCarBas>maxKetNCarBas ? maxBraNCarBas : maxKetNCarBas;
	SymmJKDerivDigest jkDigest(ginfor,maxNCarBas);

	// set up theshold value
	Double thresh   = ginfor.pickUpThresh();
	Double spThresh = ginfor.pickUpSPThresh();

	// now real work begins
	// loop over batches first
	for( UInt colBatchIndex=0; colBatchIndex<ketInfor.getNBatch(); colBatchIndex++) {
		UInt rowBeginBatchIndex = 0;
		if (braInfor == ketInfor) rowBeginBatchIndex = colBatchIndex;
		for( UInt rowBatchIndex=rowBeginBatchIndex; rowBatchIndex<braInfor.getNBatch(); rowBatchIndex++) {

			// now let's see whether we neglect this batch pair?
			Double braIntBoundary = braInfor.getBatchLowerLimitBoundary(rowBatchIndex);
			Double ketIntBoundary = ketInfor.getBatchLowerLimitBoundary(colBatchIndex);
			if (braIntBoundary*ketIntBoundary<thresh) continue;

			// initilize the index etc.
			UInt rowBeginIndex = braInfor.getStartingSigAtomSPIndex(rowBatchIndex);
			UInt colBeginIndex = ketInfor.getStartingSigAtomSPIndex(colBatchIndex);
			UInt nBraSigAtomSP = braInfor.getNSigAtomSPInBatch(rowBatchIndex);
			UInt nKetSigAtomSP = ketInfor.getNSigAtomSPInBatch(colBatchIndex);

			// let's go into the batch pair and loop over atom shell pairs
			for( UInt colAtomSPIndex=colBeginIndex; colAtomSPIndex<colBeginIndex+nKetSigAtomSP; colAtomSPIndex++) {

				//
				// for atom shell pair data for ket part
				//
				const SigAtomShellPairInfor& sigKetAtomSPInfor = ketInfor.getSigAtomShellPairInfor(colAtomSPIndex); 
				ket.init(ket1,ket2,sigKetAtomSPInfor,spThresh);

				// now loop over bra side
				UInt rowAtomSPIndex0 = rowBeginIndex;
				if (braInfor == ketInfor && rowBatchIndex == colBatchIndex) rowAtomSPIndex0 = colAtomSPIndex;
				for( UInt rowAtomSPIndex=rowAtomSPIndex0; rowAtomSPIndex<rowBeginIndex+nBraSigAtomSP; rowAtomSPIndex++) {

					//
					//get the atom shell pair data information
					//
					const SigAtomShellPairInfor& sigBraAtomSPInfor = braInfor.getSigAtomShellPairInfor(rowAtomSPIndex); 

					//
					// whether the atom shell quartet is neligible?
					//
					bool isJImportant = false;
					bool isKImportant = false;
					UInt intJob = ginfor.getIntJob();
					if (doJ(intJob)) {
						if (onlyDoJ(intJob)) {
							if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,thresh)) isJImportant = true;
						}else{
							if(denMtrxInfor.isJDerivSig(thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isJImportant = true;
						}
					}
					if (doK(intJob)) {
						if (denMtrxInfor.isKDerivSig(0,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isKImportant = true;
						if (ginfor.getNSpin() == 2) {
							if (! isKImportant && denMtrxInfor.isKDerivSig(1,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) {
								isKImportant = true;
							}
						}
					}
					if (! isJImportant && ! isKImportant) continue;

					//
					// initialize the bra shell pair
					//
					bra.init(bra1,bra2,sigBraAtomSPInfor,spThresh);

					//
					//generate integral result
					//
					UInt iNSP     = bra.getNShellPairs(); 
					UInt jNSP     = ket.getNShellPairs(); 
					for(UInt jSP=0; jSP<jNSP; jSP++) {
						UInt iSPStart = 0;
						if (braInfor == ketInfor && bra == ket) iSPStart = jSP;
						for(UInt iSP=iSPStart; iSP<iNSP; iSP++) {

							// get shell pair data
							// and begin to prepare data
							const ShellPair& braSP = bra.getShellPair(iSP);
							const ShellPair& ketSP = ket.getShellPair(jSP);
							//braSP.print(3);
							//ketSP.print(3);

							// also perform screening test using cauchy-schwarz creteria
							// for the given shell quartet
							bool isJSig = false;
							bool isKSig = false;
							Double pMax = ZERO;
							if (doJ(intJob)) {
								Double jPPairMax = denMtrxInfor.getJPPairMax(braSP,ketSP);
								if (onlyDoJ(intJob)) {
									if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,thresh)) isJSig = true;
								}else{
									if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,jPPairMax,thresh)) isJSig = true;
								}
								if (jPPairMax>pMax) pMax = jPPairMax;
							}
							if (doK(intJob)) {
								Double kPPairMax = denMtrxInfor.getKPPairMax(braSP,ketSP);
								if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,kPPairMax,thresh)) isKSig = true;
								if (kPPairMax>pMax) pMax = kPPairMax;
							}

							// if both J and K insigfinicant, we just continue
							if (! isJSig && ! isKSig) continue;

							// angular momentum
							UInt iLmin = -1;
							UInt iLmax = -1;
							UInt jLmin = -1;
							UInt jLmax = -1;
							braSP.getL(iLmin,iLmax,jLmin,jLmax);
							UInt maxL1 = braSP.getLmax();
							UInt kLmin = -1;
							UInt kLmax = -1;
							UInt lLmin = -1;
							UInt lLmax = -1;
							ketSP.getL(kLmin,kLmax,lLmin,lLmax);
							UInt maxL2 = ketSP.getLmax();
							bool switchSP = switchShellPair(iLmin,iLmax,jLmin,jLmax,kLmin,kLmax,lLmin,lLmax);
							UInt maxL  = maxL1 >= maxL2 ? maxL1 : maxL2;

							// LCode
							LInt LCode = 0;
							UInt braLCode = braSP.getLCode();
							UInt ketLCode = ketSP.getLCode();
							if (switchSP) {
								LCode = codeSPCodes(ketLCode,braLCode);
							}else{
								LCode = codeSPCodes(braLCode,ketLCode);
							}

							// primitive pairs data
							UInt inp2          = braSP.getNP2();
							UInt jnp2          = ketSP.getNP2();
							const Double* ic2  = braSP.getC2();
							const Double* ie2  = braSP.getE2();
							const Double* A    = braSP.getA();
							const Double* B    = braSP.getB();
							const Double* P    = braSP.getP();
							const Double* fbra = braSP.getFac();
							const Double* jc2  = ketSP.getC2();
							const Double* je2  = ketSP.getE2();
							const Double* C    = ketSP.getA();
							const Double* D    = ketSP.getB();
							const Double* Q    = ketSP.getP();
							const Double* fket = ketSP.getFac();
							const Double* ie2diff= braSP.getExpDiff();
							const Double* je2diff= ketSP.getExpDiff();

							// initilize the integral
							UInt iNCarBas = 0; 
							UInt jNCarBas = 0;
							braSP.getNCarBas(iNCarBas,jNCarBas);
							UInt kNCarBas = 0; 
							UInt lNCarBas = 0;
							ketSP.getNCarBas(kNCarBas,lNCarBas);
							UInt intLen = iNCarBas*jNCarBas*kNCarBas*lNCarBas;
							vset(&intArray.front(),ZERO,intLen*nElemDeriv);

							// now ERI calculation
							// right now we just set for normal ERI, so omega = 0
							Double omega = ZERO;
							if (switchSP) {
								if (ginfor.getJobOrder() == 1) {
									hgp_os_eri_d1(LCode,maxL,jnp2,inp2,pMax,omega,jc2,je2,je2diff,fket,Q,C,D,
											ic2,ie2,ie2diff,fbra,P,A,B,&intArray.front(),scr);
								}else if (ginfor.getJobOrder() == 2) {
									hgp_os_eri_d2(LCode,maxL,jnp2,inp2,pMax,omega,jc2,je2,je2diff,fket,Q,C,D,
											ic2,ie2,ie2diff,fbra,P,A,B,&intArray.front(),scr);
								}
							}else{
								if (ginfor.getJobOrder() == 1) {
									hgp_os_eri_d1(LCode,maxL,inp2,jnp2,pMax,omega,ic2,ie2,ie2diff,fbra,P,A,B,
											jc2,je2,je2diff,fket,Q,C,D,&intArray.front(),scr);
								}else if (ginfor.getJobOrder() == 2) {
									hgp_os_eri_d2(LCode,maxL,inp2,jnp2,pMax,omega,ic2,ie2,ie2diff,fbra,P,A,B,
											jc2,je2,je2diff,fket,Q,C,D,&intArray.front(),scr);
								}
							}

							// now get the integral infor and let's prepare for digestion
							const SingleIntegralInfor& intInfor = ginfor.getIntegralInfor(LCode);
							if (intInfor.getMemInfor() > 0) {
								scr.reset();
							}

							// let's judge that whether we should do digestion
							bool doDigestion = false;
							for(UInt i=0; i<intLen*nElemDeriv; i++) {
								Double val = pMax*fabs(intArray[i]);
								if (val>thresh) {
									doDigestion = true;
									break;
								}
							}

							// now let's go to the job
							if (! doDigestion) continue;

							// now initialize the jk digestion information
							// and digest raw integrals into result
							jkDigest.initDerivInfor(switchSP,braSP,ketSP,bra,ket,alphaDen,betaDen);
							jkDigest.digestJKInts(intInfor,&intArray.front(),result);
						}
					}
				}
			}
		}
	}
}

void GInts4Deriv::doJKDeriv(const MolShell& s, const Molecule& mol, DenMtrx& den, 
		Mtrx& result, bool printTiming, bool printResult) 
{
	// get the global information
	const GlobalInfor& globInfor = ginfor.getGlobalInfor();

	// firstly, we need to transform the density matrix into Cartesian order
	// this transformation is in fact a C2P transformation, but we just use
	// it's transpose form so that it looks to be P2C
	CPTransBasisNorm trans1(globInfor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW_COL);
	trans1.transform(s,s,den.getMtrx(0));
	if (ginfor.getNSpin() == 2) {
		trans1.transform(s,s,den.getMtrx(1));
	}

	// now let's form the density matrix infor
	DenMtrxInfor denMtrxInfor(globInfor,s,s,den);

	// form the significant shell pair information
	SigMolShellPairInfor braInfor(s,s,denMtrxInfor,ginfor);
	SigMolShellPairInfor ketInfor(braInfor);

	// possible timing code
	// and initilize the energy
	tick_count t0 = tick_count::now();

	if (globInfor.useMultiThreads()) {

		// set up the scheduler
		task_scheduler_init init(task_scheduler_init::deferred);
		init.initialize(globInfor.getNCPUThreads());

		// whether it's same shell pair
		bool sameShellPairs = false;
		if (braInfor == ketInfor) sameShellPairs = true;

		// now loop over batches
		for( UInt colBatchIndex=0; colBatchIndex<ketInfor.getNBatch(); colBatchIndex++) {
			UInt rowBeginBatchIndex = 0;
			if (braInfor == ketInfor) rowBeginBatchIndex = colBatchIndex;
			for( UInt rowBatchIndex=rowBeginBatchIndex; rowBatchIndex<braInfor.getNBatch(); rowBatchIndex++) {

				// now let's see whether we neglect this batch pair?
				Double braIntBoundary = braInfor.getBatchLowerLimitBoundary(rowBatchIndex);
				Double ketIntBoundary = ketInfor.getBatchLowerLimitBoundary(colBatchIndex);
				if (braIntBoundary*ketIntBoundary<ginfor.pickUpThresh()) continue;

				// initilize the length of work
				UInt nBraSigAtomSP = braInfor.getNSigAtomSPInBatch(rowBatchIndex);
				UInt nKetSigAtomSP = ketInfor.getNSigAtomSPInBatch(colBatchIndex);
				UInt len = nBraSigAtomSP*nKetSigAtomSP;
				if (sameShellPairs && rowBatchIndex == colBatchIndex) {
					len = (1+nBraSigAtomSP)*nBraSigAtomSP/2;
				}

				// do normal integral matrix calculation
				TBB_GInts4Deriv intWork(sameShellPairs,rowBatchIndex,colBatchIndex,mol,ginfor,braInfor,ketInfor,s,s,s,s,
						den.getMtrx(0),den.getMtrx(1),denMtrxInfor);
				parallel_reduce(blocked_range<UInt>(0,len),intWork);

				// let's update the result
				result.add(intWork.getResultMatrix());
			}
		}

		// possible timing code
		tick_count t1 = tick_count::now();
		Double t = (t1-t0).seconds();
		if (printTiming) {
			printf("%s  %-12.6f\n", "GInts4D derivatives work with parallel code, time in seconds  ", t);
		}

	}else{

		// do serial calculation here
		do4DAODeriv(braInfor,ketInfor,s,s,s,s,den.getMtrx(0),den.getMtrx(1),denMtrxInfor,result);

		// possible timing code
		tick_count t1 = tick_count::now();
		Double t = (t1-t0).seconds();
		if (printTiming) {
			printf("%s  %-12.6f\n", "GInts4D derivatives work with serial code, time in seconds  ", t);
		}
	}

	// now let's transorm the density matrix back, to undo the previous C2P,
	// we need the P2C data matrix
	CPTransBasisNorm trans2(globInfor,P2C_WITH_L00,UNDO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW_COL);
	trans2.transform(s,s,den.getMtrx(0));
	if (ginfor.getNSpin() == 2) {
		trans2.transform(s,s,den.getMtrx(1));
	}

	// finally let's print out the result matrix
	if (printResult) {
		result.print("result JK derivatives for gints4deriv");
	}
}

