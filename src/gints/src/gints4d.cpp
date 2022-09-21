/**
 * functions and classes associated with gints4d.h
 * working class to handle the 4D analytical integrals calculation
 * \author Fenglai Liu 
 */
#include<cmath>
#include<cstdio>
#include <boost/filesystem.hpp>
#include "tbb/tbb.h"
#include "blas.h"
#include "blas1.h"
#include "excep.h"
#include "shell.h"
#include "matrix.h"
#include "filerw.h"
#include "denmtrx.h"
#include "shellpair.h"
#include "shellprop.h"
#include "digestutil.h"
#include "spinmatrix.h"
#include "gintsinfor.h"
#include "blockmatrix.h"
#include "hgp_os_ints.h"
#include "localmemscr.h"
#include "integraljobs.h"
#include "denmtrxinfor.h"
#include "symmjkdigest.h"
#include "atomfockmtrx.h"
#include "blockmatrixlist.h"
#include "sigshellpairinfor.h"
#include "gints4d.h"
using namespace boost::filesystem;
using namespace blas;
using namespace excep;
using namespace shell;
using namespace matrix;
using namespace filerw;
using namespace denmtrx;
using namespace shellpair;
using namespace shellprop;
using namespace digestutil;
using namespace spinmatrix;
using namespace gintsinfor;
using namespace blockmatrix;
using namespace localmemscr;
using namespace integraljobs;
using namespace denmtrxinfor;
using namespace symmjkdigest;
using namespace atomfockmtrx;
using namespace blockmatrixlist;
using namespace sigshellpairinfor;
using namespace gints4d;

//
// define a global mutex just for gints4d part
//
tbb::spin_mutex GLOBAL_GINTS4D_MUTEX;

void TBB_GInts4D::operator()(const blocked_range<UInt>& r) const
{
	//
	// setup scratch data for ERI calculation
	// 1  atom shell pair (bra and ket)
	// 2  heap memory management if necessary
	// 3  integral array
	// 4  digestion class for J/K
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

	// integral array
	UInt maxBraNCarBas = braInfo.getMaxNBasisForShell(); 
	UInt maxKetNCarBas = ketInfo.getMaxNBasisForShell(); 
	UInt maxIntLen     = maxBraNCarBas*maxBraNCarBas*maxKetNCarBas*maxKetNCarBas;
	DoubleVec intArray(maxIntLen);

	// symm jk digestion
	UInt maxNCarBas = maxBraNCarBas>maxKetNCarBas ? maxBraNCarBas : maxKetNCarBas;
	SymmJKDigest jkDigest(ginfor,maxNCarBas);

	// atom fock matrix
	UInt maxBraNCarBasAtom = braInfo.getMaxNBasisForAtomShell();
	UInt maxKetNCarBasAtom = ketInfo.getMaxNBasisForAtomShell();
	UInt maxNCarBasAtom    = maxBraNCarBasAtom>=maxKetNCarBasAtom ? maxBraNCarBasAtom : maxKetNCarBasAtom;
	AtomFockMtrx atomFockMtrx(maxNCarBasAtom);

	// set up block matrix list corresponding alpha matrix and beta matrix
	BlockMtrxList alphaList(ginfor.getNAtomBlockResults(),ginfor.getNAdditionalAtomBlockResults(),
			maxNCarBasAtom,maxNCarBasAtom);
	BlockMtrxList betaList;
	if (ginfor.getNSpin() == 2) {
		betaList.initMem(ginfor.getNAtomBlockResults(),ginfor.getNAdditionalAtomBlockResults(),
				maxNCarBasAtom,maxNCarBasAtom);
	}

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
		// according to cauchy-schwarz creteria
		//
		UInt intJob = ginfor.getIntJob();
		bool isJSig = false;
		bool isAlphaKSig = false;
		bool isBetaKSig  = false;
		if (doJ(intJob)) {
			if (denMtrxInfor.isJSig(thresh,sigBraAtomSPInfor,sigKetAtomSPInfor))  isJSig = true;
		}
		if (doK(intJob)) {
			if (denMtrxInfor.isKSig(0,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isAlphaKSig = true;
			if (nSpin == 2) {
				if (denMtrxInfor.isKSig(1,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isBetaKSig = true;
			}
		}

		// now let's see whether we omit this atom shell quartet
		if (!isJSig && !isAlphaKSig && !isBetaKSig) continue;

		//
		//form atom shell pair data
		//
		bra.init(bra1,bra2,sigBraAtomSPInfor,spThresh);
		ket.init(ket1,ket2,sigKetAtomSPInfor,spThresh);

		//
		//initialize the atom fock matrix in terms of atom shell information
		//
		atomFockMtrx.initInfor(ginfor,bra,ket);

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

				// get density matrix maximum value for J/K digestion
				// also perform screening test using cauchy-schwarz creteria
				// for the given shell quartet
				Double pMax = ZERO;
				bool isJSig = false;
				bool isKSig = false;
				if (doJ(intJob)) {
					Double jPMax = denMtrxInfor.getJPMax(braSP,ketSP);
					if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,jPMax,thresh)) isJSig = true;
					if (jPMax>pMax) pMax = jPMax;
				}
				if (doK(intJob)) {
					Double kPMax = denMtrxInfor.getKPMax(braSP,ketSP);
					if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,kPMax,thresh)) isKSig = true;
					if (kPMax>pMax) pMax = kPMax;
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

				// initilize the integral
				UInt iNCarBas = 0; 
				UInt jNCarBas = 0;
				braSP.getNCarBas(iNCarBas,jNCarBas);
				UInt kNCarBas = 0; 
				UInt lNCarBas = 0;
				ketSP.getNCarBas(kNCarBas,lNCarBas);
				UInt intLen = iNCarBas*jNCarBas*kNCarBas*lNCarBas;
				vset(&intArray.front(),ZERO,intLen);

				// now ERI calculation
				Double omega = ZERO;
				if (switchSP) {
					hgp_os_eri(LCode,maxL,jnp2,inp2,pMax,omega,jc2,je2,fket,Q,C,D,ic2,ie2,fbra,
							P,A,B,&intArray.front(),scr);
				}else{
					hgp_os_eri(LCode,maxL,inp2,jnp2,pMax,omega,ic2,ie2,fbra,P,A,B,jc2,je2,fket,
							Q,C,D,&intArray.front(),scr);
				}

				// now get the integral infor and let's prepare for digestion
				const SingleIntegralInfor& intInfor = ginfor.getIntegralInfor(LCode);
				if (intInfor.getMemInfor() > 0) {
					scr.reset();
				}

				// some times we may have the result integral array
				// as insignificant, so we need to check it here
				Double maxIntVal = fabs(intArray[0]);
				for(UInt i=1;i<intLen;i++) {
					Double val = fabs(intArray[i]);
					if (val>maxIntVal) maxIntVal = val;
				}
				if (maxIntVal*pMax<thresh) continue;

				// now initialize the jk digestion information
				// and digest raw integrals into result
				jkDigest.initShellInfor(switchSP,braSP,ketSP,alphaDen,betaDen);
				if (doJKTogether(intJob)) {
					jkDigest.symmJKIntsDigest(&intArray.front());
				}else{
					if (doJ(intJob)) {
						jkDigest.symmJIntsDigest(&intArray.front());
					}
					if (doK(intJob)) {
						jkDigest.symmKIntsDigest(&intArray.front());
					}
				}

				// we need to merge the local result into the result in atom shell dimension
				atomFockMtrx.resetLinkInfor(switchSP,braSP,ketSP);
				atomFockMtrx.updateResults(ginfor,jkDigest);
			}
		}

		// now let's merge the atom dimension result into result list
		atomFockMtrx.merge(ginfor,alphaList,0);
		if(nSpin == 2) {
			atomFockMtrx.merge(ginfor,betaList,1);
		}

		// let's see whether the list reaches the length limit
		// or we just reach the end of the program
		UInt lastOne = r.end()-1;
		if (alphaList.reachLenLimit() || n == lastOne) {

			// now it's the critical section
			// to update the global results
			{
				tbb::spin_mutex::scoped_lock lock(GLOBAL_GINTS4D_MUTEX);
				alphaList.updateMatrix(alphaMatrix);
				if(nSpin == 2) {
					betaList.updateMatrix(betaMatrix);
				}
			}

			// now let's clear the data for the list
			alphaList.clear();
			if(nSpin == 2) {
				betaList.clear();
			}
		}
	}
}

void GInts4D::do4DAOMtrx(const SigMolShellPairInfor& braInfor, const SigMolShellPairInfor& ketInfor,
	const MolShell& bra1, const MolShell& bra2, 
	const MolShell& ket1, const MolShell& ket2, 
	const Mtrx& alphaDen, const Mtrx& betaDen, const DenMtrxInfor& denMtrxInfor,
	Mtrx& alphaMatrix, Mtrx& betaMatrix) const
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

	// integral array
	UInt maxBraNCarBas = braInfor.getMaxNBasisForShell(); 
	UInt maxKetNCarBas = ketInfor.getMaxNBasisForShell(); 
	UInt maxIntLen     = maxBraNCarBas*maxBraNCarBas*maxKetNCarBas*maxKetNCarBas;
	DoubleVec intArray(maxIntLen);

	// symm jk digestion
	UInt maxNCarBas = maxBraNCarBas>maxKetNCarBas ? maxBraNCarBas : maxKetNCarBas;
	SymmJKDigest jkDigest(ginfor,maxNCarBas);

	// thresh value
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
					// now test the atom shell pair significance
					//
					UInt intJob = ginfor.getIntJob();
					bool isJSig = false;
					bool isAlphaKSig = false;
					bool isBetaKSig  = false;
					if (doJ(intJob)) {
						if (denMtrxInfor.isJSig(thresh,sigBraAtomSPInfor,sigKetAtomSPInfor))  isJSig = true;
					}
					if (doK(intJob)) {
						if (denMtrxInfor.isKSig(0,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isAlphaKSig = true;
						if (ginfor.getNSpin() == 2) {
							if (denMtrxInfor.isKSig(1,thresh,sigBraAtomSPInfor,sigKetAtomSPInfor)) isBetaKSig = true;
						}
					}

					// now let's see whether we omit this atom shell quartet
					if (!isJSig && !isAlphaKSig && !isBetaKSig) continue;

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

							// get density matrix maximum value for J/K digestion
							// also perform screening test using cauchy-schwarz creteria
							// for the given shell quartet
							Double pMax = ZERO;
							bool isJSig = false;
							bool isKSig = false;
							if (doJ(intJob)) {
								Double jPMax = denMtrxInfor.getJPMax(braSP,ketSP);
								if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,jPMax,thresh)) isJSig = true;
								if (jPMax>pMax) pMax = jPMax;
							}
							if (doK(intJob)) {
								Double kPMax = denMtrxInfor.getKPMax(braSP,ketSP);
								if (! sigBraAtomSPInfor.isNegligible(iSP,jSP,sigKetAtomSPInfor,kPMax,thresh)) isKSig = true;
								if (kPMax>pMax) pMax = kPMax;
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

							// initilize the integral
							UInt iNCarBas = 0; 
							UInt jNCarBas = 0;
							braSP.getNCarBas(iNCarBas,jNCarBas);
							UInt kNCarBas = 0; 
							UInt lNCarBas = 0;
							ketSP.getNCarBas(kNCarBas,lNCarBas);
							UInt intLen = iNCarBas*jNCarBas*kNCarBas*lNCarBas;
							vset(&intArray.front(),ZERO,intLen);

							// now ERI calculation
							Double omega = ZERO;
							if (switchSP) {
								hgp_os_eri(LCode,maxL,jnp2,inp2,pMax,omega,jc2,je2,fket,Q,C,D,ic2,ie2,fbra,
										P,A,B,&intArray.front(),scr);
							}else{
								hgp_os_eri(LCode,maxL,inp2,jnp2,pMax,omega,ic2,ie2,fbra,P,A,B,jc2,je2,fket,
										Q,C,D,&intArray.front(),scr);
							}

							// now get the integral infor and let's prepare for digestion
							const SingleIntegralInfor& intInfor = ginfor.getIntegralInfor(LCode);
							if (intInfor.getMemInfor() > 0) {
								scr.reset();
							}

							// some times we may have the result integral array
							// as insignificant, so we need to check it here
							Double maxIntVal = fabs(intArray[0]);
							for(UInt i=1;i<intLen;i++) {
								Double val = fabs(intArray[i]);
								if (val>maxIntVal) maxIntVal = val;
							}
							if (maxIntVal*pMax<thresh) continue;

							// now initialize the jk digestion information
							// and digest raw integrals into result
							jkDigest.initShellInfor(switchSP,braSP,ketSP,alphaDen,betaDen);
							if (doJKTogether(intJob)) {
								jkDigest.symmJKIntsDigest(&intArray.front());
							}else{
								if (doJ(intJob)) {
									jkDigest.symmJIntsDigest(&intArray.front());
								}
								if (doK(intJob)) {
									jkDigest.symmKIntsDigest(&intArray.front());
								}
							}

							// finally update the matrix
							jkDigest.updateFockMatrix(alphaMatrix,betaMatrix);
						}
					}
				}
			}
		}
	}
}

void GInts4D::doJKMtrx(const MolShell& s, DenMtrx& den, SpinMatrix& intMtrx, 
		Double& eJK, bool printTiming, bool printMatrix, bool saveMatrix) 
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

	// for JK work, we need a blank result matrix
	// they will be in Cartesian dimension
	UInt nCarBas = s.getNCarBas();
	Mtrx alphaResult(nCarBas,nCarBas);
	Mtrx betaResult;
	if (intMtrx.getNSpin()==2) {
		betaResult.init(nCarBas,nCarBas);
	}

	// possible timing code
	// and initilize the energy
	tick_count t0 = tick_count::now();
	eJK = ZERO;

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
				TBB_GInts4D intWork(sameShellPairs,rowBatchIndex,colBatchIndex,ginfor,
						braInfor,ketInfor,s,s,s,s,den.getMtrx(0),den.getMtrx(1),denMtrxInfor,alphaResult,betaResult);
				parallel_for(blocked_range<UInt>(0,len),intWork);
			}
		}

		// possible timing code
		tick_count t1 = tick_count::now();
		Double t = (t1-t0).seconds();
		if (printTiming) {
			printf("%s  %-12.6f\n", "GInts4D work with parallel code, time in seconds  ", t);
		}

	}else{

		// do serial calculation here
		do4DAOMtrx(braInfor,ketInfor,s,s,s,s,den.getMtrx(0),den.getMtrx(1),denMtrxInfor,alphaResult,betaResult);

		// possible timing code
		tick_count t1 = tick_count::now();
		Double t = (t1-t0).seconds();
		if (printTiming) {
			printf("%s  %-12.6f\n", "GInts4D work with serial code, time in seconds  ", t);
		}
	}

	// this JK is not finished yet, we still need to do post-process
	postProcessJKMatrix(s,alphaResult);
	if(intMtrx.getNSpin() == 2) {
		postProcessJKMatrix(s,betaResult);
	}

	// we need to transform the integral matrix
	CPTransBasisNorm trans(globInfor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW_COL);
	trans.transform(s,s,alphaResult);
	if (ginfor.getNSpin() == 2) {
		trans.transform(s,s,betaResult);
	}

	// now let's transorm the density matrix back, to undo the previous C2P,
	// we need the P2C data matrix
	CPTransBasisNorm trans2(globInfor,P2C_WITH_L00,UNDO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW_COL);
	trans2.transform(s,s,den.getMtrx(0));
	if (ginfor.getNSpin() == 2) {
		trans2.transform(s,s,den.getMtrx(1));
	}

	// update result
	intMtrx.getMtrx(0).add(alphaResult);
	if (intMtrx.getNSpin() == 2) {
		intMtrx.getMtrx(1).add(betaResult);
	}

	// let's see whether we print out the result matrix?
	if (printMatrix) {
		string t0 = getJobName(ginfor.getIntJob());
		string alphaTitle = t0 + ": alpha matrix";
		string betaTitle  = t0 + ": beta  matrix";
		alphaResult.print(alphaTitle);
		if (intMtrx.getNSpin() == 2) {
			betaResult.print(betaTitle);
		}
	}

	// do you also want to save result matrix on disk?
	if (saveMatrix) {
		saveJKMtrx(den.getSec(),intMtrx.getNSpin(),globInfor,alphaResult,betaResult);
	}

	// update energy
	bool symmMatrix = true;
	const Mtrx& alphaDenMtrx = den.getMtrx(0);
	eJK += HALF*alphaResult.dotProduct(alphaDenMtrx,symmMatrix);
	if (ginfor.isCloseShell()) {
		eJK *= TWO;
	}else if (ginfor.getNSpin() == 2) {
		const Mtrx& betaDenMtrx = den.getMtrx(1);
		eJK += HALF*betaResult.dotProduct(betaDenMtrx,symmMatrix);
	}
}

void GInts4D::saveJKMtrx(UInt section, UInt nSpin, const GlobalInfor& globInfor, 
		const Mtrx& A, const Mtrx& B) const
{
	// get the data path 
	const string& scratchDir = globInfor.getScratchDir();
	path p(scratchDir.c_str());
	string sec = boost::lexical_cast<string>(section);
	path sect(sec.c_str());
	string gp = "gints";
	path gpath(gp.c_str());
	p /= sect;
	p /= gpath;

	// remove old data if there's same name folder exists
	// else create new one
	if (exists(p)) {
		remove_all(p);
	}
	create_directories(p);

	// now let's set up the name
	string t0 = getJobName(ginfor.getIntJob());

	// now let's begin to write data
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// let's get the name
		string name = t0 + "_a.bin";
		if (iSpin == 1) {
			name = t0 + "_b.bin";
		}

		// now do file writing
		FileReadWrite fileRW(p.string(),name);
		if (iSpin == 0) {
			const DoubleVec& va = A.getVec();
			fileRW.write(&va.front(),va.size());
		}else{
			const DoubleVec& vb = B.getVec();
			fileRW.write(&vb.front(),vb.size());
		}
	}
}

void GInts4D::recover(const string& dataPath, SpinMatrix& M) const
{
	// data dir location
	path p(dataPath.c_str());

	// test that whether we have the mo data for the path
	if (! exists(p)) {
		string info = "the corresponding data file does not exist: " + p.string();
		Excep excep("GInts4d","recover",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}

	// now begin to read in
	string t0 = getJobName(ginfor.getIntJob());
	for(UInt iSpin=0; iSpin<M.getNSpin(); iSpin++) {

		// let's get the name
		// this should be exactly same with the write function
		string name = t0 + "_a.bin";
		if (iSpin == 1) {
			name = t0 + "_b.bin";
		}

		// check whether the give file is over there?
		path p0(p);
		path p1(name.c_str());
		p0 /= p1;
		if (! exists(p0)) {
			string info = "the corresponding data file does not exist: " + p0.string();
			Excep excep("GInts4d","recover",EXCEPTION_FILE_MISSING,info);
			handleExcep(excep);
		}

		// now read in data
		Mtrx& Mi = M.getMtrx(iSpin);
		DoubleVec& vec = Mi.getVec();
		FileReadWrite fileRW(p.string(),name);
		fileRW.read(&vec.front(),vec.size());
	}
}

