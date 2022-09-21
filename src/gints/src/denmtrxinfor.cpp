/**
 * \file    denmtrxinfor.cpp
 * \author  Fenglai Liu 
 */
#include<cmath>
#include "blas.h"
#include "blas1.h"
#include "shell.h"
#include "excep.h"
#include "textread.h"
#include "parameterparsing.h"
#include "globalinfor.h"
#include "shellpair.h"
#include "gintsinfor.h"
#include "sigshellpairinfor.h"
#include "denmtrx.h"
#include "denmtrxinfor.h"
using namespace blas;
using namespace shell;
using namespace excep;
using namespace textread;
using namespace parameterparsing;
using namespace globalinfor;
using namespace shellpair;
using namespace gintsinfor;
using namespace sigshellpairinfor;
using namespace denmtrx;
using namespace denmtrxinfor;

void TBB_DenMtrxInfor::operator()(const blocked_range<UInt>& r) const 
{
	// the top loop is over the number of spin
	UInt nSpin = denMtrx.getNSpin();
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// now let's set the spin matrix
		const Mtrx& P = denMtrx.getMtrx(iSpin);
		Mtrx& M1      = denMtrxShellInfor.getMtrx(iSpin);
		Mtrx& M2      = denMtrxAtomInfor.getMtrx(iSpin);

		// whether the two shells are identical?
		bool sameShell = false;
		if (rs == cs) sameShell = true;

		// this is loop over the atom shell block
		for( UInt n=r.begin(); n!=r.end(); ++n ) {

			//
			//get the atom shell index
			//
			UInt rowAtomShellIndex = -1;
			UInt colAtomShellIndex = -1;
			UInt nRowAtomShells    = rs.getNAtomShells();
			if (sameShell) {
				getRowColIndexFromPackIndex(nRowAtomShells,n,rowAtomShellIndex,colAtomShellIndex); 
			}else{
				colAtomShellIndex   = n/nRowAtomShells;
				rowAtomShellIndex   = n-colAtomShellIndex*nRowAtomShells;
			}

			// get the atom shell reference
			const AtomShell& iAS = rs.getAtomShell(rowAtomShellIndex);
			const AtomShell& jAS = cs.getAtomShell(colAtomShellIndex);

			// set the value for atom-atom pair
			Double maxVal = ZERO;

			// now loop over shell-shell block
			// for same atom shell, we only goes the lower 
			// tri-angular part
			for(UInt jShell=0; jShell<jAS.getNShell(); jShell++) {
				UInt start = 0;
				if (sameShell && rowAtomShellIndex == colAtomShellIndex) start = jShell;
				for(UInt iShell=start; iShell<iAS.getNShell(); iShell++) {

					// get the shell-shell block
					const Shell& is = iAS.getShell(iShell);
					const Shell& js = jAS.getShell(jShell);

					// get the offset and dimension for the given block
					UInt iNBas = is.getNBas(type);
					UInt iOffset = is.getBasisIndex(0,type);
					UInt jNBas = js.getNBas(type);
					UInt jOffset = js.getBasisIndex(0,type);

					// now get the pMax for the given block
					Double max = ZERO;
					for(UInt j=0; j<jNBas; j++) {
						for(UInt i=0; i<iNBas; i++) {
							Double val = fabs(P.val(iOffset+i,jOffset+j));
							if (val>max) max = val;
						}
					}

					// now write it back
					UInt iShellIndex = is.getGlobalShellIndex();
					UInt jShellIndex = js.getGlobalShellIndex();
					M1(iShellIndex,jShellIndex) = max;
					if (sameShell) {
						M1(jShellIndex,iShellIndex) = max;
					}

					// let's see whether this is the largest for
					// atom-atom pair?
					if (max>maxVal) {
						maxVal = max;
					}
				}   
			}

			// now write it back for atom-atom pair
			UInt iAtom = iAS.getAtomShellIndex();
			UInt jAtom = jAS.getAtomShellIndex();
			M2(iAtom,jAtom) = maxVal;
			if (sameShell) {
				M2(jAtom,iAtom) = maxVal;
			}
		}
	}
}

DenMtrxInfor::DenMtrxInfor(const GlobalInfor& infor, const MolShell& rs, 
		const MolShell& cs, const DenMtrx& denMtrx):upperDimEnableThreads(100),type(TYPE_CART),
	maxDenMtrxPerRow(denMtrx.getNSpin()*rs.getNShell()),
	maxDenMtrxPerCol(denMtrx.getNSpin()*cs.getNShell()),
	denMtrxAtomInfor(denMtrx.getNSpin(),rs.getNAtomShells(),cs.getNAtomShells()), 
	denMtrxShellInfor(denMtrx.getNSpin(),rs.getNShell(),cs.getNShell())
{
	init(infor,rs,cs,denMtrx);
}

void DenMtrxInfor::init(const GlobalInfor& infor, const MolShell& rs, 
		const MolShell& cs, const DenMtrx& denMtrx)
{
	// do we have empty matrix data?
	bool needInit = false;
	if (maxDenMtrxPerRow.size() == 0) {
		needInit = true;
	}
	
	// if the data is not initialized, we do it here
	if (needInit) {
		denMtrxAtomInfor.init(rs.getNAtomShells(),cs.getNAtomShells()); 
		denMtrxShellInfor.init(rs.getNShell(),cs.getNShell());
		maxDenMtrxPerRow.assign(denMtrx.getNSpin()*rs.getNShell(),ZERO);
		maxDenMtrxPerCol.assign(denMtrx.getNSpin()*cs.getNShell(),ZERO);
		type = TYPE_CART;
		upperDimEnableThreads = 100;
	}

	// let's check the input data
	// the dimension of density matrix should be in accordance with
	// the input shell data
	// also set the type information
	for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {
		const Mtrx& P  = denMtrx.getMtrx(iSpin);

		// check density matrix diemension
		// they must share the same type
		UInt nRow = P.getRow();
		UInt t1 = 0;
		if(nRow == rs.getNCarBas()) t1=1;
		if(nRow == rs.getNBas())    t1=2;
		if(t1 == 0) {
			string info = "the input density matrix row dimension conflict with row shell data";
			Excep excep("DenMtrxInfor","constructor",EXCEPTION_PMAXINFOR_ERROR,info);
			handleExcep(excep);
		}	

		// col
		UInt nCol = P.getCol();
		UInt t2 = 0;
		if(nCol == cs.getNCarBas()) t2=1;
		if(nCol == cs.getNBas())    t2=2;
		if(t2 == 0) {
			string info = "the input density matrix col dimension conflict with col shell data";
			Excep excep("DenMtrxInfor","constructor",EXCEPTION_PMAXINFOR_ERROR,info);
			handleExcep(excep);
		}	

		// row shell type is not same with col shell type
		if(t1 != t2) {
			string info = "the input density matrix row and col dimension have different shell type";
			Excep excep("DenMtrxInfor","constructor",EXCEPTION_PMAXINFOR_ERROR,info);
			handleExcep(excep);
		}	

		// now set the type information
		UInt type0 = TYPE_NORM;
		if(t2 == 1) {
			type0 = TYPE_CART;
		}

		// additional check for different spin state
		if (iSpin == 1) {
			if(type0 != type) {
				string info = "the spin matrices have different shell type in input density matrix ";
				Excep excep("DenMtrxInfor","constructor",EXCEPTION_PMAXINFOR_ERROR,info);
				handleExcep(excep);
			}	
		}else{
			type = type0;
		}
	}

	// now everything good, let's process information if we define it from input
	string input = infor.getInputFile();
	UInt section = denMtrx.getSec();
	ParameterParsing pp(input,"denmtrxinfor",section);
	if (pp.hasAnyParameters()) {
		string key = "atom_number_limit_use_threads";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "atom_number_limit_use_threads can not be processed, not an integer";
				Excep excep("DenMtrxInfor","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			upperDimEnableThreads = tmp;
		}
	}

	// do we use threads in forming the data?
	bool useThreads = false;
	UInt minAtomShellNum = rs.getNAtomShells();
	if (minAtomShellNum<cs.getNAtomShells()) {
		minAtomShellNum = cs.getNAtomShells();
	}
	if (minAtomShellNum>=upperDimEnableThreads) {
		useThreads = true;
	}

	// form the data
	if (useThreads) {
		formPMaxInforInThreads(infor,rs,cs,denMtrx); 
	}else{
		directFormPMaxInfor(rs,cs,denMtrx); 
	}

	// finally form the maxDenMtrxPerRow/Col
	// the row/col dimension should be with the alpha state
	const Mtrx& M0 = denMtrxShellInfor.getMtrx(0);
	UInt nRow = M0.getRow();
	UInt nCol = M0.getCol();
	for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {
		const Mtrx& M  = denMtrxShellInfor.getMtrx(iSpin);

		// firstly let's search the column
		for(UInt iCol=0; iCol<M.getCol(); iCol++) {
			Double maxVal = maxSearch(M.getPtr(0,iCol), M.getRow());
			maxDenMtrxPerCol[iCol+nCol*iSpin] = maxVal;
		}

		// if the row shell same with col shell, thing is simper
		// on the other hand, we need to do another round of search
		bool sameShell = false;
		if (rs == cs) sameShell = true;
		if (sameShell) {
			maxDenMtrxPerRow = maxDenMtrxPerCol;
		}else{
			for(UInt iRow=0; iRow<M.getRow(); iRow++) {
				Double maxVal = ZERO;
				for(UInt iCol=0; iCol<M.getCol(); iCol++) {
					if (M.val(iRow,iCol)>maxVal) maxVal = M.val(iRow,iCol);
				}
				maxDenMtrxPerRow[iRow+nRow*iSpin] = maxVal;
			}
		}
	}
}

void DenMtrxInfor::directFormPMaxInfor(const MolShell& rs, const MolShell& cs, const DenMtrx& denMtrx) 
{
	// do we only do it for lower tri-angular part?
	// this is true for same shell case
	bool sameShell = false;
	if (rs == cs) sameShell = true;

	// loop over the spin of pMax
	for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {

		// get the input density matrix
		// and result matrix
		const Mtrx& P = denMtrx.getMtrx(iSpin);
		Mtrx& M1 = denMtrxShellInfor.getMtrx(iSpin);
		Mtrx& M2 = denMtrxAtomInfor.getMtrx(iSpin);

		// loop over the column shell
		for(UInt jAtom=0; jAtom<cs.getNAtomShells(); jAtom++) {
			const AtomShell& jAS = cs.getAtomShell(jAtom);

			// loop over row shell
			UInt begin = 0;
			if (sameShell) begin = jAtom;
			for(UInt iAtom=begin; iAtom<rs.getNAtomShells(); iAtom++) {
				const AtomShell& iAS = rs.getAtomShell(iAtom);

				// set the value for atom-atom pair
				Double maxVal = ZERO;

				// now loop over shell-shell block
				// for same atom shell, we only goes the lower 
				// tri-angular part
				for(UInt jShell=0; jShell<jAS.getNShell(); jShell++) {
					UInt start = 0;
					if (sameShell && iAtom == jAtom) start = jShell;
					for(UInt iShell=start; iShell<iAS.getNShell(); iShell++) {

						// get the shell-shell block
						const Shell& js = jAS.getShell(jShell);
						const Shell& is = iAS.getShell(iShell);

						// get the offset and dimension for the given block
						UInt iNBas = is.getNBas(type);
						UInt iOffset = is.getBasisIndex(0,type);
						UInt jNBas = js.getNBas(type);
						UInt jOffset = js.getBasisIndex(0,type);

						// now get the pMax for the given block
						Double max = ZERO;
						for(UInt j=0; j<jNBas; j++) {
							for(UInt i=0; i<iNBas; i++) {
								Double val = fabs(P.val(iOffset+i,jOffset+j));
								if (val>max) max = val;
							}
						}

						// now write it back
						UInt iShellIndex = is.getGlobalShellIndex();
						UInt jShellIndex = js.getGlobalShellIndex();
						M1(iShellIndex,jShellIndex) = max;
						if (sameShell) {
							M1(jShellIndex,iShellIndex) = max;
						}

						// let's see whether this is the largest for
						// atom-atom pair?
						if (max>maxVal) {
							maxVal = max;
						}
					}
				}   

				// now write it back for atom-atom pair
				M2(iAtom,jAtom) = maxVal;
				if (sameShell) {
					M2(jAtom,iAtom) = maxVal;
				}
			}
		}
	}
}

void DenMtrxInfor::formPMaxInforInThreads(const GlobalInfor& infor, 
		const MolShell& rs, const MolShell& cs, const DenMtrx& denMtrx)
{
	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	init.initialize(infor.getNCPUThreads());

	// possible timing code
	//tick_count t0 = tick_count::now();
	TBB_DenMtrxInfor  inforWork(type,rs,cs,denMtrx,denMtrxAtomInfor,denMtrxShellInfor); 
	UInt len = rs.getNAtomShells()*cs.getNAtomShells();
	if (rs == cs) {
		UInt num = rs.getNAtomShells(); 
		len = num*(num+1)/2;
	}
	parallel_for(blocked_range<UInt>(0,len), inforWork);

	// possible timing code
	//tick_count t1 = tick_count::now();
	//Double t = (t1-t0).seconds();
	//printf("%s  %-12.6f\n", "form density matrix infor work with TBB threads, time in seconds  ", t);
}

Double DenMtrxInfor::getJPMax(const ShellPair& braSP, const ShellPair& ketSP) const
{
	// get the shell index
	UInt bra1Index = -1;
	UInt bra2Index = -1;
	UInt ket1Index = -1;
	UInt ket2Index = -1;
	braSP.getGlobalShellIndex(bra1Index,bra2Index);
	ketSP.getGlobalShellIndex(ket1Index,ket2Index);

	// matrix data
	const Mtrx& alphaDenInfor = denMtrxShellInfor.getMtrx(0);
	const Mtrx& betaDenInfor  = denMtrxShellInfor.getMtrx(1);
	Double braPMax = alphaDenInfor.val(bra1Index,bra2Index)+betaDenInfor.val(bra1Index,bra2Index);
	Double ketPMax = alphaDenInfor.val(ket1Index,ket2Index)+betaDenInfor.val(ket1Index,ket2Index);
	Double pMax  = braPMax>ketPMax ? braPMax : ketPMax;
	return pMax;
}

Double DenMtrxInfor::getKPMax(const ShellPair& braSP, const ShellPair& ketSP) const 
{
	// get the shell index
	UInt bra1Index = -1;
	UInt bra2Index = -1;
	UInt ket1Index = -1;
	UInt ket2Index = -1;
	braSP.getGlobalShellIndex(bra1Index,bra2Index);
	ketSP.getGlobalShellIndex(ket1Index,ket2Index);

	// now get the pMax value
	Double maxDenVal = ZERO;
	for(UInt iSpin=0; iSpin<denMtrxShellInfor.getNSpin(); iSpin++) {

		// get the data
		const Mtrx& denInfor = denMtrxShellInfor.getMtrx(iSpin);
		Double pMax13 = denInfor.val(bra1Index,ket1Index);
		Double pMax14 = denInfor.val(bra1Index,ket2Index);
		Double pMax23 = denInfor.val(bra2Index,ket1Index);
		Double pMax24 = denInfor.val(bra2Index,ket2Index);

		// compare
		Double pMax1 = pMax13>pMax14 ? pMax13 : pMax14;
		Double pMax2 = pMax23>pMax24 ? pMax23 : pMax24;
		Double pMax  = pMax1>pMax2 ? pMax1 : pMax2;
		if (pMax>maxDenVal) maxDenVal = pMax;
	}

	// return
	return maxDenVal;
}

Double DenMtrxInfor::getJPPairMax(const ShellPair& braSP, const ShellPair& ketSP) const
{
	// get the shell index
	UInt bra1Index = -1;
	UInt bra2Index = -1;
	UInt ket1Index = -1;
	UInt ket2Index = -1;
	braSP.getGlobalShellIndex(bra1Index,bra2Index);
	ketSP.getGlobalShellIndex(ket1Index,ket2Index);

	// matrix data
	const Mtrx& alphaDenInfor = denMtrxShellInfor.getMtrx(0);
	const Mtrx& betaDenInfor  = denMtrxShellInfor.getMtrx(1);
	Double braPMax = alphaDenInfor.val(bra1Index,bra2Index)+betaDenInfor.val(bra1Index,bra2Index);
	Double ketPMax = alphaDenInfor.val(ket1Index,ket2Index)+betaDenInfor.val(ket1Index,ket2Index);
	return braPMax*ketPMax;
}

Double DenMtrxInfor::getKPPairMax(const ShellPair& braSP, const ShellPair& ketSP) const 
{
	// get the shell index
	UInt bra1Index = -1;
	UInt bra2Index = -1;
	UInt ket1Index = -1;
	UInt ket2Index = -1;
	braSP.getGlobalShellIndex(bra1Index,bra2Index);
	ketSP.getGlobalShellIndex(ket1Index,ket2Index);

	// now get the pMax value
	Double maxDenVal = ZERO;
	for(UInt iSpin=0; iSpin<denMtrxShellInfor.getNSpin(); iSpin++) {

		// get the data
		const Mtrx& denInfor = denMtrxShellInfor.getMtrx(iSpin);
		Double pMax13 = denInfor.val(bra1Index,ket1Index);
		Double pMax23 = denInfor.val(bra2Index,ket1Index);
		Double pMax14 = denInfor.val(bra1Index,ket2Index);
		Double pMax24 = denInfor.val(bra2Index,ket2Index);
		Double pPair1324 = pMax13*pMax24;
		Double pPair1423 = pMax14*pMax23;

		// compare
		Double pMax  = pPair1324>pPair1423 ? pPair1324 : pPair1423;
		if (pMax>maxDenVal) maxDenVal = pMax;
	}

	// return
	return maxDenVal;
}

bool DenMtrxInfor::isJSig(const Double& thresh, const SigAtomShellPairInfor& sigBraAtomSPInfor, 
		const SigAtomShellPairInfor& sigKetAtomSPInfor) const
{
	// get matrix data
	const Mtrx& alphaDenInfor = denMtrxAtomInfor.getMtrx(0);
	const Mtrx& betaDenInfor  = denMtrxAtomInfor.getMtrx(1);

	// get the max P12 and test it
	UInt atomShellBra1 = sigBraAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellBra2 = sigBraAtomSPInfor.getColAtomShellIndex(); 
	Double pMax12 = alphaDenInfor.val(atomShellBra1,atomShellBra2)+betaDenInfor.val(atomShellBra1,atomShellBra2);
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax12,thresh)) {
		return true;
	}
	// let's see the other side of coulomb digestion 
	UInt atomShellKet1 = sigKetAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellKet2 = sigKetAtomSPInfor.getColAtomShellIndex(); 
	Double pMax34 = alphaDenInfor.val(atomShellKet1,atomShellKet2)+betaDenInfor.val(atomShellKet1,atomShellKet2);
	if (! sigKetAtomSPInfor.isNegligible(sigBraAtomSPInfor,pMax34,thresh)) {
		return true;
	}

	// now we see J digestion is not significant
	return false;
}

bool DenMtrxInfor::isKSig(const UInt& iSpin, const Double& thresh, 
		const SigAtomShellPairInfor& sigBraAtomSPInfor, 
		const SigAtomShellPairInfor& sigKetAtomSPInfor) const
{
	// get the atom shell index
	UInt atomShellBra1 = sigBraAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellBra2 = sigBraAtomSPInfor.getColAtomShellIndex(); 
	UInt atomShellKet1 = sigKetAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellKet2 = sigKetAtomSPInfor.getColAtomShellIndex(); 

	// get matrix data
	const Mtrx& alphaDenInfor = denMtrxAtomInfor.getMtrx(0);
	const Mtrx& betaDenInfor  = denMtrxAtomInfor.getMtrx(1);

	// check 13 digestion
	Double pMax = alphaDenInfor.val(atomShellBra1,atomShellKet1);
	if (iSpin == 1) {
		pMax = betaDenInfor.val(atomShellBra1,atomShellKet1);
	}
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax,thresh)) {
		return true;
	}

	// 14 digestion
	pMax = alphaDenInfor.val(atomShellBra1,atomShellKet2);
	if (iSpin == 1) {
		pMax = betaDenInfor.val(atomShellBra1,atomShellKet2);
	}
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax,thresh)) {
		return true;
	}

	// 23 digestion
	pMax = alphaDenInfor.val(atomShellBra2,atomShellKet1);
	if (iSpin == 1) {
		pMax = betaDenInfor.val(atomShellBra2,atomShellKet1);
	}
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax,thresh)) {
		return true;
	}

	// finally let's see K24
	pMax = alphaDenInfor.val(atomShellBra2,atomShellKet2);
	if (iSpin == 1) {
		pMax = betaDenInfor.val(atomShellBra2,atomShellKet2);
	}
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax,thresh)) {
		return true;
	}

	// finally, if we arrive here 
	// just return false
	return false;
}

bool DenMtrxInfor::isJDerivSig(const Double& thresh, const SigAtomShellPairInfor& sigBraAtomSPInfor, 
		const SigAtomShellPairInfor& sigKetAtomSPInfor) const
{
	// get matrix data
	const Mtrx& alphaDenInfor = denMtrxAtomInfor.getMtrx(0);
	const Mtrx& betaDenInfor  = denMtrxAtomInfor.getMtrx(1);

	// get the max P12 and test it
	UInt atomShellBra1 = sigBraAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellBra2 = sigBraAtomSPInfor.getColAtomShellIndex(); 
	UInt atomShellKet1 = sigKetAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellKet2 = sigKetAtomSPInfor.getColAtomShellIndex(); 
	Double pMax12 = alphaDenInfor.val(atomShellBra1,atomShellBra2)+betaDenInfor.val(atomShellBra1,atomShellBra2);
	Double pMax34 = alphaDenInfor.val(atomShellKet1,atomShellKet2)+betaDenInfor.val(atomShellKet1,atomShellKet2);
	Double pMax2  = pMax12*pMax34;
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax2,thresh)) {
		return true;
	}

	// now we see J digestion is not significant
	return false;
}

bool DenMtrxInfor::isKDerivSig(const UInt& iSpin, const Double& thresh, 
		const SigAtomShellPairInfor& sigBraAtomSPInfor, 
		const SigAtomShellPairInfor& sigKetAtomSPInfor) const
{
	// get the atom shell index
	UInt atomShellBra1 = sigBraAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellBra2 = sigBraAtomSPInfor.getColAtomShellIndex(); 
	UInt atomShellKet1 = sigKetAtomSPInfor.getRowAtomShellIndex();
	UInt atomShellKet2 = sigKetAtomSPInfor.getColAtomShellIndex(); 

	// get matrix data
	const Mtrx& alphaDenInfor = denMtrxAtomInfor.getMtrx(0);
	const Mtrx& betaDenInfor  = denMtrxAtomInfor.getMtrx(1);

	// get the pMax values
	Double pMax13 = alphaDenInfor.val(atomShellBra1,atomShellKet1);
	Double pMax14 = alphaDenInfor.val(atomShellBra1,atomShellKet2);
	Double pMax23 = alphaDenInfor.val(atomShellBra2,atomShellKet1);
	Double pMax24 = alphaDenInfor.val(atomShellBra2,atomShellKet2);
	if (iSpin == 1) {
		pMax13 = betaDenInfor.val(atomShellBra1,atomShellKet1);
		pMax14 = betaDenInfor.val(atomShellBra1,atomShellKet2);
		pMax23 = betaDenInfor.val(atomShellBra2,atomShellKet1);
		pMax24 = betaDenInfor.val(atomShellBra2,atomShellKet2);
	}

	// now let's judge,firstly it's P13*P24
	Double pMax = pMax13*pMax24;
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax,thresh)) {
		return true;
	}

	// P23*P14
	pMax = pMax23*pMax14;
	if (! sigBraAtomSPInfor.isNegligible(sigKetAtomSPInfor,pMax,thresh)) {
		return true;
	}

	// finally, if we arrive here 
	// just return false
	return false;
}

