/**
 * CPP files corresponding to the atomfockmtrx.h
 * \brief   describing the fock matrix pieces with atom shell dimension
 * \author  Fenglai Liu
 */
#include<cstdio>
#include<iostream>
#include "excep.h"
#include "shell.h"
#include "shellpair.h"
#include "gintsinfor.h"
#include "symmjkdigest.h"
#include "integraljobs.h"
#include "blockmatrixlist.h"
#include "atomfockmtrx.h"
using namespace excep;
using namespace shell;
using namespace shellpair;
using namespace gintsinfor;
using namespace symmjkdigest;
using namespace integraljobs;
using namespace blockmatrixlist;
using namespace atomfockmtrx;

void AtomFockMtrx::initInfor(const GIntJobInfor& infor, 
		const AtomShellPair& braAtomSP,const AtomShellPair& ketAtomSP)
{
	// get the dimension etc. information for block matrix results
	// remember that for atom shell quartet (12|34), no shell pair switch
	// therefore bra and ket positions are fixed
	UInt bra1Offset = -1;
	UInt bra2Offset = -1;
	UInt ket1Offset = -1;
	UInt ket2Offset = -1;
	braAtomSP.getGlobalCarBasOffSet(bra1Offset,bra2Offset);
	ketAtomSP.getGlobalCarBasOffSet(ket1Offset,ket2Offset);
	UInt bra1NCarBas = -1;
	UInt bra2NCarBas = -1;
	UInt ket1NCarBas = -1;
	UInt ket2NCarBas = -1;
	braAtomSP.getNCarBas(bra1NCarBas,bra2NCarBas);
	ketAtomSP.getNCarBas(ket1NCarBas,ket2NCarBas);

	// now let's reset the data
	UInt job =  infor.getIntJob();
	if (doK(job)) {
		aK13.init(bra1NCarBas,ket1NCarBas,bra1Offset,ket1Offset); 
		aK14.init(bra1NCarBas,ket2NCarBas,bra1Offset,ket2Offset); 
		aK23.init(bra2NCarBas,ket1NCarBas,bra2Offset,ket1Offset);    
		aK24.init(bra2NCarBas,ket2NCarBas,bra2Offset,ket2Offset);   

		// beta
		if (infor.getNSpin() == 2) {
			bK13.init(bra1NCarBas,ket1NCarBas,bra1Offset,ket1Offset); 
			bK14.init(bra1NCarBas,ket2NCarBas,bra1Offset,ket2Offset); 
			bK23.init(bra2NCarBas,ket1NCarBas,bra2Offset,ket1Offset);    
			bK24.init(bra2NCarBas,ket2NCarBas,bra2Offset,ket2Offset);   
		}
	}
	
	// coulomb
	if (doJ(job)) {
		J12.init(bra1NCarBas,bra2NCarBas,bra1Offset,bra2Offset);     
		J34.init(ket1NCarBas,ket2NCarBas,ket1Offset,ket2Offset);    
	}
}

void AtomFockMtrx::resetLinkInfor(bool switchSP, const ShellPair& braSP, const ShellPair& ketSP)
{
	//
	// in this function, according to shell pair situation
	// we need to look up how to connect the results in shell dimension
	// with the results here in atom shell dimension
	//
	// 1,2,3,4 represents the original shell position
	//
	
	// work on J information
	if (switchSP) {
		// this case may include (34|12), (34|21), (43|12) and (43|21) 
		j34Pos     = LINKING_J12;          
		j12Pos     = LINKING_J34;          
		j34InTrans = false;      
		j12InTrans = false;      
		if (braSP.inverseShells()) j34InTrans = true;      
		if (ketSP.inverseShells()) j12InTrans = true;      
	}else{
		// this case may include (12|34), (12|43), (21|34) and (21|43) 
		j12Pos     = LINKING_J12;          
		j34Pos     = LINKING_J34;          
		j12InTrans = false;      
		j34InTrans = false;      
		if (braSP.inverseShells()) j12InTrans = true;      
		if (ketSP.inverseShells()) j34InTrans = true;      
	}

	// work on K information
	if (switchSP) {

		// if bra and ket are switched, because the 
		// exchange matrix result takes row index from
		// bra and col index from ket, then row and col
		// must all switched
		// therefore for all of results in shell dimension,
		// they must be in transpose form of the block matrix
		// here
		k13InTrans = true;      
		k14InTrans = true;      
		k23InTrans = true;      
		k24InTrans = true;      

		// determine the linking of each sub matrix
		if (braSP.inverseShells()) {
			if (ketSP.inverseShells()) {
				// (43|21)
				k13Pos     = LINKING_K24;          
				k14Pos     = LINKING_K14;          
				k23Pos     = LINKING_K23;          
				k24Pos     = LINKING_K13;          
			}else{
				// (34|21)
				k13Pos     = LINKING_K23;          
				k14Pos     = LINKING_K13;          
				k23Pos     = LINKING_K24;          
				k24Pos     = LINKING_K14;          
			}
		}else{
			if (ketSP.inverseShells()) {
				// (43|12)
				k13Pos     = LINKING_K14;          
				k14Pos     = LINKING_K24;          
				k23Pos     = LINKING_K13;          
				k24Pos     = LINKING_K23;          
			}else{
				// (34|12)
				k13Pos     = LINKING_K13;          
				k14Pos     = LINKING_K23;          
				k23Pos     = LINKING_K14;          
				k24Pos     = LINKING_K24;          
			}
		}	
	}else{

		// for similar reasons, all  of the sub matrix
		// must be in normal status and not in transpose
		k13InTrans = false;      
		k14InTrans = false;      
		k23InTrans = false;      
		k24InTrans = false;

		// work on linking information
		if (braSP.inverseShells()) {
			if (ketSP.inverseShells()) {
				// (21|43)
				k13Pos     = LINKING_K24;          
				k14Pos     = LINKING_K23;          
				k23Pos     = LINKING_K14;          
				k24Pos     = LINKING_K13;          
			}else{
				// (21|34)
				k13Pos     = LINKING_K23;          
				k14Pos     = LINKING_K24;          
				k23Pos     = LINKING_K13;          
				k24Pos     = LINKING_K14;          
			}
		}else{
			if (ketSP.inverseShells()) {
				// (12|43)
				k13Pos     = LINKING_K14;          
				k14Pos     = LINKING_K13;          
				k23Pos     = LINKING_K24;          
				k24Pos     = LINKING_K23;          
			}else{
				// (12|34)
				k13Pos     = LINKING_K13;          
				k14Pos     = LINKING_K14;          
				k23Pos     = LINKING_K23;          
				k24Pos     = LINKING_K24;          
			}
		}	
	}
}

void AtomFockMtrx::updateResults(const GIntJobInfor& infor, const SymmJKDigest& localResults)
{
	// get the integral job
	UInt job =  infor.getIntJob();

	// firstly let's do J
	if (doJ(job)) {
		const BlockMtrx& j12 = localResults.getJ12();
		if (j12Pos == LINKING_J12) {
			j12.updateData(J12,j12InTrans);
		}else{
			j12.updateData(J34,j12InTrans);
		}
		const BlockMtrx& j34 = localResults.getJ34();
		if (j34Pos == LINKING_J12) {
			j34.updateData(J12,j34InTrans);
		}else{
			j34.updateData(J34,j34InTrans);
		}
	}

	// now let's do exchange
	if (doK(job)) {

		// alpha spin
		// k13 update
		UInt iSpin=0;
		const BlockMtrx& k13 = localResults.getK13(iSpin);
		if (k13Pos == LINKING_K13) {
			k13.updateData(aK13,k13InTrans);
		}else if (k13Pos == LINKING_K14) {
			k13.updateData(aK14,k13InTrans);
		}else if (k13Pos == LINKING_K23) {
			k13.updateData(aK23,k13InTrans);
		}else{
			k13.updateData(aK24,k13InTrans);
		}

		// k14 update
		const BlockMtrx& k14 = localResults.getK14(iSpin);
		if (k14Pos == LINKING_K13) {
			k14.updateData(aK13,k14InTrans);
		}else if (k14Pos == LINKING_K14) {
			k14.updateData(aK14,k14InTrans);
		}else if (k14Pos == LINKING_K23) {
			k14.updateData(aK23,k14InTrans);
		}else{
			k14.updateData(aK24,k14InTrans);
		}

		// k23 update
		const BlockMtrx& k23 = localResults.getK23(iSpin);
		if (k23Pos == LINKING_K13) {
			k23.updateData(aK13,k23InTrans);
		}else if (k23Pos == LINKING_K14) {
			k23.updateData(aK14,k23InTrans);
		}else if (k23Pos == LINKING_K23) {
			k23.updateData(aK23,k23InTrans);
		}else{
			k23.updateData(aK24,k23InTrans);
		}

		// k24 update
		const BlockMtrx& k24 = localResults.getK24(iSpin);
		if (k24Pos == LINKING_K13) {
			k24.updateData(aK13,k24InTrans);
		}else if (k24Pos == LINKING_K14) {
			k24.updateData(aK14,k24InTrans);
		}else if (k24Pos == LINKING_K23) {
			k24.updateData(aK23,k24InTrans);
		}else{
			k24.updateData(aK24,k24InTrans);
		}

		// do we have beta spin?
		if (infor.getNSpin() == 2) {
			iSpin=1;

			// k13 update
			const BlockMtrx& k13 = localResults.getK13(iSpin);
			if (k13Pos == LINKING_K13) {
				k13.updateData(bK13,k13InTrans);
			}else if (k13Pos == LINKING_K14) {
				k13.updateData(bK14,k13InTrans);
			}else if (k13Pos == LINKING_K23) {
				k13.updateData(bK23,k13InTrans);
			}else{
				k13.updateData(bK24,k13InTrans);
			}

			// k14 update
			const BlockMtrx& k14 = localResults.getK14(iSpin);
			if (k14Pos == LINKING_K13) {
				k14.updateData(bK13,k14InTrans);
			}else if (k14Pos == LINKING_K14) {
				k14.updateData(bK14,k14InTrans);
			}else if (k14Pos == LINKING_K23) {
				k14.updateData(bK23,k14InTrans);
			}else{
				k14.updateData(bK24,k14InTrans);
			}

			// k23 update
			const BlockMtrx& k23 = localResults.getK23(iSpin);
			if (k23Pos == LINKING_K13) {
				k23.updateData(bK13,k23InTrans);
			}else if (k23Pos == LINKING_K14) {
				k23.updateData(bK14,k23InTrans);
			}else if (k23Pos == LINKING_K23) {
				k23.updateData(bK23,k23InTrans);
			}else{
				k23.updateData(bK24,k23InTrans);
			}

			// k24 update
			const BlockMtrx& k24 = localResults.getK24(iSpin);
			if (k24Pos == LINKING_K13) {
				k24.updateData(bK13,k24InTrans);
			}else if (k24Pos == LINKING_K14) {
				k24.updateData(bK14,k24InTrans);
			}else if (k24Pos == LINKING_K23) {
				k24.updateData(bK23,k24InTrans);
			}else{
				k24.updateData(bK24,k24InTrans);
			}
		}
	}
}

void AtomFockMtrx::merge(const GIntJobInfor& infor, BlockMtrxList& list, UInt iSpin) const
{
	// get the integral job
	UInt job =  infor.getIntJob();

	// firstly let's do J
	if (doJ(job)) {
		list.update(J12);
		list.update(J34);
	}

	// now let's do K
	if (doK(job)) {
		if (iSpin == 0) {
			list.update(aK13);
			list.update(aK14);
			list.update(aK23);
			list.update(aK24);
		}else{
			list.update(bK13);
			list.update(bK14);
			list.update(bK23);
			list.update(bK24);
		}
	}
}
