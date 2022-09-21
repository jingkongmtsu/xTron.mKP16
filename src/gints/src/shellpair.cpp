/**
 * CPP files corresponding to the shellpair.h
 * \brief   describing the shell pairs data for two given atom shells
 * \author  Fenglai Liu
 */
#include<cstdio>
#include<iostream>
#include "excep.h"
#include "shell.h"
#include "shellprop.h"
#include "oneintformula.h"
#include "sigshellpairinfor.h"
#include "shellpair.h"
using namespace excep;
using namespace shell;
using namespace shellprop;
using namespace oneintformula;
using namespace sigshellpairinfor; 
using namespace shellpair;

void ShellPair::init(const Shell& is, const Shell& js, const Double* ixyz, 
		const Double* jxyz, const Double& threshold)
{
	// shall we inverse the shells?
	shellInversed = false;
	if (is<js) shellInversed = true;

	// if shell are inversed, then we need to switch the data
	// also perform the screening work
	if(shellInversed) {

		// angular momentum
		pureShellI = js.isPure();
		pureShellJ = is.isPure();
		iLmin = js.getLmin();
		iLmax = js.getLmax();
		jLmin = is.getLmin();
		jLmax = is.getLmax();
		LCode = codeL(iLmin,iLmax,jLmin,jLmax);

		// basis set
		iBasOffset     = js.getLocalBasisIndex(0,TYPE_NORM);
		jBasOffset     = is.getLocalBasisIndex(0,TYPE_NORM);
		iGBasOffset    = js.getBasisIndex(0,TYPE_NORM);
		jGBasOffset    = is.getBasisIndex(0,TYPE_NORM);
		iNBas          = js.getNBas();
		jNBas          = is.getNBas();
		iCarBasOffset  = js.getLocalBasisIndex(0,TYPE_CART);
		jCarBasOffset  = is.getLocalBasisIndex(0,TYPE_CART);
		iGCarBasOffset = js.getBasisIndex(0,TYPE_CART);
		jGCarBasOffset = is.getBasisIndex(0,TYPE_CART);
		iNCarBas       = js.getNCarBas();
		jNCarBas       = is.getNCarBas();

		// shell index
		iShellIndex    = js.getLocalShellIndex();
		jShellIndex    = is.getLocalShellIndex();
		iGShellIndex   = js.getGlobalShellIndex();
		jGShellIndex   = is.getGlobalShellIndex();

		// center
		A[0] = jxyz[0];
		A[1] = jxyz[1];
		A[2] = jxyz[2];
		B[0] = ixyz[0];
		B[1] = ixyz[1];
		B[2] = ixyz[2];

		// now form the final primitives data
		sigPrimPairInfor.assign(is.getNPrim()*js.getNPrim(),1);
		formPrimPairsData(js,is,threshold);

	}else{

		// angular momentum
		pureShellI = is.isPure();
		pureShellJ = js.isPure();
		iLmin = is.getLmin();
		iLmax = is.getLmax();
		jLmin = js.getLmin();
		jLmax = js.getLmax();
		LCode = codeL(iLmin,iLmax,jLmin,jLmax);

		// basis set
		iBasOffset     = is.getLocalBasisIndex(0,TYPE_NORM);
		jBasOffset     = js.getLocalBasisIndex(0,TYPE_NORM);
		iGBasOffset    = is.getBasisIndex(0,TYPE_NORM);
		jGBasOffset    = js.getBasisIndex(0,TYPE_NORM);
		iNBas          = is.getNBas();
		jNBas          = js.getNBas();
		iCarBasOffset  = is.getLocalBasisIndex(0,TYPE_CART);
		jCarBasOffset  = js.getLocalBasisIndex(0,TYPE_CART);
		iGCarBasOffset = is.getBasisIndex(0,TYPE_CART);
		jGCarBasOffset = js.getBasisIndex(0,TYPE_CART);
		iNCarBas       = is.getNCarBas();
		jNCarBas       = js.getNCarBas();

		// shell index
		iShellIndex    = is.getLocalShellIndex();
		jShellIndex    = js.getLocalShellIndex();
		iGShellIndex   = is.getGlobalShellIndex();
		jGShellIndex   = js.getGlobalShellIndex();

		// center
		A[0] = ixyz[0];
		A[1] = ixyz[1];
		A[2] = ixyz[2];
		B[0] = jxyz[0];
		B[1] = jxyz[1];
		B[2] = jxyz[2];

		// now form the final primitives data
		sigPrimPairInfor.assign(is.getNPrim()*js.getNPrim(),1);
		formPrimPairsData(is,js,threshold);
	}

}

void ShellPair::formPrimPairsData(const Shell& is, const Shell& js, const Double& thresh) 
{

	//
	// whether this primitive pair is significant?
	//
	// this is judged by forming the overlap for the given primitive
	// pair data so that to see whether the primitive pairs is significant
	// or not
	//
	// Basically, we compute the overlap between two "S" type primitive 
	// functions even though the two primitive functions have their own
	// angular momentum. We do not compute the real overlap integrals
	// here, that is because; it's difficult to find out the specific 
	// angular momentum pair which produces the largest overlap except
	// we produce all of overlap integrals and find the largest value.
	//
	// For example, for two D-D primitive pairs, due to it's geometry
	// the largest overlap may come from the Dzz-Dzz pair, it's possible
	// that Dxx-Dxx could be the smallest. However, we may take the 
	// Dxx-Dxx as a trial test to see whether the primitive pair is 
	// significant. Then if in this case the Dzz-Dzz integral is 
	// 100 times larger than the Dxx-Dxx one? Our calculation could 
	// introduce larger error here.
	//
	// On the other hand, the S-S type of integral is only determined 
	// by distance of two primitive functions, S type of primitive 
	// is symmetrical for all of directions. Therefore it could be 
	// well expected that S-S type of integral may not be the largest 
	// integral, and it's not the smallest one, neither. Therefore,
	// S-S type of integral is better for evaluation the significance 
	// of the primitive pairs.
	//
	// We do not count in the coefficients here, since they are already
	// been normalized according to the shell angular momentum.
	//
	// finally, we note that the input shell should obey the order below:
	// L(is) >= L(js)
	//
	const Double* iexp = is.getExp();
	const Double* jexp = js.getExp();
	Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);

	//
	// now let's go to see each primitive pair
	//
	np2 = 0;
	UInt count = 0;
	for(UInt jp=0;jp<js.getNPrim();jp++) {
		for(UInt ip=0;ip<is.getNPrim();ip++) {

			// pre-calculation
			// exponents and overlap integral
			Double ia   = iexp[ip];
			Double ja   = jexp[jp];
			Double alpha= ONE/(ia+ja); 
			Double ab   = -ia*ja*alpha;
			//yw, jk. The new way exp(ab*AB2) factor is combined with 
			//contraction (c2) instead of being part of fac. This change is for hfx hole so
			//that we don't need to use 'fac'.  'fac' is still used in hgp_os for now.
			//Still needs Fenglai to check because I don't know exactly how this 'fac'
			//is used exactly in the hgp_os code.  It could be used for screening, but 
			//it seems to multiply c2 always anyway because our modification here did
			//not change results.
			//Double ov   = exp(ab*AB2)*pow(PI*alpha,1.5E0); 
			Double ov   = pow(PI*alpha,1.5E0); //yw modify
		       
			// form P point according to the 
			// Gaussian pritimive product theorem
			Double adab = ia*alpha; 
			Double bdab = ja*alpha; 
			Double Px   = A[0]*adab + B[0]*bdab;
			Double Py   = A[1]*adab + B[1]*bdab;
			Double Pz   = A[2]*adab + B[2]*bdab;

			//The primitive pair reduction is disabled because the hfx hole part needs
			//all the primitive pairs there.
			//!!! Needs to multiply exp(ab*AB2) if to be used again!
			// trial test to see the significance
			//if (fabs(ov)<thresh) {
			//	sigPrimPairInfor[count] = 0;
			//	count++;
			//		continue;
			//}

			// now everything is significant
			// we push it into result data
			e2[np2]     = alpha;
			fac[np2]    = ov;
			e2diff[np2] = ia-ja;
			P[3*np2+0]  = Px;
			P[3*np2+1]  = Py;
			P[3*np2+2]  = Pz;

			// increment the counting
			count++;
			np2++;
		} 
	} 

	// Now form coefficients array
	UInt ilmin = is.getLmin();
	UInt ilmax = is.getLmax();
	UInt jlmin = js.getLmin();
	UInt jlmax = js.getLmax();
	count = 0;
	for(UInt J=jlmin; J<=jlmax; J++){
		for(UInt I=ilmin; I<=ilmax; I++){

			// get the coefficients for each sub-shell
			const Double* ic = is.getCoe(I);
			const Double* jc = js.getCoe(J);

			UInt index=0;
			for(UInt j=0;j<js.getNPrim();j++) {
				for(UInt i=0;i<is.getNPrim();i++) {

					//The primitive pair reduction is disabled because the hfx hole part needs.
					//See comments above.
					// whether this primitive pair is significant?
					//if (sigPrimPairInfor[index] == 0) {
					//	index++;
					//	continue;
					//}
					Double ia   = iexp[i]; //yw modify
					Double ja   = jexp[j]; //yw modify
					Double alpha= ONE/(ia+ja); //yw modify
					Double ab   = -ia*ja*alpha; //yw modify
				       
					// now push in data
					c2[count] = ic[i]*jc[j]*exp(ab*AB2); //yw modify
					index++;
					count++;
				}
			}
		}
	}
}

void ShellPair::print(UInt iprint) const {

	// print out the angular momentuum information and basis sets
	string iName = getShellName(iLmin,iLmax);
	string jName = getShellName(jLmin,jLmax);
	cout << "Shell i is: " << iName << " Shell j is: " << jName << endl;
	if (shellInversed) {
		cout << "The input shells have been inversed " << endl;
	}
	cout << "Contraction degree is                           : " << getNP2() << endl;
	cout << "Normal global basis set offset for shell i is   : " << iGBasOffset << endl;
	cout << "Normal global basis set offset for shell j is   : " << jGBasOffset << endl;
	cout << "Cartesian global basis set offset for shell i is: " << iGCarBasOffset << endl;
	cout << "Cartesian global basis set offset for shell j is: " << jGCarBasOffset << endl;
	cout << "Normal local  basis set offset for shell i is   : " << iBasOffset << endl;
	cout << "Normal local  basis set offset for shell j is   : " << jBasOffset << endl;
	cout << "Cartesian local basis set offset for shell i is : " << iCarBasOffset << endl;
	cout << "Cartesian local basis set offset for shell j is : " << jCarBasOffset << endl;
	cout << "local  shell index for shell i is   : " << iShellIndex << endl;
	cout << "local  shell index for shell j is   : " << jShellIndex << endl;
	cout << "global shell index for shell i is   : " << iGShellIndex << endl;
	cout << "global shell index for shell j is   : " << jGShellIndex << endl;

	// print out center information
	printf("%-12s %-12s %-12s %-12s %-12s %-12s\n", "Ax", "Ay", "Az", "Bx", "By", "Bz");
	printf("%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n", A[0], A[1], A[2], B[0], B[1], B[2]); 

	// now primitive information
	if (iprint >= 1) {
		cout << "Gaussian primitives pairs information for this shell pair: \n"; 
		UInt np = 0;
		for(UInt J=jLmin; J<=jLmax; J++){
			for(UInt I=iLmin; I<=iLmax; I++){
				printf("Angular mom: %d, %d\n", (Int)I ,(Int)J);
				printf("%-12s %-12s %-15s %-12s %-12s %-12s\n", "Coef1*Coef2", 
						"ONE/ExpSum", "Overlap", "Px", "Py", "Pz");
				for(UInt i=0;i<getNP2();i++) {
					printf("%-12.6f %-12.6f %-15.12f %-12.6f %-12.6f %-12.6f\n", 
							c2[np+i], e2[i], fac[i], P[3*i], P[3*i+1], P[3*i+2]);
				}
				cout << endl;
				np += getNP2();
			}
		}
	}
}

void AtomShellPair::init(const MolShell& rs, const MolShell& cs, 
		const SigAtomShellPairInfor& infor, const Double& thresh)
{
	// get shell data
	rowAtomShellIndex = infor.getRowAtomShellIndex();
	colAtomShellIndex = infor.getColAtomShellIndex();
	const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
	const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);

	// dimension information
	iGBasOffset    = rowAtomShell.getBasisStartIndex(TYPE_NORM);
	jGBasOffset    = colAtomShell.getBasisStartIndex(TYPE_NORM);
	iNBas          = rowAtomShell.getNBas();
	jNBas          = colAtomShell.getNBas();
	iGCarBasOffset = rowAtomShell.getBasisStartIndex(TYPE_CART);
	jGCarBasOffset = colAtomShell.getBasisStartIndex(TYPE_CART);
	iNCarBas       = rowAtomShell.getNCarBas();
	jNCarBas       = colAtomShell.getNCarBas();

	// form shell pair data
	nShellPairs = 0;
	const Double* iXYZ = rowAtomShell.getXYZ();
	const Double* jXYZ = colAtomShell.getXYZ();
	for(UInt iSP=0;iSP<infor.getNSigShellPairs();iSP++) {
		const SigShellPairInfor& info = infor.getSigShellPairInfor(iSP);
		const Shell& iS = rowAtomShell.getShell(info.getRowShellIndex());
		const Shell& jS = colAtomShell.getShell(info.getColShellIndex());
		spList[nShellPairs].init(iS,jS,iXYZ,jXYZ,thresh);
		nShellPairs++;
	}
}

#ifdef GINTS_PHI_DEBUG

/*
 *
bool AtomShellPair::same(const UInt& maxNP2, const UInt& maxNL, const Double& thresh, 
		const MolShell& rs, const MolShell& cs, const CompactShellPairData_Phi& comASP) const
{
	//
	// this function is used for debugging purpose
	// to test that whether the atom shell pair data 
	// is same with the input comASP
	//
	// however, currently because the data inside they are 
	// not same, so if they are not same it does not mean
	// the comASP is wrong
	//
	// if fatal error occurs, we return false
	// however, on the other hand if data is not same;
	// but it's reasonable; we will print out details
	// and return true
	//

	// get shell data
	const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
	const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);
	ShellPair spTest(maxNP2,maxNL);

	// let's see whether the total number of shell pairs are same
	UInt nTolSP = 0;
	UInt nAngCode = compactshellpairdata_phi::getNAngCode(comASP);
	for(UInt iAng=0; iAng<nAngCode; iAng++) {
		nTolSP += compactshellpairdata_phi::getNSP(comASP,iAng);
	}
	UInt nSigSP = 0;
	for (UInt iSP=0; iSP<getNShellPairs(); iSP++) {
		const ShellPair& sp = getShellPair(iSP);
		if (sp.isSig()) nSigSP++;
	}
	if (nTolSP != nSigSP) {
		printf("the number of shell pair data is not same, the compact one has %d, and normal one has %d\n",
				(Int)nTolSP, (Int)nSigSP);
		if (nTolSP < nSigSP) {

			// now let's search
			printf("there are the shell pairs appear in the normal one, but missing in the compact one\n");
			for (UInt i=0; i<getNShellPairs(); i++) {
				const ShellPair& sp = getShellPair(i);
				if (! sp.isSig()) continue;
				UInt rowOffsetIndex0, colOffsetIndex0;
				sp.getLocalCarBasOffSet(rowOffsetIndex0,colOffsetIndex0);
				if (sp.inverseShells()) {
					UInt tmp = rowOffsetIndex0;
					rowOffsetIndex0 = colOffsetIndex0;
					colOffsetIndex0 = tmp;
				}
				UInt iSPIndex = getNShellPairs()+1;
				for(UInt iAng=0; iAng<nAngCode; iAng++) {
					UInt code   = compactshellpairdata_phi::getAngCode(comASP,iAng);
					UInt nSP    = compactshellpairdata_phi::getNSP(comASP,iAng);
					UInt offset = compactshellpairdata_phi::getSPOffset(comASP,iAng);
					for(UInt iSP=offset; iSP<offset+nSP; iSP++) {
						const UInt* locBasOffset = compactshellpairdata_phi::getLocCarBasOffsetArray(comASP,iSP);
						UInt rowOffsetIndex = locBasOffset[0];
						UInt colOffsetIndex = locBasOffset[1];
						const UInt* shellStatusArray = compactshellpairdata_phi::getShellStatusArray(comASP,iSP);
						UInt status = shellStatusArray[0];
						if (status == WITH_SHELL_SWITCH) {
							UInt tmp = rowOffsetIndex;
							rowOffsetIndex = colOffsetIndex;
							colOffsetIndex = tmp;
						}
						if (rowOffsetIndex == rowOffsetIndex0 && colOffsetIndex == colOffsetIndex0) {
							iSPIndex = iSP;
							break;
						}
					}
				}
				if (iSPIndex>getNShellPairs()) {
					sp.print(3);
				}
			}

		}else{
			// now let's search
			printf("there are the shell pairs appear in the compact one, but missing in the normal one\n");
			for(UInt iAng=0; iAng<nAngCode; iAng++) {
				UInt code   = compactshellpairdata_phi::getAngCode(comASP,iAng);
				UInt nSP    = compactshellpairdata_phi::getNSP(comASP,iAng);
				UInt offset = compactshellpairdata_phi::getSPOffset(comASP,iAng);
				UInt iLmin, iLmax, jLmin, jLmax;
				decodeL(code,iLmin, iLmax, jLmin, jLmax);
				const Double* c2 = NULL;
				if (iLmin != iLmax || jLmin != jLmax) {
					c2 = compactshellpairdata_phi::getC2(comASP,code);
				}
				for(UInt iSP=offset; iSP<offset+nSP; iSP++) {
					const UInt* locBasOffset = compactshellpairdata_phi::getLocCarBasOffsetArray(comASP,iSP);
					UInt rowOffsetIndex = locBasOffset[0];
					UInt colOffsetIndex = locBasOffset[1];
					const UInt* shellStatusArray = compactshellpairdata_phi::getShellStatusArray(comASP,iSP);
					UInt status = shellStatusArray[0];
					if (status == WITH_SHELL_SWITCH) {
						UInt tmp = rowOffsetIndex;
						rowOffsetIndex = colOffsetIndex;
						colOffsetIndex = tmp;
					}
					if (c2 != NULL) {
						UInt np2 = compactshellpairdata_phi::getNP2(comASP,iSP);
						UInt nC = 2;
						if (iLmin != iLmax && jLmin != jLmax) {
							nC = 4;
						}
						c2 += nC*np2;
					}
					UInt iSPIndex = getNShellPairs()+1;
					for (UInt i=0; i<getNShellPairs(); i++) {
						const ShellPair& sp = getShellPair(i);
						if (! sp.isSig()) continue;
						UInt rowOffsetIndex0, colOffsetIndex0;
						sp.getLocalCarBasOffSet(rowOffsetIndex0,colOffsetIndex0);
						if (sp.inverseShells()) {
							UInt tmp = rowOffsetIndex0;
							rowOffsetIndex0 = colOffsetIndex0;
							colOffsetIndex0 = tmp;
						}
						if (rowOffsetIndex == rowOffsetIndex0 && colOffsetIndex == colOffsetIndex0) {
							iSPIndex = iSP;
							break;
						}
					}
					if (iSPIndex>getNShellPairs()) {
						compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
					}
				}
			}
		}
		return true;
	}

	// for the comASP, let's loop over the angular momentum type
	const Double* A = compactshellpairdata_phi::getA(comASP);
	const Double* B = compactshellpairdata_phi::getB(comASP);
	for(UInt iAng=0; iAng<nAngCode; iAng++) {

		// get the angular momentum code
		// we note that the arrangement of the shell pair data
		// is same with the data in comASP
		// so let's derive the shell pair index
		UInt code   = compactshellpairdata_phi::getAngCode(comASP,iAng);
		UInt nSP    = compactshellpairdata_phi::getNSP(comASP,iAng);
		UInt offset = compactshellpairdata_phi::getSPOffset(comASP,iAng);

		// let's see whether it's composite shell pair
		UInt iLmin, iLmax, jLmin, jLmax;
		decodeL(code,iLmin, iLmax, jLmin, jLmax);
		const Double* c2 = NULL;
		if (iLmin != iLmax || jLmin != jLmax) {
			c2 = compactshellpairdata_phi::getC2(comASP,code);
		}

		// now let's loop over shell pair
		for(UInt iSP=offset; iSP<offset+nSP; iSP++) {

			// this is the data for the comASP
			// let's test the data one by one
			// local basis set offset
			// we need to restore it's original order
			const UInt* locBasOffset = compactshellpairdata_phi::getLocCarBasOffsetArray(comASP,iSP);
			UInt rowOffsetIndex = locBasOffset[0];
			UInt colOffsetIndex = locBasOffset[1];
			const UInt* shellStatusArray = compactshellpairdata_phi::getShellStatusArray(comASP,iSP);
			UInt status = shellStatusArray[0];
			if (status == WITH_SHELL_SWITCH) {
				UInt tmp = rowOffsetIndex;
				rowOffsetIndex = colOffsetIndex;
				colOffsetIndex = tmp;
			}

			// let's see whether we can find the corresponding shell pair data
			UInt iSPIndex = getNShellPairs()+1;
			for (UInt i=0; i<getNShellPairs(); i++) {
				const ShellPair& sp = getShellPair(i);

				// get the data for compare
				UInt rowOffsetIndex0, colOffsetIndex0;
				sp.getLocalCarBasOffSet(rowOffsetIndex0,colOffsetIndex0);

				// we need to consider that whether the shell pair are switched
				if (sp.inverseShells()) {
					UInt tmp = rowOffsetIndex0;
					rowOffsetIndex0 = colOffsetIndex0;
					colOffsetIndex0 = tmp;
				}

				if (rowOffsetIndex != rowOffsetIndex0 || colOffsetIndex != colOffsetIndex0) {
					continue;
				}else{
					iSPIndex = i;
					break;
				}
			}

			// if we did not get it, now print out detail for a look
			if (iSPIndex > getNShellPairs()) {
				printf("We failed to find the corresponding shell pair data in comASP for this shell pair data\n");
				printf("let me generate it for detail look\n");
				const UInt* shellStatusArray = compactshellpairdata_phi::getShellStatusArray(comASP,iSP);
				UInt status = shellStatusArray[0];
				if (status == WITH_SHELL_SWITCH) {
					const Shell& jS = rowAtomShell.getShell(colOffsetIndex);
					const Shell& iS = colAtomShell.getShell(rowOffsetIndex);
					spTest.init(iS,jS,B,A,thresh);
					spTest.print(3);
				}else{
					const Shell& iS = rowAtomShell.getShell(rowOffsetIndex);
					const Shell& jS = colAtomShell.getShell(colOffsetIndex);
					spTest.init(iS,jS,A,B,thresh);
					spTest.print(3);
				}
				return false;
			}
			const ShellPair& sp = getShellPair(iSPIndex);

			// shell status
			if (status == NO_SHELL_SWITCH && sp.inverseShells()) {
				printf("failed the comASP shell is not revised, but the sp shell is reversed\n");
				compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
				sp.print(3);
				return false;
			}else if (status == WITH_SHELL_SWITCH && ! sp.inverseShells()) {
				printf("failed the comASP shell is revised, but the sp shell is not reversed\n");
				compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
				sp.print(3);
				return false;
			}

			// test the LCode
			if (code != sp.getLCode()) {
				printf("failed the angular code is not same\n");
				compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
				sp.print(3);
				return false;
			}

			// now center A and B
			const Double* A0 = sp.getA();
			const Double* B0 = sp.getB();
			if (status == WITH_SHELL_SWITCH) {
				A0 = sp.getB();
				B0 = sp.getA();
			}
			if (fabs(A[0]-A0[0])>1.0E-12 ||
					fabs(A[1]-A0[1])>1.0E-12 ||
					fabs(A[2]-A0[2])>1.0E-12 ) {
				printf("failed the center A is not same\n");
				compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
				sp.print(3);
				return false;
			}
			if (fabs(B[0]-B0[0])>1.0E-12 ||
					fabs(B[1]-B0[1])>1.0E-12 ||
					fabs(B[2]-B0[2])>1.0E-12 ) {
				printf("failed the center B is not same\n");
				compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
				sp.print(3);
				return false;
			}

			// now np2
			UInt np2 = compactshellpairdata_phi::getNP2(comASP,iSP);
			if (np2 != sp.getNP2()) {
				printf("the number of primitive pair data does not equal with each other, "
						"compact one has %d, and normal one has %d\n", (Int)np2, (Int)sp.getNP2());
				//compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
				//sp.print(3);

				// however, we need to advance the c2 pointer so that
				// later we can print out correct data
				if (c2 != NULL) {
					UInt nC = 2;
					if (iLmin != iLmax && jLmin != jLmax) {
						nC = 4;
					}
					c2 += nC*np2;
				}
			}

			// now let's test the primitive pair data 
			if (np2 == sp.getNP2()) {

				// test the e2
				const Double* e2    = compactshellpairdata_phi::getE2(comASP,iSP);
				const Double* e2_0  = sp.getE2();
				bool omitRest = false;
				for(UInt ip2=0; ip2<np2; ip2++) {
					if (fabs(e2[ip2]-e2_0[ip2])>1.0E-12) { 
						omitRest = true;
					}
				}

				// we do not need to compare the rest
				if (omitRest) {
					printf("failed the e2 array is not same\n");
					compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
					sp.print(3);
					continue;
				}

				// test the P center
				const Double* P  = compactshellpairdata_phi::getP(comASP,iSP);
				const Double* P0 = sp.getP();
				for(UInt ip2=0; ip2<np2; ip2++) {
					if (fabs(P[3*ip2]-P0[3*ip2])>1.0E-12 ||
							fabs(P[3*ip2+1]-P0[3*ip2+1])>1.0E-12 ||
							fabs(P[3*ip2+2]-P0[3*ip2+2])>1.0E-12 ) {
						printf("failed the center P is not same\n");
						compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
						sp.print(3);
						return false;
					}
				}

				// test the c2 and fac
				if (c2 != NULL) {

					// coefficient array
					UInt nC = 2;
					if (iLmin != iLmax && jLmin != jLmax) {
						nC = 4;
					}

					// set the value
					const Double* oldC2 = sp.getC2();
					bool C2Failed = false;
					for(UInt ip2=0; ip2<np2; ip2++) {

						// set the c
						Double c_1,c_2,c_3,c_4;
						if (nC == 2) {
							c_1 = c2[2*ip2  ];
							c_2 = c2[2*ip2+1];
						}else{
							c_1 = c2[4*ip2  ];
							c_2 = c2[4*ip2+1];
							c_3 = c2[4*ip2+2];
							c_4 = c2[4*ip2+3];
						}

						// get the old C
						Double oldC_1,oldC_2,oldC_3,oldC_4;
						if (nC == 2) {
							oldC_1 = oldC2[ip2];
							oldC_2 = oldC2[ip2+np2];
						}else{
							oldC_1 = oldC2[ip2];
							oldC_2 = oldC2[ip2+np2];
							oldC_3 = oldC2[ip2+2*np2];
							oldC_4 = oldC2[ip2+3*np2];
						}

						// now let's see the c2
						if (nC == 2) {
							if (fabs(c_1-oldC_1)>1.0E-12 || fabs(c_2-oldC_2)>1.0E-12) { 
								C2Failed = true;
							}
						}else{
							if (fabs(c_1-oldC_1)>1.0E-12 || fabs(c_2-oldC_2)>1.0E-12 ||
									fabs(c_3-oldC_3)>1.0E-12 || fabs(c_4-oldC_4)>1.0E-12 ) { 
								C2Failed = true;
							}
						}
					}

					// print out c2 array
					if (C2Failed) {
						printf("failed the coefficient data is not same\n");
						compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
						sp.print(3);
						return false;
					}

					// now increase the c2
					c2 += nC*np2;

					// factor array
					const Double* f    = compactshellpairdata_phi::getFac(comASP,iSP);
					const Double* f_0  = sp.getFac();
					for(UInt ip2=0; ip2<np2; ip2++) {
						if (fabs(f[ip2]-f_0[ip2])>1.0E-12) { 
							printf("failed the fac array is not same\n");
							compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
							sp.print(3);
							return false;
						}
					}
				}else{

					// factor array, this case we need to multiply the c2
					const Double* f    = compactshellpairdata_phi::getFac(comASP,iSP);
					const Double* f_0  = sp.getFac();
					const Double* c2_0 = sp.getC2();
					for(UInt ip2=0; ip2<np2; ip2++) {
						Double newFac = f_0[ip2]*c2_0[ip2];
						if (fabs(f[ip2]-newFac)>1.0E-12) { 
							printf("failed the fac array is not same\n");
							compactshellpairdata_phi::printSPData(code,iSP,c2,comASP);
							sp.print(3);
							return false;
						}
					}
				}
			}
		}
	}

	// now let's see the return status
	return true;
}

*/

#endif

void AtomShellPair::print(UInt iprint) const
{
	cout << endl;
	cout << "*******************************************" << endl;
	cout << "*      AtomShellPair Object Print         *" << endl;
	cout << "*******************************************" << endl;
	cout << "total shell pairs(including insignificant ones): " << getNShellPairs() << endl;
	cout << "row atom is                        : " << rowAtomShellIndex << endl;
	cout << "col atom is                        : " << colAtomShellIndex << endl;
	cout << "row atom normal basis offset       : " << iGBasOffset << endl;
	cout << "col atom normal basis offset       : " << jGBasOffset << endl;
	cout << "row atom Cartesian basis offset    : " << iGCarBasOffset << endl;
	cout << "col atom Cartesian basis offset    : " << jGCarBasOffset << endl;
	cout << "number of row atom norm. basis set : " << iNBas << endl;
	cout << "number of col atom norm. basis set : " << jNBas << endl;
	cout << "number of row atom Cart. basis set : " << iNCarBas << endl;
	cout << "number of col atom Cart. basis set : " << jNCarBas << endl;
	for (UInt iSP=0; iSP<getNShellPairs(); iSP++) {
		const ShellPair& sp = getShellPair(iSP);
		if (sp.isSig()) {
			sp.print(iprint);
		}
	}
	cout << endl;
}
