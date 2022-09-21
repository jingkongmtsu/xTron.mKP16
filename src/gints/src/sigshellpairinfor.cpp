/**
 * CPP files corresponding to the sigshellpairinfor.h
 * \author fenglai liu 
 */
#include<cmath>
#include<map>
#include<cstdio>
#include<iostream>
#include<algorithm>
#include "digestutil.h"
#include "excep.h"
#include "blas.h"
#include "blas1.h"
#include "localmemscr.h"
#include "shell.h"
#include "shellprop.h"
#include "integraljobs.h"
#include "shellsize.h"
#include "denmtrxinfor.h"
#include "hgp_os_ints.h"
#include "gintsinfor.h"
#include "shellpair.h"
#ifdef GINTS_PHI_DEBUG
#include "compactsigmolshellpairinfor.h"
#endif
#include "sigshellpairinfor.h"
using namespace digestutil;
using namespace blas;
using namespace excep;
using namespace localmemscr;
using namespace shell;
using namespace shellprop;
using namespace shellsize;
using namespace gintsinfor;
using namespace shellpair;
using namespace denmtrxinfor;
using namespace integraljobs;
#ifdef GINTS_PHI_DEBUG
using namespace compactsigmolshellpairinfor;
#endif
using namespace sigshellpairinfor;

SigAtomShellPairInfor::SigAtomShellPairInfor(const AtomShell& rs, const AtomShell& cs, 
		const AtomShellSize& rsize, const AtomShellSize& csize,
		const Double& thresh, bool sameMolShell):rowAtomShellIndex(rs.getAtomShellIndex()),
	colAtomShellIndex(cs.getAtomShellIndex()),nSigSP(0),intBoundary(ONE),lenC2(0),lenNP2(0), 
	nAngCodeInfor(0)
{
	// reserve the number of possible shell pairs
	UInt nSP = rs.getNShell()*cs.getNShell();
	sigSPList.reserve(nSP);   

	// prepare the input vectors to evaluate significant shell pairs
	UInt iMaxNP = rs.getMaxNP();
	UInt jMaxNP = cs.getMaxNP();
	UInt iNL    = rs.getMaxCompDegree();
	UInt jNL    = cs.getMaxCompDegree();
	UInt maxNP2 = iMaxNP*jMaxNP;
	DoubleVec iexp2(maxNP2);
	DoubleVec fbra(maxNP2);
	DoubleVec P(3*maxNP2);
	DoubleVec icoe2(maxNP2*iNL*jNL);

	// prepare the basis set scaling vectors as well as 
	// OV result for significant shell pairs
	// be careful for the SP shell case
	// therefore if the maxL is P, we made it into SP
	// in case that P shell is actually SP shell 
	UInt iMaxL  = rs.getMaxL(); 
	UInt jMaxL  = cs.getMaxL(); 
	UInt iNCBas = getCartBas(iMaxL,iMaxL);
	if (iMaxL==1) iNCBas = getCartBas(0,iMaxL);
	UInt jNCBas = getCartBas(jMaxL,jMaxL);
	if (jMaxL==1) jNCBas = getCartBas(0,jMaxL);
	UInt maxNBas= iNCBas>jNCBas ? iNCBas:jNCBas;
	DoubleVec ov(iNCBas*jNCBas);
	DoubleVec rowNormVec(maxNBas);
	DoubleVec colNormVec(maxNBas);

	// let's pre-calculate the distance between two atoms
	const Double* A = rs.getXYZ();
	const Double* B = cs.getXYZ();
	Double distance = rs.getDistance(B);

	// now evaluate the significant shell pairs
	// for each shell pair data inside the atom shell pair,
	// we sort the angular momentum according to the order
	// defined in the constant array
	for(UInt iOrder=0; iOrder<MAX_SHELL_PAIR_NUMBER; iOrder++) {

		// now get the corresponding shell code
		UInt angCode   = SHELL_PAIR_ORDER_ARRAY[iOrder];
		UInt oldOffset = nSigSP;

		// now let's loop over the shell pair data
		for(UInt jShell=0;jShell<cs.getNShell();jShell++) {
			const Shell& jS = cs.getShell(jShell);
			for(UInt iShell=0;iShell<rs.getNShell();iShell++) {

				// check that whether the angular momentum match the case
				const Shell& iS = rs.getShell(iShell);
				UInt iLmin = iS.getLmin();
				UInt iLmax = iS.getLmax();
				UInt jLmin = jS.getLmin();
				UInt jLmax = jS.getLmax();
				UInt code  = codeL(iLmin,iLmax,jLmin,jLmax);
				if (code != angCode) continue;

				//
				// firstly, let's deal with a special case:
				// the two input atom shell are on the same center
				// this happens no matter whether the rs and cs are same molecule shell
				//
				// for this case we can not use the normal overlap integrals
				// to evaluate the significance of shell pair. Due to the 
				// symmetry of shells, the result overlap may be zero (for 
				// example, the (P|S) overlap on the same atom is zero). However,
				// such shell pair is definitely significant. Therefore, we 
				// need to treat it as a special case.
				//
				// On the other hand, the shell pairs on the same atom will be 
				// always significant in terms of the overlap integral. This 
				// is because the AB2 is always 0, so the prefactor in the overlap
				// is always 1 (see the val of pref in formPrimPairData). The 
				// prefactor is the dominant factor to determine the significance
				// of the shell pair
				//
				//
				Double distanceThresh = 0.000001E0;
				if (distance<distanceThresh) {

					// for the same atom shell, we only take the lower tri-angular part
					// for the shell pairs
					if (sameMolShell & (iShell < jShell)) continue;

					// consider all of shell pairs are significant
					SigShellPairInfor sigShellPair(iS,jS);
					sigSPList.push_back(sigShellPair);
					nSigSP = nSigSP + 1;

					// update the lenC2 etc.
					UInt np2 = iS.getNPrim()*jS.getNPrim();
					lenNP2  += np2;
					if (iLmin != iLmax || jLmin != jLmax) {
						lenC2 += (iLmax-iLmin+1)*(jLmax-jLmin+1)*np2;
					}
					continue;
				}

				// now let's determine that whether two shells are out of shell radius?
				if (rsize.outofShellRadius(iShell,jShell,csize,distance)) continue;

				// do we need to switch the shell in pair?
				bool needSwitch = switchShell(iLmin,iLmax,jLmin,jLmax);

				//
				// now let's do the screening test
				// 
				// The sparsity between two shells are evaluated simply by the overlap
				// integral. We will form the overlap integrals in terms of the Cartesian
				// type of basis sets (even if the basis set data is in pure form). The ov
				// will be normalized in terms of the Cartesian type of basis set.
				//
				// For the pure form of overlap integral, it's clear to see that they are
				// roughly in the same magnitude with the Cartesian type of overlap 
				// integrals (see the function of getC2PFromLxLyLz in the purecart.cpp 
				// in shell folder, this function records the Cartesian to Pure transformation
				// coefficients). However, for L >6, the difference is getting bigger; 
				// so use it with caution for high L case.
				//
				// Since the overlap between two same atoms has been handled above, therefore
				// here we only deal with the shell pairs whose centers are different. From 
				// the overlap integral direct calculation formula (see oneintformula.cpp),
				// it's clear that some given lmn combination  of the integral could be 
				// zero (that is because the PAi and PBi is zero and i=x,y,z, see the function
				// of Fk in the oneintformula.cpp). However; since the two centers are 
				// different then there must be one dimension (i=x,y,z) of PAi and PBi is 
				// not zero, accordingly there's at least one overlap integrals in the 
				// <x|y> is not zero. That is the reason why we can use the overlap 
				// integrals to judge the significance: if the overlap integral between 
				// the two given shell pair is zero for all of integral elements, then
				// the two shells form an insignificant shell pair.  
				//
				// On the other hand, since in overlap integrals there's a requirement that
				// L(is) >= L(js), therefore the input shell is required to be in correct
				// order, that's the reason we judge whether to switch shell
				//
				// Here we do not consider the redundancy of primitive pairs. The redundancy
				// will be considered when we form the real shell pair data
				//

				// form the prim pair data 
				if (needSwitch) {
					formPrimPairData(jS,iS,B,A,iexp2,icoe2,fbra,P);
				}else{
					formPrimPairData(iS,jS,A,B,iexp2,icoe2,fbra,P);
				}

				// clear the result vector
				UInt len = iS.getNCarBas()*jS.getNCarBas();
				ov.assign(len,ZERO);

				// now let's do the work
				LInt LCode= static_cast<LInt>(code);
				UInt np2  = iS.getNPrim()*jS.getNPrim();
				if (needSwitch) {
					hgp_os_twobodyoverlap(LCode,np2,&icoe2.front(),&iexp2.front(),
							&fbra.front(),&P.front(),B,A,&ov.front());
				}else{
					hgp_os_twobodyoverlap(LCode,np2,&icoe2.front(),&iexp2.front(),
							&fbra.front(),&P.front(),A,B,&ov.front());
				}


				// post work - to norm integrals 
				if (needSwitch) {
					normOV(jS,iS,rowNormVec,colNormVec,ov);
				}else{
					normOV(iS,jS,rowNormVec,colNormVec,ov);
				}

				//
				// debug input
				// see ov
				//
				/*
					UInt iBas = iS.getBasisIndex(0,TYPE_NORM);
					UInt jBas = jS.getBasisIndex(0,TYPE_NORM);
					if (needSwitch) {
					UInt tmp = iBas;
					iBas = jBas;
					jBas = tmp;
					}
					if (iBas ==  && jBas == ) {
					cout << "result in sigshellpair infor " << endl;
					UInt nRow = iS.getNCarBas();
					UInt nCol = jS.getNCarBas();
					if (needSwitch) {
					UInt tmp = nRow;
					nRow = nCol;
					nCol = tmp;
					}
					for(UInt i=0; i<nRow; i++) {
					for(UInt j=0; j<nCol; j++) {
					printf("%-16.13f  ", ov[i+j*nRow]);
					}
					printf("\n");
					}
					printf("\n");
					}
					*/

				// finally, let's find the maximum integral value
				Double maxIntegral = maxSearch(&ov.front(),len);
				if (maxIntegral<thresh) continue;

				// update the significance data
				SigShellPairInfor sigShellPair(iS,jS);
				sigSPList.push_back(sigShellPair);
				nSigSP = nSigSP + 1;

				// update the lenC2 etc.
				lenNP2 += np2;
				if (iLmin != iLmax || jLmin != jLmax) {
					lenC2 += (iLmax-iLmin+1)*(jLmax-jLmin+1)*np2;
				}
			}
		}

		// now let's see the infor
		// no shell pair data for this angular momentum
		if (oldOffset == nSigSP) continue;

		// now set the data
		UInt index                 = nAngCodeInfor;
		angCodeInfor[index]        = angCode;
		angCodeSPOffset[2*index  ] = oldOffset;
		angCodeSPOffset[2*index+1] = nSigSP - oldOffset;
		nAngCodeInfor = nAngCodeInfor + 1;
	}
}

void SigAtomShellPairInfor::formPrimPairData(const Shell& is, const Shell& js,
		const Double* A, const Double* B, DoubleVec& iexp2,
		DoubleVec& icoe2, DoubleVec& fbra, DoubleVec& P) const
{
	//
	// this function is used to form the input data
	// for calling integral routines
	// before filling the data, the input vector 
	// must be assigned with enough memory
	//
	// make sure that L(is) >= L(js)
	//

	UInt inp = is.getNPrim();
	UInt jnp = js.getNPrim();
	const Double* iexp = is.getExp();
	const Double* jexp = js.getExp();
	Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
	UInt count = 0;
	for(UInt jp=0; jp<jnp; jp++) {
		for(UInt ip=0; ip<inp; ip++) {

			// prefactors etc.
			Double ia   = iexp[ip];
			Double ja   = jexp[jp];
			Double alpla= ia+ja; 
			Double ab   = -ia*ja/alpla;
			Double pref = exp(ab*AB2)*pow(PI/alpla,1.5E0);
			iexp2[count]= ONE/alpla;
			fbra[count] = pref;

			// form P point according to the 
			// Gaussian pritimive product theorem
			Double adab = ia/alpla; 
			Double bdab = ja/alpla; 
			Double Px   = A[0]*adab + B[0]*bdab;
			Double Py   = A[1]*adab + B[1]*bdab;
			Double Pz   = A[2]*adab + B[2]*bdab;
			P[3*count+0]= Px;
			P[3*count+1]= Py;
			P[3*count+2]= Pz;
			count++;
		}
	}

	//
	// form coefficients data
	// here we need to be careful
	// since the two shells may be switched 
	// due to fixed L pair in integral calculation
	//
	count = 0;
	UInt iLmax  = is.getLmax();
	UInt iLmin  = is.getLmin();
	UInt jLmax  = js.getLmax();
	UInt jLmin  = js.getLmin();
	UInt nLBra1 = iLmax-iLmin+1;
	UInt nLBra2 = jLmax-jLmin+1;
	UInt nBra1P = is.getNPrim();
	UInt nBra2P = js.getNPrim();
	const Double* icoe = is.getCoe();
	const Double* jcoe = js.getCoe();
	for(UInt j=0; j<nLBra2; j++) {
		for(UInt i=0; i<nLBra1; i++) {
			for(UInt jp=0; jp<nBra2P; jp++) {
				for(UInt ip=0; ip<nBra1P; ip++) {
					Double ic = icoe[ip+i*nBra1P];
					Double jc = jcoe[jp+j*nBra2P];
					icoe2[count] = ic*jc;
					count++;
				}
			}
		}
	}

	//
	// debug print out
	// active it if you have any question 
	// for the overlap input parameters calculation
	//
	/*
	UInt np  = 0;
	UInt np2 = inp*jnp;
	for(UInt J=jLmin; J<=jLmax; J++){
		for(UInt I=iLmin; I<=iLmax; I++){
			printf("Angular mom: %d, %d\n", (Int)I ,(Int)J);
			for(UInt i=0;i<np2;i++) {
				printf("%-12s %-12s %-15s %-12s %-12s %-12s\n", "Coef1*Coef2", 
						"ONE/ExpSum", "Overlap", "Px", "Py", "Pz");
				printf("%-12.6f %-12.6f %-15.9f %-12.6f %-12.6f %-12.6f\n", 
						icoe2[np+i], iexp2[i], fbra[i], P[3*i], P[3*i+1], P[3*i+2]);
			}
			np += np2;
		}
	}
	*/
}

void SigAtomShellPairInfor::normOV(const Shell& is, const Shell& js, 
		DoubleVec& rowScaleNormVec, DoubleVec& colScaleNormVec, 
		DoubleVec& ov) const
{
	//
	// make sure that L(is) >= L(js)
	//
	UInt iLmax  = is.getLmax();
	UInt jLmax  = js.getLmax();
	if (iLmax>=2 || jLmax>=2) {

		// do we switch shells in integral calculation?
		UInt row = is.getNCarBas();
		UInt col = js.getNCarBas();

		// firstly do the row shell scale vector
		if (iLmax>=2) {
			getNormScaleVector(iLmax,&rowScaleNormVec.front()); 
		}

		// secondly do the column shell scale vector
		if (jLmax>=2) {
			getNormScaleVector(jLmax,&colScaleNormVec.front()); 
		}

		// do the work
		for(UInt iCol=0; iCol<col; iCol++) {

			// scale column basis set
			if (jLmax>=2) {
				Double N = colScaleNormVec[iCol];
				vscal(&ov[iCol*row],N,row);
			}

			// scale row basis set
			if (iLmax>=2) {
				vmul(&ov[iCol*row],&rowScaleNormVec.front(),&ov[iCol*row],row);
			}
		}
	}
}

void SigAtomShellPairInfor::copySigData(const SigAtomShellPairInfor& infor, const Double& thresh) 
{
	rowAtomShellIndex = infor.rowAtomShellIndex; 
	colAtomShellIndex = infor.colAtomShellIndex;          
	nSigSP            = 0;
	nAngCodeInfor     = 0;
	intBoundary       = infor.intBoundary;
	lenNP2            = infor.lenNP2;
	lenC2             = infor.lenC2;
	sigSPList.clear();
	for(UInt iAng=0; iAng<infor.getNAngCode(); iAng++) {
		UInt angCode   = infor.getAngCode(iAng);
		UInt oldOffset = nSigSP;
		UInt beginIndex,nSP;
		infor.getSPOffsetData(iAng,beginIndex,nSP);
		for(UInt iSP=beginIndex; iSP<beginIndex+nSP; iSP++) {
			const SigShellPairInfor& sp = infor.getSigShellPairInfor(iSP);
			Double boundary = sp.getIntBoundaryVal();
			if (boundary<thresh) continue;
			sigSPList.push_back(sp);
			nSigSP = nSigSP + 1;
		}
		if (nSigSP>oldOffset) {
			UInt index                 = nAngCodeInfor;
			angCodeInfor[index]        = angCode;
			angCodeSPOffset[2*index  ] = oldOffset;
			angCodeSPOffset[2*index+1] = nSigSP - oldOffset;
			nAngCodeInfor = nAngCodeInfor + 1;
		}
	}
}

#ifdef GINTS_PHI_DEBUG
void SigAtomShellPairInfor::init(const CompactSigAtomShellPairInfor& aspInfor, 
					const CompactSigShellPairInfor* spInfor)
{
	// set data
	rowAtomShellIndex = compactsigmolshellpairinfor::getRowAtomShellIndex(aspInfor); 
	colAtomShellIndex = compactsigmolshellpairinfor::getColAtomShellIndex(aspInfor);          
	nSigSP            = compactsigmolshellpairinfor::getNSigSP(aspInfor);
	intBoundary       = compactsigmolshellpairinfor::getAtomSPIntBoundaryVal(aspInfor);
	nAngCodeInfor     = compactsigmolshellpairinfor::getNAngCode(aspInfor);
	lenC2             = 0;
	lenNP2            = 0;
	for(UInt iSP=0; iSP<nSigSP; iSP++) {
		UInt rowShlIndex   = compactsigmolshellpairinfor::getRowShellIndex(spInfor[iSP]);
		UInt colShlIndex   = compactsigmolshellpairinfor::getColShellIndex(spInfor[iSP]);
		Double intBoundary = compactsigmolshellpairinfor::getSPIntBoundaryVal(spInfor[iSP]);
		SigShellPairInfor newSP(rowShlIndex, colShlIndex, intBoundary);
		sigSPList[iSP] = newSP;
	}
	for(UInt iAng=0; iAng<nAngCodeInfor; iAng++) {
		UInt angCode = compactsigmolshellpairinfor::getAngCode(iAng,aspInfor);
		UInt offset  = compactsigmolshellpairinfor::getSPOffsetData(iAng,aspInfor);
		UInt nSP     = compactsigmolshellpairinfor::getNSPForAngCode(iAng,aspInfor);
		angCodeInfor[iAng]        = angCode;
		angCodeSPOffset[2*iAng  ] = offset;
		angCodeSPOffset[2*iAng+1] = nSP;
	}
}
#endif

bool SigAtomShellPairInfor::operator==(const SigAtomShellPairInfor& asp) const 
{
	// here we omit the two data, lenC2 and lenNP2
	// they are in the sigatomshellpairinfor because we need to know the maximum value
	// for all of atom shell pairs, it's only the maximum value be used 
	// however, to derive the maximum value of lenC2 and lenNP2 we need to keep track of 
	// the two values for all of atom shell pairs. So we have lenC2 and lenNP2 at SigAtomShellPairInfor,
	// but they are not available in the compactsigshellpairinfor class.
	//
	// therefore if the input asp is converted from the CompactSigShellPairInfor,
	// we will not have the information over there

	// test the foundamental data
	if (rowAtomShellIndex != asp.rowAtomShellIndex ||
			colAtomShellIndex != asp.colAtomShellIndex ||         
			nSigSP != asp.nSigSP) return false;           

	// test integral boundary
	if (fabs(intBoundary-asp.intBoundary)>THRESHOLD_MATH) return false;

	// test that whehter the angular momentum code information are same
	if (nAngCodeInfor != asp.nAngCodeInfor) return false;
	for(UInt i=0; i<nAngCodeInfor; i++) {
		if (angCodeInfor[i] != asp.angCodeInfor[i]) return false;
		if (angCodeSPOffset[2*i  ] != asp.angCodeSPOffset[2*i  ]) return false;
		if (angCodeSPOffset[2*i+1] != asp.angCodeSPOffset[2*i+1]) return false;
	}

	// test each shell pair
	for(UInt i=0; i<nSigSP; i++) {
		const SigShellPairInfor& sigSP0 = asp.getSigShellPairInfor(i);
		const SigShellPairInfor& sigSP1 = getSigShellPairInfor(i);
		UInt index1_0 = sigSP0.getRowShellIndex();
		UInt index2_0 = sigSP0.getColShellIndex();
		Double bound_0= sigSP0.getIntBoundaryVal();
		UInt index1_1 = sigSP1.getRowShellIndex();
		UInt index2_1 = sigSP1.getColShellIndex();
		Double bound_1= sigSP1.getIntBoundaryVal();
		if (index1_0 != index1_1 || index2_0 != index2_1) return false;
		if (fabs(bound_0-bound_1)>THRESHOLD_MATH) return false;
	}

	// now two are same
	return true;
}

void SigAtomShellPairInfor::print() const {
	printf("%s  %-7d\n", "row atom shell index:", (Int)rowAtomShellIndex);
	printf("%s  %-7d\n", "col atom shell index:", (Int)colAtomShellIndex);
	printf("%s  %-7d\n", "Num. of Sig. ShlPair:", (Int)getNSigShellPairs());
	printf("%s  %-7d\n", "num. of angular code:", (Int)nAngCodeInfor);
	printf("%s  %-7d\n", "num. of prim.  pairs:", (Int)lenNP2);
	printf("%s  %-7d\n", "num. of coeff. pairs:", (Int)lenC2);
	printf("%s  %-15.12f\n", "maximum integral boundary value:", intBoundary);
	printf("%-15s   %-15s   %-15s   %-15s\n", "angCode", "row shell index", "col shell index", "int. boundary");
	for(UInt iAng=0; iAng<getNAngCode(); iAng++) {
		UInt angCode   = getAngCode(iAng);
		UInt beginIndex,nSP;
		getSPOffsetData(iAng,beginIndex,nSP);
		for(UInt iSP=beginIndex; iSP<beginIndex+nSP; iSP++) {
			const SigShellPairInfor& sigSP = getSigShellPairInfor(iSP);
			UInt index1 = sigSP.getRowShellIndex();
			UInt index2 = sigSP.getColShellIndex();
			Double bound= sigSP.getIntBoundaryVal();
			printf("%-15d   %-15d   %-15d   %-15.12f\n", (Int)angCode, (Int)index1, (Int)index2, bound);
		}
	}
}

SigMolShellPairInfor::SigMolShellPairInfor(const MolShell& rs, const MolShell& cs, 
		const GIntJobInfor& ginfor):rowShellCode(rs.getCode()),colShellCode(cs.getCode()),
	maxL(0),maxNShellPairs(0),maxNP2(0),maxLContraction(1),maxASPNP2(0),maxASPLenC2(0),
	maxNBasShell(0),maxNBasAtom(0)
{
	//
	// we can form some data only with shell
	//
	maxNP2 = rs.getMaxNP()*cs.getMaxNP();
	UInt maxRowLContraction = rs.getMaxLContraction();
	UInt maxColLContraction = cs.getMaxLContraction();
	maxLContraction = maxRowLContraction>maxColLContraction ? maxRowLContraction : maxColLContraction;

	// form maxL
	UInt maxL1 = rs.getMaxL();
	UInt maxL2 = cs.getMaxL();
	maxL = maxL1>maxL2 ? maxL1 : maxL2;

	// form max number basis set for shell
	// if maxL is 1, we need to be careful since we may
	// have SP shell
	maxNBasShell = getCartBas(maxL,maxL); 
	if (maxL == 1) maxNBasShell = getCartBas(0,maxL);

	// now let's form the nbas based on atom shell
	UInt maxRowNBasAtom = rs.getMaxBasis(TYPE_CART);
	UInt maxColNBasAtom = cs.getMaxBasis(TYPE_CART);
	maxNBasAtom = maxRowNBasAtom>maxColNBasAtom ? maxRowNBasAtom : maxColNBasAtom;

	// let's build the significant shell pair information data
	// for this version it does not use the density matrix information
	// therefore we directly build the shell pairs data
	formSigAtomSPList(rs,cs,ginfor,sigAtomSPList);
}

SigMolShellPairInfor::SigMolShellPairInfor(const MolShell& rs, const MolShell& cs, 
		const DenMtrxInfor& denMtrxInfor, const GIntJobInfor& ginfor):rowShellCode(rs.getCode()),
	colShellCode(cs.getCode()),maxL(0),maxNShellPairs(0),maxNP2(0),maxLContraction(1),
	maxASPNP2(0),maxASPLenC2(0),maxNBasShell(0),maxNBasAtom(0)
{
	//
	// we can form some data only with shell
	//
	maxNP2 = rs.getMaxNP()*cs.getMaxNP();
	UInt maxRowLContraction = rs.getMaxLContraction();
	UInt maxColLContraction = cs.getMaxLContraction();
	maxLContraction = maxRowLContraction>maxColLContraction ? maxRowLContraction : maxColLContraction;

	// form maxL
	UInt maxL1 = rs.getMaxL();
	UInt maxL2 = cs.getMaxL();
	maxL = maxL1>maxL2 ? maxL1 : maxL2;

	// form max number basis set for shell
	// if maxL is 1, we need to be careful since we may
	// have SP shell
	maxNBasShell = getCartBas(maxL,maxL); 
	if (maxL == 1) maxNBasShell = getCartBas(0,maxL);

	// now let's form the nbas based on atom shell
	UInt maxRowNBasAtom = rs.getMaxBasis(TYPE_CART);
	UInt maxColNBasAtom = cs.getMaxBasis(TYPE_CART);
	maxNBasAtom = maxRowNBasAtom>maxColNBasAtom ? maxRowNBasAtom : maxColNBasAtom;

	// for this version of shell pair information, it's used for integral jobs
	// with digesting density matrix data. we will check it here
	UInt jobName = ginfor.getIntJob();
	if (! digestWithDensityMatrix(jobName)) {
		string job = getJobName(jobName);
		printf("this integral job does not need density matrix infor for forming sig shell pair%s\n",
				job.c_str());
		string infor = "the given job should use another constructor without density matrix information";
		Excep excep("SigMolShellPairInfor","constructor",EXCEPTION_GINTS_SHELL_PAIR_ERROR,infor);
		handleExcep(excep);
	}
	if (! doERI(jobName)) {
		string job = getJobName(jobName);
		printf("this integral job can not use Cauchy-Schwarz inequality to set up integral boundary%s\n",
				job.c_str());
		string infor = "it's only ERI type of job can use this constructor to set up sig shell pair infor";
		Excep excep("SigMolShellPairInfor","constructor",EXCEPTION_GINTS_SHELL_PAIR_ERROR,infor);
		handleExcep(excep);
	}

	// let's build the significant shell pair information data
	// for this version, it will work together with the density 
	// matrix information.
	//
	// We will form the si atom shell pair list first,
	// then set up the integral boundary with possible density matrix infor;
	// next to sort the shell pairs according to the integral boundary
	//
	// currently because we only apply the cauchy-schwarz inequality 
	// for the ERI type of integrals, so we inpose the if condition
	// as below

	// form the tmp list
	SigAtomShellPairInforVec list;
	formSigAtomSPList(rs,cs,ginfor,list);

	// further more do the integral boundary check
	setIntBoundary(ginfor,rs,cs,denMtrxInfor,list);

	// squeeze the shell pair list to form the final one
	squeezeSPList(list,ginfor);

	// let's sort the sig atom shell pair information
	// according to the result integral boundary values
	std::stable_sort(sigAtomSPList.begin(),sigAtomSPList.end());

	// set up the batch information for future use
	setupBatchSigAtomShellPairInfor(ginfor); 
}

void SigMolShellPairInfor::init(const MolShell& rs, const MolShell& cs, 
		const GIntJobInfor& ginfor)
{
	// initialize the data
	rowShellCode    = rs.getCode();
	colShellCode    = cs.getCode();
	maxL            = 0;
	maxNShellPairs  = 0;
	maxNP2          = 0;
	maxLContraction = 1;
	maxNBasShell    = 0;
	maxNBasAtom     = 0;
	maxASPNP2       = 0;
	maxASPLenC2     = 0;

	//
	// we can form some data only with shell
	//
	maxNP2 = rs.getMaxNP()*cs.getMaxNP();
	UInt maxRowLContraction = rs.getMaxLContraction();
	UInt maxColLContraction = cs.getMaxLContraction();
	maxLContraction = maxRowLContraction>maxColLContraction ? maxRowLContraction : maxColLContraction;

	// form maxL
	UInt maxL1 = rs.getMaxL();
	UInt maxL2 = cs.getMaxL();
	maxL = maxL1>maxL2 ? maxL1 : maxL2;

	// form max number basis set for shell
	// if maxL is 1, we need to be careful since we may
	// have SP shell
	maxNBasShell = getCartBas(maxL,maxL); 
	if (maxL == 1) maxNBasShell = getCartBas(0,maxL);

	// now let's form the nbas based on atom shell
	UInt maxRowNBasAtom = rs.getMaxBasis(TYPE_CART);
	UInt maxColNBasAtom = cs.getMaxBasis(TYPE_CART);
	maxNBasAtom = maxRowNBasAtom>maxColNBasAtom ? maxRowNBasAtom : maxColNBasAtom;

	// let's build the significant shell pair information data
	// for this version it does not use the density matrix information
	// therefore we directly build the shell pairs data
	formSigAtomSPList(rs,cs,ginfor,sigAtomSPList);
}

void SigMolShellPairInfor::init(const MolShell& rs, const MolShell& cs, 
		const DenMtrxInfor& denMtrxInfor, const GIntJobInfor& ginfor)
{
	// initialize the data
	rowShellCode    = rs.getCode();
	colShellCode    = cs.getCode();
	maxL            = 0;
	maxNShellPairs  = 0;
	maxNP2          = 0;
	maxLContraction = 1;
	maxASPNP2       = 0;
	maxASPLenC2     = 0;
	maxNBasShell    = 0;
	maxNBasAtom     = 0;

	//
	// we can form some data only with shell
	//
	maxNP2 = rs.getMaxNP()*cs.getMaxNP();
	UInt maxRowLContraction = rs.getMaxLContraction();
	UInt maxColLContraction = cs.getMaxLContraction();
	maxLContraction = maxRowLContraction>maxColLContraction ? maxRowLContraction : maxColLContraction;

	// form maxL
	UInt maxL1 = rs.getMaxL();
	UInt maxL2 = cs.getMaxL();
	maxL = maxL1>maxL2 ? maxL1 : maxL2;

	// form max number basis set for shell
	// if maxL is 1, we need to be careful since we may
	// have SP shell
	maxNBasShell = getCartBas(maxL,maxL); 
	if (maxL == 1) maxNBasShell = getCartBas(0,maxL);

	// now let's form the nbas based on atom shell
	UInt maxRowNBasAtom = rs.getMaxBasis(TYPE_CART);
	UInt maxColNBasAtom = cs.getMaxBasis(TYPE_CART);
	maxNBasAtom = maxRowNBasAtom>maxColNBasAtom ? maxRowNBasAtom : maxColNBasAtom;

	// for this version of shell pair information, it's used for integral jobs
	// with digesting density matrix data. we will check it here
	UInt jobName = ginfor.getIntJob();
	if (! digestWithDensityMatrix(jobName)) {
		string job = getJobName(jobName);
		printf("this integral job does not need density matrix infor for forming sig shell pair%s\n",
				job.c_str());
		string infor = "the given job should use another constructor without density matrix information";
		Excep excep("SigMolShellPairInfor","constructor",EXCEPTION_GINTS_SHELL_PAIR_ERROR,infor);
		handleExcep(excep);
	}
	if (! doERI(jobName)) {
		string job = getJobName(jobName);
		printf("this integral job can not use Cauchy-Schwarz inequality to set up integral boundary%s\n",
				job.c_str());
		string infor = "it's only ERI type of job can use this constructor to set up sig shell pair infor";
		Excep excep("SigMolShellPairInfor","constructor",EXCEPTION_GINTS_SHELL_PAIR_ERROR,infor);
		handleExcep(excep);
	}

	// let's build the significant shell pair information data
	// for this version, it will work together with the density 
	// matrix information.
	//
	// We will form the si atom shell pair list first,
	// then set up the integral boundary with possible density matrix infor;
	// next to sort the shell pairs according to the integral boundary
	//
	// currently because we only apply the cauchy-schwarz inequality 
	// for the ERI type of integrals, so we inpose the if condition
	// as below

	// form the tmp list
	SigAtomShellPairInforVec list;
	formSigAtomSPList(rs,cs,ginfor,list);

	// further more do the integral boundary check
	setIntBoundary(ginfor,rs,cs,denMtrxInfor,list);

	// squeeze the shell pair list to form the final one
	squeezeSPList(list,ginfor);

	// let's sort the sig atom shell pair information
	// according to the result integral boundary values
	std::stable_sort(sigAtomSPList.begin(),sigAtomSPList.end());

	// set up the batch information for future use
	setupBatchSigAtomShellPairInfor(ginfor); 
}

void SigMolShellPairInfor::formSigAtomSPList(const MolShell& rs, const MolShell& cs, 
		const GIntJobInfor& infor, SigAtomShellPairInforVec& list)
{
	// get the threshold for shell pair
	Double threshold = infor.pickUpSPThresh();

	// reset the maximum number of shell pairs
	maxNShellPairs = 0;
	maxASPNP2      = 0;
	maxASPLenC2    = 0;

	// form the shell size data
	MolShellSize rsize(rs,threshold);
	MolShellSize csize(cs,threshold);

	// are the two shells same?
	bool sameMolShell = false;
	if (rs == cs) sameMolShell = true;

	// reserve the space for forming atom shell pairs
	UInt nAtomShellPairs = 0;
	UInt nRowAtomShell   = 0;
	UInt nColAtomShell   = 0;
	if(sameMolShell) {
		UInt nAtomShell = rs.getNAtomShells();
		nRowAtomShell   = nAtomShell;
		nColAtomShell   = nAtomShell;
		nAtomShellPairs = nAtomShell*(1+nAtomShell)/2;
	}else{
		nRowAtomShell   = rs.getNAtomShells();
		nColAtomShell   = cs.getNAtomShells();
		nAtomShellPairs = nRowAtomShell*nColAtomShell;
	}
	list.reserve(nAtomShellPairs);

	// now let's loop over atoms in the shell
	for(UInt colIndex=0; colIndex<nColAtomShell; colIndex++) {
		for(UInt rowIndex=0; rowIndex<nRowAtomShell; rowIndex++) {

			// now get the atom shell data
			const AtomShell& rowAtomShell = rs.getAtomShell(rowIndex);
			const AtomShell& colAtomShell = cs.getAtomShell(colIndex);

			// if the two molshell are same, 
			// we only pick up data in the 
			// lower-triangular plus diagonal part
			UInt rowAtomIndex = rowAtomShell.getAtomShellIndex();
			UInt colAtomIndex = colAtomShell.getAtomShellIndex();
			if (sameMolShell && rowAtomIndex<colAtomIndex) continue;

			// whether the two atom shell out of atom size?
			UInt atom1 = rowAtomShell.getAtomic();
			const AtomShellSize& rowAtomShellSize = rsize.getAtomShellSizeInfor(atom1);
			UInt atom2 = colAtomShell.getAtomic();
			const AtomShellSize& colAtomShellSize = csize.getAtomShellSizeInfor(atom2);
			Double distance = rowAtomShell.getDistance(colAtomShell.getXYZ());
			if(rowAtomShellSize.outofAtomRadius(colAtomShellSize,distance)) continue;

			// now let's try to form the sig atom shell pair
			SigAtomShellPairInfor sigPair(rowAtomShell,colAtomShell,rowAtomShellSize,
					colAtomShellSize,threshold,sameMolShell);	
			if (sigPair.isInsig()) continue;
			list.push_back(sigPair);

			// see the maxNShellPair for this sp
			if (sigPair.getNSigShellPairs() >maxNShellPairs) {
				maxNShellPairs = sigPair.getNSigShellPairs();
			}

			// also update the lenC2 etc.
			if (sigPair.getLenC2()>maxASPLenC2) maxASPLenC2 = sigPair.getLenC2();
			if (sigPair.getLenNP2()>maxASPNP2 ) maxASPNP2   = sigPair.getLenNP2();
		}
	}
}

void SigMolShellPairInfor::squeezeSPList(const SigAtomShellPairInforVec& list, const GIntJobInfor& ginfor) 
{
	// now get the length of list
	UInt nSigSP = 0;
	Double thresh = ginfor.pickUpThresh();
	for(UInt i=0; i<list.size(); i++) {
		const SigAtomShellPairInfor& infor = list[i];
		Double boundary = infor.getIntBoundaryVal();
		if (infor.isInsig() || boundary<thresh) continue;
		nSigSP++;
	}
	sigAtomSPList.reserve(nSigSP);

	// now push data in
	SigAtomShellPairInfor newSPInfor;
	for(UInt i=0; i<list.size(); i++) {
		const SigAtomShellPairInfor& infor = list[i];
		Double boundary = infor.getIntBoundaryVal();
		if (infor.isInsig() || boundary<thresh) continue;
		newSPInfor.copySigData(infor,thresh);
		sigAtomSPList.push_back(newSPInfor);
	}
}

void SigMolShellPairInfor::setIntBoundary(const GIntJobInfor& ginfor, const MolShell& bra1, 
		const MolShell& bra2, const DenMtrxInfor& denInfor, 
		SigAtomShellPairInforVec& sigAtomShellPairList) const
{
	// now it's time to set up the TBB object 
	TBB_CSIntBoundary cs(ginfor,denInfor,bra1,bra2,*this,sigAtomShellPairList);

	// get the global information
	const GlobalInfor& globInfor = ginfor.getGlobalInfor();

	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	if (globInfor.useMultiThreads()) {
		init.initialize(globInfor.getNCPUThreads());
	}else{
		init.initialize(1);
	}
	
	// possible timing code
	//tick_count t0 = tick_count::now();

	// this is the parallel body
	UInt len = sigAtomShellPairList.size();
	parallel_for(blocked_range<UInt>(0,len), cs);

	// possible timing code
	//tick_count t1 = tick_count::now();
	//Double t = (t1-t0).seconds();
	//if (globInfor.printTiming()) {
	//printf("%s  %-12.6f\n", "cauchy-schwarz(CS) integral boundary with TBB threads, time in seconds ", t);
	//}
}

void SigMolShellPairInfor::setupBatchSigAtomShellPairInfor(const GIntJobInfor& ginfor) 
{
	// create the number of levels, the first batch is from 0.1 to infinity
	// then [0.1,  +infinity)  - level 0
	// then [0.01, 0.1)        - level 1
	// then [0.001,0.01)       - level 2
	// etc.
	// set all of boundary in default -1
	// we may have incontinous level,
	// for example; [0.01, 0.1) may be followd by
	// interval of [0.00001, 0.0001) 
	Double thresh = ginfor.pickUpThresh();
	UInt nLevels = static_cast<UInt>(ceil(-log10(thresh)));
	UIntVec batchArray(2*nLevels,-1);
	DoubleVec limitVals(nLevels,ZERO);

	// now let's count the batches, this is half-closed interval
	// we take [lower, upper)
	for (UInt iSP=0; iSP<getNSigAtomShellPairs(); iSP++) {
		const SigAtomShellPairInfor& sp = getSigAtomShellPairInfor(iSP);
		Double val = sp.getIntBoundaryVal();

		// determine which batch it belongs
		// we note, that if val>=0.1 it's the first batch
		UInt iBatch  = 0;
		if (val<0.1E0) {
			if (val<thresh) {
				UInt order = static_cast<UInt>(ceil(-log10(thresh)));
				iBatch = order-1;
			}else{
				UInt order = static_cast<UInt>(ceil(-log10(val)));
				iBatch = order-1;
			}
		}

		// check the iBatch
		// it should not exceed the range of batchArray
		if (iBatch>=nLevels) {
				string infor = "the batch index calculation is wrong, seems the boundary is too small";
				Excep excep("SigMolShellPairInfor","setupBatchSigAtomShellPairInfor",
						EXCEPTION_SIG_SHELL_PAIR_INT_BOUNDARY_ERROR,infor);
				handleExcep(excep);
		}

		// for this batch whether the lower boundary is set?
		UInt lowerBoundIndex = batchArray[2*iBatch+0];
		if (lowerBoundIndex == static_cast<UInt>(-1)) {
			batchArray[2*iBatch+0] = iSP;

			// for the new batch we also update the lower limit value
			Double power = 0.1E0;
			if (iBatch>0) {
				power = pow(0.1E0,iBatch+1);
			}
			limitVals[iBatch] = power;
		}

		// also set the upper bound
		// this is half interval, so upper is not reached
		batchArray[2*iBatch+1] = iSP+1;
	}

	// now let's copy the batch information
	batchSigAtomShellPairInfor.reserve(2*nLevels);
	batchThreshLimits.reserve(nLevels);    
	for(UInt iBatch=0; iBatch<batchArray.size()/2; iBatch++) {
		UInt lowerBoundIndex = batchArray[2*iBatch+0];
		UInt upperBoundIndex = batchArray[2*iBatch+1];
		if (lowerBoundIndex == static_cast<UInt>(-1) || 
				upperBoundIndex == static_cast<UInt>(-1)) {
			continue;
		}
		batchSigAtomShellPairInfor.push_back(lowerBoundIndex);
		batchSigAtomShellPairInfor.push_back(upperBoundIndex);
		Double boundaryVal = limitVals[iBatch];
		batchThreshLimits.push_back(boundaryVal);
	}
}

void SigMolShellPairInfor::print(UInt level) const 
{
	cout << endl;
	cout << "*******************************************" << endl;
	cout << "*    SigMolShellPairInfor Object Print    *" << endl;
	cout << "*******************************************" << endl;
	cout << "max number of shell pairs: " << getMaxNSP() << endl;
	cout << "max L contraction        : " << getMaxNL() << endl;
	cout << "max number of prim. pairs: " << getMaxNP2() << endl;
	cout << "max prim. pairs for asp  : " << getMaxAtomSPNP2() << endl;
	cout << "max len of c2 for asp    : " << getMaxAtomSPLenC2() << endl;

	// batch information
	cout << endl;
	if (level>=1) {
		if (batchSigAtomShellPairInfor.size() > 0) {
			cout << "print out the batch significant atom shell pair information " << endl;
			for(UInt iBatch=0; iBatch<getNBatch(); iBatch++) {
				UInt lowerBoundIndex = getStartingSigAtomSPIndex(iBatch); 
				UInt totalNAtomSP    = getNSigAtomSPInBatch(iBatch);
				UInt upperBoundIndex = lowerBoundIndex+totalNAtomSP-1;
				Double lower = getBatchLowerLimitBoundary(iBatch);
				Double upper = lower*10.0E0;
				if (iBatch == 0) {
					printf("batch ranges from [%-E to +infinity), index ranges from %-12d to %-12d, totally %-12d\n", 
							lower, (Int)lowerBoundIndex, (Int)upperBoundIndex, (Int)(totalNAtomSP));
				}else{
					printf("batch ranges from [%-E to %-E), index ranges from %-12d to %-12d, totally %-12d\n", 
							lower, upper, (Int)lowerBoundIndex, (Int)upperBoundIndex, (Int)(totalNAtomSP));
				}
			}
		}
	}

	// print out the detailed information
	if (level>=2) {
		cout << endl;
		for (UInt iSP=0; iSP<getNSigAtomShellPairs(); iSP++) {
			const SigAtomShellPairInfor& sp = getSigAtomShellPairInfor(iSP);
			cout << "significant atom shell pair infor index " << iSP << endl;
			sp.print();
		}
	}
}

void SigMolShellPairInfor::statPrint(const Double& thresh, const MolShell& rs, const MolShell& cs) const
{
	cout << endl;
	cout << "*******************************************" << endl;
	cout << "* Statistic Report for Integral Boundary  *" << endl;
	cout << "*******************************************" << endl;

	//
	// set up number of levels
	// we do not check the number here, must be < 1 and > 0
	//
	UInt nLevels = static_cast<UInt>(-log10(thresh));

	//
	// now let's set up the vector of map
	// so to record for each threshold, the L pair arrangement
	//
	map<UInt,UInt> emptyTemplate;
	vector<map<UInt,UInt> >angInfor(nLevels,emptyTemplate);

	//
	// for example, if thresh is 1.0^-8,
	// then we will see the integral result ranging from 0.1 to 1.0^-8
	// so totally 8 levels
	// also the statistic report is based on both atom shell pair and 
	// shell pair counting
	//
	UIntVec atomShellPair(nLevels,0);
	UIntVec shellPair(nLevels,0);

	//
	// now we loop over the information to get statistics
	//
	for (UInt iSP=0; iSP<getNSigAtomShellPairs(); iSP++) {
		const SigAtomShellPairInfor& sp = getSigAtomShellPairInfor(iSP);
		UInt rowAtomShellIndex = sp.getRowAtomShellIndex();
		UInt colAtomShellIndex = sp.getColAtomShellIndex();
		const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
		const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);

		// atom shell pair level
		Double val = sp.getIntBoundaryVal();
		if (val>0.1E0) {
			atomShellPair[0] += 1;
		}else if (val<thresh) {
			atomShellPair[nLevels-1] += 1;
		}else {
			UInt pos = (UInt)floor(-log10(val));
			atomShellPair[pos] += 1;
		}

		// now we loop over all of shell pairs for this atom pair
		for (UInt iP2=0; iP2<sp.getNSigShellPairs(); iP2++) {
			const SigShellPairInfor& p2 = sp.getSigShellPairInfor(iP2);

			// this is for simple statistic
			UInt pos = 0;
			Double val= p2.getIntBoundaryVal();
			if (val>0.1E0) {
				shellPair[0] += 1;
			}else if (val<thresh) {
				pos = nLevels-1;
				shellPair[pos] += 1;
			}else {
				pos = (UInt)floor(-log10(val));
				shellPair[pos] += 1;
			}

			// now check the angular momentum information
			map<UInt,UInt>& LInfor = angInfor[pos];
			const Shell& iS = rowAtomShell.getShell(p2.getRowShellIndex());
			const Shell& jS = colAtomShell.getShell(p2.getColShellIndex());
			UInt iLmin = iS.getLmin();
			UInt iLmax = iS.getLmax();
			UInt jLmin = jS.getLmin();
			UInt jLmax = jS.getLmax();
			UInt LCode = codeL(iLmin,iLmax,jLmin,jLmax);
			if (LInfor.find(LCode) == LInfor.end()) {
				LInfor.insert(std::pair<UInt,UInt>(LCode,1));
			}else{
				LInfor.at(LCode) += 1; 
			}
		}
	}

	//
	// now print out all of information
	//
	for(UInt i=0; i<2; i++) {
		cout << endl << endl;
		cout << "***********************************" << endl;
		if (i==0) {
			cout << "statistic information for atom pair" << endl;
		}else{
			cout << "statistic information for shell pair" << endl;
		}
		cout << "***********************************" << endl;
		for(UInt iLevel=0; iLevel<nLevels; iLevel++) {

			// get number
			UInt n = atomShellPair[iLevel];
			if (i==1) n = shellPair[iLevel];

			// now set up the print statement
			string infor = "for range between 1(>1) to 0.1: ";
			if (iLevel==1) {
				infor = "for range between 1.0E-1 to 1.0E-2: ";
			}else if (iLevel==2) {
				infor = "for range between 1.0E-2 to 1.0E-3: ";
			}else if (iLevel==3) {
				infor = "for range between 1.0E-3 to 1.0E-4: ";
			}else if (iLevel==4) {
				infor = "for range between 1.0E-4 to 1.0E-5: ";
			}else if (iLevel==5) {
				infor = "for range between 1.0E-5 to 1.0E-6: ";
			}else if (iLevel==6) {
				infor = "for range between 1.0E-6 to 1.0E-7: ";
			}else if (iLevel==7) {
				infor = "for range between 1.0E-7 to 1.0E-8: ";
			}else if (iLevel==8) {
				infor = "for range between 1.0E-8 to 1.0E-9: ";
			}else if (iLevel==9) {
				infor = "for range between 1.0E-9 to 1.0E-10: ";
			}else if (iLevel==10){
				infor = "for range between 1.0E-10 to 1.0E-11: ";
			}else if (iLevel==11){
				infor = "for range between 1.0E-11 to 1.0E-12: ";
			}else if (iLevel==12){
				infor = "for range between 1.0E-12 to 1.0E-13: ";
			}else if (iLevel==13){
				infor = "for range between 1.0E-13 to 1.0E-14: ";
			}else if (iLevel>=14){
				infor = "for range between 1.0E-14 to lower: ";
			}

			// now print out
			printf("%-45s  %d\n", infor.c_str(), (Int)n);

			// for this level, print out the angular momentum infor
			if (i == 1) {
				cout << "***********************************" << endl;
				cout << "statistic information for L pairs  " << endl;
				cout << "***********************************" << endl;
				const map<UInt,UInt>& LInfor = angInfor[iLevel];
				for (map<UInt,UInt>::const_iterator it=LInfor.begin(); it!=LInfor.end(); ++it) {
					UInt key = it->first;
					UInt val = it->second;

					// decode L code
					UInt iLmin,iLmax,jLmin,jLmax;
					decodeL(key,iLmin,iLmax,jLmin,jLmax);
					printf("%-4d  %-4d  %-4d  %-4d:  %-d\n", (Int)iLmin,(Int)iLmax,(Int)jLmin,
							(Int)jLmax,(Int)val);
				}
				cout << endl;
			}
		}
	}
}

UInt SigMolShellPairInfor::getSigSPVecLen(const MolShell& rs, const MolShell& cs, const UInt& type) const
{
	UInt len = 0;
	for (UInt iSP=0; iSP<getNSigAtomShellPairs(); iSP++) {
		const SigAtomShellPairInfor& sp = getSigAtomShellPairInfor(iSP);
		UInt rowAtomShellIndex = sp.getRowAtomShellIndex();
		UInt colAtomShellIndex = sp.getColAtomShellIndex();
		const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
		const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);
		UInt basPair = rowAtomShell.getNBas(type)*colAtomShell.getNBas(type);
		len += basPair;
	}
	return len;
}

UInt SigMolShellPairInfor::getS2VecOffset(const MolShell& rs, const MolShell& cs,
		const UInt& rowShellIndex, const UInt& colShellIndex, const UInt& type) const
{
	// let's observe the offset
	UInt len = 0;
	for (UInt iSP=0; iSP<getNSigAtomShellPairs(); iSP++) {
		const SigAtomShellPairInfor& sp = getSigAtomShellPairInfor(iSP);
		UInt rowAtomShellIndex = sp.getRowAtomShellIndex();
		UInt colAtomShellIndex = sp.getColAtomShellIndex();
		const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
		const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);
		if (rowAtomShellIndex == rowShellIndex && colShellIndex == colAtomShellIndex) {
			return len;
		}else{
			UInt basPair = rowAtomShell.getNBas(type)*colAtomShell.getNBas(type);
			len += basPair;
		}
	}

	// if we reach here, something wrong happens
	// so issue an error
	string infor = "for the input row/col atom shell index, we did not get the corresponding offset";
	Excep excep("SigMolShellPairInfor","getS2VecOffset",EXCEPTION_GINTS_SHELL_PAIR_ERROR,infor);
	handleExcep(excep);
	return -1;
}

void SigMolShellPairInfor::convertSigSPVec(bool toVec, const GlobalInfor& globInfor, 
		const MolShell& rs, const MolShell& cs, const UInt& transWork, const UInt& scaleWork, 
		const UInt& matrixStatus, Mtrx& M, DoubleVec& v) const
{
	// let's double check if the transform C2P/P2C work is needed
	if (transWork != NO_TRANS) {

		// set the source data and target data basis set type
		// the code is copied from CPTransAtomShell constructor
		UInt sourceType = TYPE_NORM;
		UInt targetType = TYPE_NORM;
		if (transWork == C2P_WITH_L00 || transWork == C2P_WITH_XYZ) {
			sourceType = TYPE_CART;
		}else if (transWork == P2C_WITH_L00 || transWork == P2C_WITH_XYZ) {
			targetType = TYPE_CART;
		}
		if (matrixStatus == WITH_MATRIX_TRANSPOSE) {
			UInt tmp = sourceType;
			sourceType = targetType;
			targetType = tmp;
		}

		// set the matrix type and vec type
		UInt matrixType = sourceType;
		UInt vecType    = targetType; 
		if (! toVec) {
			matrixType = targetType;
			vecType    = sourceType; 
		}

		// we need to check on the matrix dimension
		// if it does not fit then we can not do the job
		if (M.getRow() != rs.getNBas(matrixType) || M.getCol() != cs.getNBas(matrixType)) {
			string infor = "input matrix dimension does not equal to the given input basis function dimension";
			Excep excep("SigMolShellPairInfor","convertSigSPVec",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}

		// check the s2 vector type
		UInt len = getSigSPVecLen(rs,cs,vecType);
		if (len != v.size()) {
			string infor = "input s2 vector dimension is conflict to the dimension determined by the basis set type";
			Excep excep("SigMolShellPairInfor","convertSigSPVec",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}
	}

	// now let's determine the basType
	// it's used when the transWork is NO_TRANS
	// in this case both of the s2 vector and matrix should have same type
	UInt basType = -1;
	if (transWork == NO_TRANS) {

		// get the final type of the target matrix
		// so that set the offset
		UInt nRow = M.getRow();
		UInt t1 = 0;
		if(nRow == rs.getNCarBas()) t1=1;
		if(nRow == rs.getNBas())    t1=2;
		if(t1 == 0) {
			string info = "the input matrix row dimension conflict with row shell data";
			Excep excep("SigMolShellPairInfor","convertSigSPVec",EXCEPTION_BLOCK_MATRIX_CPTRANS_ERROR,info);
			handleExcep(excep);
		}	

		// col
		UInt nCol = M.getCol();
		UInt t2 = 0;
		if(nCol == cs.getNCarBas()) t2=1;
		if(nCol == cs.getNBas())    t2=2;
		if(t2 == 0) {
			string info = "the input matrix col dimension conflict with col shell data";
			Excep excep("SigMolShellPairInfor","convertSigSPVec",EXCEPTION_BLOCK_MATRIX_CPTRANS_ERROR,info);
			handleExcep(excep);
		}	

		// row shell type is not same with col shell type
		if(t1 != t2) {
			string info = "the input matrix row and col dimension have different shell type";
			Excep excep("SigMolShellPairInfor","convertSigSPVec",EXCEPTION_BLOCK_MATRIX_CPTRANS_ERROR,info);
			handleExcep(excep);
		}	

		// now set the type information
		UInt type0 = TYPE_NORM;
		if(t2 == 1) {
			type0 = TYPE_CART;
		}

		// check the s2 vector type
		UInt len = getSigSPVecLen(rs,cs,type0);
		if (len != v.size()) {
			string infor = "input s2 vector dimension is conflict to the dimension determined by the basis set type";
			Excep excep("SigMolShellPairInfor","convertSigSPVec",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}

		// finally let's set the type
		basType = type0;
	}

	// now ready for the work
	// now it's time to set up the TBB object 
	TBB_MtrxToS2Vec s2vec(toVec,transWork,scaleWork,matrixStatus,basType,rs,cs,*this,M,v);

	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	if (globInfor.useMultiThreads()) {
		init.initialize(globInfor.getNCPUThreads());
	}else{
		init.initialize(1);
	}
	
	// possible timing code
	//tick_count t0 = tick_count::now();

	// this is the parallel body
	UInt len = sigAtomSPList.size();
	parallel_for(blocked_range<UInt>(0,len), s2vec);

	// possible timing code
	//tick_count t1 = tick_count::now();
	//Double t = (t1-t0).seconds();
	//printf("%s  %-12.6f\n", "convert s2 form vectors with TBB threads, time in seconds ", t);
}

#ifdef GINTS_PHI_DEBUG

bool SigMolShellPairInfor::testCompactSigMolSPInfor(const CompactSigMolShellPairInfor& spInfor) const 
{
	// calcualte the number of sig shell pairs
	UInt nSigSP = 0;
	for (UInt iSP=0; iSP<getNSigAtomShellPairs(); iSP++) {
		const SigAtomShellPairInfor& sp = getSigAtomShellPairInfor(iSP);
		nSigSP += sp.getNSigShellPairs();
	}

	// firstly test all of non-array information
	if (compactsigmolshellpairinfor::getNSigShellPairs(spInfor) != nSigSP) {
		cout << "number of total shell pairs are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getNSigAtomShellPairs(spInfor) != getNSigAtomShellPairs()) {
		cout << "number of atom shell pairs are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxNP2(spInfor) != getMaxNP2()) {
		cout << "maxNP2 are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxAtomSPNP2(spInfor) != getMaxAtomSPNP2()) {
		cout << "maxASPNP2 are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxAtomSPLenC2(spInfor) != getMaxAtomSPLenC2()) {
		cout << "maxASPLenC2 are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxNL(spInfor) != getMaxNL()) {
		cout << "maxNL are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxNSP(spInfor) != getMaxNSP()) {
		cout << "maxNSP are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxL(spInfor) != getMaxL()) {
		cout << "maxL are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxNBasisForShell(spInfor) != getMaxNBasisForShell()) {
		cout << "maxNBasisForShell are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getMaxNBasisForAtomShell(spInfor) != getMaxNBasisForAtomShell()) {
		cout << "maxNBasisForAtomShell are not same" << endl;
		return false;
	}
	if (compactsigmolshellpairinfor::getNBatch(spInfor) != getNBatch()) {
		cout << "nBatches are not same" << endl;
		return false;
	}

	// now let's test the atom shell pair infor array and shell pair infor array
	// we will re-generated the atom shell pair infor data
	SigAtomShellPairInfor newASP(compactsigmolshellpairinfor::getMaxNSP(spInfor));
	for(UInt i=0; i<compactsigmolshellpairinfor::getNSigAtomShellPairs(spInfor); i++) {
		const CompactSigAtomShellPairInfor& aspInfor = compactsigmolshellpairinfor::getCompactSigAtomShellPairInfor(i,spInfor);
		UInt beginIndex = compactsigmolshellpairinfor::getBeginIndexSigSPIndex(aspInfor);
		const CompactSigShellPairInfor* spArray = compactsigmolshellpairinfor::getCompactSigShellPairInforArray(beginIndex,
				spInfor);
		newASP.init(aspInfor,spArray);
		const SigAtomShellPairInfor& oldSP = getSigAtomShellPairInfor(i);
		if (newASP == oldSP) {
			continue;
		}else{
			cout << "the " << i << " atom shell pair infor data are not same" << endl;
			newASP.print();
			oldSP.print();
			return false;
		}
	}

	// finally if we have batch information let's test it here
	for(UInt i=0; i<compactsigmolshellpairinfor::getNBatch(spInfor); i++) {
		UInt beginIndex0 = getStartingSigAtomSPIndex(i);
		UInt nSP0        = getNSigAtomSPInBatch(i);
		Double limit0    = getBatchLowerLimitBoundary(i);
		UInt beginIndex  = compactsigmolshellpairinfor::getStartingSigAtomSPIndex(spInfor,i);
		UInt nSP         = compactsigmolshellpairinfor::getNSigAtomSPInBatch(spInfor,i);
		Double limit     = compactsigmolshellpairinfor::getBatchLowerLimitBoundary(spInfor,i);
		if (nSP0 != nSP) {
			cout << "batch testing nSP does not match" << endl;
			return false;
		}
		if (beginIndex0 != beginIndex) {
			cout << "batch testing beginIndex does not match" << endl;
			return false;
		}
		if (fabs(limit0-limit)>1.0E-15) {
			cout << "batch testing limit does not match" << endl;
			return false;
		}
	}

	// finally, everything is good
	// let's return
	return true;
}

#endif


TBB_CSIntBoundary::TBB_CSIntBoundary(const GIntJobInfor& ginfor, const DenMtrxInfor& denInfor0,
		const MolShell& b1, const MolShell& b2, const SigMolShellPairInfor& infor0, 
		SigAtomShellPairInforVec& list):thresh(ginfor.pickUpSPThresh()),
	intThresh(ginfor.pickUpThresh()),jobInfor(ginfor),denInfor(denInfor0),bra1(b1),bra2(b2),
	infor(infor0),spList(list),scale(infor0.getMaxL())
{
	// we only scale the integrals
	scale.initScaleData();

	// actually for S, P, SP scale that's all constant one
	scale_S[0]  = ONE;
	scale_P[0]  = ONE;
	scale_P[1]  = ONE;
	scale_P[2]  = ONE;
	scale_SP[0] = ONE; 
	scale_SP[1] = ONE; 
	scale_SP[2] = ONE; 
	scale_SP[3] = ONE; 
}

void TBB_CSIntBoundary::operator()(const blocked_range<UInt>& r) const
{
	//
	// setup scratch data for calculation
	// 1  atom shell pair
	// 2  heap memory management if necessary
	// 3  integral array
	//
	
	// shell pair data
	UInt maxNP2 = infor.getMaxNP2();
	UInt maxNL  = infor.getMaxNL();
	UInt maxNSP = infor.getMaxNSP();
	AtomShellPair sp(maxNSP,maxNP2,maxNL);

	// heap memory tool
	UInt maxL = infor.getMaxL();
	UInt heapMemLen = 0;
	if (jobInfor.useHeapMem(maxL)) {
		heapMemLen = jobInfor.setMemSize(maxL);
	}
	LocalMemScr scr(heapMemLen);

	// integral array
	// we do not count in the redundant ones
	UInt maxNCarBas  = infor.getMaxNBasisForShell(); 
	UInt maxIntLen   = maxNCarBas*maxNCarBas*maxNCarBas*maxNCarBas;
	DoubleVec intArray(maxIntLen);

	// now real working code
	for( UInt n=r.begin(); n!=r.end(); ++n ) {

		//
		//get shell pair infor
		//
		SigAtomShellPairInfor& sigAtomShellPairInfor = spList[n];

		//
		//form atom shell pair
		//
		sp.init(bra1,bra2,sigAtomShellPairInfor,thresh);

		//
		//generate integral boundary
		//
		Double atomShellPairIntBound = ZERO;
		UInt nSP = sp.getNShellPairs(); 
		for(UInt iSP=0; iSP<nSP; iSP++) {

			// whether the real shell data is significant?
			// in this case the integral boundary is set to zero
			const ShellPair& spData = sp.getShellPair(iSP);
			if (! spData.isSig()) {
				sigAtomShellPairInfor.setIntBoundaryVal(iSP,ZERO);
				continue;
			}

			// now let's do ERI for (uv|r12|uv)
			UInt np2          = spData.getNP2();
			UInt Lmax         = spData.getLmax();
			const Double* c2  = spData.getC2();
			const Double* e2  = spData.getE2();
			const Double* A   = spData.getA();
			const Double* B   = spData.getB();
			const Double* P   = spData.getP();
			const Double* fac = spData.getFac();

			// form LCode
			UInt code = spData.getLCode();
			LInt LCode= codeSPCodes(code,code);

			// initilize the integral
			UInt iNCarBas, jNCarBas;
			spData.getNCarBas(iNCarBas,jNCarBas);
			UInt nCarBas2 = iNCarBas*jNCarBas;
			UInt intLen = nCarBas2*nCarBas2;
			vset(&intArray.front(),ZERO,intLen);

			//
			// now do ERI, generate the raw ERI
			// we note that here we do not really check the significance inside
			// since (AB|r12|AB) is always important if shell pair (AB| is important
			// so we set PMAX = 1
			//
			// also the omega value is set to 0, so we check the real ERI
			// even though this may be erf(omega*r12)/r12 or erfc(omega*r12)/r12
			//
			// because their absolute value all less than the normal ERI
			//
			Double Pmax  = ONE;
			Double omega = ZERO;
			hgp_os_eri(LCode,Lmax,np2,np2,Pmax,omega,c2,e2,fac,P,A,B,c2,e2,fac,
					P,A,B,&intArray.front(),scr);
			
			// now get the integral infor and reset scr
			const SingleIntegralInfor& intInfor = jobInfor.getIntegralInfor(LCode);
			if (intInfor.getMemInfor() > 0) {
				scr.reset();
			}

			// scale the raw ERI so to generate the normalized ERI
			UInt l1Min,l1Max,l2Min,l2Max;
			spData.getL(l1Min,l1Max,l2Min,l2Max);
			const Double* bra1Scale = getScaleVec(l1Min,l1Max);
			const Double* bra2Scale = getScaleVec(l2Min,l2Max);
			scale4BodyInts(l1Min,l1Max,l2Min,l2Max,l1Min,l1Max,l2Min,l2Max, 
					bra1Scale,bra2Scale,bra1Scale,bra2Scale,&intArray[0]); 

			// now ERI should be in Cartesian format and it's fully normalized
			// let's see the maximum value for the integrals
			Double maxVal = fabs(maxSearch(&intArray.front(),intLen));
			Double val = sqrt(maxVal);

			// now let's see that whether this shell pair is significant or not
			// also generate the corresponding integral boundary
			if (jobInfor.getJobOrder() == 0) {

				// now this is for energy calculation
				// intBoundary = sqrt(max(uv|r12|uv)) 
				// later for the analytical integrals, we will see that
				// whether we will use 
				// (u,v|r12|lam,eta) <= sqrt(max(uv|r12|uv))*sqrt(max(lam,eta|r12|lam,eta))
				// as the judgement for significance of shell quartet (u,v|r12|lam,eta)
				//
				// here for the job involves only K, we can even screening the 
				// shell pairs in considering with P_(u,max) and P_(v,max)
				// information in terms of (uv|r12|uv). P_(u,max) is the maximum
				// density matrix element for row u.
				//
				// this is because for (u,v|r12|lam,eta), the digestion always
				// takes the P_(u,*) and P_(v,*), therefore it can see that 
				// P_(u,lam)*(u,v|r12|lam,eta) (or P_(u,eta))<=
				// P_(u,max)*sqrt(max(uv|r12|uv))*sqrt(max(lam,eta|r12|lam,eta)) =>
				//
				// P_(u,lam)*(u,v|r12|lam,eta) <= P_(u,max)*sqrt(max(uv|r12|uv))
				//
				Double intBoundaryVal = val;

				// firstly check whether the integral itself is small enough
				if (intBoundaryVal<intThresh) {
					sigAtomShellPairInfor.setIntBoundaryVal(iSP,ZERO);
					continue;
				}

				// now let's see whether it the job only with K
				UInt job = jobInfor.getIntJob(); 
				if (onlyDoK(job)) {

					// get the shell index
					UInt rowShellIndex = -1;
					UInt colShellIndex = -1;
					spData.getGlobalShellIndex(rowShellIndex,colShellIndex);

					// let's get the P Max value for both shell for u 
					// and shell for v
					// here we note that it does not matter that two
					// shells are switched
					Double aPMax1 = denInfor.getPmaxOnRowShell(rowShellIndex,0);
					Double aPMax2 = denInfor.getPmaxOnRowShell(colShellIndex,0);
					Double pMax   = aPMax1 > aPMax2 ? aPMax1 : aPMax2;
					if (denInfor.getNSpin() == 2) {
						Double bPMax1 = denInfor.getPmaxOnRowShell(rowShellIndex,1);
						Double bPMax2 = denInfor.getPmaxOnRowShell(colShellIndex,1);
						Double bPMax  = bPMax1 > bPMax2 ? bPMax1 : bPMax2;
						if (bPMax>pMax) pMax = bPMax;
					}

					// let's see whether we need to get rid of this shell pair?
					Double screening = intBoundaryVal*pMax;
					if (screening<intThresh) {
						sigAtomShellPairInfor.setIntBoundaryVal(iSP,ZERO);
						continue;
					}
				}

				// now this shell pair is sig, set up the integral boundary value
				// we will use this value to judge the significance of shell quartet
				// in real ERI calculation
				sigAtomShellPairInfor.setIntBoundaryVal(iSP,intBoundaryVal);
				if (atomShellPairIntBound<intBoundaryVal) atomShellPairIntBound = intBoundaryVal;
			}else{

				// now it's the derivatives calculation
				UInt job = jobInfor.getIntJob(); 
				Double intBoundaryVal = val;

				// let's see the J
				// the jPmax actually returns the Pmax for u, v shell pair block
				Double jPmax = denInfor.getJPMax(spData,spData);
				Double JBoundaryVal = intBoundaryVal*jPmax;
				if (onlyDoJ(job)) {

					// let's see whether the boundary value is significant enough
					if (JBoundaryVal<intThresh) {
						sigAtomShellPairInfor.setIntBoundaryVal(iSP,ZERO);
						continue;
					}

					// well it's sig enough, set the boundary value
					sigAtomShellPairInfor.setIntBoundaryVal(iSP,JBoundaryVal);
					if (atomShellPairIntBound<JBoundaryVal) atomShellPairIntBound = JBoundaryVal;

					// now let's continue
					continue;
				}

				// now let's see the exchange
				// get the shell index
				UInt rowShellIndex = -1;
				UInt colShellIndex = -1;
				spData.getGlobalShellIndex(rowShellIndex,colShellIndex);

				// let's get the P Max pair value for both shell for u 
				// and shell for v, this is used for K screening
				// this is because that for the derivatives, 
				// it's always P_u,max*P_v,max*integral
				// therefore we can detect it before real calculation
				//
				// however, for the K integral boundary value, because it's 
				// undermined, therefore we can not use any density matrix
				// information
				Double aPMax1 = denInfor.getPmaxOnRowShell(rowShellIndex,0);
				Double aPMax2 = denInfor.getPmaxOnRowShell(colShellIndex,0);
				Double pMaxPair = aPMax1*aPMax2; 
				if (denInfor.getNSpin() == 2) {
					Double bPMax1 = denInfor.getPmaxOnRowShell(rowShellIndex,1);
					Double bPMax2 = denInfor.getPmaxOnRowShell(colShellIndex,1);
					Double bPMaxPair  = bPMax1*bPMax2; 
					if (bPMaxPair>pMaxPair) pMaxPair = bPMaxPair;
				}
				Double KScreening   = intBoundaryVal*pMaxPair;
				Double KBoundaryVal = intBoundaryVal;

				// whether this is only K job?
				if (onlyDoK(job)) {

					// let's see whether the boundary value is significant enough
					if (KScreening<intThresh) {
						sigAtomShellPairInfor.setIntBoundaryVal(iSP,ZERO);
						continue;
					}

					// well it's sig enough, set the boundary value
					sigAtomShellPairInfor.setIntBoundaryVal(iSP,KBoundaryVal);
					if (atomShellPairIntBound<KBoundaryVal) atomShellPairIntBound = KBoundaryVal;

					// now let's continue
					continue;
				}

				// now the J type work is done with together with K type work
				// let's see that whether we can get rid of this shell pair?
				bool isJSig = false;
				if (doJ(job) && JBoundaryVal>intThresh) isJSig = true;
				bool isKSig = false;
				if (doK(job) && KScreening>intThresh) isKSig = true;
				if (!isJSig && !isKSig) {
					sigAtomShellPairInfor.setIntBoundaryVal(iSP,ZERO);
					continue;
				}

				// next set the boundary value
				sigAtomShellPairInfor.setIntBoundaryVal(iSP,intBoundaryVal);
				if (atomShellPairIntBound<intBoundaryVal) atomShellPairIntBound = intBoundaryVal;
			}
		}

		//
		// finally, set integral boundary for atom shell pair
		//
		sigAtomShellPairInfor.setIntBoundaryVal(atomShellPairIntBound);
	}
}

void TBB_MtrxToS2Vec::operator()(const blocked_range<UInt>& r) const 
{
	// the source data matrix 
	Mtrx source(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());

	// the target data matrix 
	Mtrx target(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());

	// scratch matrix used in CP transformation
	Mtrx tmp(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());

	// C2P/P2C transformation and basis set scaling 
	UInt maxL = infor.getMaxL();
	NormAtomShellData scale(maxL,infor.getMaxNBasisForAtomShell(),scaleStatus);
	CPTransAtomShell trans(maxL,transStatus,matrixStatus);

	///////////////////////////////////////////////////////////////////////
	// real working code - loop over atom shell pair data                //
	///////////////////////////////////////////////////////////////////////
	for( UInt n=r.begin(); n!=r.end(); ++n ) {

		//
		//get shell pair infor
		//
		const SigAtomShellPairInfor& sigAtomShellPairInfor = infor.getSigAtomShellPairInfor(n); 

		// get the atom shell data
		UInt rowAtomShellIndex = sigAtomShellPairInfor.getRowAtomShellIndex();
		UInt colAtomShellIndex = sigAtomShellPairInfor.getColAtomShellIndex();
		const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
		const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);

		// we need to determine the s2 vector the basis set data type
		UInt vecBasType = basType;
		if (transStatus != NO_TRANS) {
			if (doMtrxToVec) {
				vecBasType = trans.getTargetType();
			}else{
				vecBasType = trans.getSourceType();
			}
		}

		// get the offset for the vector
		UInt vecOffset = infor.getS2VecOffset(rs,cs,rowAtomShellIndex,colAtomShellIndex,vecBasType);

		// case 1
		// we do not do scaling or normalization work, then we just copy data
		// from source to target
		if (transStatus == NO_TRANS && scaleStatus == NO_SCALE) {

			// set the dimension
			UInt nRow = rowAtomShell.getNBas(basType);
			UInt nCol = colAtomShell.getNBas(basType);
			UInt rowOffset = rowAtomShell.getBasisStartIndex(basType);
			UInt colOffset = colAtomShell.getBasisStartIndex(basType);

			// now let's move the data
			if (doMtrxToVec) {
				M.updateToVector(rowOffset,colOffset,nRow,nCol,ONE,&v[vecOffset]);
			}else{
				M.updateFromVector(rowOffset,colOffset,nRow,nCol,ONE,&v[vecOffset]);
			}

			// now move on to the next one
			continue;
		}

		// let's do C2P or P2C transformation
		// if no such transformation needed, we just copy the data
		// since scaled work must be performed
		if (transStatus != NO_TRANS) {

			// set the source/target matrix
			UInt nSourceRow = rowAtomShell.getNBas(trans.getSourceType());
			UInt nSourceCol = colAtomShell.getNBas(trans.getSourceType());
			UInt nTargetRow = rowAtomShell.getNBas(trans.getTargetType());
			UInt nTargetCol = colAtomShell.getNBas(trans.getTargetType());
			UInt rowOffset  = rowAtomShell.getBasisStartIndex(trans.getSourceType());
			UInt colOffset  = colAtomShell.getBasisStartIndex(trans.getSourceType());

			// initialize the input matrix
			source.init(nSourceRow,nSourceCol);
			target.init(nTargetRow,nTargetCol);
			if (doMtrxToVec) {
				source.copyFromMatrix(rowOffset,colOffset,0,0,nSourceRow,nSourceCol,M);
			}else{
				source.copyFromVector(0,0,nSourceRow,nSourceCol,&v[vecOffset]);
			}

			// we will perform job based on atom shell pair
			trans.cpTransformOnAtomShellPair(rowAtomShell,colAtomShell,source,tmp,target); 
		}else{
			// set the dimension
			UInt nTargetRow = rowAtomShell.getNBas(basType);
			UInt nTargetCol = colAtomShell.getNBas(basType);

			// now directly copy the data
			target.init(nTargetRow,nTargetCol);
			if (doMtrxToVec) {
				UInt rowOffset  = rowAtomShell.getBasisStartIndex(basType);
				UInt colOffset  = colAtomShell.getBasisStartIndex(basType);
				target.copyFromMatrix(rowOffset,colOffset,0,0,nTargetRow,nTargetCol,M);
			}else{
				target.copyFromVector(0,0,nTargetRow,nTargetCol,&v[vecOffset]);
			}
		}

		// let's do scale work finally
		// it's only on target matrix
		if (scaleStatus != NO_SCALE) {
			scale.normRowData(rowAtomShell,target); 
			scale.normColData(colAtomShell,target); 
		}

		// get the final type of the target matrix
		UInt type0 = basType;
		if (transStatus != NO_TRANS) {
			type0 = trans.getTargetType();
		}

		// set the offset
		UInt nRow       = target.getRow();
		UInt nCol       = target.getCol();
		UInt rowOffset  = rowAtomShell.getBasisStartIndex(type0);
		UInt colOffset  = colAtomShell.getBasisStartIndex(type0);

		// let's copy the target matrix to the destination
		// now let's move the data
		if (doMtrxToVec) {
			target.updateToVector(0,0,nRow,nCol,ONE,&v[vecOffset]);
		}else{
			target.updateToMatrix(0,0,rowOffset,colOffset,nRow,nCol,ONE,M);
		}
	}
}

