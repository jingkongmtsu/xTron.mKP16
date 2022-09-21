/**
 * \file    xdm.cpp
 * \brief   class used to calculate VDW by XDM model
 * \author  Fenglai Liu and Jing Kong
 * \note
 *
 * XDM references:
 * 
 * E. Johnson and A. Becke, J. Chem. Phys. 123, 024101 (2005).
 *
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 122, 154104 (2005).
 *
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 123, 154101 (2005).
 * 
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 124, 014104 (2006).
 * 
 * E. R. Johnson and A. D. Becke, J. Chem. Phys. 124, 174104 (2006).
 * 
 * E. R. Johnson and A. D. Becke, Chem. Phys. Lett. 432, 600 (2006).
 * 
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 124108 (2007).
 * 
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).
 * 
 * F. O. Kannemann and A. D. Becke, J. Chem. Theory Comput. 5, 719 (2009).
 * 
 * F. O. Kannemann and A. D. Becke, J. Chem. Theory Comput. 6, 1081 (2010).
 * 
 * A. D. Becke, A. Arabi, and F. O. Kannemann, Can. J. Chem. 88, 1057 (2010).
 * 
 * E. Johnson, J. Chem. Phys. 135, 234109 (2011).
 * 
 * A. Arabi and A. Becke, J. Chem. Phys. 137, 014104 (2012).
 * 
 * A. Otero-de-la Roza and E. R. Johnson, J. Chem. Phys. 138, 054103 (2013).
 *
 */
#include <iostream>
#include "globalinfor.h"
#include "excep.h"
#include "blas.h"
#include "blas1.h"
#include "element.h"
#include "molecule.h"
#include "shell.h"
#include "denmtrx.h"
#include "atomgrids.h"
#include "batchgrid.h"
#include "sigatombasis.h"
#include "dftmatrix.h"
#include "batchbasis.h"
#include "batchvar.h"
#include "batchfunc.h"
#include "atomdenmtrx.h"
#include "hirshfeld.h"
#include "xcenergyinfor.h"
#include "vdwinfor.h"
#include "xdm.h"
using namespace globalinfor;
using namespace excep;
using namespace blas;
using namespace element;
using namespace molecule;
using namespace denmtrx;
using namespace atomgrids;
using namespace batchgrid;
using namespace sigatombasis;
using namespace dftmatrix;
using namespace batchbasis;
using namespace batchvar;
using namespace batchfunc;
using namespace atomdenmtrx;
using namespace xcenergyinfor;
using namespace hirshfeld;
using namespace vdwinfor;
using namespace xdm;
using namespace std;

//
// this piece of code is actually copied from xcints.cpp
//
void TBB_XDM::operator()(const blocked_range<UInt>& r) 
{
	// real working loop
	for( UInt iBatch=r.begin(); iBatch!=r.end(); ++iBatch ) {

		// get the atomic number and the batch grids
		// ibatch is the local batch index
		// iAtom is the atom index to fetch geometry data
		// Z is the corresponding atomic number
		UInt ibatch,Z,iAtom;
		molGrids.getBatchInfor(iBatch,iAtom,Z,ibatch);

		// get the atom grid and shell data
		const AtomShell& as = ms.getAtomShell(iAtom);
		const AtomGrids& ag = molGrids.getAtomGrid(Z);

		// batch grid 
		BatchGrid bg(ag,as,ibatch);

		// sig atom basis
		SigAtomBasis sigList(bg,infor,ms,molShellSize);
		if (sigList.empty()) continue;

		// do partition weights, this is necessary
		// for the xdm part
		bg.BeckeWeights0(infor,ms,sigList);

		// batch basis set
		// for xdm we do not have exrho variable
		BatchBasis bbasis(xcvar,bg,ms,sigList);
		bbasis.setupPhi(bg,ms,sigList);

		// batch var
		BatchVar bVar(sigList,bbasis,xcvar,denMtrx);

		// for the batch functional calculation,
		// we need to make an empty xc energy infor
		XCEnergyInfor xcProfile(infor,xcfunc);

		// batch functional
		BatchFunc bFunc(bVar,denMtrx,xcfunc,xcvar,infor);
		bFunc.doFuncDeriv(bVar,xcfunc,xcvar,xcProfile); 

		// form Hirshfeld weights and free atom vol
		// if you want to do debug printing for hirshfeld weights,
		// please un-comment the lines in the hirshfeld class
		Hirshfeld hirshfeld(infor,sigList,bg);
		hirshfeld.formWtsFreeAtomVol(ms,sigList,bbasis,bg,atomDenMatrices);
		updateFreeAtomVol(sigList,hirshfeld);

		// now finally bring all of things together to form the result
		formSquareMultipleAndAtomVol(ms,sigList,bg,bVar,bFunc,hirshfeld);
	}
}

void TBB_XDM::updateFreeAtomVol(const SigAtomBasis& sigList, const Hirshfeld& hirshfeld) 
{
	const UIntVec& atomIndexArray = sigList.getSigAtoms();
	for(UInt iAtom=0; iAtom<sigList.getNSigAtoms(); iAtom++) {
		UInt atomIndex = atomIndexArray[iAtom];
		freeAtomVol[atomIndex] += hirshfeld.getFreeAtomVol(iAtom);
	}
}

void TBB_XDM::formSquareMultipleAndAtomVol(const MolShell& ms, const SigAtomBasis& sigList, 
		const BatchGrid& bg, const BatchVar& bvar, const BatchFunc& bfunc, 
		const Hirshfeld& hirfeld)
{
	// set the factor - since it's spin resolved result
	// therefore for single spin state we multiply 2
	UInt nSpin  = infor.getNSpin();
	Double fac = ONE;
	if (nSpin == 1) fac = TWO;
	
	// now to calculate the multipole for each significant atom
	UInt nGrids = bg.getNGrids();
	DoubleVec tmp(nGrids);
	const Double* wts = bg.getGridWts();
	const UIntVec& atomIndexArray = sigList.getSigAtoms();
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// get the density and exchange hole values
		UInt var = RA;
		if (iSpin == 1) var = RB;
		UInt pos = xcvar.getVarPos(var);
		const Double* rho = bvar.getVar(pos);
		const Double* dx2 = bfunc.getFunVals() + iSpin*nGrids;

		// now go to each order and each sig atom
		UInt highestOrder = squareMultiple.getCol();
		for(UInt iCol=0; iCol<highestOrder; iCol++) {
			UInt order = iCol+1;
			for(UInt iAtom=0; iAtom<sigList.getNSigAtoms(); iAtom++) {

				// form rdx vector according to the order
				// this is corresponding to equation 10
				UInt atomIndex      = atomIndexArray[iAtom];
				const AtomShell& as = ms.getAtomShell(atomIndex);
				const Double* center= as.getXYZ();
				const Double* coord = bg.getGridCoord();
				for(UInt iG=0; iG<nGrids; iG++) {
					Double x2   = (center[0]-coord[3*iG  ])*(center[0]-coord[3*iG  ]);
					Double y2   = (center[1]-coord[3*iG+1])*(center[1]-coord[3*iG+1]);
					Double z2   = (center[2]-coord[3*iG+2])*(center[2]-coord[3*iG+2]);
					Double r    = sqrt(x2+y2+z2);
					Double dx   = sqrt(dx2[iG]); 
					Double rl   = pow(r,order);
					Double rdxl = pow((r - dx),order);
					tmp[iG]     = (rl-rdxl)*(rl-rdxl);
				}

				// form square multiple
				// this is corresponding to equation 10
				const Double* hirshfeldWts = hirfeld.getWts(iAtom);
				Double contribution = vdot4(hirshfeldWts,wts,rho,&tmp.front(),nGrids);
				squareMultiple(atomIndex,iCol) += fac*contribution;

				/*
				if (order == 3) {
					printf("%-4s  %-14s  %-14s  %-14s  %-14s\n", "Index", "Hir Wts", "Wts", 
							"rhoA", "rdx");
					for(Int iG=0; iG<nGrids; iG++) {
						printf("%-4d  %-14.7f  %-14.7f  %-14.7f  %-14.7f\n", 
								iG+1, hirshfeldWts[iG], wts[iG], rho[iG], tmp[iG]);
					}
					cout << "M2l order " << order << " atom " << atomIndex << " " 
						<< squareMultiple(atomIndex,iCol) << endl; 
				}
				*/
			}
		}
	}

	// now it's for free/effective volumn data
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// get the density 
		UInt var = RA;
		if (iSpin == 1) var = RB;
		UInt pos = xcvar.getVarPos(var);
		const Double* rho = bvar.getVar(pos);

		// loop over sig atoms
		for(UInt iAtom=0; iAtom<sigList.getNSigAtoms(); iAtom++) {

			// form r3 vector 
			UInt atomIndex      = atomIndexArray[iAtom];
			const AtomShell& as = ms.getAtomShell(atomIndex);
			const Double* center= as.getXYZ();
			const Double* coord = bg.getGridCoord();
			for(UInt iG=0; iG<nGrids; iG++) {
				Double x2   = (center[0]-coord[3*iG  ])*(center[0]-coord[3*iG  ]);
				Double y2   = (center[1]-coord[3*iG+1])*(center[1]-coord[3*iG+1]);
				Double z2   = (center[2]-coord[3*iG+2])*(center[2]-coord[3*iG+2]);
				Double r    = sqrt(x2+y2+z2);
				tmp[iG]     = r*r*r;
			}

			// \int r^{3} \rho w_{i} dr - effective volumn
			// this is equation 17
			const Double* hirshfeldWts = hirfeld.getWts(iAtom);
			Double contribution = vdot4(hirshfeldWts,wts,rho,&tmp.front(),nGrids);
			effectiveAtomVol[atomIndex] += fac*contribution;
		}
	}
}

XDM::XDM(const VDWInfor& vdwInfor,const MolShell& ms,const XCIntJobInfor& xcinfor):infor(xcinfor),
	xcfunc(vdwInfor.getEXHoleFuncName()),xcvar(xcfunc,xcinfor.getNSpin(),GROUND_STATE_DFT,0),
	size(ms,xcinfor.getThresh()),molGrids(ms,size,xcinfor),E6(ZERO),E8(ZERO),E10(ZERO),
	nonAdditiveE3(ZERO),E(ZERO)
{

}

void XDM::makeCijForIJAtomPair(const UInt& iAtom, const UInt& jAtom, const Mtrx& squareMultiple,
		const DoubleVec& effectiveAtomPol, Double& C6ij, Double& C8ij, Double& C10ij) const
{
	// calculate the C6ij, C8ij and C10ij
	// order 0: R^6 terms
	// the definition is from equation 13 
	// we note all of the Cij share the same denominator
	UInt order = 0;
	Double iAlpha = effectiveAtomPol[iAtom];
	Double jAlpha = effectiveAtomPol[jAtom];
	Double iM1    = squareMultiple.val(iAtom,order);
	Double jM1    = squareMultiple.val(jAtom,order);

	// evaluating the denominator
	Double denominator = iM1*jAlpha+jM1*iAlpha;
	Double thresh = infor.getThresh();
	if (denominator<thresh) {
		C6ij  = ZERO;
		C8ij  = ZERO;
		C10ij = ZERO;
		cout << "denominator value: " << denominator << endl;
		string infor = "for atom pairs " + boost::lexical_cast<string>(iAtom) + " " + boost::lexical_cast<string>(jAtom) 
			+ " the denominator is less than threshold, we just bypass this atom pair";
		Excep excep("XDM","makeCijForIJAtomPair",EXCEPTION_DIVIDE_BY_SMALL_NUMBER_BYPASS,infor);
		handleExcep(excep);
		return;
	}
	denominator = ONE/denominator;

	// now finish order 0 
	C6ij = iAlpha*jAlpha*iM1*jM1*denominator;

	// order 1: R^8 terms
	// the definition is from equation 14 
	order = 1;
	Double iM2    = squareMultiple.val(iAtom,order);
	Double jM2    = squareMultiple.val(jAtom,order);
	C8ij   = 1.5E0*iAlpha*jAlpha*(iM1*jM2+iM2*jM1)*denominator;

	// order 2: R^10 terms
	// the definition is from equation 15 
	order = 2;
	Double iM3    = squareMultiple.val(iAtom,order);
	Double jM3    = squareMultiple.val(jAtom,order);
	C10ij  = (10.0E0/5.0E0)*iAlpha*jAlpha*(iM1*jM3+iM3*jM1)*denominator;
	C10ij += (21.0E0/5.0E0)*iAlpha*jAlpha*iM2*jM2*denominator;
}

Double XDM::makeBeckeDampingFac(const VDWInfor& xdmInfor, const Double& C6ij, 
		const Double& C8ij, const Double& C10ij) const
{
	// obtain the fitting parameters for XDM model
	Double a1 = xdmInfor.getXDMPar1();
	Double a2 = xdmInfor.getXDMPar2();

	// now it's the damping factor
	Double RCij   = pow(C8ij/C6ij,0.5E0)+pow(C10ij/C6ij,0.25E0)+pow(C10ij/C8ij,0.5E0); 
	RCij = RCij/3.0E0;
	Double Rvdwij = a1*RCij+a2;
	return Rvdwij;
}

void XDM::xdmEnergy(const MolShell& ms, const VDWInfor& xdmInfor, 
		const Mtrx& squareMultiple, const DoubleVec& effectiveAtomVol, const DoubleVec& freeAtomVol) 
{
	// determine the effective pol data
	UInt nAtoms = ms.getNAtomShells();
	DoubleVec effectiveAtomPol(nAtoms);
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomShell& as = ms.getAtomShell(iAtom);
		Double freeAtomPol = getAtomPolarizability(as.getAtomic());
		effectiveAtomPol[iAtom] = effectiveAtomVol[iAtom]*freeAtomPol/freeAtomVol[iAtom];
	}

	// calculate the E6, E8 and E10
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomShell& iAtomShell = ms.getAtomShell(iAtom);
		for(UInt jAtom=iAtom+1; jAtom<nAtoms; jAtom++) {
			const AtomShell& jAtomShell = ms.getAtomShell(jAtom);

			// get Cij
			Double C6ij  = ZERO;
			Double C8ij  = ZERO;
			Double C10ij = ZERO;
			makeCijForIJAtomPair(iAtom,jAtom,squareMultiple,effectiveAtomPol,C6ij,C8ij,C10ij);

			// now it's the damping factor Rvdw
			Double Rvdwij   = makeBeckeDampingFac(xdmInfor,C6ij,C8ij,C10ij); 
			Double Rvdwij6  = pow(Rvdwij,6.0E0);
			Double Rvdwij8  = Rvdwij6*Rvdwij*Rvdwij;
			Double Rvdwij10 = Rvdwij8*Rvdwij*Rvdwij;

			// get the R
			const Double* ixyz = iAtomShell.getXYZ();
			const Double* jxyz = jAtomShell.getXYZ();
			Double x2ij  = (ixyz[0]-jxyz[0])*(ixyz[0]-jxyz[0]);
			Double y2ij  = (ixyz[1]-jxyz[1])*(ixyz[1]-jxyz[1]);
			Double z2ij  = (ixyz[2]-jxyz[2])*(ixyz[2]-jxyz[2]);
			Double Rij   = sqrt(x2ij+y2ij+z2ij);
			Double Rij6  = pow(Rij,6.0E0);
			Double Rij8  = Rij6*Rij*Rij;
			Double Rij10 = Rij8*Rij*Rij;

			// now assemble all of pieces into the energy
			E6  += C6ij/(Rij6+Rvdwij6);
			E8  += C8ij/(Rij8+Rvdwij8);
			E10 += C10ij/(Rij10+Rvdwij10);
			//cout << "C6 " << C6ij << " " << E6  << endl;
			//cout << "C8 " << C8ij << " " << E8 << endl;
			//cout << "C10 " << C10ij << " " << E10 << endl;

			//
			// calculate the non-addative 3 body term
			// this implementation is following the paper below:
			// "A density-functional model of the three-body dispersion interaction 
			// based on exchange dipole moment"
			//
			if (xdmInfor.hasNonAddative3BodyTerm()) {
				for(UInt kAtom=jAtom+1; kAtom<nAtoms; kAtom++) {
					const AtomShell& kAtomShell = ms.getAtomShell(kAtom);

					// get the eta term: equation 14
					UInt order = 0;
					Double iAlpha = effectiveAtomPol[iAtom];
					Double jAlpha = effectiveAtomPol[jAtom];
					Double kAlpha = effectiveAtomPol[kAtom];
					Double iM1    = squareMultiple.val(iAtom,order);
					Double jM1    = squareMultiple.val(jAtom,order);
					Double kM1    = squareMultiple.val(kAtom,order);
					Double etai   = ZERO;
					Double etaj   = ZERO;
					Double etak   = ZERO;
					Double thresh = infor.getThresh();
					if(iAlpha>thresh) etai = (TWO/THREE)*iM1/iAlpha;
					if(jAlpha>thresh) etaj = (TWO/THREE)*jM1/jAlpha;
					if(kAlpha>thresh) etak = (TWO/THREE)*kM1/kAlpha;

					// C9ijk: equation 12
					Double etaij = etai+etaj;
					Double etajk = etaj+etak;
					Double etaik = etai+etak;
					Double C9ijk = (THREE/TWO)*iAlpha*jAlpha*kAlpha*etai*etaj*etak*(etai+etaj+etak)/(etaij*etajk*etaik);

					// distance between atom i and k
					const Double* kxyz = kAtomShell.getXYZ();
					Double x2ik  = (ixyz[0]-kxyz[0])*(ixyz[0]-kxyz[0]);
					Double y2ik  = (ixyz[1]-kxyz[1])*(ixyz[1]-kxyz[1]);
					Double z2ik  = (ixyz[2]-kxyz[2])*(ixyz[2]-kxyz[2]);
					Double Rik   = sqrt(x2ik+y2ik+z2ik);

					// distance between atom j and k
					Double x2jk  = (jxyz[0]-kxyz[0])*(jxyz[0]-kxyz[0]);
					Double y2jk  = (jxyz[1]-kxyz[1])*(jxyz[1]-kxyz[1]);
					Double z2jk  = (jxyz[2]-kxyz[2])*(jxyz[2]-kxyz[2]);
					Double Rjk   = sqrt(x2jk+y2jk+z2jk);

					// cosine value of triangle ABC
					Double R2ij = Rij*Rij;
					Double R2ik = Rik*Rik;
					Double R2jk = Rjk*Rjk;
					Double cosi = -0.5E0*(R2jk - R2ij - R2ik)/(Rij*Rik);
					Double cosj = -0.5E0*(R2ik - R2ij - R2jk)/(Rij*Rjk);
					Double cosk = -0.5E0*(R2ij - R2ik - R2jk)/(Rik*Rjk);

					// how calculate damping factor for ik and jk pair
					Double C6ik   = ZERO;
					Double C8ik   = ZERO;
					Double C10ik  = ZERO;
					makeCijForIJAtomPair(iAtom,kAtom,squareMultiple,effectiveAtomPol,C6ik,C8ik,C10ik);
					Double Rvdwik = makeBeckeDampingFac(xdmInfor,C6ik,C8ik,C10ik); 
					Double C6jk   = ZERO;
					Double C8jk   = ZERO;
					Double C10jk  = ZERO;
					makeCijForIJAtomPair(jAtom,kAtom,squareMultiple,effectiveAtomPol,C6jk,C8jk,C10jk);
					Double Rvdwjk = makeBeckeDampingFac(xdmInfor,C6jk,C8jk,C10jk); 

					// now let's get the final energy term: equation 20
					Double dnorm1 = pow(Rvdwik,3.0E0) + pow(Rik,3.0E0);
					Double dnorm2 = pow(Rvdwjk,3.0E0) + pow(Rjk,3.0E0);
					Double dnorm3 = pow(Rvdwij,3.0E0) + pow(Rij,3.0E0);
					Double dnorm  = dnorm1*dnorm2*dnorm3;
					Double Eijk   = C9ijk*(THREE*cosi*cosj*cosk + ONE)/dnorm;
					nonAdditiveE3 += Eijk;
				}
			}
		}
	}

	// finally write it into the exchange energy
	UInt vdwOrder = xdmInfor.getVDWOrder();
	if (vdwOrder == 1) E = MINUS_ONE*E6;
	if (vdwOrder == 2) E = MINUS_ONE*(E6+E8);
	if (vdwOrder == 3) E = MINUS_ONE*(E6+E8+E10);
	if (xdmInfor.hasNonAddative3BodyTerm()) {
		E += nonAdditiveE3;
	}
}

void XDM::energy(const VDWInfor& xdmInfor, const MolShell& ms, const Molecule& mol, 
		const DenMtrx& denMtrx, bool printTiming)
{
	// firstly, make sure that the input ms data
	// and mol are consistent with each other
	if (ms.getNAtomShells() != mol.getNAtoms()) {
		string infor = "nAtoms in molecule is not same with nAtomShell of molShell";
		Excep excep("XDM","energy",EXCEPTION_RAW_SHELL_DATA_CONFLICT_MOL_DATA,infor);
		handleExcep(excep);
	}

	// calculate the atom density matrix data
	const GlobalInfor& ginfor = infor.getGlobalInfor();
	AtomDenMtrx atomDenMtrx(ginfor,mol);
	atomDenMtrx.formAtomDenMtrxInSCF(mol,ms);

	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	if (ginfor.useMultiThreads()) {
		init.initialize(ginfor.getNCPUThreads());
	}else{
		init.initialize(1);
	}

	// possible timing code
	tick_count t0 = tick_count::now();

	// do normal integral matrix calculation
	TBB_XDM tbb_xdm(ms.getNAtomShells(),ms,molGrids,size,xcfunc,xcvar,infor,denMtrx,atomDenMtrx); 
	UInt len = molGrids.getNTotalBatches();
	parallel_reduce(blocked_range<UInt>(0,len), tbb_xdm);

	// possible timing code
	tick_count t1 = tick_count::now();
	Double t = (t1-t0).seconds();
	if (printTiming) {
		printf("%s  %-12.6f\n", "XDM work with TBB threads, time in seconds  ", t);
	}

	// finally let's do xdm energy
	xdmEnergy(ms,xdmInfor,tbb_xdm.getSquareMultiple(),
			tbb_xdm.getEffectiveAtomVol(),tbb_xdm.getFreeAtomVol()); 
}

