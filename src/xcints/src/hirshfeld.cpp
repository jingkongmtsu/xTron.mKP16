/**
 * \file    hirshfeld.cpp
 * \brief   classes used to describe the Hirshfeld partition scheme
 * \author  Fenglai Liu and Jing Kong
 * \note
 *
 * Hirshfeld references:
 * F. L. Hirshfeld, Theor. Chim. Acta Vol 44, Page 129, 1977
 */
#include "globalinfor.h"
#include "blas.h"
#include "blas1.h"
#include "shell.h"
#include "atomdenmtrx.h"
#include "xcintsinfor.h"
#include "batchgrid.h"
#include "sigatombasis.h"
#include "dftmatrix.h"
#include "batchbasis.h"
#include "batchvar.h"
#include "xcvar.h"
#include "hirshfeld.h"
using namespace globalinfor;
using namespace blas;
using namespace shell;
using namespace atomdenmtrx;
using namespace xcintsinfor;
using namespace batchgrid;
using namespace sigatombasis;
using namespace dftmatrix;
using namespace batchbasis;
using namespace batchvar;
using namespace xcvar;
using namespace hirshfeld;

Hirshfeld::Hirshfeld(const XCIntJobInfor& infor, const SigAtomBasis& sigList, 
		const BatchGrid& bg):thresh(infor.getThresh()),xcfunc("SLATER"),
	weights(bg.getNGrids(),sigList.getNSigAtoms()),
	freeAtomVol(sigList.getNSigAtoms(),ZERO) { }

void Hirshfeld::formWtsFreeAtomVol(const MolShell& ms, 
		const SigAtomBasis& sigList, const BatchBasis& basis, const BatchGrid& bg, 
		const AtomDenMtrx& atomDenMatrices)
{
	// loop over the sig atom to produce the free atom rho etc.
	// just store it tmp in weights
	UInt nGrids = bg.getNGrids();
	DoubleVec tmp(nGrids,ZERO);
	for(UInt iAtom=0; iAtom<sigList.getNSigAtoms(); iAtom++) {

		// get the atom shell data
		const UIntVec& sigAtomList = sigList.getSigAtoms();
		UInt atomIndex = sigAtomList[iAtom];
		const AtomShell& as = ms.getAtomShell(atomIndex);

		// retrieve the free atom density matrix 
		// convert the data into sig order
		UInt atomic = as.getAtomic();
		const DenMtrx& denMtrx = atomDenMatrices.getAtomDenMtrx(atomic);
		DFTMatrix dftDen(sigList,sigList,iAtom,iAtom,denMtrx.getNSpin());
		for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {
			dftDen.intoSigOrder(iSpin,denMtrx.getMtrx(iSpin),as,as,sigList,sigList,iAtom,iAtom);
		}

		// form the new xcvar for the given free atom
		// right now we only have order 0 hirshfeld calculation
		UInt job = GROUND_STATE_DFT;
		UInt order = 0;
		XCVar xcvar(xcfunc,denMtrx.getNSpin(),job,order);
				
		// do the batch basis
		UInt nSigBas = sigList.getNSigBasis(iAtom);
		BatchBasis atomBasis(denMtrx.getNSpin(),basis.getMaxBasisDerivOrder(),basis.getNGrids(),nSigBas);
		atomBasis.batchBasisSplit(basis,sigList,iAtom);

		// now form batchVar
		BatchVar bVar(atomBasis,xcvar,dftDen,denMtrx);

		// now form spin resolved rho
		UInt pos = xcvar.getVarPos(RA);
		const Double* rhoA = bVar.getVar(pos);
		vaxpy(rhoA,weights.getPtr(0,iAtom),ONE,atomBasis.getNGrids());
		pos = xcvar.getVarPos(RB);
		const Double* rhoB = bVar.getVar(pos);
		vaxpy(rhoB,weights.getPtr(0,iAtom),ONE,atomBasis.getNGrids());

		// produce the point vector r^3: r = |point-sig atom|
      const Double* coord = bg.getGridCoord();
      for(UInt iG=0; iG<nGrids; iG++) {
         Double r  = as.getDistance(&coord[3*iG]);
         tmp[iG]   = r*r*r;
      }   

      // now turn to the free atom volumn result for 
      // this sig atom : \int r^3 \rho_{atom}(r) dr
		// this \rho_{atom}(r) is free atom density
      const Double* wts = bg.getGridWts();
      freeAtomVol[iAtom] = vdot3(&tmp.front(),weights.getPtr(0,iAtom),wts,nGrids);
	}

	// form the pro-molecule density
	tmp.assign(nGrids,ZERO);
	for(UInt iAtom=0; iAtom<sigList.getNSigAtoms(); iAtom++) {
		vaxpy(weights.getPtr(0,iAtom),&tmp.front(),ONE,nGrids);
	}

	// now let's finally form the weights
	// now make tmp vector as spin resolved total rho
	// here we need to pay attention to the case that divided by small number
	for(UInt iAtom=0; iAtom<sigList.getNSigAtoms(); iAtom++) {
		for(UInt iGrid=0; iGrid<nGrids; iGrid++) {
			Double molDen  = tmp[iGrid];
			Double atomDen = weights(iGrid,iAtom);
			if (molDen<thresh) {
				weights(iGrid,iAtom) = ZERO;
			}else{
				weights(iGrid,iAtom) = atomDen/molDen;
			}
		}
	}

	// add some debug printing
	// if you want to see the weights
	// comment below
	//weights.print("hirshfeld weights in batch");
}


