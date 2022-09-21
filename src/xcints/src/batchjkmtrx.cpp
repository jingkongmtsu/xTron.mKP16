/**
 * CPP files corresponding to the batchjkmtrx.h
 * \author  Fenglai Liu 
 */
#include "shell.h"
#include "batchbasis.h"
#include "batchgrid.h"
#include "blas.h"
#include "xcintsinfor.h"
#include "sigatombasis.h"
#include "halfjkrho.h"
#include "batchjkmtrx.h"
using namespace shell;
using namespace batchbasis;
using namespace batchgrid;
using namespace blas;
using namespace xcintsinfor;
using namespace sigatombasis;
using namespace halfjkrho;
using namespace batchjkmtrx;

BatchJKMtrx::BatchJKMtrx(const MolShell& ms, const SigAtomBasis& sigList, const HalfJKRho& halfJKRho,
		const BatchBasis& basis, const BatchGrid& grid, const XCIntJobInfor& infor):intJob(infor.getIntJob()),
	JMtrx(sigList,sigList,1),KMtrx(ms,sigList,infor.getNSpin())
{
	//
	// here we need to note that for both J and K matrix we will do add transpose
	// later in xcints. This is because we want to retain the hermetian property
	// of the fock matrix (M^{T} = M), usually this property will be damaged by
	// numerical calculation in xc part.
	//

	// firstly let's do K matrix calcualtion
	if (doK(intJob)) {

		// set up data
		UInt nSigBas= basis.getNSigBasis();
		UInt nGrids = grid.getNGrids();
		Mtrx halfFxc(nGrids,nSigBas);
		const SpinMatrix& halfExRho = halfJKRho.getHalfExRho(); 

		// now let's do each spin state
		for(UInt iSpin=0; iSpin<KMtrx.getNSpin(); iSpin++) {

			// scale the basis set value
			// we note that the scaling of HALF is because we want to 
			// cancel the effects of "add transpose" later in xcints
			// see the comment in BatchEXRhoMtrx constructor, because
			// of the partial derivatives, the 1/2 in front of energy
			// is cancelled for fock matrix building.
			// therefore, here the HALF cancells the "add transpose"
			// later in xcints, and it leads to correct result
			const Mtrx& phi = basis.getPhi(0);
			const Double* wts = grid.getGridWts();
			for(UInt i=0; i<nSigBas; i++) {
				vmul(wts,phi.getPtr(0,i),halfFxc.getPtr(0,i),nGrids);
				vscal(halfFxc.getPtr(0,i),HALF,nGrids);
			}

			// now we need to combine the halfFxc and 
			// the halfExRho together
			const Mtrx& hExRho = halfExRho.getMtrx(iSpin);
			Mtrx& Fxc = KMtrx.getMtrx(iSpin);
			Fxc.mult(hExRho,halfFxc,true,false,ONE,ZERO);
		}
	}

	// now let's do J matrix calculation
	if (doJ(intJob)) {

		// set up data
		UInt nSigBas= basis.getNSigBasis();
		UInt nGrids = grid.getNGrids();
		Mtrx halfFxc(nGrids,nSigBas);

		// scale the basis set value with the half Coulomb vector
		// we note that there's 1/2 in front of Coulomb energy
		// however, because we need to do add transpose later
		// in xcints for final j matrx, therefore we scale
		// the whole thing with (1/2)/2 = 0.25
		const Mtrx& phi = basis.getPhi(0);
		const Double* wts = grid.getGridWts();
		const DoubleVec& halfJRho = halfJKRho.getHalfCouRhoVec();
		for(UInt i=0; i<nSigBas; i++) {
			vmul3add(wts,phi.getPtr(0,i),&halfJRho[0],halfFxc.getPtr(0,i),nGrids);
			vscal(halfFxc.getPtr(0,i),0.25E0,nGrids);
		}

		// now let's combine with phi value again
		Mtrx& Fxc = JMtrx.getMtrx(0);
		Fxc.mult(phi,halfFxc,true,false,ONE,ZERO);
	}

	// finally, do we do JK together?
	// let's merge the J matrix into the K matrix
	if (doJKTogether(intJob)) {
		const Mtrx& J = JMtrx.getMtrx(0);
		UInt nSigBasis= basis.getNSigBasis();
		const UIntVec& sigBasisIndex = sigList.getGlobalSigBasisIndex();
		for(UInt iSpin=0; iSpin<KMtrx.getNSpin(); iSpin++) {
			Mtrx& Fxc = KMtrx.getMtrx(iSpin);
			for(UInt j=0; j<Fxc.getCol(); j++) {
				for(UInt i=0; i<nSigBasis; i++) {
					UInt rowIndex = sigBasisIndex[i]; 
					Fxc(rowIndex,j) += J.val(i,j);
				}
			}
		}
	}
}

void BatchJKMtrx::updateJKMtrx(const XCIntJobInfor& infor, 
		const SigAtomBasis& sigList, SpinMatrix& M, Mtrx& jMtrx) const 
{
	// update K matrix, or JK together
	if (doK(intJob)) {
		for(UInt iSpin=0; iSpin<M.getNSpin(); iSpin++) {
			KMtrx.intoNormOrder(sigList,iSpin,M.getMtrx(iSpin));
		}
	}

	// update J matrix, only with the situation that doing J only
	if (onlyDoJ(intJob)) {
		JMtrx.intoNormOrder(sigList,sigList,0,jMtrx);
	}
}
