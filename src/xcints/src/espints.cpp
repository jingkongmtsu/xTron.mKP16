/**
 * functions and classes associated with espints.h
 * \author Fenglai Liu 
 */
#include<cstdio>
#include<cstdlib>
#include "blas.h"
#include "blas1.h"
#include "excep.h"
#include "shell.h"
#include "matrix.h"
#include "cptrans.h"
#include "shellpair.h"
#include "gintsinfor.h"
#include "xcintsinfor.h"
#include "spinmatrix.h"
#include "hgp_os_ints.h"
#include "integraljobs.h"
#include "sigshellpairinfor.h"
#include "espints.h"
#include "sphereint.h"

using namespace blas;
using namespace excep;
using namespace shell;
using namespace matrix;
using namespace cptrans;
using namespace shellpair;
using namespace gintsinfor;
using namespace xcintsinfor;
using namespace spinmatrix;
using namespace integraljobs;
using namespace sigshellpairinfor;
using namespace espints;
using namespace sphereint;

void ESPInts::doESPIntWork(const DoubleVec& grids, const MolShell& ms, const SpinMatrix& denPhi, 
		const DoubleVec& cartDen, const SigMolShellPairInfor& infor, const DoubleVec& wRho, 
		SpinMatrix& halfExRho, DoubleVec& halfJRhoVec, DoubleVec& halfJRhoVecMtrx) const
{
	///////////////////////////////////////////////////////////////////////
	// setup scratch data for calculation                                //
	///////////////////////////////////////////////////////////////////////
	UInt nPts = grids.size()/3;

	// shell pair data
	UInt maxNP2 = infor.getMaxNP2();
	UInt maxNL  = infor.getMaxNL();
	UInt maxNSP = infor.getMaxNSP();
	AtomShellPair sp(maxNSP,maxNP2,maxNL);

	// integral array
	DoubleVec esp(infor.getMaxNBasisForShell()*infor.getMaxNBasisForShell()*nPts);

	// set up the scratch for reading the density matrix vector in
	// also set up scratch for the result of shell pair vector form J matrix 
	Mtrx cartAtomBlockDen;
	Mtrx cartAtomBlockJMtrx;
	Mtrx cartShellBlockDen;
	Mtrx cartShellBlockJMtrx;
	UInt intJob = ginfor.getIntJob();
	if (doJ(intJob)) {
		cartAtomBlockDen.init(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());
		cartAtomBlockJMtrx.init(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());
		cartShellBlockDen.init(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());
		cartShellBlockJMtrx.init(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());
	}

	// set up block form data for exchange energy density forming
	// for exchange energy density, the bra and ket shell pair will
	// produce different result, therefore we have denphi and halfExrho
	// result in bra/ket difference
	// also different spin result will also be different, a here
	// is alpha; b is for beta
	Mtrx aBraAtomBlockDenPhi;
	Mtrx aBraAtomBlockExRho;
	Mtrx bBraAtomBlockDenPhi;
	Mtrx bBraAtomBlockExRho;
	Mtrx aKetAtomBlockDenPhi;
	Mtrx aKetAtomBlockExRho;
	Mtrx bKetAtomBlockDenPhi;
	Mtrx bKetAtomBlockExRho;
	DoubleVec shellBlockDenPhi;
	DoubleVec shellBlockExRho;
	if (doK(intJob)) {

		// initialize the atom block input and output
		aBraAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
		aBraAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
		aKetAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
		aKetAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
		if (denPhi.getNSpin() == 2) {
			bBraAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
			bBraAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
			bKetAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
			bKetAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
		}

		// initialize the shell block input and output
		shellBlockDenPhi.assign(infor.getMaxNBasisForShell(),ZERO);
		shellBlockExRho.assign(infor.getMaxNBasisForShell(),ZERO);
	}

	//
	// we comment out the use of pMaxInfor
	//
	// we also need to keep data of pMaxInfor per each atom shell pair block
	// both rowShellPmaxInfor and colShellPmaxInfor are scratch space used
	// to extract data for the batch of grids
	// because we compute the results together for all of grid points,
	// therefore we need the pMaxInfor to store the final form result
	//UInt nMaxShell = infor.getMaxNShellForAtomShell();
	//Mtrx rowShellPmaxInfor(nMaxShell,nPts);
	//Mtrx colShellPmaxInfor(nMaxShell,nPts);
	//DoubleVec pMaxArray(nPts);

	///////////////////////////////////////////////////////////////////////
	// real working code                                                 //
	///////////////////////////////////////////////////////////////////////
	UInt cartDenOffset = 0;
	for( UInt iAtomShellPair=0; iAtomShellPair<infor.getNSigAtomShellPairs(); ++iAtomShellPair ) {

		//
		// get shell pair infor
		//
		const SigAtomShellPairInfor& sigAtomShellPairInfor = infor.getSigAtomShellPairInfor(iAtomShellPair); 

		//
		// form atom shell pair
		//
		sp.init(ms,ms,sigAtomShellPairInfor,ginfor.pickUpSPThresh());

		// set the dimension
		UInt rowNBas,colNBas;
		sp.getNCarBas(rowNBas,colNBas);

		// Coulomb part
		if (doJ(intJob)) {
			cartAtomBlockDen.reset(rowNBas,colNBas,false);
			cartAtomBlockJMtrx.reset(rowNBas,colNBas,true);
			vcopy(&cartDen[cartDenOffset],cartAtomBlockDen.getPtr(),rowNBas*colNBas);
		}

		// exchange part
		if (doK(intJob)) {

			// we need the offset data
			UInt rowBasOffset,colBasOffset;
			sp.getGlobalCarBasOffSet(rowBasOffset,colBasOffset);

			// for the grid part, it starts from 0
			UInt start = 0;

			// now let's initialize data and copy the denphi
			const Mtrx& aDenPhi = denPhi.getMtrx(0);
			aBraAtomBlockDenPhi.reset(rowNBas,nPts,false);
			aKetAtomBlockDenPhi.reset(colNBas,nPts,false);
			aBraAtomBlockDenPhi.copyFromMatrix(rowBasOffset,start,0,0,rowNBas,nPts,aDenPhi);
			aKetAtomBlockDenPhi.copyFromMatrix(colBasOffset,start,0,0,colNBas,nPts,aDenPhi);
			aBraAtomBlockExRho.reset(rowNBas,nPts,true);
			aKetAtomBlockExRho.reset(colNBas,nPts,true);

			// do it for possible beta spin
			if (denPhi.getNSpin() == 2) {
				const Mtrx& bDenPhi = denPhi.getMtrx(1);
				bBraAtomBlockDenPhi.reset(rowNBas,nPts,false);
				bKetAtomBlockDenPhi.reset(colNBas,nPts,false);
				bBraAtomBlockDenPhi.copyFromMatrix(rowBasOffset,start,0,0,rowNBas,nPts,bDenPhi);
				bKetAtomBlockDenPhi.copyFromMatrix(colBasOffset,start,0,0,colNBas,nPts,bDenPhi);
				bBraAtomBlockExRho.reset(rowNBas,nPts,true);
				bKetAtomBlockExRho.reset(colNBas,nPts,true);
			}
		}

		//
		//form Pmaxinfor for row atom shell and column atom shell
		//
		/*
			UInt rowShellIndex = -1;
			UInt colShellIndex = -1;
			UInt nRowShell     = -1;
			UInt nColShell     = -1;
			sp.getGlobalShellOffSet(rowShellIndex,colShellIndex);
			sp.getNShells(nRowShell,nColShell);
			bool reform = false;
			rowShellPmaxInfor.reset(nRowShell,nPts,reform);
			colShellPmaxInfor.reset(nColShell,nPts,reform);
			for(UInt iG=0; iG<nPts; iG++) {
			UInt gridOffset   = start+iG;
			if (rowShellIndex<colShellIndex) {
			const Double* p0  = pMaxInfor.getPtr(rowShellIndex,gridOffset);
			Double*       t0  = rowShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nRowShell; iShell++) {
			t0[iShell] = p0[iShell];
			}
			const Double* p1  = pMaxInfor.getPtr(colShellIndex,gridOffset);
			Double*       t1  = colShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nColShell; iShell++) {
			t1[iShell] = p1[iShell];
			}
			}else{
			const Double* p1  = pMaxInfor.getPtr(colShellIndex,gridOffset);
			Double*       t1  = colShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nColShell; iShell++) {
			t1[iShell] = p1[iShell];
			}
			const Double* p0  = pMaxInfor.getPtr(rowShellIndex,gridOffset);
			Double*       t0  = rowShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nRowShell; iShell++) {
			t0[iShell] = p0[iShell];
			}
			}
			}
			*/

		//
		//generate integral result
		//
		UInt nSP = sp.getNShellPairs(); 
		for(UInt iSP=0; iSP<nSP; iSP++) {

			// whether the real shell data is significant?
			const ShellPair& spData = sp.getShellPair(iSP);
			if (! spData.isSig()) {
				continue;
			}

			// now let's do integral job - get primitive pairs data
			LInt code         = static_cast<LInt>(spData.getLCode());
			UInt np2          = spData.getNP2();
			const Double* c2  = spData.getC2();
			const Double* e2  = spData.getE2();
			const Double* A   = spData.getA();
			const Double* B   = spData.getB();
			const Double* P   = spData.getP();
			const Double* fac = spData.getFac();

			// now form the pMaxArray
			// if shells are inversed we need to switch so that
			// to match the atom shell
			/*
				UInt iShellOffset = 0; 
				UInt jShellOffset = 0;
				spData.getLocalShellIndex(iShellOffset,jShellOffset);
				if (spData.inverseShells()) {
				UInt tmp     = iShellOffset;
				iShellOffset = jShellOffset;
				jShellOffset = tmp;
				}	

			// now let's derive the pMax value per each grid point
			// we just get the maximum value
			for(UInt iG=0; iG<nPts; iG++) {
			Double rowPMaxVal = rowShellPmaxInfor.val(iShellOffset,iG);
			Double colPMaxVal = colShellPmaxInfor.val(jShellOffset,iG);
			Double pMaxVal    = rowPMaxVal;
			if (rowPMaxVal<colPMaxVal) pMaxVal = colPMaxVal;
			pMaxArray[iG]     = pMaxVal;
			}
			*/

			// do the raw integral
			UInt iNCarBas = 0; 
			UInt jNCarBas = 0;
			spData.getNCarBas(iNCarBas,jNCarBas);
			UInt intLen = iNCarBas*jNCarBas*nPts;
			vset(&esp.front(),ZERO,intLen);
			hgp_os_esp(code,np2,nPts,c2,e2,fac,P,A,B,&grids[0],&esp.front()); 

			// get the global offset
			UInt iBasOffset = 0; 
			UInt jBasOffset = 0;
			spData.getGlobalCarBasOffSet(iBasOffset,jBasOffset);

			// set the local shell pair offset
			UInt iLocOffset,jLocOffset;  
			spData.getLocalCarBasOffSet(iLocOffset,jLocOffset);

			// if we do Coulomb part, now we need to get the block density matrix out
			// we need to consider the case that whether the two shells are switched
			// the copy is different if two shells are switched
			// exchange part different on each grid, so we will do it later
			if (doJ(intJob)) {
				UInt rowNBas = iNCarBas;
				UInt colNBas = jNCarBas;
				cartShellBlockDen.reset(rowNBas,colNBas,false);
				cartShellBlockJMtrx.reset(rowNBas,colNBas,true);
				if (spData.inverseShells()) {
					for(UInt i=0; i<rowNBas; i++) {
						UInt colIndex = iLocOffset + i;
						for(UInt j=0; j<colNBas; j++) {
							UInt rowIndex = jLocOffset + j;
							cartShellBlockDen(i,j) = cartAtomBlockDen.val(rowIndex,colIndex);
						}
					}
				}else{
					cartShellBlockDen.copyFromMatrix(iLocOffset,jLocOffset,0,0,rowNBas,colNBas,cartAtomBlockDen);
				}
			}

			// digestion
			for(UInt iPt=0; iPt<nPts; iPt++) {

				// get the esp data
				const Double* espData = &esp[iNCarBas*jNCarBas*iPt];

				// now let's do digestion
				// firstly it's exchange part
				if (doK(intJob)) {

					//
					// this is alpha spin part
					//

					// firstly let's do the bra side denphi digestion
					// this will produce the ket side result
					// iLocOffset and iNCarBas are all associated with the bra shell block index
					// get the denphi
					if (spData.inverseShells()) {
						vcopy(aKetAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
					}else{
						vcopy(aBraAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
					}

					// digest data to form the result, result is in ket side shell block
					// jLocOffset and jNCarBas are all associated with the ket shell block index
					for(UInt j=0; j<jNCarBas; j++) {
						shellBlockExRho[j] = MINUS_ONE*vdot(&espData[j*iNCarBas],&shellBlockDenPhi[0],iNCarBas);
					}

					// let's put the result back to atom shell block
					// the result is in ket shell block
					if (spData.inverseShells()) {
						vaxpy(&shellBlockExRho[0],aBraAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
					}else{
						vaxpy(&shellBlockExRho[0],aKetAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
					}

					// if the two shells are different, we need to do the digestion for the ket side denphi
					if (jBasOffset != iBasOffset) {

						// read in denphi
						// this is the ket side denphi digestion
						// jLocOffset and jNCarBas are all associated with the ket shell block index
						if (spData.inverseShells()) {
							vcopy(aBraAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
						}else{
							vcopy(aKetAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
						}

						// form results
						// result is in bra side shell block
						// iLocOffset and iNCarBas are all associated with the ket shell block index
						for(UInt i=0; i<iNCarBas; i++) shellBlockExRho[i] = ZERO;
						for(UInt i=0; i<iNCarBas; i++) {
							for(UInt j=0; j<jNCarBas; j++) {
								shellBlockExRho[i] += MINUS_ONE*espData[i+j*iNCarBas]*shellBlockDenPhi[j];
							}
						}

						// let's put the result into the atom block result
						// remember that the result is in bra side
						// iLocOffset and iNCarBas are all associated with the bra shell block index
						if (spData.inverseShells()) {
							vaxpy(&shellBlockExRho[0],aKetAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
						}else{
							vaxpy(&shellBlockExRho[0],aBraAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
						}
					}

					//
					// this is beta spin part
					//
					if (denPhi.getNSpin() == 2) {

						// firstly let's do the bra side denphi digestion
						// this will produce the ket side result
						// iLocOffset and iNCarBas are all associated with the bra shell block index
						// get the denphi
						if (spData.inverseShells()) {
							vcopy(bKetAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
						}else{
							vcopy(bBraAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
						}

						// digest data to form the result, result is in ket side shell block
						// jLocOffset and jNCarBas are all associated with the ket shell block index
						for(UInt j=0; j<jNCarBas; j++) {
							shellBlockExRho[j] = MINUS_ONE*vdot(&espData[j*iNCarBas],&shellBlockDenPhi[0],iNCarBas);
						}

						// let's put the result back to atom shell block
						// the result is in ket shell block
						if (spData.inverseShells()) {
							vaxpy(&shellBlockExRho[0],bBraAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
						}else{
							vaxpy(&shellBlockExRho[0],bKetAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
						}

						// if the two shells are different, we need to do the digestion for the ket side denphi
						if (jBasOffset != iBasOffset) {

							// read in denphi
							// this is the ket side denphi digestion
							// jLocOffset and jNCarBas are all associated with the ket shell block index
							if (spData.inverseShells()) {
								vcopy(bBraAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
							}else{
								vcopy(bKetAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
							}

							// form results
							// result is in bra side shell block
							// iLocOffset and iNCarBas are all associated with the ket shell block index
							for(UInt i=0; i<iNCarBas; i++) shellBlockExRho[i] = ZERO;
							for(UInt i=0; i<iNCarBas; i++) {
								for(UInt j=0; j<jNCarBas; j++) {
									shellBlockExRho[i] += MINUS_ONE*espData[i+j*iNCarBas]*shellBlockDenPhi[j];
								}
							}

							// let's put the result into the atom block result
							// remember that the result is in bra side
							// iLocOffset and iNCarBas are all associated with the bra shell block index
							if (spData.inverseShells()) {
								vaxpy(&shellBlockExRho[0],bKetAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
							}else{
								vaxpy(&shellBlockExRho[0],bBraAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
							}
						}
					}
				}

				// now let's do Coulomb
				if (doJ(intJob)) {

					// this is just a dot product
					Double jVal  = vdot(&espData[0],cartShellBlockDen.getPtr(),iNCarBas*jNCarBas);

					// because for the coulomb part the mu and nu loop over all of index, therefore
					// there will be a factor of 2 for the mu != nu
					if (jBasOffset != iBasOffset) jVal = jVal*TWO;

					// now it's result
					halfJRhoVec[iPt] += jVal;

					// now let's do the matrix part
					Double wrho = wRho[iPt];
					vaxpy(&espData[0],cartShellBlockJMtrx.getPtr(0,0),wrho,iNCarBas*jNCarBas);
				}
			}

			// at the end of shell pair block, we need to update the 
			// shell block J matrix to atom block
			if (doJ(intJob)) {

				// firstly, let's see whether the shell block data is significant enough
				/*
				bool isSig = false;
				const DoubleVec& shellBlockVec = cartShellBlockJMtrx.getVec();
				for (UInt i=0; i<iNCarBas*jNCarBas; i++) {
					if (fabs(shellBlockVec[i])>ginfor.pickUpThresh()) {
						isSig = true;
						break;
					}
				}
				*/

				// we only copy the data when it's significant
				if (spData.inverseShells()) {
					UInt rowNBas = iNCarBas;
					UInt colNBas = jNCarBas;
					for(UInt i=0; i<rowNBas; i++) {
						UInt colIndex = iLocOffset + i;
						for(UInt j=0; j<colNBas; j++) {
							UInt rowIndex = jLocOffset + j;
							cartAtomBlockJMtrx(rowIndex,colIndex) += cartShellBlockJMtrx(i,j);
						}
					}
				}else{
					cartShellBlockJMtrx.updateToMatrix(0,0,iLocOffset,jLocOffset,iNCarBas,jNCarBas,
							ONE,cartAtomBlockJMtrx);
				}
			}
		}

		// debug j matrix
		//cartAtomBlockJMtrx.print("atom block matrix");

		// now after the calculation on atom block
		// we need to write the result back to the global ones
		if (doJ(intJob)) {
			// for Coulomb, only the J matrix part need to be updated here
			UInt len = cartAtomBlockJMtrx.getRow()*cartAtomBlockJMtrx.getCol();
			vcopy(cartAtomBlockJMtrx.getPtr(),&halfJRhoVecMtrx[cartDenOffset],len);
			cartDenOffset += len;
		}

		// for exchange, we just copy them to the result
		if (doK(intJob)) {

			// set the offset and dimension
			UInt rowBasOffset,colBasOffset;
			sp.getGlobalCarBasOffSet(rowBasOffset,colBasOffset);
			UInt rowNBas,colNBas;
			sp.getNCarBas(rowNBas,colNBas);

			// the grid set starts from 0
			UInt start = 0;

			// now let's copy the data
			Mtrx& aHalfExRho = halfExRho.getMtrx(0);
			aBraAtomBlockExRho.updateToMatrix(0,0,rowBasOffset,start,rowNBas,nPts,ONE,aHalfExRho);
			aKetAtomBlockExRho.updateToMatrix(0,0,colBasOffset,start,colNBas,nPts,ONE,aHalfExRho);

			// do it for possible beta spin
			if (halfExRho.getNSpin() == 2) {
				Mtrx& bHalfExRho = halfExRho.getMtrx(1);
				bBraAtomBlockExRho.updateToMatrix(0,0,rowBasOffset,start,rowNBas,nPts,ONE,bHalfExRho);
				bKetAtomBlockExRho.updateToMatrix(0,0,colBasOffset,start,colNBas,nPts,ONE,bHalfExRho);
			}
		}
	}
}

void ESPInts::doMtrx(const XCIntJobInfor& xcInfor, const MolShell& ms, const SigMolShellPairInfor& spInfor, 
		const DoubleVec& grid, const SpinMatrix& denPhi, const DoubleVec& cartDen, const DoubleVec& wRho,
		SpinMatrix& halfExRhoResult, DoubleVec& halfJRhoVec, DoubleVec& halfJRhoVecMtrx) const
{
	// we need a temp result to hold the exchange rho
	// this is because we loop over grid points, it's better
	// to set the grid dimension on the column dimension
	UInt nGrids = grid.size()/3;
	UInt intJob = ginfor.getIntJob();
	SpinMatrix halfExRho(denPhi.getNSpin());
	if (doK(intJob)) {
		halfExRho.init(ms.getNCarBas(),nGrids);
	}

	// possible timing code
	tick_count t0 = tick_count::now();

	// now do the job
	doESPIntWork(grid,ms,denPhi,cartDen,spInfor,wRho,halfExRho,halfJRhoVec,halfJRhoVecMtrx);

	// possible timing code
	tick_count t1 = tick_count::now();
	Double t = (t1-t0).seconds();
	//printf("espints work with number of grids %d, time in seconds %-12.6f\n ", (Int)nGrids, t);

	// we need to do final transformation to the exchange rho
	if (doK(intJob)) {

		// we need to scale the column dimension
		const GlobalInfor& globInfor = ginfor.getGlobalInfor();
		CPTransBasisNorm trans(globInfor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_COL);

		// however, the final result is (nGrid,nBas)
		// therefore we need to do transpose
		// then add in data
		for(UInt iSpin=0; iSpin<halfExRho.getNSpin(); iSpin++) {

			// transpose it, the nBas dimension 
			// should be in second dimension
			Mtrx& M = halfExRho.getMtrx(iSpin);
			bool inPlace = false;
			M.transpose(inPlace);

			// transform the col data
			trans.transform(ms,ms,M);

			Mtrx& T = halfExRhoResult.getMtrx(iSpin);
			T.add(M);
		}
	}
}

void ESPInts::doESPIntWorkXHole(const DoubleVec& grids, const MolShell& ms, const SpinMatrix& denPhi, 
		const DoubleVec& cartDen, const SigMolShellPairInfor& infor, const DoubleVec& wRho, 
			   SpinMatrix& halfExRho, DoubleVec& halfJRhoVec, DoubleVec& halfJRhoVecMtrx, Double sValue) const//yw
{
	///////////////////////////////////////////////////////////////////////
	// setup scratch data for calculation                                //
	///////////////////////////////////////////////////////////////////////
	UInt nPts = grids.size()/3;

	// shell pair data
	UInt maxNP2 = infor.getMaxNP2();
	UInt maxNL  = infor.getMaxNL();
	UInt maxNSP = infor.getMaxNSP();
	AtomShellPair sp(maxNSP,maxNP2,maxNL);

	// integral array
	DoubleVec esp(infor.getMaxNBasisForShell()*infor.getMaxNBasisForShell()*nPts);

	// set up the scratch for reading the density matrix vector in
	// also set up scratch for the result of shell pair vector form J matrix 
	Mtrx cartAtomBlockDen;
	Mtrx cartAtomBlockJMtrx;
	Mtrx cartShellBlockDen;
	Mtrx cartShellBlockJMtrx;
	UInt intJob = ginfor.getIntJob();
	if (doJ(intJob)) {
		cartAtomBlockDen.init(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());
		cartAtomBlockJMtrx.init(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());
		cartShellBlockDen.init(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());
		cartShellBlockJMtrx.init(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());
	}

	// set up block form data for exchange energy density forming
	// for exchange energy density, the bra and ket shell pair will
	// produce different result, therefore we have denphi and halfExrho
	// result in bra/ket difference
	// also different spin result will also be different, a here
	// is alpha; b is for beta
	Mtrx aBraAtomBlockDenPhi;
	Mtrx aBraAtomBlockExRho;
	Mtrx bBraAtomBlockDenPhi;
	Mtrx bBraAtomBlockExRho;
	Mtrx aKetAtomBlockDenPhi;
	Mtrx aKetAtomBlockExRho;
	Mtrx bKetAtomBlockDenPhi;
	Mtrx bKetAtomBlockExRho;
	DoubleVec shellBlockDenPhi;
	DoubleVec shellBlockExRho;
	if (doK(intJob)) {

		// initialize the atom block input and output
		aBraAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
		aBraAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
		aKetAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
		aKetAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
		if (denPhi.getNSpin() == 2) {
			bBraAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
			bBraAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
			bKetAtomBlockDenPhi.init(infor.getMaxNBasisForAtomShell(),nPts);
			bKetAtomBlockExRho.init(infor.getMaxNBasisForAtomShell(),nPts);
		}

		// initialize the shell block input and output
		shellBlockDenPhi.assign(infor.getMaxNBasisForShell(),ZERO);
		shellBlockExRho.assign(infor.getMaxNBasisForShell(),ZERO);
	}

	//
	// we comment out the use of pMaxInfor
	//
	// we also need to keep data of pMaxInfor per each atom shell pair block
	// both rowShellPmaxInfor and colShellPmaxInfor are scratch space used
	// to extract data for the batch of grids
	// because we compute the results together for all of grid points,
	// therefore we need the pMaxInfor to store the final form result
	//UInt nMaxShell = infor.getMaxNShellForAtomShell();
	//Mtrx rowShellPmaxInfor(nMaxShell,nPts);
	//Mtrx colShellPmaxInfor(nMaxShell,nPts);
	//DoubleVec pMaxArray(nPts);

	///////////////////////////////////////////////////////////////////////
	// real working code                                                 //
	///////////////////////////////////////////////////////////////////////
	UInt cartDenOffset = 0;
	for( UInt iAtomShellPair=0; iAtomShellPair<infor.getNSigAtomShellPairs(); ++iAtomShellPair ) {

		//
		// get shell pair infor
		//
		const SigAtomShellPairInfor& sigAtomShellPairInfor = infor.getSigAtomShellPairInfor(iAtomShellPair); 

		//
		// form atom shell pair
		//
		sp.init(ms,ms,sigAtomShellPairInfor,ginfor.pickUpSPThresh());

		// set the dimension
		UInt rowNBas,colNBas;
		sp.getNCarBas(rowNBas,colNBas);

		// Coulomb part
		if (doJ(intJob)) {
			cartAtomBlockDen.reset(rowNBas,colNBas,false);
			cartAtomBlockJMtrx.reset(rowNBas,colNBas,true);
			vcopy(&cartDen[cartDenOffset],cartAtomBlockDen.getPtr(),rowNBas*colNBas);
		}

		// exchange part
		if (doK(intJob)) {

			// we need the offset data
			UInt rowBasOffset,colBasOffset;
			sp.getGlobalCarBasOffSet(rowBasOffset,colBasOffset);

			// for the grid part, it starts from 0
			UInt start = 0;

			// now let's initialize data and copy the denphi
			const Mtrx& aDenPhi = denPhi.getMtrx(0);
			aBraAtomBlockDenPhi.reset(rowNBas,nPts,false);
			aKetAtomBlockDenPhi.reset(colNBas,nPts,false);
			aBraAtomBlockDenPhi.copyFromMatrix(rowBasOffset,start,0,0,rowNBas,nPts,aDenPhi);
			aKetAtomBlockDenPhi.copyFromMatrix(colBasOffset,start,0,0,colNBas,nPts,aDenPhi);
			aBraAtomBlockExRho.reset(rowNBas,nPts,true);
			aKetAtomBlockExRho.reset(colNBas,nPts,true);

			// do it for possible beta spin
			if (denPhi.getNSpin() == 2) {
				const Mtrx& bDenPhi = denPhi.getMtrx(1);
				bBraAtomBlockDenPhi.reset(rowNBas,nPts,false);
				bKetAtomBlockDenPhi.reset(colNBas,nPts,false);
				bBraAtomBlockDenPhi.copyFromMatrix(rowBasOffset,start,0,0,rowNBas,nPts,bDenPhi);
				bKetAtomBlockDenPhi.copyFromMatrix(colBasOffset,start,0,0,colNBas,nPts,bDenPhi);
				bBraAtomBlockExRho.reset(rowNBas,nPts,true);
				bKetAtomBlockExRho.reset(colNBas,nPts,true);
			}
		}

		//
		//form Pmaxinfor for row atom shell and column atom shell
		//
		/*
			UInt rowShellIndex = -1;
			UInt colShellIndex = -1;
			UInt nRowShell     = -1;
			UInt nColShell     = -1;
			sp.getGlobalShellOffSet(rowShellIndex,colShellIndex);
			sp.getNShells(nRowShell,nColShell);
			bool reform = false;
			rowShellPmaxInfor.reset(nRowShell,nPts,reform);
			colShellPmaxInfor.reset(nColShell,nPts,reform);
			for(UInt iG=0; iG<nPts; iG++) {
			UInt gridOffset   = start+iG;
			if (rowShellIndex<colShellIndex) {
			const Double* p0  = pMaxInfor.getPtr(rowShellIndex,gridOffset);
			Double*       t0  = rowShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nRowShell; iShell++) {
			t0[iShell] = p0[iShell];
			}
			const Double* p1  = pMaxInfor.getPtr(colShellIndex,gridOffset);
			Double*       t1  = colShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nColShell; iShell++) {
			t1[iShell] = p1[iShell];
			}
			}else{
			const Double* p1  = pMaxInfor.getPtr(colShellIndex,gridOffset);
			Double*       t1  = colShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nColShell; iShell++) {
			t1[iShell] = p1[iShell];
			}
			const Double* p0  = pMaxInfor.getPtr(rowShellIndex,gridOffset);
			Double*       t0  = rowShellPmaxInfor.getPtr(0,iG);
			for(UInt iShell=0; iShell<nRowShell; iShell++) {
			t0[iShell] = p0[iShell];
			}
			}
			}
			*/

		//
		//generate integral result
		UInt nSP = sp.getNShellPairs(); 
		for(UInt iSP=0; iSP<nSP; iSP++) {

			// whether the real shell data is significant?
			const ShellPair& spData = sp.getShellPair(iSP);
			if (! spData.isSig()) {
				continue;
			}

			// now let's do integral job - get primitive pairs data
			LInt code         = static_cast<LInt>(spData.getLCode());
			UInt np2          = spData.getNP2();
			const Double* c2  = spData.getC2();
			const Double* e2  = spData.getE2();
			const Double* e2diff = spData.getExpDiff();//yw
			const Double* A   = spData.getA();
			const Double* B   = spData.getB();
			const Double* P   = spData.getP();
			const Double* fac = spData.getFac();

			// now form the pMaxArray
			// if shells are inversed we need to switch so that
			// to match the atom shell
			/*
				UInt iShellOffset = 0; 
				UInt jShellOffset = 0;
				spData.getLocalShellIndex(iShellOffset,jShellOffset);
				if (spData.inverseShells()) {
				UInt tmp     = iShellOffset;
				iShellOffset = jShellOffset;
				jShellOffset = tmp;
				}	

			// now let's derive the pMax value per each grid point
			// we just get the maximum value
			for(UInt iG=0; iG<nPts; iG++) {
			Double rowPMaxVal = rowShellPmaxInfor.val(iShellOffset,iG);
			Double colPMaxVal = colShellPmaxInfor.val(jShellOffset,iG);
			Double pMaxVal    = rowPMaxVal;
			if (rowPMaxVal<colPMaxVal) pMaxVal = colPMaxVal;
			pMaxArray[iG]     = pMaxVal;
			}
			*/

			// do the raw integral
			UInt iNCarBas = 0; 
			UInt jNCarBas = 0;
			spData.getNCarBas(iNCarBas,jNCarBas);
			UInt intLen = iNCarBas*jNCarBas*nPts;
			vset(&esp.front(),ZERO,intLen);
			
                        //hgp_os_esp(code,np2,nPts,c2,e2,fac,P,A,B,&grids[0],&esp.front()); 
			sphereInt(code,np2,nPts,c2,e2,P,A,B,&grids[0],&esp.front(),sValue); //yw
			

			// get the global offset
			UInt iBasOffset = 0; 
			UInt jBasOffset = 0;
			spData.getGlobalCarBasOffSet(iBasOffset,jBasOffset);

			// set the local shell pair offset
			UInt iLocOffset,jLocOffset;  
			spData.getLocalCarBasOffSet(iLocOffset,jLocOffset);

			// if we do Coulomb part, now we need to get the block density matrix out
			// we need to consider the case that whether the two shells are switched
			// the copy is different if two shells are switched
			// exchange part different on each grid, so we will do it later
			if (doJ(intJob)) {
				UInt rowNBas = iNCarBas;
				UInt colNBas = jNCarBas;
				cartShellBlockDen.reset(rowNBas,colNBas,false);
				cartShellBlockJMtrx.reset(rowNBas,colNBas,true);
				if (spData.inverseShells()) {
					for(UInt i=0; i<rowNBas; i++) {
						UInt colIndex = iLocOffset + i;
						for(UInt j=0; j<colNBas; j++) {
							UInt rowIndex = jLocOffset + j;
							cartShellBlockDen(i,j) = cartAtomBlockDen.val(rowIndex,colIndex);
						}
					}
				}else{
					cartShellBlockDen.copyFromMatrix(iLocOffset,jLocOffset,0,0,rowNBas,colNBas,cartAtomBlockDen);
				}
			}

			// digestion
			for(UInt iPt=0; iPt<nPts; iPt++) {

				// get the esp data
				const Double* espData = &esp[iNCarBas*jNCarBas*iPt];

				// now let's do digestion
				// firstly it's exchange part
				if (doK(intJob)) {

					//
					// this is alpha spin part
					//

					// firstly let's do the bra side denphi digestion
					// this will produce the ket side result
					// iLocOffset and iNCarBas are all associated with the bra shell block index
					// get the denphi
					if (spData.inverseShells()) {
						vcopy(aKetAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
					}else{
						vcopy(aBraAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
					}

					// digest data to form the result, result is in ket side shell block
					// jLocOffset and jNCarBas are all associated with the ket shell block index
					for(UInt j=0; j<jNCarBas; j++) {
						shellBlockExRho[j] = MINUS_ONE*vdot(&espData[j*iNCarBas],&shellBlockDenPhi[0],iNCarBas);
					}

					// let's put the result back to atom shell block
					// the result is in ket shell block
					if (spData.inverseShells()) {
						vaxpy(&shellBlockExRho[0],aBraAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
					}else{
						vaxpy(&shellBlockExRho[0],aKetAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
					}

					// if the two shells are different, we need to do the digestion for the ket side denphi
					if (jBasOffset != iBasOffset) {

						// read in denphi
						// this is the ket side denphi digestion
						// jLocOffset and jNCarBas are all associated with the ket shell block index
						if (spData.inverseShells()) {
							vcopy(aBraAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
						}else{
							vcopy(aKetAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
						}

						// form results
						// result is in bra side shell block
						// iLocOffset and iNCarBas are all associated with the ket shell block index
						for(UInt i=0; i<iNCarBas; i++) shellBlockExRho[i] = ZERO;
						for(UInt i=0; i<iNCarBas; i++) {
							for(UInt j=0; j<jNCarBas; j++) {
								shellBlockExRho[i] += MINUS_ONE*espData[i+j*iNCarBas]*shellBlockDenPhi[j];
							}
						}

						// let's put the result into the atom block result
						// remember that the result is in bra side
						// iLocOffset and iNCarBas are all associated with the bra shell block index
						if (spData.inverseShells()) {
							vaxpy(&shellBlockExRho[0],aKetAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
						}else{
							vaxpy(&shellBlockExRho[0],aBraAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
						}
					}

					//
					// this is beta spin part
					//
					if (denPhi.getNSpin() == 2) {

						// firstly let's do the bra side denphi digestion
						// this will produce the ket side result
						// iLocOffset and iNCarBas are all associated with the bra shell block index
						// get the denphi
						if (spData.inverseShells()) {
							vcopy(bKetAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
						}else{
							vcopy(bBraAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
						}

						// digest data to form the result, result is in ket side shell block
						// jLocOffset and jNCarBas are all associated with the ket shell block index
						for(UInt j=0; j<jNCarBas; j++) {
							shellBlockExRho[j] = MINUS_ONE*vdot(&espData[j*iNCarBas],&shellBlockDenPhi[0],iNCarBas);
						}

						// let's put the result back to atom shell block
						// the result is in ket shell block
						if (spData.inverseShells()) {
							vaxpy(&shellBlockExRho[0],bBraAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
						}else{
							vaxpy(&shellBlockExRho[0],bKetAtomBlockExRho.getPtr(jLocOffset,iPt),ONE,jNCarBas);
						}

						// if the two shells are different, we need to do the digestion for the ket side denphi
						if (jBasOffset != iBasOffset) {

							// read in denphi
							// this is the ket side denphi digestion
							// jLocOffset and jNCarBas are all associated with the ket shell block index
							if (spData.inverseShells()) {
								vcopy(bBraAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
							}else{
								vcopy(bKetAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
							}

							// form results
							// result is in bra side shell block
							// iLocOffset and iNCarBas are all associated with the ket shell block index
							for(UInt i=0; i<iNCarBas; i++) shellBlockExRho[i] = ZERO;
							for(UInt i=0; i<iNCarBas; i++) {
								for(UInt j=0; j<jNCarBas; j++) {
									shellBlockExRho[i] += MINUS_ONE*espData[i+j*iNCarBas]*shellBlockDenPhi[j];
								}
							}

							// let's put the result into the atom block result
							// remember that the result is in bra side
							// iLocOffset and iNCarBas are all associated with the bra shell block index
							if (spData.inverseShells()) {
								vaxpy(&shellBlockExRho[0],bKetAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
							}else{
								vaxpy(&shellBlockExRho[0],bBraAtomBlockExRho.getPtr(iLocOffset,iPt),ONE,iNCarBas);
							}
						}
					}
				}

				// now let's do Coulomb
				if (doJ(intJob)) {

					// this is just a dot product
					Double jVal  = vdot(&espData[0],cartShellBlockDen.getPtr(),iNCarBas*jNCarBas);

					// because for the coulomb part the mu and nu loop over all of index, therefore
					// there will be a factor of 2 for the mu != nu
					if (jBasOffset != iBasOffset) jVal = jVal*TWO;

					// now it's result
					halfJRhoVec[iPt] += jVal;

					// now let's do the matrix part
					Double wrho = wRho[iPt];
					vaxpy(&espData[0],cartShellBlockJMtrx.getPtr(0,0),wrho,iNCarBas*jNCarBas);
				}
			}

			// at the end of shell pair block, we need to update the 
			// shell block J matrix to atom block
			if (doJ(intJob)) {

				// firstly, let's see whether the shell block data is significant enough
				/*
				bool isSig = false;
				const DoubleVec& shellBlockVec = cartShellBlockJMtrx.getVec();
				for (UInt i=0; i<iNCarBas*jNCarBas; i++) {
					if (fabs(shellBlockVec[i])>ginfor.pickUpThresh()) {
						isSig = true;
						break;
					}
				}
				*/

				// we only copy the data when it's significant
				if (spData.inverseShells()) {
					UInt rowNBas = iNCarBas;
					UInt colNBas = jNCarBas;
					for(UInt i=0; i<rowNBas; i++) {
						UInt colIndex = iLocOffset + i;
						for(UInt j=0; j<colNBas; j++) {
							UInt rowIndex = jLocOffset + j;
							cartAtomBlockJMtrx(rowIndex,colIndex) += cartShellBlockJMtrx(i,j);
						}
					}
				}else{
					cartShellBlockJMtrx.updateToMatrix(0,0,iLocOffset,jLocOffset,iNCarBas,jNCarBas,
							ONE,cartAtomBlockJMtrx);
				}
			}
		}

		// debug j matrix
		//cartAtomBlockJMtrx.print("atom block matrix");

		// now after the calculation on atom block
		// we need to write the result back to the global ones
		if (doJ(intJob)) {
			// for Coulomb, only the J matrix part need to be updated here
			UInt len = cartAtomBlockJMtrx.getRow()*cartAtomBlockJMtrx.getCol();
			vcopy(cartAtomBlockJMtrx.getPtr(),&halfJRhoVecMtrx[cartDenOffset],len);
			cartDenOffset += len;
		}

		// for exchange, we just copy them to the result
		if (doK(intJob)) {

			// set the offset and dimension
			UInt rowBasOffset,colBasOffset;
			sp.getGlobalCarBasOffSet(rowBasOffset,colBasOffset);
			UInt rowNBas,colNBas;
			sp.getNCarBas(rowNBas,colNBas);

			// the grid set starts from 0
			UInt start = 0;

			// now let's copy the data
			Mtrx& aHalfExRho = halfExRho.getMtrx(0);
			aBraAtomBlockExRho.updateToMatrix(0,0,rowBasOffset,start,rowNBas,nPts,ONE,aHalfExRho);
			aKetAtomBlockExRho.updateToMatrix(0,0,colBasOffset,start,colNBas,nPts,ONE,aHalfExRho);

			// do it for possible beta spin
			if (halfExRho.getNSpin() == 2) {
				Mtrx& bHalfExRho = halfExRho.getMtrx(1);
				bBraAtomBlockExRho.updateToMatrix(0,0,rowBasOffset,start,rowNBas,nPts,ONE,bHalfExRho);
				bKetAtomBlockExRho.updateToMatrix(0,0,colBasOffset,start,colNBas,nPts,ONE,bHalfExRho);
			}
		}
	}
}

void ESPInts::doMtrxXHole(const XCIntJobInfor& xcInfor, const MolShell& ms, const SigMolShellPairInfor& spInfor, 
		const DoubleVec& grid, const SpinMatrix& denPhi, const DoubleVec& cartDen, const DoubleVec& wRho,
		     SpinMatrix& halfExRhoResult, DoubleVec& halfJRhoVec, DoubleVec& halfJRhoVecMtrx, Double sValue) const//yw
{
	// we need a temp result to hold the exchange rho
	// this is because we loop over grid points, it's better
	// to set the grid dimension on the column dimension
	UInt nGrids = grid.size()/3;
	UInt intJob = ginfor.getIntJob();
	SpinMatrix halfExRho(denPhi.getNSpin());
	if (doK(intJob)) {
		halfExRho.init(ms.getNCarBas(),nGrids);
	}

	// possible timing code
	tick_count t0 = tick_count::now();

	// now do the job
	doESPIntWorkXHole(grid,ms,denPhi,cartDen,spInfor,wRho,halfExRho,halfJRhoVec,halfJRhoVecMtrx, sValue);//yw

	// possible timing code
	tick_count t1 = tick_count::now();
	Double t = (t1-t0).seconds();
	//printf("espints work with number of grids %d, time in seconds %-12.6f\n ", (Int)nGrids, t);

	// we need to do final transformation to the exchange rho
	if (doK(intJob)) {

		// we need to scale the column dimension
		const GlobalInfor& globInfor = ginfor.getGlobalInfor();
		CPTransBasisNorm trans(globInfor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_COL);

		// however, the final result is (nGrid,nBas)
		// therefore we need to do transpose
		// then add in data
		for(UInt iSpin=0; iSpin<halfExRho.getNSpin(); iSpin++) {

			// transpose it, the nBas dimension 
			// should be in second dimension
			Mtrx& M = halfExRho.getMtrx(iSpin);
			bool inPlace = false;
			M.transpose(inPlace);

			// transform the col data
			trans.transform(ms,ms,M);

			Mtrx& T = halfExRhoResult.getMtrx(iSpin);
			T.add(M);
		}
	}
}
