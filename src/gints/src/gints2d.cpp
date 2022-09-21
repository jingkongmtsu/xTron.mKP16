/**
 * functions and classes associated with gints2d.h
 * \author Fenglai Liu 
 */
#include<cstdio>
#include "blas1.h"
#include "blas.h"
#include "excep.h"
#include "shell.h"
#include "matrix.h"
#include "cptrans.h"
#include "molecule.h"
#include "shellpair.h"
#include "gintsinfor.h"
#include "blockmatrix.h"
#include "integraljobs.h"
#include "hgp_os_ints.h"
#include "sigshellpairinfor.h"
#include "gints2d.h"
using namespace blas;
using namespace excep;
using namespace shell;
using namespace matrix;
using namespace cptrans;
using namespace molecule;
using namespace shellpair;
using namespace gintsinfor;
using namespace blockmatrix;
using namespace integraljobs;
using namespace sigshellpairinfor;
using namespace gints2d;

TBB_GInts2DIntMatrix::TBB_GInts2DIntMatrix(const GIntJobInfor& ginfor0, 
		const SigMolShellPairInfor& shellPairInfor, const MolShell& rowShell, 
		const MolShell& colShell, const Molecule& mol0, Mtrx& result):mol(mol0),
	rs(rowShell),cs(colShell),infor(shellPairInfor),
			ginfor(ginfor0),intMtrx(result),coord(3*mol.getNAtoms()),Z(mol.getNAtoms()) 
{ 
	// stuff the Z and coord array
	UInt natoms = mol.getNAtoms();
	for(UInt iAtom=0; iAtom<natoms; iAtom++) {
		Z[iAtom] = mol.getAtomic(iAtom);
		const Double* atomXYZ = mol.getXYZ(iAtom);
		coord[0+iAtom*3] = atomXYZ[0];
		coord[1+iAtom*3] = atomXYZ[1];
		coord[2+iAtom*3] = atomXYZ[2];
	}

	// finally let's check the job - whether it matches here
	UInt job = ginfor.getIntJob();
	if (job != TWO_BODY_OVERLAP && job != KINETIC && job != NUCLEAR_ATTRACTION) {
		string infor = "The integral job can not be processed here";
		Excep excep("TBB_GInts2DIntMatrix","constructor",EXCEPTION_GINTS_INVALID_INTEGRAL_JOB,infor);
		handleExcep(excep);
	}
}

TBB_GInts2DMultiIntMatrix::TBB_GInts2DMultiIntMatrix(const GIntJobInfor& ginfor0, 
		const SigMolShellPairInfor& shellPairInfor, const MolShell& rowShell, 
		const MolShell& colShell, MtrxVec& result):rs(rowShell),cs(colShell),infor(shellPairInfor),
	ginfor(ginfor0),intMtrxVec(result)
{ 
	// for this constructor, the center of C is all zero
	C[0] = ZERO;
	C[1] = ZERO;
	C[2] = ZERO;

	// finally let's check the job - whether it matches here
	UInt job = ginfor.getIntJob();
	if (job != MOM_P) {
		string infor = "The integral job can not be processed here";
		Excep excep("TBB_GInts2DMultiIntMatrix","constructor",EXCEPTION_GINTS_INVALID_INTEGRAL_JOB,infor);
		handleExcep(excep);
	}
}

TBB_GInts2DMultiIntMatrix::TBB_GInts2DMultiIntMatrix(const GIntJobInfor& ginfor0, 
		const SigMolShellPairInfor& shellPairInfor, const MolShell& rowShell, 
		const MolShell& colShell, const Double* angCenter, 
		MtrxVec& result):rs(rowShell),cs(colShell),infor(shellPairInfor),ginfor(ginfor0),intMtrxVec(result)
{ 
	// for this constructor, the center of C will be copied 
	// from the angCenter
	C[0] = angCenter[0];
	C[1] = angCenter[1];
	C[2] = angCenter[2];

	// finally let's check the job - whether it matches here
	UInt job = ginfor.getIntJob();
	if (job != MOM_P) {
		string infor = "The integral job can not be processed here";
		Excep excep("TBB_GInts2DMultiIntMatrix","constructor",EXCEPTION_GINTS_INVALID_INTEGRAL_JOB,infor);
		handleExcep(excep);
	}
}

void TBB_GInts2DIntMatrix::operator()(const blocked_range<UInt>& r) const 
{
	///////////////////////////////////////////////////////////////////////
	// setup scratch data for calculation                                //
	// 1  atom shell pair                                                //
	// 2  intermediate result based on Cart. shell dimension             //
	// 3  intermediate result based on Cart. atom shell dimension        //
	// 4  intermediate result based on norm. atom shell dimension        //
	// 5  C2P transformation plus possible basis set normalization       //
	///////////////////////////////////////////////////////////////////////

	// shell pair data
	UInt maxNP2 = infor.getMaxNP2();
	UInt maxNL  = infor.getMaxNL();
	UInt maxNSP = infor.getMaxNSP();
	AtomShellPair sp(maxNSP,maxNP2,maxNL);

	// block matrix which is used to store the data on shell dimension
	BlockMtrx shellResultMtrx(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());

	// block matrix which is used to store the data on cartesian atom shell dimension
	BlockMtrx cartAtomShellMtrx(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());

	// block matrix which is used to store the data on cartesian atom shell dimension
	BlockMtrx normAtomShellMtrx(infor.getMaxNBasisForAtomShell(),infor.getMaxNBasisForAtomShell());

	// C2P transformation and basis set scaling 
	// we perform C2P transformation (cartesian to pure)
	// as well as normal basis set scaling work for each shell pair block
	UInt maxL = infor.getMaxL();
	NormAtomShellData scale(maxL,infor.getMaxNBasisForAtomShell(),DO_SCALE);
	CPTransAtomShell trans(maxL,C2P_WITH_L00,WITH_MATRIX_ITSELF);

	// possible scratch matrix used in CP transformation
	Mtrx tmp(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());

	///////////////////////////////////////////////////////////////////////
	// real working code - loop over atom shell pair data                //
	///////////////////////////////////////////////////////////////////////
	for( UInt n=r.begin(); n!=r.end(); ++n ) {

		//
		//get shell pair infor
		//
		const SigAtomShellPairInfor& sigAtomShellPairInfor = infor.getSigAtomShellPairInfor(n); 

		//
		//form atom shell pair
		//
		sp.init(rs,cs,sigAtomShellPairInfor,ginfor.pickUpSPThresh());

		//
		//form atom shell pair result matrix
		//
		UInt nRowAtomBasis     = 0;
		UInt nColAtomBasis     = 0;
		UInt nRowAtomBasOffset = 0;
		UInt nColAtomBasOffset = 0;
		sp.getNCarBas(nRowAtomBasis,nColAtomBasis);
		sp.getGlobalCarBasOffSet(nRowAtomBasOffset,nColAtomBasOffset);
		cartAtomShellMtrx.init(nRowAtomBasis,nColAtomBasis,nRowAtomBasOffset,nColAtomBasOffset);

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
			const Double* diff= spData.getExpDiff();

			//
			// now set up integral results
			//
			UInt nRowBasisSP, nColBasisSP;
			spData.getNCarBas(nRowBasisSP,nColBasisSP);
			UInt nRowBasisOffsetSP, nColBasisOffsetSP;
			spData.getLocalCarBasOffSet(nRowBasisOffsetSP,nColBasisOffsetSP);
			shellResultMtrx.init(nRowBasisSP,nColBasisSP,nRowBasisOffsetSP,nColBasisOffsetSP);

			// do we need to perform updation with transpose of shellResultMtrx?
			bool inTrans = spData.inverseShells();

			// get the job type, and get raw integral
			UInt jobType = ginfor.getIntJob();
			switch(jobType) {
				case TWO_BODY_OVERLAP:
					hgp_os_twobodyoverlap(code,np2,c2,e2,fac,P,A,B,shellResultMtrx.getPtr());
					shellResultMtrx.updateData(cartAtomShellMtrx,inTrans);
					break;
				case KINETIC:
					hgp_os_kinetic(code,np2,c2,e2,diff,fac,P,A,B,shellResultMtrx.getPtr());
					shellResultMtrx.updateData(cartAtomShellMtrx,inTrans);
					break;
				case NUCLEAR_ATTRACTION:
					hgp_os_nai(code,np2,Z.size(),c2,e2,fac,P,A,B,&coord.front(),&Z.front(),
							shellResultMtrx.getPtr());
					shellResultMtrx.updateData(cartAtomShellMtrx,inTrans);
					break;
			}

			//
			// debug the original integral matrix
			//
			//if (nRowAtomBasOffset == 9 && nColAtomBasOffset == 0) {
			//	if (nRowBasisOffsetSP == 3 && nColBasisOffsetSP == 0) {
			//		spData.print();
			//		shellResultMtrx.blockPrint("initial integral block");
			//	}
			//	cartAtomShellMtrx.blockPrint("initial integral atom block");
			//}
		}

		//
		// now let's process the result matrix on atom shell pair
		// if both rs and cs are in Cartesian format,
		// then the cartAtomShellMtrx will just reflect the final data 
		// block in the result. In this case, we do not need normAtomShellMtrx
		//
		if (!rs.allCart() || !cs.allCart()) {

			UInt nRowAtomBasis     = 0;
			UInt nColAtomBasis     = 0;
			UInt nRowAtomBasOffset = 0;
			UInt nColAtomBasOffset = 0;
			sp.getNBas(nRowAtomBasis,nColAtomBasis);
			sp.getGlobalBasOffSet(nRowAtomBasOffset,nColAtomBasOffset);
			normAtomShellMtrx.init(nRowAtomBasis,nColAtomBasis,nRowAtomBasOffset,nColAtomBasOffset);

			// get the atom shell data
			UInt rowAtomShellIndex = sigAtomShellPairInfor.getRowAtomShellIndex();
			UInt colAtomShellIndex = sigAtomShellPairInfor.getColAtomShellIndex();
			const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
			const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);

			// we will perform job based on atom shell pair
			trans.cpTransformOnAtomShellPair(rowAtomShell,colAtomShell, 
					cartAtomShellMtrx, tmp, normAtomShellMtrx); 

			// it's possible that we also need scale it
			// this is true if the atom basis is mixed with Cart. basis set 
			// and spherical basis sets
			scale.normRowData(rowAtomShell,normAtomShellMtrx); 
			scale.normColData(colAtomShell,normAtomShellMtrx); 

			// now update the result
			// because the matrix order is correct, 
			// we do it without transpose
			bool withTranspose = false;
			normAtomShellMtrx.updateData(intMtrx,withTranspose);
		}else{

			// in this case let's do scale work
			UInt rowAtomShellIndex = sigAtomShellPairInfor.getRowAtomShellIndex();
			UInt colAtomShellIndex = sigAtomShellPairInfor.getColAtomShellIndex();
			const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
			const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);
			scale.normRowData(rowAtomShell,cartAtomShellMtrx); 
			scale.normColData(colAtomShell,cartAtomShellMtrx); 

			// now update the result
			// because the matrix order is correct, 
			// we do it without transpose
			bool withTranspose = false;
			cartAtomShellMtrx.updateData(intMtrx,withTranspose);
		}
	}
}

void TBB_GInts2DMultiIntMatrix::operator()(const blocked_range<UInt>& r) const 
{
	///////////////////////////////////////////////////////////////////////
	// setup scratch data for calculation                                //
	// 1  atom shell pair                                                //
	// 2  intermediate result based on Cart. shell dimension             //
	// 3  intermediate result based on norm. shell dimension             //
	// 4  C2P transformation plus possible basis set normalization       //
	// 5  temp array to hold the integral results on Cartesian dimension //
	///////////////////////////////////////////////////////////////////////

	// shell pair data
	UInt maxNP2 = infor.getMaxNP2();
	UInt maxNL  = infor.getMaxNL();
	UInt maxNSP = infor.getMaxNSP();
	AtomShellPair sp(maxNSP,maxNP2,maxNL);

	// integral result array, here we need to consider the multiple matrix dimension
	DoubleVec intResultArray(infor.getMaxNBasisForShell()*infor.getMaxNBasisForShell()*intMtrxVec.size());

	// block matrix which is used to store the data on cartesian atom shell dimension
	BlockMtrx cartShellMtrx(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());

	// block matrix which is used to store the data on cartesian atom shell dimension
	BlockMtrx normShellMtrx(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());

	// C2P transformation and basis set scaling data
	// we perform C2P transformation (cartesian to pure)
	// as well as normal basis set scaling work for each shell pair block
	UInt maxL = infor.getMaxL();
	CPTransData cpWork(maxL,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF);

	// possible scratch matrix used in CP transformation
	Mtrx tmp(infor.getMaxNBasisForShell(),infor.getMaxNBasisForShell());

	///////////////////////////////////////////////////////////////////////
	// real working code - loop over atom shell pair data                //
	///////////////////////////////////////////////////////////////////////
	for( UInt n=r.begin(); n!=r.end(); ++n ) {

		//
		//get shell pair infor
		//
		const SigAtomShellPairInfor& sigAtomShellPairInfor = infor.getSigAtomShellPairInfor(n); 

		//
		//form atom shell pair
		//
		sp.init(rs,cs,sigAtomShellPairInfor,ginfor.pickUpSPThresh());

		// get the atom shell data
		UInt rowAtomShellIndex = sigAtomShellPairInfor.getRowAtomShellIndex();
		UInt colAtomShellIndex = sigAtomShellPairInfor.getColAtomShellIndex();
		const AtomShell& rowAtomShell = rs.getAtomShell(rowAtomShellIndex);
		const AtomShell& colAtomShell = cs.getAtomShell(colAtomShellIndex);

		//
		//generate integral result for each shell pair
		//
		UInt nSP = sp.getNShellPairs(); 
		for(UInt iSP=0; iSP<nSP; iSP++) {

			// whether the real shell data is significant?
			const ShellPair& spData = sp.getShellPair(iSP);
			if (! spData.isSig()) {
				continue;
			}

			// for the angular moment code
			UInt iLmin = -1;
			UInt iLmax = -1;
			UInt jLmin = -1;
			UInt jLmax = -1;
			spData.getL(iLmin,iLmax,jLmin,jLmax);
			UInt kLmin = -1;
			UInt kLmax = -1;
			UInt jobType = ginfor.getIntJob();
			if (jobType == MOM_P) {
				kLmin = 1;
				kLmax = 1;
			}
			LInt code = codeL(iLmin,iLmax,jLmin,jLmax,kLmin,kLmax);

			// now let's do integral job - get primitive pairs data
			UInt np2          = spData.getNP2();
			const Double* c2  = spData.getC2();
			const Double* e2  = spData.getE2();
			const Double* A   = spData.getA();
			const Double* B   = spData.getB();
			const Double* P   = spData.getP();
			const Double* fac = spData.getFac();

			// clear the result array
			UInt nRowCarBas,nColCarBas;
			spData.getNCarBas(nRowCarBas,nColCarBas);
			vset(&intResultArray[0],ZERO,nRowCarBas*nColCarBas*intMtrxVec.size());

			// get the job type, and get raw integral
			switch(jobType) {
				case MOM_P:
					hgp_os_mom_p(code,np2,c2,e2,fac,P,A,B,C,&intResultArray.front());
					break;
			}

			// whether the shell pair data is in inverse order?
			bool inTrans = spData.inverseShells();

			//
			// now set up Cartesian integral matrix block
			// the matrix block will be in same status as it's mother matrix
			// this will make thing goes easier
			//
			UInt nRowCarBasOffsetSP, nColCarBasOffsetSP;
			spData.getGlobalCarBasOffSet(nRowCarBasOffsetSP,nColCarBasOffsetSP);
			if (inTrans) {
				UInt tmp    = nRowCarBas;
				nRowCarBas  = nColCarBas;
				nColCarBas  = tmp;
				tmp                 = nRowCarBasOffsetSP;
				nRowCarBasOffsetSP  = nColCarBasOffsetSP;
				nColCarBasOffsetSP  = tmp;
			}
			cartShellMtrx.init(nRowCarBas,nColCarBas,nRowCarBasOffsetSP,nColCarBasOffsetSP);

			//
			// now set up normal integral matrix block
			// the matrix block will be in same status as it's mother matrix
			// this will make thing goes easier
			//
			UInt nRowBas,nColBas;
			spData.getNBas(nRowBas,nColBas);
			UInt nRowBasOffsetSP, nColBasOffsetSP;
			spData.getGlobalBasOffSet(nRowBasOffsetSP,nColBasOffsetSP);
			if (inTrans) {
				UInt tmp = nRowBas;
				nRowBas  = nColBas;
				nColBas  = tmp;
				tmp              = nRowBasOffsetSP;
				nRowBasOffsetSP  = nColBasOffsetSP;
				nColBasOffsetSP  = tmp;
			}
			normShellMtrx.init(nRowBas,nColBas,nRowBasOffsetSP,nColBasOffsetSP);

			// get the shell data
			UInt rowShellIndex, colShellIndex;
			if (inTrans) {
				spData.getLocalShellIndex(colShellIndex,rowShellIndex);
			}else{
				spData.getLocalShellIndex(rowShellIndex,colShellIndex);
			}
			const Shell& rowShell = rowAtomShell.getShell(rowShellIndex);
			const Shell& colShell = colAtomShell.getShell(colShellIndex);

			// now let's digest the data
			for(UInt iMtrx = 0; iMtrx<intMtrxVec.size(); iMtrx++) {

				//
				// copy the result integral results to the shell block part
				// we will consider the transpose status
				//
				const Double* rawIntegralResult = &intResultArray[nRowCarBas*nColCarBas*iMtrx];
				Double* cartShellMtrxPtr = cartShellMtrx.getPtr();
				if (inTrans) {
					for(UInt iCol=0; iCol<nColCarBas; iCol++) {
						for(UInt iRow=0; iRow<nRowCarBas; iRow++) {
							cartShellMtrxPtr[iRow+iCol*nRowCarBas] = rawIntegralResult[iCol+iRow*nColCarBas];
						}
					}
				}else{
					for(UInt i=0; i<nRowCarBas*nColCarBas; i++) {
						cartShellMtrxPtr[i] = rawIntegralResult[i];
					}
				}

				//
				// now let's do transform work
				//
				if (rowShell.isPure()) {

					// let's multiple the transform matrix
					// we note here the Lmax must be >=2
					// also lmin should equal lmax
					if (rowShell.getLmax() != rowShell.getLmin()) {
						string infor = "in the C2P transformation process the shell angular momentum is not supported";
						Excep excep("TBB_GInts2DMultiIntMatrix","operator()",
								EXCEPTION_UNSUPPORTED_ANGULAR_MOMENTUM,infor);
						handleExcep(excep);
					}

					// c2p transformation result = C^{T}*M
					const Mtrx& c2p = cpWork.getConvertMatrix(rowShell.getLmax());
					tmp.reset(nRowBas,nColCarBas,true);
					tmp.mult(c2p,cartShellMtrx,true,false,ONE,ZERO);

					// now let's see the col shell
					if (colShell.isPure()) {

						// let's multiple the transform matrix
						// we note here the Lmax must be >=2
						// also lmin should equal lmax
						if (colShell.getLmax() != colShell.getLmin()) {
							string infor = "in the C2P transformation process the shell angular momentum is not supported";
							Excep excep("TBB_GInts2DMultiIntMatrix","operator()",
									EXCEPTION_UNSUPPORTED_ANGULAR_MOMENTUM,infor);
							handleExcep(excep);
						}

						// result = M*C
						const Mtrx& newC2p = cpWork.getConvertMatrix(colShell.getLmax());
						normShellMtrx.mult(tmp,newC2p,false,false,ONE,ZERO);
					}else{

						// do we need to scale the col shell?
						if (colShell.getLmax() >= 2) {
							const DoubleVec& convert = cpWork.getConvertVec(colShell.getLmax());
							for(UInt i=0; i<tmp.getCol(); i++) {
								vscal(tmp.getPtr(0,i),convert[i],tmp.getRow());
							}
						}

						// now let's copy the data to the norm shell mtrx 
						normShellMtrx.copyMatrix(tmp);
					}
				}else{

					// do we need to scale the row shell?
					if (rowShell.getLmax() >= 2) {
						const DoubleVec& convert = cpWork.getConvertVec(rowShell.getLmax());
						for(UInt i=0; i<cartShellMtrx.getCol(); i++) {
							vmul(&convert[0],cartShellMtrx.getPtr(0,i),cartShellMtrx.getPtr(0,i),cartShellMtrx.getRow());
						}
					}

					// now let's see the col shell
					if (colShell.isPure()) {

						// let's multiple the transform matrix
						// we note here the Lmax must be >=2
						// also lmin should equal lmax
						if (colShell.getLmax() != colShell.getLmin()) {
							string infor = "in the C2P transformation process the shell angular momentum is not supported";
							Excep excep("TBB_GInts2DMultiIntMatrix","operator()",
									EXCEPTION_UNSUPPORTED_ANGULAR_MOMENTUM,infor);
							handleExcep(excep);
						}

						// result = M*C
						const Mtrx& c2p = cpWork.getConvertMatrix(colShell.getLmax());
						normShellMtrx.mult(cartShellMtrx,c2p,false,false,ONE,ZERO);
					}else{
						// do we need to scale the col shell?
						if (colShell.getLmax() >= 2) {
							const DoubleVec& convert = cpWork.getConvertVec(colShell.getLmax());
							for(UInt i=0; i<cartShellMtrx.getCol(); i++) {
								vscal(cartShellMtrx.getPtr(0,i),convert[i],cartShellMtrx.getRow());
							}
						}

						// now let's copy the data to the norm shell mtrx 
						// remember we can not copy the index, so we use copyMatrix
						// to only copy it's content
						normShellMtrx.copyMatrix(cartShellMtrx);
					}
				}

				// now let's update the result
				Mtrx& result = intMtrxVec[iMtrx];
				normShellMtrx.updateData(result,false);
			}
		}
	}
}

void GInts2D::doMtrx(const MolShell& rs, const MolShell& cs, const Molecule& mol, 
					Mtrx& intMtrx, bool printTiming) const
{
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
	tick_count t0 = tick_count::now();

	// do normal integral matrix calculation
	if (ginfor.getJobOrder() == 0) {
		TBB_GInts2DIntMatrix  intWork(ginfor,infor,rs,cs,mol,intMtrx); 
		UInt len = infor.getNSigAtomShellPairs();
		parallel_for(blocked_range<UInt>(0,len), intWork);
	}

	// possible timing code
	tick_count t1 = tick_count::now();
	Double t = (t1-t0).seconds();
	if (printTiming) {
		printf("%s  %-12.6f\n", "GInts2D work with TBB threads, time in seconds  ", t);
	}
}

void GInts2D::doMultiMtrx(const MolShell& rs, const MolShell& cs,
					MtrxVec& intMtrxList, bool printTiming) const
{
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
	tick_count t0 = tick_count::now();

	// do normal integral matrix calculation
	if (ginfor.getJobOrder() == 0) {
		TBB_GInts2DMultiIntMatrix  intWork(ginfor,infor,rs,cs,intMtrxList); 
		UInt len = infor.getNSigAtomShellPairs();
		parallel_for(blocked_range<UInt>(0,len), intWork);
	}

	// possible timing code
	tick_count t1 = tick_count::now();
	Double t = (t1-t0).seconds();
	if (printTiming) {
		printf("%s  %-12.6f\n", "GInts2D work for multiple result matrices with TBB threads, time in seconds  ", t);
	}
}
