/**
 * \file    xcints.cpp
 * \brief   top interface class for doing DFT work
 * \author  Fenglai Liu and Jing Kong
 */
#include<cstdio>
#include <boost/filesystem.hpp>
#include "tbb/tbb.h"
#include "globalinfor.h"
#include "atomdenmtrx.h"
#include "gintsinfor.h"
#include "excep.h"
#include "shell.h"
#include "matrix.h"
#include "filerw.h"
#include "denmtrx.h"
#include "atomgrids.h"
#include "batchgrid.h"
#include "halfjkrho.h"
#include "cptrans.h"
#include "sigatombasis.h"
#include "sigshellpairinfor.h"
#include "dftmatrix.h"
#include "batchbasis.h"
#include "batchvar.h"
#include "batchfunc.h"
#include "batchxcmtrx.h"
#include "batchjkmtrx.h"
#include "hirshfeld.h"
#include "oddelec.h"
#include "xcenergyinfor.h"
#include "xcints.h"
#include "xcintsprops.h"

using namespace boost::filesystem;
using namespace globalinfor;
using namespace atomdenmtrx;
using namespace gintsinfor;
using namespace excep;
using namespace shell;
using namespace matrix;
using namespace filerw;
using namespace denmtrx;
using namespace atomgrids;
using namespace batchgrid;
using namespace halfjkrho;
using namespace cptrans;
using namespace sigatombasis;
using namespace sigshellpairinfor;
using namespace dftmatrix;
using namespace batchbasis;
using namespace batchvar;
using namespace batchfunc;
using namespace batchxcmtrx;
using namespace batchjkmtrx;
using namespace hirshfeld;
using namespace oddelec;
using namespace xcenergyinfor;
using namespace xcints;
using namespace xcintsprops;

//
// define a global mutex just for xcints part
//
tbb::spin_mutex GLOBAL_XCINTS_MUTEX;

void TBB_XCIntsMatrix::operator()(const blocked_range<UInt>& r) const
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
		if (infor.doDebug("sigatombasis")) {
			sigList.print();
		}

		// do partition weights?
		if (infor.doParWeights()) {
			bg.BeckeWeights0(infor,ms,sigList);
		}
		if (infor.doDebug("batchgrid")) {
			bg.print();
		}

		// batch basis set
		BatchBasis bbasis(xcvar,bg,ms,sigList);
		bbasis.setupPhi(bg,ms,sigList);
		if (infor.doDebug("batchbasis")) {
			bbasis.print();
		}

		// form half jk rho
		HalfJKRho halfJKRho(infor,xcvar,bg,ms,sigList,cartDen);
		halfJKRho.doHalfJKRho(infor,xcvar,ms,bg,sigList,bbasis,denMtrx,cartDen,spData); 

		// if this is not a XC job, we do not go over the functional and 
		// variable calculation stuff; so we do it here
		UInt intJob = infor.getIntJob();
		if (! isDFTJob(intJob)) {

			// based on the batch basis and half JK rho, we just proceed
			// to form the JK matrix
			BatchJKMtrx batchJKMtrx(ms,sigList,halfJKRho,bbasis,bg,infor);

			// now it's the critical section
			// to update the global results
			{
				tbb::spin_mutex::scoped_lock lock(GLOBAL_XCINTS_MUTEX);
				batchJKMtrx.updateJKMtrx(infor,sigList,xcMtrx,jMtrx);
				if (doJ(intJob)) {
					halfJKRho.updateHalfCouRhoVecMatrix(s2VecJMtrx);
				}
			}

			// now go to another batch
			continue;
		}

		// batch var
		BatchVar bVar(sigList,bbasis,xcvar,denMtrx);
		if (xcvar.hasExRho()) {
			bVar.buildExRho(ms,sigList,denMtrx,bbasis,halfJKRho,xcvar); 
		}
		if (infor.doDebug("exrho")) {

			// let's calculate the alpha and beta exchange energy
			Double alphaEx0 = bVar.getEX(0,xcvar,bg);
			Double betaEx0  = bVar.getEX(1,xcvar,bg);

			// now let's put a mutex on and write the result into global one
			{
				tbb::spin_mutex::scoped_lock lock(GLOBAL_XCINTS_MUTEX);
				alphaEx += alphaEx0;
				betaEx  += betaEx0;
			}

			// get rid of the rest of code
			continue;
		}
		if (infor.doDebug("batchvar")) {
			bVar.print(xcvar);
		}

		// hirshfeld weights possibly used
		Hirshfeld hirshfeld(infor,sigList,bg);
		if (xcfunc.useHirshfeldWeights()) {
			hirshfeld.formWtsFreeAtomVol(ms,sigList,bbasis,bg,atomDenMtrx);
		}

		// if we do XC energy profiling,
		// we need to initialize it before pass
		// it to batch func
		XCEnergyInfor xcEProfile0(infor,xcfunc);
		if (xcEProfile0.doXCEnergyProfiling()) {
			xcEProfile0.init(bg);
		}

		// batch functional
		BatchFunc bFunc(bVar,denMtrx,xcfunc,xcvar,infor);
		if (xcfunc.useHirshfeldWeights()) {
			bFunc.formBR89EXHoleWeights(ms,sigList,hirshfeld); 
		}
		bFunc.doFuncDeriv(bVar,xcfunc,xcvar,xcEProfile0); 
		if (infor.doDebug("batchfunc")) {
			bFunc.print();
		}

		// also collect the batch energy result against
		// the weights
		if (xcEProfile0.doXCEnergyProfiling()) {
			xcEProfile0.formBatchXCEnergy(bg);
		}

		// possible odd electron calculation
		ODDElec oddElec0(ms,infor);
		if (infor.doOddElec()) {
			oddElec0.doOddElecPopulation(infor,bg,bVar,xcvar,denMtrx); 
		}

		// update the enegy
		Double batchExc = ZERO;
		bFunc.getBatchExc(bg,batchExc);

		// batch xc matrix and update the result matrix
		BatchXCMtrx Fxc(sigList,bbasis,bFunc,bg,xcvar,infor);
		if (infor.doDebug("batchxcmtrx")) {
			Fxc.print();
		}

		// let's calculate the exrho part contribution
		BatchEXRhoMtrx exRhoFxc(ms,sigList,bbasis,halfJKRho,bFunc,bg,xcvar,infor);

		// now it's the critical section
		// to update the global results
		{
			tbb::spin_mutex::scoped_lock lock(GLOBAL_XCINTS_MUTEX);
			exc += batchExc;
			Fxc.updateXCMtrx(sigList,xcMtrx);
			if (xcvar.hasExRho()) {
				exRhoFxc.updateXCMtrx(sigList,xcMtrx);
			}
			if (infor.doOddElec()) {
				oddElec.add(infor,oddElec0);
			}
			if (xcEProfile0.doXCEnergyProfiling()) {
				xcEnergyProfile.add(xcEProfile0);
			}
		}
	}
}

XCInts::XCInts(const GlobalInfor& gInfor, const XCFunc& xcfunc0, const MolShell& ms, 
		const XCIntsInfor& xcinfor, const GIntsInfor& gintsInfor, const UInt& job, 
		const UInt& jobOrder):exc(ZERO),infor(gInfor,gintsInfor,xcinfor,job,jobOrder),
	xcfunc(xcfunc0),xcvar(xcfunc0,infor.getNSpin(),job,jobOrder),
	size(ms,infor.getThresh()),molGrids(ms,size,infor),xcMtrx(infor.getNSpin())
{
	// now let's set up the result 
	if (jobOrder == 0) {

		// if this is DFT job, then the result matrix is xc fock matrix
		if (isDFTJob(job)) {
			xcMtrx.init(ms.getNBas(),ms.getNBas());
		}else{
			// if we have the numerical K matrix, xcMtrx will be used to store K matrix
			if (doK(job)) {
				xcMtrx.init(ms.getNBas(),ms.getNBas());
			}

			// if only J calculation is required, then we will initialize jMtrx
			if (onlyDoJ(job)) {
				jMtrx.init(ms.getNBas(),ms.getNBas());
			}
		}
	}else{
		string infor = "currently XCInts does not support gradient/Hessian calculation";
		Excep excep("XCInts","Constructor",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
		handleExcep(excep);
	}

	// debug code
	if (infor.doDebug("xcvar")) {
		xcvar.print();
	}
	if (infor.doDebug("molegrids")) {
		molGrids.print(xcinfor.debugLevelMoleGrids());
	}
}

void XCInts::doXCMtrx(const MolShell& ms, const DenMtrx& denMtrx, SpinMatrix& result, 
		bool printTiming, bool printMatrix)
{
	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	if (infor.useMultiThreads()) {
		init.initialize(infor.getNCPUThreads());
	}else{
		init.initialize(1);
	}

	// some functionals replies on atom density matrix data
	// so here let's do it
	// we assume that it should be already existing in the hard disk
	// so only we only reload it
	const GlobalInfor& gInfor = infor.getGlobalInfor(); 
	AtomDenMtrx atomDenMtrx(gInfor,denMtrx.getSec());
	if (xcfunc.useAtomDenMtrxData()) {
		atomDenMtrx.recover(ms);
	}

	// if the odd electron population is used, let's set it up here
	ODDElec oddElec(ms,infor);

	// if we do energy profiling, we need to set up the class
	XCEnergyInfor xcEnergyProfile(infor,xcfunc);

	// set up the shell pair data
	// currently for XC job only the exrho variable will 
	// trigger the gints job and needs the shell pair data
	SigMolShellPairInfor spData;
	if (xcvar.hasExRho()) {
		GIntJobInfor jobInfor(infor.getGIntsInfor(),EX_DENSITY);
		spData.init(ms,ms,jobInfor);
	}

	// make an empty Cartesian Density, this is used for Coulomb
	// will never be used here
	DoubleVec cartDen;

	// make an empty s2 form J matrix
	// it is also not used here
	DoubleVec s2VecJMtrx;

	// for debugging the exrho
	Double alphaEx = ZERO;
	Double betaEx  = ZERO;

	// possible timing code
	tick_count t0 = tick_count::now();
	// Change here to turn on hole calculation.
  if ( false )
  {
		XCInput xcInp(xcfunc,ms,denMtrx);
		XCIterInput xcItInp(spData,molGrids,size,xcvar,infor,cartDen,atomDenMtrx);
		//Change here 2. Choose a type of hole to do and the number of s values.  
		//HFXHole xcHoleAvg(denMtrx.getNSpin(), 100);  
		MetaGGAHole xcHoleAvg(denMtrx.getNSpin(), 25); 
		XCIntsOper xcOper(xcInp,xcItInp, xcHoleAvg);
		UInt len = molGrids.getNTotalBatches();
		//Change here 3. Uncomment for one-point. 
		//len = 1; 
		parallel_for(blocked_range<UInt>(0,len), xcOper);
		xcHoleAvg.dump();
		exit(0);
  }
	// do normal integral matrix calculation
	TBB_XCIntsMatrix  tbb_xcints(ms,spData,molGrids,size,xcfunc,xcvar,infor,denMtrx,cartDen,
			atomDenMtrx,exc,alphaEx,betaEx,xcMtrx,jMtrx,s2VecJMtrx,oddElec,xcEnergyProfile); 
	UInt len = molGrids.getNTotalBatches();
	parallel_for(blocked_range<UInt>(0,len), tbb_xcints);

	// do we do exrho debugging
	if (infor.doDebug("exrho")) {
		printf("in xcints, alpha exchange energy: %-20.14f\n", alphaEx);
		printf("in xcints, beta  exchange energy: %-20.14f\n", betaEx);
		return;
	}

	// update result
	// here we need to do add transpose to symmetrize the final result
	for(UInt iSpin=0; iSpin<result.getNSpin(); iSpin++) {
		Mtrx& M  = result.getMtrx(iSpin);
		Mtrx& M0 = xcMtrx.getMtrx(iSpin);
		M0.addTranspose();
		M.add(M0);
		if (printMatrix) {
			string title = "alpha XC Matrix";
			if (iSpin == 1) {
				title = "beta XC Matrix";
			}
			M0.print(title);
		}
	}

	// if we do odd electron population, update the results
	if (infor.doOddElec()) {
		oddElec.printResults(ms);
	}

	// if we do xc energy profile, we can print it here
	if (xcEnergyProfile.doXCEnergyProfiling()) {
		xcEnergyProfile.print();
	}

	// possible timing code
	tick_count t1 = tick_count::now();
	Double t = (t1-t0).seconds();
	if (printTiming) {
		printf("%s  %-12.6f\n", "XCInts work with threads, time in seconds  ", t);
	}
}

void XCInts::doJKMtrx(const MolShell& ms, DenMtrx& denMtrx, SpinMatrix& result, 
		bool printTiming, bool printMatrix)
{
	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	if (infor.useMultiThreads()) {
		init.initialize(infor.getNCPUThreads());
	}else{
		init.initialize(1);
	}

	// we just build an empty one
	const GlobalInfor& gInfor = infor.getGlobalInfor(); 
	AtomDenMtrx atomDenMtrx(gInfor,denMtrx.getSec());

	// build an empty one
	ODDElec oddElec(ms,infor);

	// build an empty xcenergyinfor
	XCEnergyInfor xcEnergyProfile(infor,xcfunc);

	// cartDen and s2VecJMtrx are both s2 form vector
	// if the Coulomb job is performed, we need to set them up
	DoubleVec cartDen;
	DoubleVec s2VecJMtrx;

	// build an shell pair infor
	UInt intJob = infor.getIntJob();
	GIntJobInfor gintJobInfor(infor.getGIntsInfor(),intJob);
	SigMolShellPairInfor spInfor(ms,ms,gintJobInfor);  
	if (doJ(intJob)) {

		// get the length
		UInt type = TYPE_CART;
		UInt len  = spInfor.getSigSPVecLen(ms,ms,type);
		cartDen.assign(len,ZERO);
		s2VecJMtrx.assign(len,ZERO);

		// now let's convert the density matrix into s2 form
		// we will do alpha + beta form
		// this is to generate the cartDen
		UInt scaleWork  = DO_SCALE;
		UInt transWork  = C2P_WITH_L00;
		UInt mtrxStatus = WITH_MATRIX_TRANSPOSE;
		bool toVec      = true;
		for(UInt i=0; i<denMtrx.getNSpin(); i++) {
			Mtrx& M = denMtrx.getMtrx(i);
			spInfor.convertSigSPVec(toVec,gInfor,ms,ms,transWork,scaleWork,mtrxStatus,M,cartDen);
		}

		// do debug to see that whether the cartDen is correct
		/*
			bool doDebugCartDen = false;
			if (doDebugCartDen) {

		// set the matrix 
		Mtrx tmp(ms.getNCarBas(),ms.getNCarBas());
		tmp.reset(ms.getNBas(),ms.getNBas(), false);
		tmp.copyMatrix(denMtrx.getMtrx(0));

		// now do the transformation
		CPTransBasisNorm trans1(gInfor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW_COL);
		trans1.transform(ms,ms,tmp);

		// set up the new cart den
		DoubleVec newCartDen(len);
		spInfor.convertSigSPVec(toVec,gInfor,ms,ms,NO_TRANS,NO_SCALE,mtrxStatus,tmp,newCartDen);

		// now testing the two version of cart den
		for(UInt i=0; i<len; i++) {
		Double diff = fabs(cartDen[i]-newCartDen[i]);
		if (diff>1.0E-12) {
		printf("for element %d the diff is %-16.14f\n", Int(i), diff);
		}
		}
		}
		*/

		// make it into alpha + beta part
		if (infor.closeShell() && infor.getNSpin()==1) {
			vscal(&cartDen[0],TWO,len); 
		}
	}

	// set empty alpha EX etc.
	Double alphaEx = ZERO;
	Double betaEx  = ZERO;

	// possible timing code
	tick_count t0 = tick_count::now();

	// do normal integral matrix calculation
	TBB_XCIntsMatrix  tbb_xcints(ms,spInfor,molGrids,size,xcfunc,xcvar,infor,denMtrx,cartDen,
			atomDenMtrx,exc,alphaEx,betaEx,xcMtrx,jMtrx,s2VecJMtrx,oddElec,xcEnergyProfile); 
	UInt len = molGrids.getNTotalBatches();
	parallel_for(blocked_range<UInt>(0,len), tbb_xcints);

	// possible timing code
	tick_count t1 = tick_count::now();
	Double t = (t1-t0).seconds();
	if (printTiming) {
		printf("%s  %-12.6f\n", "XCInts work with threads, time in seconds  ", t);
	}

	// update result for K job
	if (doK(intJob)) {
		for(UInt iSpin=0; iSpin<result.getNSpin(); iSpin++) {
			Mtrx& M  = result.getMtrx(iSpin);
			Mtrx& M0 = xcMtrx.getMtrx(iSpin);
			M0.addTranspose();
			M.add(M0);
		}
	}

	// if only has J, we will use JMtrx
	if (onlyDoJ(intJob)) {
		jMtrx.addTranspose();
		for(UInt iSpin=0; iSpin<result.getNSpin(); iSpin++) {
			Mtrx& M  = result.getMtrx(iSpin);
			M.add(jMtrx);
		}
	}

	// additionally, if we do J matrix, we need to count in
	// the s2 vector J matrix
	if (doJ(intJob)) {
		UInt scaleWork  = DO_SCALE;
		UInt transWork  = C2P_WITH_L00;
		UInt mtrxStatus = WITH_MATRIX_ITSELF;
		bool toVec      = false;
		for(UInt iSpin=0; iSpin<result.getNSpin(); iSpin++) {
			Mtrx& M  = result.getMtrx(iSpin);
			spInfor.convertSigSPVec(toVec,gInfor,ms,ms,transWork,scaleWork,mtrxStatus,M,s2VecJMtrx);
		}
	}

	// finally let's see whether we do print matrix
	if (printMatrix) {
		for(UInt iSpin=0; iSpin<result.getNSpin(); iSpin++) {
			string title = "alpha JK Matrix derived from xcints";
			if (iSpin == 1) {
				title = "beta JK Matrix derived from xcints";
			}
			const Mtrx& M  = result.getMtrx(iSpin);
			M.print(title);
		}
	}
}

void XCInts::saveToDisk(const string& dataPath, const SpinMatrix& result) const
{
	// data dir location
	path p(dataPath.c_str());

	// test that whether we have the mo data for the path
	if (exists(p)) {
		remove_all(p);
	}
	create_directories(p);

	// now begin to read in
	for(UInt iSpin=0; iSpin<result.getNSpin(); iSpin++) {

		// let's get the name
		// this should be exactly same with the write function
		string name = "a.bin";
		if (iSpin == 1) {
			name = "b.bin";
		}

		// now read in data
		const Mtrx& Mi       = result.getMtrx(iSpin);
		const DoubleVec& vec = Mi.getVec();
		FileReadWrite fileRW(p.string(),name);
		fileRW.write(&vec.front(),vec.size());
	}
}

