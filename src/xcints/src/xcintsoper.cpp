/*
 * \file    xcints.cpp
 * \brief   top interface class for doing DFT work
 * \author  Fenglai Liu and Jing Kong
 */
#include<cstdio>
#include <iostream>
#include <boost/filesystem.hpp>
#include "tbb/tbb.h"
#include "tbb/compat/thread"
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

using namespace std;
using namespace std::this_thread;
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
tbb::spin_mutex GLOBAL_XCINTS_MUTEX1;

void XCIntsOper::operator()(const blocked_range<UInt>& r) const
{
	// real working loop
	//bibd should be declared here.
	for( UInt iBatch=r.begin(); iBatch!=r.end(); ++iBatch ) {

		// get the atomic number and the batch grids
		// ibatch is the local batch index
		// iAtom is the atom index to fetch geometry data
		// Z is the corresponding atomic number
		UInt ibatch,Z,iAtom;
		xcItInp.molGrids.getBatchInfor(iBatch,iAtom,Z,ibatch);

		// get the atom grid and shell data
		const AtomShell& as = xcInp.ms.getAtomShell(iAtom);
		const AtomGrids& ag = xcItInp.molGrids.getAtomGrid(Z);

		// batch grid 
		//Change here, 'true' for one point.
		BatchGrid bg(ag,as,ibatch, false);

		// sig atom basis
		SigAtomBasis sigList(bg,xcItInp.infor,xcInp.ms,xcItInp.molShellSize);
		if (sigList.empty()) continue;
		if (xcItInp.infor.doDebug("sigatombasis")) {
			sigList.print();
		}

		// do partition weights?
		if (xcItInp.infor.doParWeights()) {
			bg.BeckeWeights0(xcItInp.infor,xcInp.ms,sigList);
		}
		if (xcItInp.infor.doDebug("batchgrid")) {
			bg.print();
		}

		// batch basis set
		BatchBasis bbasis(xcItInp.xcvar,bg,xcInp.ms,sigList);
		bbasis.setupPhi(bg,xcInp.ms,sigList);
		if (xcItInp.infor.doDebug("batchbasis")) {
			bbasis.print();
		}

		// form half jk rho
		HalfJKRho halfJKRho(xcItInp.infor,xcItInp.xcvar,bg,xcInp.ms,sigList,
		                    xcItInp.cartDen);
		halfJKRho.doHalfJKRho(xcItInp.infor,xcItInp.xcvar,xcInp.ms,bg,sigList,bbasis,
		                      xcInp.denMtrx,xcItInp.cartDen,xcItInp.spData);

		// if this is not a XC job, we do not go over the functional and 
		// variable calculation stuff; so we do it here

		// batch var
		BatchVar bVar(sigList,bbasis,xcItInp.xcvar,xcInp.denMtrx);
		if (xcItInp.xcvar.hasExRho()) {
			bVar.buildExRho(xcInp.ms,sigList,xcInp.denMtrx,bbasis,halfJKRho,xcItInp.xcvar);
 		}
		if (xcItInp.infor.doDebug("batchvar")) {
			bVar.print(xcItInp.xcvar);
		}

		// hirshfeld weights possibly used
		Hirshfeld hirshfeld(xcItInp.infor,sigList,bg);
		if (xcInp.xcfunc.useHirshfeldWeights()) {
			hirshfeld.formWtsFreeAtomVol(xcInp.ms,sigList,bbasis,bg,xcItInp.atomDenMtrx);
		}

		{
			//Reform. This line is not necessary if bg, sigList, etc are data
			//members of BatchIterBasicData, which in turn is included in BatchIter.
			BatchIterBasicData* bibd = new BatchIterBasicData(bg,sigList,bbasis,bVar);

			//This will be part of something, either XCInts (which then needs an op()
			//for tbb), or  batch iteractor class.

			//The following line will become BatchProps bProps(...). The make* calls
			//will be inside the constructor.
			//bHoleAvg in XCHoleAvg should be made as unique_ptr so that it does not
			//have to be a member of XCHoleAvg and I dont have to worry about its
			//deletion. It can be passed to accumulate. 
			//cout << "tbbxcmtrx 1  thread::id = " << get_id() << endl;
			BatchProp* bProp = xcProp.makeBatchProp(xcInp, xcItInp, *bibd);
			bProp->integrate();  //This line should be go after doFuncDeriv call.
			{
				tbb::spin_mutex::scoped_lock lock(GLOBAL_XCINTS_MUTEX1);
				xcProp.accumulate(bProp);
			}	
			delete bibd;  //This line will eventually become unecessary.
			//This line will go to the destructor of BatchProps or unecessary for unique_ptr.
			delete bProp;  
			continue;
		}
	}
}

