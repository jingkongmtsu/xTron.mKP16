/*
 * \file    hfxhole.cpp
 * \brief   Compute the spherically averated HF exchange hole function.
 * \author  Yiting Wang and Jing Kong
 */
#include <iostream>
#include<cstdio>
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
//#include "xcenergyinfor.h"
#include "xcintsprops.h"
#include "xcints.h"

using namespace boost::filesystem;
//using namespace globalinfor;
//using namespace atomdenmtrx;
//using namespace gintsinfor;
using namespace std;
using namespace std::this_thread;
using namespace excep;
using namespace shell;
using namespace matrix;
using namespace filerw;
using namespace denmtrx;
//using namespace atomgrids;
using namespace batchgrid;
using namespace halfjkrho;
//using namespace cptrans;
using namespace sigatombasis;
using namespace sigshellpairinfor;
//using namespace dftmatrix;
using namespace batchbasis;
using namespace batchvar;
using namespace batchfunc;
//using namespace batchxcmtrx;
//using namespace batchjkmtrx;
//using namespace xcintsprops;
using namespace xcints;
using namespace xcintsprops;

//The output of this function is to a file, for now.
BatchProp* XCHoleAvg::makeBatchProp(const XCInput& xcInp, 
           const XCIterInput& xcItInp, const BatchIterBasicData& bibd) const
{
	//cout << "makebatchprop 1: thread::id = " << get_id() << endl;
	//Reform: this must be made a unique_ptr.
	BatchHoleAvg* bHoleAvg = new BatchHoleAvg(*this, xcInp, xcItInp, bibd);
	makeBatchHole(*bHoleAvg, xcInp, xcItInp, bibd);
	return bHoleAvg;
}

void HFXHole::makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp, 
                        const XCIterInput& xcItInp, 
                        const BatchIterBasicData& bibd) const
{
	//This is not necessary because there is only one hfx hole.  This just
	//serves as an example for other types of holes.
	(*this.*cpHFXHole)(bHoleAvg, xcInp, xcItInp, bibd);
}

void HFXHole::hfxHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp, 
                        const XCIterInput& xcItInp, 
                        const BatchIterBasicData& bibd) const
{
  //Double normXHole = 0;
  //Double espXHole = 0;
	//Old.  To be deleted.  Change here.  Uncomment for one-point
  //FILE* myfile = fopen("sValue.txt","w");
  for ( UInt is = 0; is < sValues.size(); is++ )
  {
   	HalfJKRho halfJKRho(xcItInp.infor,xcItInp.xcvar,bibd.bg,xcInp.ms,
			bibd.sigList,xcItInp.cartDen);
   	Double sValue = sValues[is];
   	//Double sValueWeight = sValueWeightVec[is];
   	vector<DoubleVec>& xhole = bHoleAvg.holes[is]; 
   	halfJKRho.doHalfJKRhoXHole(xcItInp.infor,xcItInp.xcvar,xcInp.ms,
		              bibd.bg,bibd.sigList,bibd.bbasis,xcInp.denMtrx,
		              xcItInp.cartDen,xcItInp.spData,sValue);
   	BatchVar bVar(bibd.sigList,bibd.bbasis,xcItInp.xcvar,xcInp.denMtrx);
   	bVar.buildExRho(xcInp.ms,bibd.sigList,xcInp.denMtrx,bibd.bbasis,
		                halfJKRho,xcItInp.xcvar); 
   	bVar.xhole(xcItInp.xcvar, xhole);
		{
			//Change here.  Uncomment the following 2 lines for one-point hole.
			//Old.  To be deleted.  Change here.  Uncomment for one-point
    	//fprintf(myfile,"xholeVal[%u] =  %12.30lf\n", 0, xhole[0][0]);
    	//fprintf(myfile,"sValue = %12.30lf\n",sValue);

    	//normXHole += sValueWeight*pow(sValue,2)*xholeVal[0];
    	//espXHole += sValueWeight*sValue*xholeVal[0];
  	}
	}
	//Change here.  Uncomment for one-point
	//Old.  To be deleted.  Change here.  Uncomment for one-point
  //fclose(myfile);
}

void XCHoleAvg::accumulate(const BatchProp* bProp)
{
	const BatchHoleAvg* bHoleAvg = (const BatchHoleAvg*)bProp;
	for (UInt ispin = 0; ispin < holeAvg.size(); ispin++)
	{
		for (UInt isval = 0; isval < holeAvg[ispin].size(); isval++)
		{
			holeAvg[ispin][isval] += bHoleAvg->holeAvg[ispin][isval];
		}
	}
}

