/**
 * CPP files corresponding to the xcenergyinfor.h
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include<cstdio>
#include<cmath>
#include "blas1.h"
#include "excep.h"
#include "xcfunc.h"
#include "xcintsinfor.h"
#include "batchgrid.h"
#include "batchvar.h"
#include "xcints.h"
#include "xcintsprops.h"
#include "denmtrx.h"
#include "functionallist.h"

using namespace excep;
using namespace blas;
using namespace xcfunc;
using namespace xcintsinfor;
using namespace batchgrid;
using namespace batchvar;
using namespace xcints;
using namespace xcintsprops;
using namespace denmtrx;
//using namespace xcenergyinfor;
using namespace std;

BatchHoleAvg::BatchHoleAvg
  (const XCHoleAvg& hAvg, const XCInput& xi, const XCIterInput& xii, 
   const BatchIterBasicData& bibd) 
: BatchProp(xi,xii,bibd), xcHoleAvg(hAvg), 
  holes(hAvg.sValues.size(), vector<DoubleVec>(xi.denMtrx.getNSpin(), 
    DoubleVec(bibd.bg.getNGrids()))), 
  holeAvg(xi.denMtrx.getNSpin(), DoubleVec(hAvg.holeAvg.size())) {}

//fout name should be read in from an infor class from the caller routine
//in the future with a field.  See 'xtron desirables' on group docs.
XCHoleAvg::XCHoleAvg(UInt nSpin, UInt nSValue)
: holeAvg(nSpin, DoubleVec(nSValue)), sValues(nSValue), fout("sValueAvg.txt")
{
  Double radian = Double(1.50);
  for(UInt irad=0; irad<nSValue; irad++) {
    Double t = irad+1;
    Double u = nSValue-irad;
    Double tdu = t/u;
    sValues[irad] = radian*pow(tdu,TWO);
    //sValueWeightVec[irad] = TWO*pow(radian,THREE)*(nSValue+1)*
		//pow(tdu,FIVE)/pow(u,TWO);
  }
}

void XCHoleAvg::dump() const
{
	FILE* myfile = fopen(fout.c_str(),"w");
  //for (UInt ispin = 0; ispin < holeAvg.size(); ispin++)
	UInt ispin = 0;
  {
  	for (UInt isval = 0; isval < holeAvg[ispin].size(); isval++)
    {
    	fprintf(myfile,"xholeVal[%u] =  %12.30lf\n",isval, holeAvg[ispin][isval]);
    	fprintf(myfile,"sValue = %12.30lf\n",sValues[isval]);
    }
  }
  fclose(myfile);
}

void BatchHoleAvg::integrate()
{
	UInt nSpin = xcInp.denMtrx.getNSpin();
	UInt nGrid = bibd.bg.getNGrids();
	UInt nSValue = holes.size();
	vector<UInt> posDen(nSpin);
	posDen[0] = xcItInp.xcvar.getVarPos(RA);
	posDen[nSpin-1] = xcItInp.xcvar.getVarPos(RB);
	for (UInt isval = 0; isval < nSValue; isval++)
	{
		//for (UInt ispin = 0; ispin < nSpin; ispin++)
  // Looping over ispin leads to problem, maybe mismatch of variable dimensions
		UInt ispin = 0;
		{
			holeAvg[ispin][isval] = 0;
			const Double* wts = bibd.bg.getGridWts();
			const Double* den = bibd.bvar.getVar(posDen[ispin]);
			for (UInt ig = 0; ig < nGrid; ig++)
			{
				//Change here 2. one-point uses the second line.
				holeAvg[ispin][isval] += wts[ig]*den[ig]*holes[isval][ispin][ig];
				//holeAvg[ispin][isval]  = holes[isval][ispin][ig];
			}
		}
	}
}

MetaGGAHole::MetaGGAHole(UInt nSpin, UInt nSValue)
  : XCHoleAvg(nSpin, nSValue)
{
	//Change here for different MetaGGA holes.
	cpMetaGGAHole = &br89xhole_s;
}

void MetaGGAHole::makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp,
     const XCIterInput& xcItInp, const BatchIterBasicData& bibd) const
{
  const Double *rhoA(0), *rhoB(0), *DRA(0), *DRB(0), *TA(0), *TB(0), 
	             *LA(0), *LB(0), *EXRA(0), *EXRB(0);
  bibd.bvar.varForFort(rhoA, rhoB, DRA, DRB, TA, TB, LA, LB, EXRA, EXRB, 
	                     xcItInp.xcvar);
  Int ng_     = static_cast<Int>(bibd.bg.getNGrids());
  Int nDen_   = static_cast<Int>(xcItInp.infor.getNDen());

	for ( UInt is = 0; is < sValues.size(); is++ )
  {
		Double sValue = sValues[is];
		vector<DoubleVec>& xhole = bHoleAvg.holes[is];
		(*cpMetaGGAHole)(&xhole[0].front(), &sValue, rhoA, rhoB, DRA, DRB, LA, LB,
      TA, TB, &ng_, &nDen_);
	}
}
