/**
 * CPP files corresponding to the atomgrids.h
 * \author  Fenglai liu and Jing kong
 */
#include <string>
#include <cstdio>
#include <iostream>
#include<boost/lexical_cast.hpp>
#include "excep.h"
#include "element.h"
#include "xcintsinfor.h"
#include "gridinfor.h"
#include "atomgrids.h"
using namespace element;
using namespace excep;
using namespace xcintsinfor;
using namespace gridinfor;
using namespace atomgrids;

///////////////////////////////////////////////////////////////////
//            batch grid information for one batch               //
///////////////////////////////////////////////////////////////////
BatchGridInfor::BatchGridInfor(const AtomGrids& atomGrids, const UInt& iBatch, const UInt& initBatchSize)
{
	// first step, evaluate that how much radial points we are going to account
	initRadIndex      = getRadBegin(atomGrids,iBatch);
	UInt initAngIndex = getAngBegin(atomGrids,iBatch);
	UInt nRad         = atomGrids.getNRadPts();   // total number of radial points in atom grid
	UInt batchLen     = initBatchSize;            // initial batch size
	UInt nTotalPts    = 0;                        // used to record the point number
	nRadPts           = 0;                        // number of radial points involved in this batch
	bool batchEnd     = false;
	for (UInt iRad = initRadIndex; iRad<nRad; iRad++) {
		nRadPts++;
		Double radius = atomGrids.getRadius(iRad);
		UInt region = atomGrids.getRegion(radius);
		UInt nAng   = atomGrids.getNAngForRegion(region);
		UInt begin  = 0;
		if (iRad == initRadIndex) begin = initAngIndex;
		for(UInt iAng = begin; iAng<nAng; iAng++) {
			nTotalPts++;
			if (nTotalPts==batchLen) {
				batchEnd = true;
				break;
			}
		}
		if (batchEnd) break;
	}
	batchSize = nTotalPts;

	// now let's reserve data for the result
	angGridOffset.reserve(nRadPts); 
	angGridIndex.reserve(2*nRadPts); 

	// now go to see the angular part information
	nTotalPts = 0;
	lastPointAngGridInParts = true;
	batchEnd = false;
	for (UInt iRad = 0; iRad<nRadPts; iRad++) {
		UInt rad    = initRadIndex + iRad;
		Double radius = atomGrids.getRadius(rad);
		UInt region = atomGrids.getRegion(radius);
		UInt nAng   = atomGrids.getNAngForRegion(region);
		UInt begin  = 0;
		if (iRad == 0) begin = initAngIndex;
		UInt offset = atomGrids.getAngGridOffset(region);
		angGridOffset.push_back(offset);
		angGridIndex.push_back(begin);
		UInt iAng;
		for(iAng = begin; iAng<nAng; iAng++) {
			nTotalPts++;
			if (nTotalPts==batchSize) {
				if (iAng==nAng-1) lastPointAngGridInParts = false;
				batchEnd = true;
				break;
			}
		}

		// push into the end index
		// if it goes over all of angular point, the last index 
		// is actually nAng, but we need nAng-1
		if (iAng<nAng) {
			angGridIndex.push_back(iAng);
		}else{
			angGridIndex.push_back(nAng-1);
		}
		if (batchEnd) break;
	}
}

UInt BatchGridInfor::getRadBegin(const AtomGrids& atomGrids, const UInt& iBatch) const {
	if (iBatch == 0) return 0;
	UInt previousBatch        = iBatch - 1;
	const BatchGridInfor& pb  = atomGrids.getBatchGridInfor(previousBatch);
	// this should be global radial point index in radPts array in atomGrids
	UInt lastRadPointIndex    = pb.initRadIndex+pb.nRadPts-1; 

	// does the last point in previous batch really end?
	if (pb.lastPointAngGridInParts) return lastRadPointIndex;
	return lastRadPointIndex + 1;
}

UInt BatchGridInfor::getAngBegin(const AtomGrids& atomGrids, const UInt& iBatch) const {
	if (iBatch == 0) return 0;

	// calculate the ang index if we have previous batch
	UInt previousBatch  = iBatch - 1;
	const BatchGridInfor& pb = atomGrids.getBatchGridInfor(previousBatch);
	if (pb.lastPointAngGridInParts) {
		UInt lastRadPoint = pb.nRadPts-1;  // this is the last point index inside the batch
		UInt angEndIndex  = pb.angGridIndex[2*lastRadPoint+1];
		return angEndIndex + 1;
	}

	// if the final radial point in last batch has all of its angular
	// points, then this new batch's ang part should start from 0
	return 0;
}

///////////////////////////////////////////////////////////////////
//                        atom grid class                        //
///////////////////////////////////////////////////////////////////
AtomGrids::AtomGrids(const UInt& atomic, const Double& atomSize0, 
		const XCIntJobInfor& xcInfor):pruneMethod(xcInfor.getPruneGridMethod()),
	choice(xcInfor.getGridChoice()),nRegions(0),Z(atomic),atomSize(atomSize0), 
	atomRadii(ZERO), nTotalPts(0), nRadPts(xcInfor.getNRad()),nAngPts(xcInfor.getNAng()), 
	initBatchSize(xcInfor.getInitBatchSize())
{
	atomRadii = element::getAtomRadii(Z);
	preCheck();
	formQuadraturePtsWts(xcInfor);
	formBatch(xcInfor);
}

void AtomGrids::preCheck() const
{
	// check the atomic number
	if (Z == 0 || Z > MAX_ATOMS_TYPE+1) {
		string infor = "Improper Z is: " + boost::lexical_cast<string>(Z);
		Excep excep("AtomGrids","preCheck",INVALID_ATOMIC_NUMBER,infor);
		handleExcep(excep);
	}

	// check the grid information
	if (! usePruneGrids()) {

		// check lebdev grids
		bool doingGood = false;
		for(UInt i=0; i<NUMBER_VALID_LEBEDEV; i++) {
			if (nAngPts == VALID_LEBEDEV_GRID_NUMBER[i]) {
				doingGood = true;
				break;
			}
		}

		// now print the information
		if (!doingGood) {
			cout << "Improper number of angular points: " << nAngPts << endl;
			cout << "In Lebedev method, valid numbers of grid are: " << endl;

			// eight numbers per line
			UInt n = NUMBER_VALID_LEBEDEV%8;
			UInt m = 0;
			UInt i = 0;
			while(m < n) {
				cout << VALID_LEBEDEV_GRID_NUMBER[i] << "  " 
					<< VALID_LEBEDEV_GRID_NUMBER[i+1] << "  "
					<< VALID_LEBEDEV_GRID_NUMBER[i+2] << "  "
					<< VALID_LEBEDEV_GRID_NUMBER[i+3] << "  "
					<< VALID_LEBEDEV_GRID_NUMBER[i+4] << "  "
					<< VALID_LEBEDEV_GRID_NUMBER[i+5] << "  "
					<< VALID_LEBEDEV_GRID_NUMBER[i+6] << "  "
					<< VALID_LEBEDEV_GRID_NUMBER[i+7] << endl;
				m++;
				i = i + 8;
			}

			// now issue error infor
			string infor = "invalid Lebdev point is given";
			Excep excep("AtomGrids","preCheck",EXCEPTION_XCINTS_INVALID_GRID_CHOICE,infor);
			handleExcep(excep);
		}	
	}
}

void AtomGrids::print(UInt level) const {

	if (level >= 0) {
		cout << "Atom's atomic number " << Z << endl;
		cout << "Atom size is " << atomSize << " Atom radii is " << atomRadii << endl;
		cout << "Total number of grid points " << nTotalPts << endl;
		cout << "Total number of radial points " << nRadPts << endl;
		if (usePruneGrids()) {
			cout << "Grid is pruned " << endl;
		}else{
			cout << "Grid is not pruned " << endl;
		}
		cout << endl << endl;
	}

	if (level >= 1) {
		cout << "Grid Information Given in Batch" << endl;
		for(UInt i=0; i<getNBatch(); i++) {

			// get batch grid data
			const BatchGridInfor& bg = getBatchGridInfor(i);
			UInt radStart = bg.initRadIndex;
			UInt nRad = bg.nRadPts;

			// printing index information for each batch
			cout << "For current batch: " << i+1 << " Number of points " << getBatchSize(i) << endl;
			printf ("%-12s  %-12s  %-12s  %-12s  %-12s  %-12s  %-12s\n", "Num.", "Rad index", 
					"Rad value", "NAng", "Ang Start", "Ang End", "Ang Offset");
			UInt count = 1;
			for(UInt iRad=0; iRad<nRad; iRad++) {
				UInt rindex = radStart + iRad;
				UInt begin, end;
				getAngGridIndex(i,iRad,begin,end);
				Double r = getRadius(rindex);
				UInt region = getRegion(r);
				UInt nAng   = getNAngForRegion(region);
				UInt offset = getAngGridOffset(region);
				printf ("%-12d  %-12d  %-12.6f  %-12d  %-12d  %-12d  %-12d\n", (Int)count, (Int)rindex, r,
						(Int)nAng, (Int)begin, (Int)end, (Int)offset);
				count++;
			}
			cout << endl << endl;
		}
	}

	if (level >= 2) {

		// printing radial information
		printf ("%-12s  %-16s  %-16s\n", "Rad index", "Rad value", "Rad Wts");
		for(UInt iRad=0; iRad<nRadPts; iRad++) {
			Double r = radPts[iRad];
			Double w = radWts[iRad];
			printf ("%-12d  %-16.12f  %-16.12f\n", (Int)iRad, r, w);
		}
		cout << endl << endl;

		// printing angular information 
		cout << "Angular information" << endl;
		printf ("%12s  %16s  %16s  %16s  %16s\n", "Num.", 
				"Ang X", "Ang Y", "Ang Z", "Ang weight");
		UInt count = 1;
		UInt offset = 0;
		for(UInt iRegion=0; iRegion<nRegions; iRegion++) {
			UInt nAng = getNAngForRegion(iRegion);
			cout << "Region " << iRegion << " with number of " << nAng << " points" << endl;
			for(UInt iAng=0; iAng<nAng; iAng++) {
				Double X = angPts[3*(offset+iAng)+0];
				Double Y = angPts[3*(offset+iAng)+1];
				Double Z = angPts[3*(offset+iAng)+2];
				Double W = angWts[offset+iAng];
				printf ("%12d  %16.12f  %16.12f  %16.12f  %16.12f\n", (Int)count, X, Y, Z, W);
				count++;
			}
			offset += nAng;
			cout << endl << endl;
		}
		cout << endl << endl;
	}
}

///////////////////////////////////////////////////////
// !!! functions related to the grid information     //
///////////////////////////////////////////////////////
void AtomGrids::formGridInfor()
{
	if (usePruneGrids()) {

		//
		// here in this section we record the prune grid information
		// all of information here is hard-coding
		//
		if (pruneMethod == SG1_GRID) {

			// scale factor for each region
			nRegions = 5;
			scaleFac.reserve(nRegions-1);
			if (Z <= 2) {
				scaleFac.push_back(0.25E0);  
				scaleFac.push_back(0.5E0);  
				scaleFac.push_back(1.0E0);  
				scaleFac.push_back(4.5E0);  
			}else if (Z > 2 && Z <= 10) {
				scaleFac.push_back(0.1667E0);  
				scaleFac.push_back(0.5E0);  
				scaleFac.push_back(0.9E0);  
				scaleFac.push_back(3.5E0);  
			}else if (Z > 10 && Z <= 18) {
				scaleFac.push_back(0.1E0);  
				scaleFac.push_back(0.4E0);  
				scaleFac.push_back(0.8E0);  
				scaleFac.push_back(2.5E0);  
			}else {
				cout << "Atomic number is " << Z << endl;
				string infor = "currently SG1 only support pruned grids for Z<=18"; 
				Excep excep("AtomGrids","formGridInfor",EXCEPTION_XCINTS_INVALID_GRID_CHOICE,infor);
				handleExcep(excep);
			}

			// angular information
			angInRegion.assign(nRegions,0);
			angInRegion[0] = 6;
			angInRegion[1] = 38;
			angInRegion[2] = 86;
			angInRegion[3] = 194;
			angInRegion[4] = 86;

			// finally, we will set the radial points
			// SG1 is based on (50,194)
			nRadPts = 50;

		}else if (pruneMethod == BAKER_GRID) {

			// have an eye on the choice
			if (choice != STANDARD_GRID && choice != COARSE_GRID && choice != FINE_GRID) {
				string infor = "invalid choice for grid calibration when using BAKER grids"; 
				Excep excep("AtomGrids","formGridInfor",EXCEPTION_XCINTS_INVALID_GRID_CHOICE,infor);
				handleExcep(excep);
			}

			// now form the angular grid
			if (choice == STANDARD_GRID) {

				// number of regions
				nRegions = 4;
				scaleFac.reserve(nRegions-1);
				angInRegion.reserve(nRegions);

				// H and He
				if (Z <= 2) {

					// scaleFac data
					scaleFac.push_back(0.5E0);  
					scaleFac.push_back(1.0E0);       
					scaleFac.push_back(7.5E0);

					// number of angular points
					angInRegion.push_back(14);
					angInRegion.push_back(50);
					angInRegion.push_back(194);
					angInRegion.push_back(50);
				}else{

					// scaleFac data
					scaleFac.push_back(0.25E0);  
					scaleFac.push_back(1.0E0);       
					scaleFac.push_back(7.5E0);

					// number of angular points
					angInRegion.push_back(26);
					angInRegion.push_back(110);
					angInRegion.push_back(302);
					angInRegion.push_back(110);
				}

			}else if (choice == COARSE_GRID) {

				// number of regions
				nRegions = 2;
				if (Z > 18) nRegions = 4;
				scaleFac.reserve(nRegions-1);
				angInRegion.reserve(nRegions);

				if (Z <= 2) {

					// scaleFac data
					scaleFac.push_back(0.7E0);  

					// number of angular points
					angInRegion.push_back(14);
					angInRegion.push_back(50);

				}else if (Z >= 3 && Z<= 18){

					// scaleFac data
					scaleFac.push_back(0.35E0);  

					// number of angular points
					angInRegion.push_back(26);
					angInRegion.push_back(110);

				}else {

					// scaleFac data
					scaleFac.push_back(0.25E0);  
					scaleFac.push_back(0.75E0);  
					scaleFac.push_back(5.0E0);  

					// number of angular points
					angInRegion.push_back(14);
					angInRegion.push_back(50);
					angInRegion.push_back(194);
					angInRegion.push_back(50);

				}	

			}else {

				// number of regions
				nRegions = 3;
				scaleFac.reserve(nRegions-1);
				angInRegion.reserve(nRegions);

				// this is fine grid
				if (Z <= 2) {

					// scaleFac data
					scaleFac.push_back(0.25E0);  
					scaleFac.push_back(0.75E0);  

					// number of angular points
					angInRegion.push_back(26);
					angInRegion.push_back(110);
					angInRegion.push_back(302);

				}else {

					// scaleFac data
					scaleFac.push_back(0.175E0);  
					scaleFac.push_back(0.625E0);  

					// number of angular points
					angInRegion.push_back(50);
					angInRegion.push_back(194);
					angInRegion.push_back(434);
				}	
			}

			// finally, we will set the radial points
			// overwrite the original value
			if (Z<=2) {
				nRadPts = 30;
			}else if (Z>3 && Z<=10) {
				nRadPts = 50;
			}else if (Z>11 && Z<=18) {
				nRadPts = 70;
			}else if (Z>19 && Z<=36) {
				nRadPts = 90;
			}else{
				nRadPts = 120;
			}
		}
	}else{
		// no scale factor need to be defined
		// we only need the region information
		nRegions = 1;
		angInRegion.assign(nRegions,0);
		angInRegion[0] = nAngPts;
	}
}

UInt AtomGrids::getRegion(const Double& radius) const
{
	// for pruned grid, we need to search it's region
	// else we just return 0
	if (usePruneGrids()) {
		for(UInt i=0; i<nRegions; i++) {

			// determine the lower limit and upper limit
			// for the given region. We note that
			// for the upper limit, we just give a very 
			// large distance, which could be 
			// never reached for any atom
			Double lowerLimit = ZERO;
			if (i>0) lowerLimit = scaleFac[i-1]*atomRadii;
			Double upperLimit = 100000E0; 
			if (i<nRegions-1) upperLimit = scaleFac[i]*atomRadii;

			// now we determine the region
			if (radius >= lowerLimit && radius < upperLimit) return i;
		}
	}
	return 0;
}

UInt AtomGrids::getTotalAngPtsInType() const 
{
	UInt np = 0;
	for (UInt iRegion=0; iRegion<nRegions; iRegion++) {
		UInt nAng  = angInRegion[iRegion];
		np += nAng;
	}
	return np;
}

UInt AtomGrids::getAngGridOffset(const UInt& region) const
{
	UInt offset = 0;
	for (UInt iRegion=0; iRegion<nRegions; iRegion++) {
		if (iRegion == region) return offset;
		UInt nAng  = angInRegion[iRegion];
		offset += nAng;
	}
	return offset;
}

void AtomGrids::eulmac(const UInt& nRad, Double* Pts, Double* Wts) const {

	//
	// create the points by using Eular-Maclaurine formula
	// see paper "A standard grid for density functional calculations"
	// Peter M. W. Gill and Benny G. Johnson and John A. Pople
	// Chemical Physics Letters, 1993, Vol. 209, issue 5-6, Page 506-512
	//
	Double r = atomRadii;
	for(UInt i=0; i<nRad; i++) {
		Double t = i+1;
		Double u = nRad-i;
		Double tdu = t/u;
		Pts[i] = r*pow(tdu,TWO);
		Wts[i] = TWO*pow(r,THREE)*(nRad+1)*pow(tdu,FIVE)/pow(u,TWO);
	}
}

void AtomGrids::formQuadraturePtsWts(const XCIntJobInfor& infor)
{

	// form the grid information
	formGridInfor();

	// forming the radial points
	DoubleVec radialPts(nRadPts,ZERO);
	DoubleVec radialWts(nRadPts,ZERO);
	eulmac(nRadPts,&radialPts.front(),&radialWts.front());

	// now consider the atom size we have to cut points off 
	// also if we do insignificant grid check
	// then we do it here
	radPts.reserve(nRadPts);
	radWts.reserve(nRadPts);
	UInt count = 0;
	bool insigGridsCheck = infor.doSigGridCheck();
	Double gridCheckThresh = infor.getInsigGridCheckThresh();
	for(UInt i=0; i<nRadPts; i++) {

		// insig grid check
		// threshold value is set near to the double limit
		if (insigGridsCheck) {
			if (radialWts[i] < gridCheckThresh) continue;
		}

		// grid cut off?
		if (radialPts[i] > atomSize) break;
		radPts.push_back(radialPts[i]);
		radWts.push_back(radialWts[i]);
		count++;
	}
	nRadPts = count; // this is the real number of radical points

	// set the tmp array for doing ang pts calculation
	UInt nTmpAng = 0;
	for(UInt i=0; i<nRegions; i++) {
		UInt nAng = angInRegion[i];
		if (nAng>nTmpAng) nTmpAng = nAng;
	}
	DoublePrecVec tmpAngPts(3*nTmpAng,ZERO);
	DoublePrecVec tmpAngWts(nTmpAng,ZERO);

	// forming the angular points
	// we note the process is same between pruned or un-pruned grids
	UInt nAngGrids = getTotalAngPtsInType();
	angPts.assign(3*nAngGrids,ZERO);
	angWts.assign(nAngGrids,ZERO);
	UInt offset = 0; 
	for(UInt i=0; i<nRegions; i++) {

		// make up some tmp array to hold the data
		Int nAng = static_cast<Int>(angInRegion[i]);

		// make up the angular data
		lebedev(&tmpAngPts[0],&tmpAngPts[nAng],&tmpAngPts[2*nAng],&tmpAngWts[0],&nAng);

		// copy the tmp data into the real data, and we scale the weight with 4*PI
		UInt pos = 3*offset;
		for(Int i=0; i<nAng; i++) {
			angPts[pos+3*i  ]  = tmpAngPts[i];
			angPts[pos+3*i+1]  = tmpAngPts[i+nAng];
			angPts[pos+3*i+2]  = tmpAngPts[i+2*nAng];
			angWts[offset+i ]  = tmpAngWts[i]*FOUR*PI; 
		}

		// increment the offset
		offset += nAng;
	}
}

///////////////////////////////////////////////////////
// !!! functions related to the batch information    //
///////////////////////////////////////////////////////
void AtomGrids::formBatch(const XCIntJobInfor& xcInfor)
{
	// firstly, we need to calculate the total number of points
	UInt np = 0;
	for(UInt iRad=0; iRad<nRadPts; iRad++) {
		Double radius  = getRadius(iRad);
		UInt region    = getRegion(radius);
		UInt nAng      = getNAngForRegion(region);
		np += nAng;
	}
	nTotalPts = np;

	// form the batch 
	//
	// there are many ways to form the batch for this selected atom grids.
	//
	// One of the method to form the batch, is to make all of available 
	// threads share "evenly" on the work load. For example, if this atom
	// has 4000 grids and we have 20 threads; then if we let all of threads
	// working on this atom simutaniously then each thread should work on
	// 200 grid pints.
	//
	// Comparing this way, another way is to stick to the input batch size
	// which is given from the user. For example, this atom still has
	// 4000 grids; and the user wants each thread work on 500 grids then
	// the workload for this atom will be shared by 8 threads; the rest
	// of the 12 threads will be assigned to the following atoms.
	//
	// there are also other ways to form the batch. For example, each angular
	// surface (for example, for 128 302 grid set each 302 set we can make a 
	// batch) can be wrapped into a batch. However, we think the first and 
	// second methods are best for selection. Because they are considering 
	// the relation between the workload and working threads, this should 
	// be the main concern here.
	//
	// We prefer to the second way. In the long run the machine may be quite
	// different on running the jobs. For example, the normal Xeon CPU may
	// have 20-40 threads, but the Knight Landing CPU has about 140 working 
	// threads (each core has 2 working threads). The user should decide
	// the setting for the batch for best select, not the program itself.
	//
	
	// OK, first let's see how many batches we may have
	// if the lest grids larger than half of the batch size,
	// we consider it to be an independent batch
	UInt nBatches  = 1;
	if (initBatchSize<nTotalPts) {
		UInt leftGrids = nTotalPts%initBatchSize;
		nBatches       = (nTotalPts-leftGrids)/initBatchSize;
		UInt halfBatchSize = static_cast<UInt>(initBatchSize*HALF);
		if (leftGrids>=halfBatchSize) {
			nBatches += 1;
		}
	}
	infor.reserve(nBatches);

	// now let's form the batch
	// according to this algorithm, only the last batch may have
	// different batch size than the initial batch size
	UInt leftPts = nTotalPts;
	for(UInt iBatch = 0; iBatch<nBatches; iBatch++) {
		UInt batchSize = initBatchSize;
		if (iBatch == nBatches - 1) batchSize = leftPts;
		leftPts -= batchSize;
		BatchGridInfor batch(*this,iBatch,batchSize);
		infor.push_back(batch);
	}
}

