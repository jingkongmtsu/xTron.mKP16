/**
 * CPP files corresponding to the batch grid 
 * \author  Fenglai liu and Jing kong
 *
 *  Becke weight(this is also called partition weights) reference:
 *
 *  Becke A D. 
 *  A multicenter numerical integration scheme for polyatomic molecules[J]. 
 *  The Journal of Chemical Physics, 1988, 88: 2547.
 *
 */
#include<iostream>
#include<cstdio>
#include "atom.h"
#include "shell.h"
#include "excep.h"
#include "atomgrids.h"
#include "xcintsinfor.h"
#include "sigatombasis.h"
#include "batchgrid.h"
using namespace atom;
using namespace shell;
using namespace excep;
using namespace xcintsinfor;
using namespace sigatombasis;
using namespace atomgrids;
using namespace batchgrid;
using namespace std;

BatchGrid::BatchGrid(const AtomGrids& atomG, const Atom& atom, 
		const UInt& iBatch, const bool& onePoint) :
		inRad(ZERO),outRad(ZERO)
{
	if ( onePoint )
	{
		init1Point(atomG, atom, iBatch);
	}
	else
	{
		init(atomG, atom, iBatch);
	}
}

void BatchGrid::init1Point(const AtomGrids& atomG, const Atom& atom, const UInt& iBatch)
{
	nGrids = 1;
	coord = DoubleVec(3);
	coord[0] = 0; coord[1] = 0; coord[2] = 0.7180959;  //Change here.
	inRad = pow(atom.getXYZ()[0] - coord[0], 2) 
		+ pow(atom.getXYZ()[1] - coord[1], 2) 
		+ pow(atom.getXYZ()[2] - coord[2], 2);
	inRad = sqrt(inRad);
	outRad = inRad;
	batchAtom = atom.getAtomIndex();
	wts = DoubleVec(nGrids,ZERO);
	wts[0] = 1;
}
	

void BatchGrid::init(const AtomGrids& atomG, const Atom& atom, const UInt& iBatch)
{
	nGrids = atomG.getBatchSize(iBatch);
	batchAtom = atom.getAtomIndex();
	coord = DoubleVec(3*nGrids,ZERO);
	wts = DoubleVec(nGrids,ZERO);
	// data array from atom grids
	const Double* rad = atomG.getRadPts(iBatch);
	const Double* radWts = atomG.getRadWts(iBatch);
	UInt nRad = atomG.getNRadInBatch(iBatch);

	// the batch grid is actually located in the inner ring
	// and outer ring - it's possible that the inner ring radius
	// is same with the outer ring, however; they may be different
	inRad = rad[0];
	outRad= rad[nRad-1];

	// now let's generate the qudrature weights and coordinates
	// partition weights is later added in
	UInt i=0;
	for(UInt iRad=0; iRad<nRad; iRad++) {

		// radial point value
		Double radP = rad[iRad];
		Double radW = radWts[iRad];

		// angular information
		const Double* ang = atomG.getAngPts(iBatch,iRad);
		const Double* angWts = atomG.getAngWts(iBatch,iRad);
		UInt angBegin,angEnd;
		atomG.getAngGridIndex(iBatch,iRad,angBegin,angEnd);
		const Double* xyz = atom.getXYZ();
		for(UInt iAng=angBegin; iAng<=angEnd; iAng++) {
			coord[3*i  ] = radP*ang[3*iAng  ] + xyz[0]; // x
			coord[3*i+1] = radP*ang[3*iAng+1] + xyz[1]; // y
			coord[3*i+2] = radP*ang[3*iAng+2] + xyz[2]; // z
			wts[i]       = radW*angWts[iAng];           // weight
			i++;
		}
	}
}

void BatchGrid::BeckeWeights0(const XCIntJobInfor& xcinfor, const MolShell& ms, const SigAtomBasis& sigList) 
{
	// firstly, construct the neighbors for the batch grid
	// We use the sigfinicant atoms data
	const UIntVec& neigbors = sigList.getSigAtoms();

	// tmp data - based on the neighbor atoms
	UInt nAtoms = neigbors.size();
	DoubleVec rAG(nAtoms);
	DoubleVec sij(nAtoms*nAtoms);

	// we need to identify which atom is the mother 
	UInt mother = nAtoms+1;
	for(UInt i=0; i<nAtoms; i++) {
		if (neigbors[i] == batchAtom) {
			mother = i;
			break;
		}
	}
	if (mother == nAtoms+1) {
		string infor = "we did not find the batch grid mother atom from the input sigList";
		Excep excep("BatchGrid","BeckeWeights0",EXCEPTION_XCINTS_FAIL_TO_LOCATE_BATCH_ATOM_IN_SIGATOMLIST,
				infor);
		handleExcep(excep);
	}

	// for saving time, we also need a vector to store the 
	// distance between atom pairs
	DoubleVec atomPairDis(nAtoms*nAtoms,ZERO);
	for(UInt i=0; i<nAtoms; i++) {
		UInt iAtom = neigbors[i];
		const AtomShell& iA = ms.getAtomShell(iAtom);
		const Double* xyz = iA.getXYZ();
		for(UInt j=i; j<nAtoms; j++) {
			if (j==i) continue;
			UInt jAtom = neigbors[j];
			const AtomShell& jA = ms.getAtomShell(jAtom);
			Double dis = jA.getDistance(xyz);
			atomPairDis[j+i*nAtoms] = dis;
			atomPairDis[i+j*nAtoms] = dis;
		}
	}

	// now it's the real job
	for(UInt iG=0; iG<nGrids; iG++) {

		// for a given grid point, get the distance matrix
		rAG.assign(nAtoms,ZERO);
		for(UInt i=0; i<nAtoms; i++) {
			UInt iAtom = neigbors[i];
			const AtomShell& iA = ms.getAtomShell(iAtom);
			rAG[i] = iA.getDistance(&coord[3*iG]);
		}

		//
		// now for the given point, let's calculate the step function 
		// of s(i,j).
		// in this array, i is the reference point and j will loop
		// over all of its neigbor atoms. Furthermore, as j == i
		// the step function we set to one so to make it not
		// affect the cell function(cell function is multiplication of 
		// step functions).
		// furthermore, we note that the function defined in Becke's original
		// paper is actually an odd function. Therefore, if we
		// have mu_{ij} -> mu_{ji}, then the final function 
		// just change the sign-it's easy to see, that s_{ij}+s_{ji} = 1
		//
		// furthermore, for how to deriving the step function; we have multiple
		// ways. One way is Becke's original way, the other way is given by
		// R. Eric Stratmann and Gustavo E. Scuseria and Michael J. Frisch
		// "Achieving linear scaling in exchange-correlation density functional 
		// quadratures"
		// Chemical Physics Letters, 257, 213 - 223, 1996
		//
		// Becke's original paper is here:
		// "A multicenter numerical integration scheme for polyatomic molecules"
		// The Journal of Chemical Physics, 88, 2547-2553, 1988
		//
		sij.assign(nAtoms*nAtoms,ZERO);
		if (xcinfor.getStepFunc() == BECKE_ORI_STEP_FUNCTION) {
			for(UInt i=0; i<nAtoms; i++) {
				sij[i+i*nAtoms] = ONE;
				Double ri = rAG[i];
				for(UInt j=i+1; j<nAtoms; j++) {

					// get the distance
					Double rij = atomPairDis[j+i*nAtoms];

					// perform real calculation
					Double mu  = (ri-rAG[j])/rij;
					Double f1  = HALF*mu*(THREE-mu*mu);
					Double f2  = HALF*f1*(THREE-f1*f1);
					Double f3  = HALF*f2*(THREE-f2*f2);
					Double s   = (ONE-f3)*HALF;
					sij[j+i*nAtoms] = s;
					sij[i+j*nAtoms] = ONE - s;
				}
			}
		}else if (xcinfor.getStepFunc() == STRATMANN_STEP_FUNCTION) {

			// currently we only have two ways. this is CPL, 1996
			// they set the empirical parameter of a
			Double a = 0.64E0;

			// now compute the sij
			for(UInt i=0; i<nAtoms; i++) {
				sij[i+i*nAtoms] = ONE;
				Double ri = rAG[i];
				for(UInt j=i+1; j<nAtoms; j++) {

					// get the distance
					Double rij = atomPairDis[j+i*nAtoms];

					// perform real calculation
					// this is a piece wise function
					Double mu  = (ri-rAG[j])/rij;
					Double f   = ZERO;
					if (mu<=-a) {
						f = -1.0E0;
					}else if (mu>=a) {
						f = 1.0E0;
					}else {
						Double mua  = mu/a;
						Double mua3 = mua*mua*mua;
						Double mua5 = mua3*mua*mua;
						Double mua7 = mua5*mua*mua;
						f = 0.0625E0*(35.0E0*mua-35.0E0*mua3+21.0E0*mua5-5.0E0*mua7);
					}
					Double s   = (ONE-f)*HALF;
					sij[j+i*nAtoms] = s;
					sij[i+j*nAtoms] = ONE - s;
				}
			}
		}

		// calculate cell function 
		Double pi = ONE;    // cell function for the parent atom
		Double pt = ZERO;   // sum of cell function - appear in the denominator
		for(UInt j=0; j<nAtoms; j++) {
			pi *= sij[j+mother*nAtoms];  // cell function for mother atom
			Double pj = ONE;             // cell function for the atom j
			for(UInt k=0; k<nAtoms; k++) {
				pj *= sij[k+j*nAtoms];
			}
			pt += pj;
		}
		wts[iG] *= pi/pt;
		//printf("becke weights %d  %-12.8f\n", iG+1, pi/pt);
	}
}

void BatchGrid::print() const {

	//print out the batch points
	cout << "*******************************************" << endl;
	cout << "*           BatchGrid Result              *" << endl;
	cout << "*******************************************" << endl;
	printf("%-5s  %-14s  %-14s  %-14s  %-14s\n", "Num", "X ", "Y ", "Z ", "Weight");
	cout << "BatchGrid nGrids = " << nGrids << endl;
	for(UInt i=0; i<nGrids; i++) {
		printf("%-5d  % -14.7f  % -14.7f  % -14.7f  % -14.7f\n", (Int)i+1, coord[3*i  ], 
				coord[3*i+1], coord[3*i+2], wts[i]);
	}
	cout << endl << endl;
}
