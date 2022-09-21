/**
 * CPP files corresponding to the molecule grids
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include<cstdio>
#include<boost/lexical_cast.hpp>
#include "element.h"
#include "shellsize.h"
#include "xcintsinfor.h"
#include "shell.h"
#include "excep.h"
#include "atom.h"
#include "molegrids.h"
using namespace atom;
using namespace shell;
using namespace excep;
using namespace element;
using namespace shellsize;
using namespace xcintsinfor;
using namespace molegrids;

const AtomGrids& MoleGrids::getAtomGrid(const UInt& Z) const {

	// search the one that match the atomic number
	for(UInt i=0; i<getNAtomGrids(); i++) {
		const AtomGrids& ag = molGrids[i];
		if (Z == ag.getAtomic()) return ag;
	}

	// something wrong here
	string infor = "improper Z is: " + boost::lexical_cast<string>(Z);
	Excep excep("MoleGrids","getAtomGrid",EXCEPTION_XCINTS_FAIL_TO_FIND_ATOMGRIDS,infor);
	handleExcep(excep);
	return molGrids[0];
}

MoleGrids::MoleGrids(const MolShell& ms, const MolShellSize& shellSize, const XCIntJobInfor& infor)
{
	// creating molecule grids based on atom types
	UInt nAtoms = shellSize.getNAtomTypes(); 
	UIntVec atomTypes;
	atomTypes.reserve(nAtoms);
	shellSize.getAtomTypeInfor(atomTypes);

	// we may need to take a look that wether we have included
	// the ghost atom
	UInt nAtomGrids = nAtoms;
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		UInt Z = atomTypes[iAtom];
		if (isGhostAtom(Z)) {
			nAtomGrids = nAtomGrids - 1;
			break;
		}
	}

	// now let's construct atom grids
	molGrids.reserve(nAtomGrids);
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		UInt Z = atomTypes[iAtom];
		if (isGhostAtom(Z)) continue;
		Double size = shellSize.getAtomSize(Z);
		AtomGrids atomGrids(Z,size,infor);
		molGrids.push_back(atomGrids);
	}

	// now let's form the batch information based
	// on the molGrids
	// firstly let's evaluate that how many batches
	// we have
	// we will get the geometry information from the atom shell
	// we will need the atom information
	nAtoms = ms.getNAtomShells();
	UInt nBatches = 0;
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const Atom& atom = ms.getAtomShell(iAtom);
		UInt atomic = atom.getAtomic();
		if (isGhostAtom(atomic)) continue;
		const AtomGrids& ag = getAtomGrid(atomic);
		nBatches += ag.getNBatch();
	}

	// ok, let's store the batch information right now
	batchInfor.reserve(3*nBatches);
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const Atom& atom = ms.getAtomShell(iAtom);
		UInt atomic = atom.getAtomic();
		if (isGhostAtom(atomic)) continue;
		const AtomGrids& ag = getAtomGrid(atomic);
		for(UInt iBatch=0; iBatch<ag.getNBatch(); iBatch++) {
			batchInfor.push_back(iAtom);
			batchInfor.push_back(atomic);
			batchInfor.push_back(iBatch);
		}
	}
}

UInt MoleGrids::getNMaxBatchGrids() const
{
	// please remember, that here nAtoms is actually
	// the number of atom types
	UInt nAtoms = getNAtomGrids();
	UInt maxBatchSize = 0;
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomGrids& ag = molGrids[iAtom];
		for(UInt iBatch=0; iBatch<ag.getNBatch(); iBatch++) {
			UInt size = ag.getBatchSize(iBatch);
			if (maxBatchSize<size) maxBatchSize = size;
		}
	}

	// now let's return 
	return maxBatchSize;
}

void MoleGrids::print(UInt level) const
{
	cout << "*******************************************" << endl;
	cout << "*              Grids Results              *" << endl;
	cout << "*******************************************" << endl;
	// print out the batch information
	UInt nBatches = getNTotalBatches();
	printf("print out the batch information for the whole molecule\n");
	printf("maximum batch grid size is %d\n", getNMaxBatchGrids());
	printf("%-12s %-6s %-6s %-17s\n", "batch index", "atom", "atomic", "local batch index");
	for(UInt iBatch=0; iBatch<nBatches; iBatch++) {
		UInt Z, index, iAtom;
		getBatchInfor(iBatch,iAtom,Z,index);
		printf("%-12d %-6d %-6d %-6d\n", (Int)iBatch, (Int)iAtom,(Int)Z,(Int)index);
	}
	
	// print out the grid information
	UInt nAtoms = getNAtomGrids();
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		cout << "Atom " << iAtom << endl;
		molGrids[iAtom].print(level);
		cout << endl << endl;
	}
}
