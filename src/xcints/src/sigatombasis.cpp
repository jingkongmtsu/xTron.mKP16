/**
 * CPP files corresponding to the sigatombasis.h
 * \author fenglai liu and jing kong
 */
#include<boost/lexical_cast.hpp>
#include<iostream>
#include<algorithm>
#include "batchgrid.h"
#include "excep.h"
#include "shell.h"
#include "shellsize.h"
#include "xcintsinfor.h"
#include "sigatombasis.h"
using namespace batchgrid;
using namespace shell;
using namespace excep;
using namespace xcintsinfor;
using namespace shellsize;
using namespace sigatombasis;
using namespace std;

void SigAtomBasis::formSigDataWithAtomCenterApprox(const BatchGrid& grid, const MolShell& ms)
{
	// get the batch atom infor
	UInt batchAtom = grid.AtomInCurrentBatch();
	const AtomShell& atomShl = ms.getAtomShell(batchAtom);
	UInt nSigBasis = atomShl.getNBas();
	UInt nSigShell = atomShl.getNShell();

	// now ready to push data
	// we only have one significant atom
	sigAtoms.assign(1,0);
	sigAtoms[0] = batchAtom;
	atomBasisMapping.assign(2,0);
	atomBasisMapping[0] = 0;
	atomBasisMapping[1] = nSigBasis;
	atomShellMapping.assign(2,0);
	atomShellMapping[0] = 0;
	atomShellMapping[1] = nSigShell;

	// produce the sig shell infor
	sigBas.reserve(nSigBasis); 
	sigShell.reserve(nSigShell); 
	for(UInt iShell=0; iShell<atomShl.getNShell(); iShell++) {
		const Shell& s = atomShl.getShell(iShell);
		UInt basOffset = s.getBasisIndex(0,TYPE_NORM);
		for(UInt i=0; i<s.getNBas(); i++) {
			UInt gIndex = basOffset + i;
			sigBas.push_back(gIndex);
		}
		sigShell.push_back(s.getLocalShellIndex());
	}
}

void SigAtomBasis::formSigData(const BatchGrid& grid, const MolShell& ms, 
		const MolShellSize& shellSize) 
{
	// get the batch grid size
	Double inSize, outSize;
	grid.getBatchGridRadius(inSize,outSize);

	// firstly, let's try to approximate the significant atoms
	UInt nSigAtoms = 0;               
	UInt nSigBasis = 0;    
	UInt nSigShell = 0;    
	UInt nG        = grid.getNGrids();
	DoubleVec atomGridDistance(nG);
	UIntVec neighbors;
	neighbors.reserve(ms.getNAtomShells());
	UInt parentAtom = grid.AtomInCurrentBatch();
	const AtomShell& parentAtomShl = ms.getAtomShell(parentAtom);
	for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {

		//
		// two situations:
		// 1 if the atom and the parentAtom are very far away from
		// each other, then their distance is larger than the sum
		// of atom's radius plus the given grid outer shell radius;
		// then in this case it means the atom can not cover the 
		// grid, this is insignificant atom (r > outSize+iAtomSize)
		//
		// 2 if both two atoms are close with each other, then the 
		// grid is far away from the parentAtom; in this case;
		// the inner shell radius of the grid is even larger than
		// the sum of atom distance and the atom radius, then
		// this atom is also insignificant (r + iAtomSize < inSize)
		//
		// Here we note that this coarse estimation we have an
		// assumption, that all of batch grids are forming an
		// entire concentric. 
		//
		//
		const AtomShell& atomShl = ms.getAtomShell(iAtom);
		const Double* xyz = atomShl.getXYZ();
		Double r = parentAtomShl.getDistance(xyz);
		UInt atomic = atomShl.getAtomic();
		Double iAtomSize = shellSize.getAtomSize(atomic);
		if ((r + iAtomSize < inSize) || (r > outSize+iAtomSize)) {
			continue;
		}

		//
		// now let's step into more accurate estimation
		// compute the distance between the grid and the atom
		//
		const Double* coord = grid.getGridCoord();
		for(UInt iG=0; iG<nG; iG++) {
			atomGridDistance[iG] = atomShl.getDistance(&coord[3*iG]);
		}

		// now we need accurately estimate based on the atom size
		bool isSigAtom = false;
		for(UInt iG=0; iG<nG; iG++) {
			if (iAtomSize >= atomGridDistance[iG]) {
				isSigAtom = true;
				break;
			}
		}

		// now we can have a conclustion
		if (isSigAtom) {
			nSigAtoms++;
			nSigBasis += atomShl.getNBas();
			nSigShell += atomShl.getNShell();
			neighbors.push_back(iAtom);
		}
	}

	// now initilize the data
	// we note that this is only a rough estimation
	sigAtoms.reserve(nSigAtoms);
	sigBas.reserve(nSigBasis); 
	sigShell.reserve(nSigShell); 
	atomBasisMapping.reserve(nSigAtoms+1);
	atomShellMapping.reserve(nSigAtoms+1);

	// here "0" means the first index in the sigBas and sigShell
	atomBasisMapping.push_back(0);
	atomShellMapping.push_back(0);

	// real working part
	// we starts from the neighbor list to a fine Calibration
	nSigBasis = 0;
	nSigShell = 0;
	for(UInt i=0; i<neighbors.size(); i++) {

		// get the atom index
		// the atom in the neighbor list are all significant atoms
		// based on the procedure we made above
		UInt iAtom = neighbors[i];
		sigAtoms.push_back(iAtom);

		// now let's calculate the atom-grid distance
		// again, for evaluating the shell
		const AtomShell& atomShl = ms.getAtomShell(iAtom);
		const Double* coord = grid.getGridCoord();
		for(UInt iG=0; iG<nG; iG++) {
			atomGridDistance[iG] = atomShl.getDistance(&coord[3*iG]);
		}

		// counting shells within this atom
		UInt atomic = atomShl.getAtomic();
		const AtomShellSize& atomShellSize = shellSize.getAtomShellSizeInfor(atomic);
		for(UInt iShell=0; iShell<atomShl.getNShell(); iShell++) {

			// is it a sig shell
			Double shlSize = atomShellSize.getShellRadius(iShell);
			bool isSigShell = false;
			for(UInt iG=0; iG<nG; iG++) {
				if (shlSize >= atomGridDistance[iG]) {
					isSigShell = true;
					break;
				}
			}

			// now add in data
			if (isSigShell) {
				const Shell& s = atomShl.getShell(iShell);

				// form basis set information
				UInt basOffset = s.getBasisIndex(0,TYPE_NORM);
				for(UInt i=0; i<s.getNBas(); i++) {
					UInt gIndex = basOffset + i;
					sigBas.push_back(gIndex);
					nSigBasis++;
				}

				// form shell information
				sigShell.push_back(s.getLocalShellIndex());
				nSigShell++;
			}
		}

		// now deal with the mapping array
		atomBasisMapping.push_back(nSigBasis);
		atomShellMapping.push_back(nSigShell);
	}
}

SigAtomBasis::SigAtomBasis(const BatchGrid& grid, const XCIntJobInfor& infor,
		const MolShell& ms, const MolShellSize& atomShlSize): 
	withAtomApprox(false)
{
	if (withAtomApprox) {
		formSigDataWithAtomCenterApprox(grid,ms);
	}else{
		formSigData(grid,ms,atomShlSize); 
	}
}

void SigAtomBasis::print() const {

	cout << "*******************************************" << endl;
	cout << "*        SigAtomBasis Result              *" << endl;
	cout << "*******************************************" << endl;

	// print out the significant array information
	cout << "The significant atoms and shells index" << endl;
	cout << "The basis set offset is in normal type order" << endl;
	for(UInt iAtom=0; iAtom<getNSigAtoms(); iAtom++) {
		cout << "Atom " << sigAtoms[iAtom] << " is significant; " << endl;
		UInt startIndex = getSigShellBeginIndex(iAtom);
		for(UInt iShell=0; iShell<getNSigShell(iAtom); iShell++) {
			UInt index = startIndex + iShell;
			cout << "for this atom the local shell " << sigShell[index] << " is significant; " << endl;
		}
	}
	cout << endl;
	for(UInt iBas=0; iBas<getNSigBasis(); iBas++) {
		cout << "Basis set " << sigBas[iBas] << " is significant; " << endl;
	}
	cout << endl;
}

bool SigAtomBasis::isSigAtoms(const UInt& iAtom) const {
	UIntVec::const_iterator it;
	it = find(sigAtoms.begin(), sigAtoms.end(), iAtom);
	if (it == sigAtoms.end()) return false;
	return true;
}

bool SigAtomBasis::isSigShell(const UInt& basOffset) const {
	UIntVec::const_iterator it;
	it = find(sigBas.begin(), sigBas.end(), basOffset);
	if (it == sigBas.end()) {
		return false;
	}
	return true;
}

UInt SigAtomBasis::getSigBasisIndex(const UInt& iBasis) const {
	UInt iSigBasis = getNSigBasis() + 1;
	for(UInt iBas=0; iBas<getNSigBasis(); iBas++) {
		if (sigBas[iBas] == iBasis) {
			iSigBasis = iBas;
			break;
		}
	}
	if (iSigBasis>getNSigBasis()) {
		string infor = "fail to find the corresponding sigbasis index for the given normal basis set: " + 
			boost::lexical_cast<string>(iBasis);
		Excep excep("SigAtomBasis","getSigBasisIndex",EXCEPTION_XCINTS_INVALID_SIGBASIS_INDEX,infor);
		handleExcep(excep);
	}
	return iSigBasis;
}

UInt SigAtomBasis::getLocalAtomBasisIndex(const AtomShell& as, const UInt& iSigAtom, 
		const UInt& iSigBasis) const
{

	// we need to get the basis set offset for the whole atom
	UInt basOffset = as.getBasisStartIndex(TYPE_NORM);

	// now let's try to get the global basis set index for the 
	// given sig atom and sig basis 
	UInt sigBasIndex = getGlobalBasisIndex(iSigAtom,iSigBasis);

	// additional check 
	// make sure that the given sigBasIndex is within the given atom shell
	if (sigBasIndex < basOffset || sigBasIndex >= basOffset + as.getNBas()) {
		string infor;
		if (sigBasIndex < basOffset) {
			string info1 = "sigBasisIndex < basis offset: ";
			string info2 = boost::lexical_cast<string>(sigBasIndex);
			string info3 = boost::lexical_cast<string>(basOffset);
			infor = info1 + info2 + " " + info3;
		}else{
			string info1 = "sigBasisIndex >= basis offset+as.getNBas(): ";
			string info2 = boost::lexical_cast<string>(sigBasIndex);
			string info3 = boost::lexical_cast<string>(basOffset+as.getNBas());
			infor = info1 + info2 + " " + info3;
		}
		Excep excep("SigAtomBasis","getLocalAtomBasisIndex",EXCEPTION_XCINTS_INVALID_SIGBASIS_INDEX,infor);
		handleExcep(excep);
	}

	// the local index would be the difference between them
	return sigBasIndex - basOffset;
}

