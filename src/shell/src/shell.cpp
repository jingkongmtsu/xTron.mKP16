/**
 * CPP files corresponding to the shell.h
 * \author  Fenglai liu and Jing kong
 */
//#include<math.h>
#include<iostream>
#include<algorithm>
#include<cstdio>
#include "excep.h"
#include "shellprop.h"
#include "parseshell.h"
#include "element.h"
#include "molecule.h"
#include "shell.h"
#include "boost/lexical_cast.hpp"
#ifdef SHELL_PHI_DEBUG
#include "atomshelltypedatainfor.h"
#endif
using namespace excep;
using namespace shellprop;
using namespace parseshell;
using namespace element;
using namespace molecule;
#ifdef SHELL_PHI_DEBUG
using namespace atomshelltypedatainfor; 
#endif
using namespace shell;

/***********************************************************************************
 *                      functions related to the shell class                       *
 ***********************************************************************************/
Shell::Shell(const RawShell& sr):isDummy(false),localShellIndex(-1),globalShellIndex(-1),
	normBasStart(-1),cartBasStart(-1),localNormBasStart(-1),localCartBasStart(-1),
	nBas(-1),nCarBas(-1),lmin(-1),lmax(-1),nPrim(0)
{
	// determine the shell type etc.
	string type = sr.shellType;
	getL(type,lmin,lmax);
	nCarBas   = getCartBas(lmin,lmax);
	nBas      = nCarBas;
	bool pure = sr.isPure;
	UInt nPurBas   = getPureBas(lmin,lmax);
	if (pure) nBas = nPurBas;

	// stuff the exponents and coefficients for this shell
	nPrim = sr.ngau;
	e.assign(nPrim,ZERO);
	for(UInt i=0;i<nPrim;i++){
		Double exp = sr.exp[i];
		Double s   = pow(sr.scale, TWO);// usually it's 1
		exp        = exp*s; 
		e[i]       = exp;
	}

	// for the coefficient array, we have to do more treatment. The
	// treatment below is that each Gaussian primitives also need
	// normalized! the final normalization step will be done later
	c.assign(nPrim*(lmax-lmin+1),ZERO);
	for(UInt iShell=lmin;iShell<=lmax;iShell++) {

		// get the correct index for the coe array
		UInt np = (iShell-lmin)*nPrim; 

		// stuff the whole coe array
		// also proceed with normalization for primitive function
		UInt offset = nPrim*(iShell-lmin);
		for(UInt i=0; i<nPrim; i++){
			Double coe  = sr.coe[np+i];
			Double s    = pow(e[i],(iShell+THREE/TWO)/TWO);
			c[offset+i] = coe*s;
		}
	}

	// normalization for basis set
	// here we assume that all of basis sets are in type of x^Ly^0z^0
	// later in the transformation calculation we will scale it back
	// see shell transform and scale part of code
	for(UInt iShell=lmin;iShell<=lmax;iShell++) {
		UInt np = (iShell-lmin)*nPrim; 
		Double norm = normCartBasisSets(iShell, 0, 0, nPrim, &c[np], &e[0]);
		for(UInt i=0; i<nPrim; i++){
			c[np+i] = c[np+i]*norm;
		}
	}
}	

void Shell::formDummyShell()
{
	// set the data
	// dummy shell is still a shell: S shell!!
	isDummy           = true;
	localShellIndex   = 0;
	globalShellIndex  = 0;
	normBasStart      = 0;
	cartBasStart      = 0;
	localNormBasStart = 0;
	localCartBasStart = 0;
	nBas              = 1;
	nCarBas           = 1;
	lmin              = 0;
	lmax              = 0;
	nPrim             = 1;

	// in default, the coefficient is 1 and the exponent is 0
	c.assign(1,ONE);
	e.assign(1,ZERO);
}

UInt Shell::getBasisIndex(const UInt& iBas, const UInt& type) const 
{
	switch (type) {
		case TYPE_CART:
			return cartBasStart+iBas;
			break;
		case TYPE_NORM:
			return normBasStart+iBas;
			break;
		default:
#ifdef DEBUG
			string infor = " "; 
			Excep excep("Shell","getBasisIndex",
					EXCEPTION_SHELL_BASIS_INDEX_TYPE_INVALID,infor);
			handleExcep(excep);
#endif
			return -1; 
			break;
	}
}

UInt Shell::getLocalBasisIndex(const UInt& iBas, const UInt& type) const 
{
	switch (type) {
		case TYPE_CART:
			return localCartBasStart+iBas;
			break;
		case TYPE_NORM:
			return localNormBasStart+iBas;
			break;
		default:
#ifdef DEBUG
			string infor = " "; 
			Excep excep("Shell","getLocalBasisIndex",
					EXCEPTION_SHELL_BASIS_INDEX_TYPE_INVALID,infor);
			handleExcep(excep);
#endif
			return -1; 
			break;
	}
}

#ifdef SHELL_PHI_DEBUG
void Shell::resetShellData(const UInt& lmin0, const UInt& lmax0, const UInt& shellStatus, 
		const UInt& np, const Double* coef, const Double* primExp)
{
	// fill in the simple data
	isDummy = false;    
	lmin    = lmin0;             
	lmax    = lmax0;             
	nPrim   = np;            
	nCarBas = getCartBas(lmin,lmax);
	UInt nPurBas = getPureBas(lmin,lmax);
	nBas    = nPurBas;
	if (shellStatus == CART_SHELL_TYPE) nBas = nCarBas;

	// now fill in array
	UInt lenCoe = nPrim*(lmax-lmin+1);
	for(UInt i=0; i<lenCoe; i++) {
		c[i] = coef[i];
	}		
	for(UInt i=0; i<nPrim; i++) {
		e[i] = primExp[i]; 
	}
}
#endif

void Shell::print(UInt iprint) const {

	bool printmore = false;
	if (iprint >=1) {
		printmore = true;
	}

	// for dummy shell, we do not print
	if (isDummyShell()) {
		cout << "This is dummy shell so we do not print anything" << endl;
	}

	string type = getShellName(lmin,lmax);
	const char* usePure = "CART";
	if (lmax>1 && isPure()) usePure = "PURE";
	printf ("%-3s  %-4s\n", type.c_str(), usePure);
	cout << "local  shell index is : " << localShellIndex << endl; 
	cout << "global shell index is : " << globalShellIndex << endl; 
	cout << "Normal    basis set starting from: " << normBasStart << " to " << 
		normBasStart + nBas - 1 << endl;
	cout << "Cartesian basis set starting from: " << cartBasStart << " to " << 
		cartBasStart + nCarBas - 1<< endl;
	cout << "local normal    basis set starting from: " << localNormBasStart << " to " << 
		localNormBasStart + nBas - 1 << endl;
	cout << "local cartesian basis set starting from: " << localCartBasStart << " to " << 
		localCartBasStart + nCarBas - 1<< endl;
	if (printmore) {
		UInt nPrim = getNPrim();
		for (UInt k=0;k<nPrim;k++) {
			UInt np = nPrim;
			if (type == "SP"){
				printf ("%-14.8f  %-14.8f  %-14.8f  \n", e[k], c[k], c[np+k]);
			}else{
				printf ("%-14.8f  %-14.8f  \n", e[k], c[k]);
			}
		}
		cout << "****" << endl;
	}
}

bool Shell::sameShellData(const Shell& s) const
{
	// test all of data fields
	if( isDummy != s.isDummy) return false;          
	if( localShellIndex != s.localShellIndex) return false;  
	if( globalShellIndex != s.globalShellIndex) return false; 
	if( normBasStart != s.normBasStart) return false;     
	if( cartBasStart != s.cartBasStart) return false;     
	if( localNormBasStart != s.localNormBasStart) return false;
	if( localCartBasStart != s.localCartBasStart) return false;
	if( nBas != s.nBas) return false;             
	if( nCarBas != s.nCarBas) return false;          
	if( lmin != s.lmin) return false;             
	if( lmax != s.lmax) return false;             
	if( nPrim != s.nPrim) return false;            

	// finally test the array
	for (UInt i=0; i<nPrim; i++) {
		if (fabs(e[i]-s.e[i])>THRESHOLD_MATH) return false;
	}
	for (UInt i=0; i<nPrim*(lmax-lmin+1); i++) {
		if (fabs(c[i]-s.c[i])>THRESHOLD_MATH) return false;
	}

	// now everything is good
	return true;
}

/***********************************************************************************
 *                      functions related to the shells class                      *
 ***********************************************************************************/
AtomShell::AtomShell(const RawAtomShell& shls, const Atom& atom, bool splitSPShell):Atom(atom),
	allCartesian(true),allSpherical(true),nBas(0),nCarBas(0),
	normBasStart(-1),cartBasStart(-1),maxL(0),maxNP(0),shlIndexStart(-1),nShells(0)
{
	// let's estimate the number of shells may need
	UInt ns = 0;
	for(UInt iShell=0; iShell<shls.getNShell(); iShell++) {
		const RawShell& sr = shls.getShell(iShell);
		if (sr.shellType == "SP" && splitSPShell) {
			ns = ns + 2;
		}else{
			ns = ns + 1;
		}
	}
	shells.reserve(ns);

	// fill the shell data
	for(UInt iShell=0; iShell<shls.getNShell(); iShell++) {
		const RawShell& sr = shls.getShell(iShell);

		// do we need to split the SP shell data?
		if (sr.shellType == "SP" && splitSPShell) {
			for(UInt i=0; i<2; i++) {

				// set a new raw shell data
				string shelltype = "S";
				if (i == 1) shelltype = "P";
				RawShell newSR(sr.isPure,sr.scale,1,sr.ngau,shelltype);
				newSR.splitSPShellData(sr);

				// now read in this new rawshell data
				Shell s(newSR);
				nBas     += s.getNBas();
				nCarBas  += s.getNCarBas();
				nShells  += 1;
				shells.push_back(s);
				if (s.getLmax() > maxL) maxL = s.getLmax();
				if (s.getNPrim() > maxNP) maxNP = s.getNPrim();
			}
		}else{
			Shell s(sr);
			nBas     += s.getNBas();
			nCarBas  += s.getNCarBas();
			nShells  += 1;
			shells.push_back(s);
			if (s.getLmax() > maxL) maxL = s.getLmax();
			if (s.getNPrim() > maxNP) maxNP = s.getNPrim();
			// here we have to note, that S,P shell is considered as
			// cart. So all pure or all cart judgement should
			// be made for shells which is larger or equal D
			if (s.getLmax() >=2 && s.isPure()) {
				allCartesian = false;
			}
			if (s.getLmax() >=2 && ! s.isPure()) {
				allSpherical = false;
			}
		}
	}

	// finally we are able to form the local index/offset data for atom shell
	formLocalShellOffset(); 
}

void AtomShell::formDummyAtomShell()
{
	// however, we need to keep on eye 
	// that whether the shell data is 
	// empty
	if (! empty()) {
		string infor = "  "; 
		Excep excep("AtomShell","formDummyShell",EXCEPTION_SHELL_NOT_EMPTY,infor);
		handleExcep(excep);
	}

	// append a dummy shell to this atom shell
	// this dummy shell, is a S shell
	// we also note, that here we need to update
	// atom index information - it's one atom,
	// which contains a dummy shell!
	allCartesian = true;
	allSpherical = false;
	nBas         = 1;
	nCarBas      = 1;
	normBasStart = 0;
	cartBasStart = 0;
	maxL         = 0;
	maxNP        = 1;
	nShells      = 1;
	shlIndexStart= 0;
	updateAtomIndex(0);

	// now let's append a dummy shell
	Shell dummy;
	dummy.formDummyShell();
	shells.reserve(1);
	shells.push_back(dummy);
}

void AtomShell::formShellOffset(UInt basOffset, UInt carBasOffset) 
{
	normBasStart  = basOffset; 
	cartBasStart  = carBasOffset; 
	for(UInt iShell=0; iShell<getNShell(); iShell++) {
		shells[iShell].formShellOffset(basOffset,carBasOffset);
		basOffset    += shells[iShell].getNBas();
		carBasOffset += shells[iShell].getNCarBas();
	}
}

void AtomShell::formLocalShellOffset() 
{
	UInt localNormBasStart = 0;
	UInt localCartBasStart = 0;
	for(UInt iShell=0; iShell<getNShell(); iShell++) {
		shells[iShell].formLocalShellOffset(localNormBasStart,localCartBasStart);
		localNormBasStart += shells[iShell].getNBas();
		localCartBasStart += shells[iShell].getNCarBas();
	}
}

void AtomShell::formShellIndex() 
{
	// check the shell starting index
	// whether it has been properly initilized
	if (shlIndexStart == static_cast<UInt>(-1)) {
		string infor = "the starting global shell index is not initialized"; 
		Excep excep("MolShell","formShellIndex",EXCEPTION_INVALID_SHELL_DATA,infor);
		handleExcep(excep);
	}

	// now do index forming
	for(UInt iShell=0; iShell<getNShell(); iShell++) {
		shells[iShell].formShellIndex(shlIndexStart,iShell);
	}
}

bool AtomShell::hasHighL() const
{
	for(UInt iShell=0; iShell<getNShell(); iShell++) {
		const Shell& s = getShell(iShell);
		if (s.getLmax() >= 2) return true;
	}
	return false;
}

UInt AtomShell::getBasisStartIndex(const UInt& type) const 
{
	switch (type) {
		case TYPE_CART:
			return cartBasStart;
			break;
		case TYPE_NORM:
			return normBasStart;
			break;
		default:
#ifdef DEBUG
			string infor = " "; 
			Excep excep("AtomShell","getBasisStartIndex",
					EXCEPTION_SHELL_BASIS_INDEX_TYPE_INVALID,infor);
			handleExcep(excep);
#endif
			return -1; 
			break;
	}
}

UInt AtomShell::getMaxCompDegree() const
{
	UInt degree = 1;
	UInt nShell = getNShell();
	for(UInt iShell=0; iShell<nShell; iShell++) {
		const Shell& s = getShell(iShell);
		UInt lmin = s.getLmin();
		UInt lmax = s.getLmax();
		UInt comp = lmax-lmin+1;
		if (comp>degree) degree = comp;
	}
	return degree;
}

#ifdef SHELL_PHI_DEBUG
void AtomShell::formAtomShellFromCompactForm(const UInt& atomShellIndex, const Double* coord, 
		const AtomShellTypeDataInfor& compactInfor, UInt cartBasOffset, UInt normBasOffset, UInt shellOffset)
{
	// fill the atom information
	atomIndex = atomShellIndex;
	xyz[0]    = coord[0]; 
	xyz[1]    = coord[1]; 
	xyz[2]    = coord[2]; 
	atomic    = atomshelltypedatainfor::getAtomic(compactInfor);

	// reset all of data
	allCartesian = true; 
	allSpherical = true;
	nBas         = 0;
	nCarBas      = 0;
	maxL         = 0;
	maxNP        = 0;
	nShells      = atomshelltypedatainfor::getNShell(compactInfor);

	// now let's filling in the shell information
	for(UInt iShell=0; iShell<nShells; iShell++) {

		// form the given shell
		UInt lmin0, lmax0;
		getLInfor(compactInfor,iShell,lmin0,lmax0);
		UInt shellStatus   = getShellStatus(compactInfor,iShell);
		UInt np            = getNPrim(compactInfor,iShell);
		const Double* coef = getCoe(compactInfor,iShell);
		const Double* pexp = getExp(compactInfor,iShell);
		Shell& s = shells[iShell];
		s.resetShellData(lmin0,lmax0,shellStatus,np,coef,pexp);

		// update the data
		nBas     += s.getNBas();
		nCarBas  += s.getNCarBas();
		if (s.getLmax() > maxL) maxL = s.getLmax();
		if (s.getNPrim() > maxNP) maxNP = s.getNPrim();
		// here we have to note, that S,P shell is considered as
		// cart. So all pure or all cart judgement should
		// be made for shells which is larger or equal D
		if (s.getLmax() >=2 && s.isPure()) {
			allCartesian = false;
		}
		if (s.getLmax() >=2 && ! s.isPure()) {
			allSpherical = false;
		}
	}

	// form the local index/offset data for atom shell
	formLocalShellOffset(); 

	// now form the global basis set offset data
	// set the entry data for this atom shell
	UInt basOffset    = normBasOffset; 
	UInt carBasOffset = cartBasOffset;
	formShellOffset(basOffset,carBasOffset);	

	// finally make the shell offset
	shlIndexStart = shellOffset;
	formShellIndex();
}
#endif

void AtomShell::print(UInt iprint) const 
{
	if (isDummyAtomShell()) {
		cout << "This is dummy atom shell, nothing printed here " << endl;
	}

	string sym = getAtomicSymbol(atomic);
	printf ("%-3s  %-14.7f  %-14.7f  %-14.7f  \n",sym.c_str(), 
			xyz[0], xyz[1], xyz[2]);
	cout << "global shell starting index " << shlIndexStart << endl;
	for (UInt i=0;i<getNShell();i++) {
		const Shell& s = getShell(i);
		s.print(iprint);
	}	
}

bool AtomShell::sameAtomShell(const AtomShell& as) const
{
	// test the simple data
	if(allCartesian  != as.allCartesian) return false; 
	if(allSpherical  != as.allSpherical) return false; 
	if(nBas          != as.nBas) return false;         
	if(nCarBas       != as.nCarBas) return false;      
	if(normBasStart  != as.normBasStart) return false;
	if(cartBasStart  != as.cartBasStart) return false; 
	if(maxL          != as.maxL) return false;         
	if(maxNP         != as.maxNP) return false;        
	if(shlIndexStart != as.shlIndexStart) return false;
	if(nShells       != as.nShells) return false;      

	// now let me test each shell
	for (UInt i=0;i<getNShell();i++) {
		const Shell& s = as.getShell(i);
		if (! shells[i].sameShellData(s)) return false;
	}

	// now everything is good
	return true;
}

/***********************************************************************************
 *                    functions related to the MolShell class                      *
 ***********************************************************************************/
void MolShell::formMolShellData(const Molecule& mol, const RawMolShell& rms)
{
	//
	// first step, we need to check whether mol and rms
	// are consistent with each other
	// that is to say, they must have same number of atoms
	//
	UInt nAtomTypes = mol.getNAtomTypes();
	if (nAtomTypes != rms.getNAtomShell()) {
		string infor = "the RawMolShell returns different atom number from input molecule"; 
		Excep excep("MolShell","formMolShellData",
				EXCEPTION_RAW_SHELL_DATA_CONFLICT_MOL_DATA,infor);
		handleExcep(excep);
	}

	// 
	// also let's check whether we need to split she SP shells?
	//
	bool needSplitSPShell = false;
	if (rms.hasHighAngShellData() && rms.hasSPShellData()) {
		needSplitSPShell = true;
	}

	// now let's fill in data
	UInt nAtoms = mol.getNAtoms();
	atomShells.reserve(nAtoms);	
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const UInt atomic = mol.getAtomic(iAtom);
		const Atom& atom  = mol.getAtom(iAtom);
		if (isGhostAtom(atomic)) continue;
		AtomShell atomShl(rms.getAtomShell(atomic),atom,needSplitSPShell);
		nBas     += atomShl.getNBas();
		nCarBas  += atomShl.getNCarBas();
		nShell   += atomShl.getNShell();
		atomShells.push_back(atomShl);
		if (atomShl.getMaxL() > maxL) maxL = atomShl.getMaxL();
		if (atomShl.getMaxNP() > maxNP) maxNP = atomShl.getMaxNP();
		if (! atomShl.allPure()) {
			allSpherical = false;
		}
		if (! atomShl.allCart()) {
			allCartesian = false;
		}
	}
}

void MolShell::formShellIndex() {

	// form the global basis set index
	// and the shell index
	UInt normBasCount = 0;
	UInt cartBasCount = 0;
	UInt iShell = 0;
	for(UInt iAtom=0; iAtom<getNAtomShells(); iAtom++) {
		atomShells[iAtom].formShellOffset(normBasCount,cartBasCount);
		atomShells[iAtom].updateShellStartIndex(iShell);
		atomShells[iAtom].formShellIndex();
		normBasCount += atomShells[iAtom].getNBas();
		cartBasCount += atomShells[iAtom].getNCarBas();
		iShell += atomShells[iAtom].getNShell();
	}

	// in safe, when call this function we always do the check
	if (normBasCount != nBas) {
		string infor = "Invalid normal basis set shell index forming in MolShell"; 
		Excep excep("MolShell","formShellIndex",EXCEPTION_SHELL_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
	if (cartBasCount != nCarBas) {
		string infor = "Invalid cartesian basis set shell index forming in MolShell"; 
		Excep excep("MolShell","formShellIndex",EXCEPTION_SHELL_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}
}

void MolShell::appendAtomShell(const AtomShell& atomShl)
{
	// appending data, and refresh the data too
	atomShells.push_back(atomShl);
	nBas     += atomShl.getNBas();
	nCarBas  += atomShl.getNCarBas();
	nShell   += atomShl.getNShell();
	if (atomShl.getMaxL() > maxL) maxL = atomShl.getMaxL();
	if (atomShl.getMaxNP() > maxNP) maxNP = atomShl.getMaxNP();
	if (! atomShl.allPure()) {
		allSpherical = false;
	}
	if (! atomShl.allCart()) {
		allCartesian = false;
	}

	// we need to see whether the angular momentum is supported in this program
	if (maxL > MAX_L) {
		printf("the maximum angular momentum supported in this program is %d\n", (Int)MAX_L);
		string infor = "in the shell data you used unsupported angular momentum shell data"; 
		Excep excep("MolShell","appendAtomShell",EXCEPTION_INVALID_SHELL_DATA,infor);
		handleExcep(excep);
	}
}

MolShell::MolShell(const string& input, const Molecule& mol, UInt part, bool debugRawShellData):
	allCartesian(true), allSpherical(true),nShell(0),nBas(0),nCarBas(0),
	maxL(0),maxNP(0),code(-1)
{
	// form the shell code and form raw shell data
	UInt section = mol.getSec();
	code = codeShell(section,part);
	RawMolShell rms(input,mol,part,debugRawShellData);

	// form the atom shell data
	formMolShellData(mol,rms);

	// form the shell index
	formShellIndex();

	// we need to see whether the angular momentum is supported in this program
	if (maxL > MAX_L) {
		printf("the maximum angular momentum supported in this program is %d\n", (Int)MAX_L);
		string infor = "in the shell data you used unsupported angular momentum shell data"; 
		Excep excep("MolShell","constructor",EXCEPTION_INVALID_SHELL_DATA,infor);
		handleExcep(excep);
	}
}

MolShell::MolShell(const AtomShell& as):allCartesian(true), allSpherical(true),
	nShell(0),nBas(0),nCarBas(0),maxL(0),maxNP(0),code(-1)
{
	// append this atom shell
	atomShells.reserve(1);
	appendAtomShell(as);

	// we need to correct the 
	// geometry data for this one atom shell
	// it should be on origin
	AtomShell& atomShell = atomShells[0];
	Double X = ZERO;
	Double Y = ZERO;
	Double Z = ZERO;
	atomShell.updateXYZ(X,Y,Z);

	// also update the atom index for the atom shell
	// it should be 0
	atomShell.updateAtomIndex(0);

	// form the shell index
	formShellIndex();
}

void MolShell::formDummyShell()
{
	// however, we need to keep on eye 
	// that whether the shell data is 
	// empty
	if (! empty()) {
		string infor = "  "; 
		Excep excep("MolShell","formDummyShell",EXCEPTION_SHELL_NOT_EMPTY,infor);
		handleExcep(excep);
	}

	// fill in data
	AtomShell as;
	as.formDummyAtomShell();
	atomShells.push_back(as);
	nShell = 1;
	nBas   = 1;
	nCarBas= 1;
	maxL   = 0;
	maxNP  = 1;
	code   = -1;
}

UInt MolShell::getNAtomShellTypes() const
{
	UInt atomTypesData[MAX_ATOMS_TYPE];
	for(UInt iAtom=0; iAtom<MAX_ATOMS_TYPE; iAtom++) {
		atomTypesData[iAtom] = 0;
	}

	UInt nAtoms = getNAtomShells();
	UInt nTypes = 0;
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomShell& atomShell = getAtomShell(iAtom);
		UInt atomic = atomShell.getAtomic();
		if (atomTypesData[atomic] > 0) continue;
		nTypes++;
		atomTypesData[atomic] = 1;
	}   
	return nTypes;
}

bool MolShell::hasAng(const UInt& lmin, const UInt& lmax) const 
{

   // firstly, we need to examine that whether we have the shell of (lmin,lmax)
   // in the given molshell?
   if (maxL < lmax) return false;

   // test for dummy shells
   if (dummy()) {
      if (lmax == lmin &&  lmin == 0) return true;
      return false;
   }   

   // then let's do the real job
	// initilize the atom type data
	// this is used to make the whole process faster
	UInt atomTypesData[MAX_ATOMS_TYPE];
	for(UInt iAtom=0; iAtom<MAX_ATOMS_TYPE; iAtom++) {
		atomTypesData[iAtom] = 0;
	}

	// now let's go to see
	UInt nAtoms = getNAtomShells();
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomShell& atomShell = getAtomShell(iAtom);
		UInt atomic = atomShell.getAtomic();
		if (atomTypesData[atomic] > 0) continue;
		for(UInt iShell=0; iShell<atomShell.getNShell(); iShell++) {
			const Shell& s = atomShell.getShell(iShell);
			if (s.getLmin() == lmin && s.getLmax() == lmax) return true;
		}   
		atomTypesData[atomic] = 1;
	}   
   return false;
}

UInt MolShell::getMaxBasis(const UInt& type) const 
{
	// debug check
#ifdef DEBUG
	if (type != TYPE_CART && type != TYPE_NORM) {
		string infor = "only TYPE_CART, TYPE_NORM and TYPE_PURE are allowed"; 
		Excep excep("MolShell","getMaxBasis",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}
#endif

   // then let's do the real job
	// initilize the atom type data
	UInt atomTypesData[MAX_ATOMS_TYPE];
	for(UInt iAtom=0; iAtom<MAX_ATOMS_TYPE; iAtom++) {
		atomTypesData[iAtom] = 0;
	}

	// now let's go to see
	UInt nAtoms = getNAtomShells();
	UInt nBasis = 0;
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomShell& atomShell = getAtomShell(iAtom);
		UInt atomic = atomShell.getAtomic();
		if (atomTypesData[atomic] > 0) continue;
		if (atomShell.getNBas(type)>nBasis) {
			nBasis = atomShell.getNBas(type);
		}
		atomTypesData[atomic] = 1;
	}   
   return nBasis;
}

UInt MolShell::getMaxNShell() const
{
   // then let's do the real job
	// initilize the atom type data
	UInt atomTypesData[MAX_ATOMS_TYPE];
	for(UInt iAtom=0; iAtom<MAX_ATOMS_TYPE; iAtom++) {
		atomTypesData[iAtom] = 0;
	}

	// now let's go to see
	UInt nAtoms  = getNAtomShells();
	UInt nMaxShells = 0;
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomShell& atomShell = getAtomShell(iAtom);
		UInt atomic = atomShell.getAtomic();
		if (atomTypesData[atomic] > 0) continue;
		if (atomShell.getNShell()>nMaxShells) {
			nMaxShells = atomShell.getNShell();
		}
		atomTypesData[atomic] = 1;
	}   
   return nMaxShells;
}

UInt MolShell::getMaxLContraction() const 
{

   // then let's do the real job
	// initilize the atom type data
	UInt atomTypesData[MAX_ATOMS_TYPE];
	for(UInt iAtom=0; iAtom<MAX_ATOMS_TYPE; iAtom++) {
		atomTypesData[iAtom] = 0;
	}

	// now let's go to see
	UInt nAtoms = getNAtomShells();
	UInt maxNL  = 1;
	for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {
		const AtomShell& atomShell = getAtomShell(iAtom);
		UInt atomic = atomShell.getAtomic();
		if (atomTypesData[atomic] > 0) continue;
		UInt nL = atomShell.getMaxCompDegree();
		if (nL>maxNL) maxNL = nL;
		atomTypesData[atomic] = 1;
	}   
   return maxNL;
}

void MolShell::print(UInt iprint) const {

	cout << endl;
	cout << "*******************************************" << endl;
	cout << "*          MolShell Object Print          *" << endl;
	cout << "*******************************************" << endl;
	if (dummy()) {
		cout << "This is dummy shell, no printing here" << endl;
	}

	// print out the basic information
	cout << "Number of shells is: " << getNShell() << endl;
	cout << "Number of basis set functions is: " << getNBas() << endl;
	cout << "Number of Cartesian basis set functions is: " << getNCarBas() << endl;
	cout << "Maximum angular momentum is: " << maxL << "  Maximum contraction degree is: "
		<< maxNP << endl;
	cout << "Maximum Cartesian basis set for all of atoms: " << getMaxBasis(TYPE_CART) << endl;
	cout << "Maximum L contraction for all of atoms: " << getMaxLContraction() << endl;
	if (iprint >=1) {
		cout << "Detailed Shell Infor" << endl;
		for (UInt i=0;i<getNAtomShells();i++) {
			const AtomShell& as = getAtomShell(i);
			as.print(iprint);
		}	
	}
	cout << endl;
}


