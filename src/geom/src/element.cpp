#include "element.h"
#include "excep.h"
#include "textread.h"
#include <boost/lexical_cast.hpp>
using namespace textread;
using namespace excep;
using namespace element;

string element::getAtomicSymbol(UInt n) 
{
#ifdef DEBUG
	if (n>MAX_ATOMS_TYPE) {
		string infor = "Atomic number is " + boost::lexical_cast<string>(n);
		Excep excep("element","getAtomicSymbol",INVALID_ATOMIC_NUMBER,infor);
		handleExcep(excep);
	}
#endif
	string s = SYMBOLS[n];
	return s;
}

UInt element::getAtomicNumber(string s)
{
	WordConvert w;
	for (UInt i=0;i<MAX_ATOMS_TYPE;i++) {
		if (w.compare(SYMBOLS[i],s)) {
			return i;
		}
	}
	// in the abnormal case, we return -1
	// this is the maximum value of size_t
	// then it will cause error in the following step
	return -1;
}

bool element::isGhostAtom(const UInt& Z) {
	if (Z == 0) return true;
	return false;
}

bool element::isAtomicSymbol(string& s){
	WordConvert w;
	for (UInt i=0; i<MAX_ATOMS_TYPE; i++) {
		if (w.compare(SYMBOLS[i],s)) {
			s = SYMBOLS[i]; 
			return true;
		}
	}
	return false;
}

Double element::getAtomRadii(const UInt& n){
#ifdef DEBUG
	if (n>MAX_ATOMS_TYPE) {
		string infor = "Atomic number is " + boost::lexical_cast<string>(n);
		Excep excep("element","getAtomRadii",INVALID_ATOMIC_NUMBER,infor);
		handleExcep(excep);
	}
#endif
	return ATOM_RADII[n]*ANGSTROM_TO_BOHR;
}

UInt element::getGroundStateAtomMult(const UInt& n){
#ifdef DEBUG
	if (n>MAX_ATOMS_TYPE) {
		string infor = "Atomic number is " + boost::lexical_cast<string>(n);
		Excep excep("element","getGroundStateAtomMult",INVALID_ATOMIC_NUMBER,infor);
		handleExcep(excep);
	}
#endif
	return GROUND_STATE_ATOM_MULTIPLICITY[n]; 
}

Double element::getAtomPolarizability(const UInt& n) {
#ifdef DEBUG
	if (n>MAX_ATOMS_TYPE) {
		string infor = "Atomic number is " + boost::lexical_cast<string>(n);
		Excep excep("element","getAtomPolarizability",INVALID_ATOMIC_NUMBER,infor);
		handleExcep(excep);
	}
#endif
	Double data = GROUND_STATE_ATOM_POLARIZABILITIES[n];
	if (data < 0) {
		string infor = "Atomic number is " + boost::lexical_cast<string>(n);
		Excep excep("element","getAtomPolarizability",EXCEPTION_POLARIZABILITIES,infor);
		handleExcep(excep);
		return ZERO;
	}else{
		return data*ANGSTROM_TO_BOHR*ANGSTROM_TO_BOHR*ANGSTROM_TO_BOHR;
	}
}


