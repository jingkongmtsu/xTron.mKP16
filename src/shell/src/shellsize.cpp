#include<cstdio>
#include<string>
#include<cmath>
#include<iostream>
#include "excep.h"
#include "shell.h"
#include "element.h"
#include "boost/lexical_cast.hpp"
#include "shellsize.h"
using namespace std;
using namespace excep;
using namespace element;
using namespace shell;
using namespace shellsize;

AtomShellSize::AtomShellSize(const AtomShell& s, const Double& thresh):atomic(s.getAtomic()),
  atomRadius(ZERO)	
{
	UInt nShell = s.getNShell();
	Double atomRad = ZERO;
	shellRadius.reserve(nShell);
	for(UInt iShell=0; iShell<nShell; iShell++) {
		const Shell& is = s.getShell(iShell);
		Double size = setShellRadius(is,thresh);
		shellRadius.push_back(size);
		if (size > atomRad) atomRad = size;
	}
	atomRadius = atomRad;
}

Double AtomShellSize::setShellRadius(const Shell& is, const Double& thresh) const 
{

	//
	// dummy shell actually has infinite shell size
	// since it's exponent part is zero, therefore it's a constant line
	// here we will return a very large number
	//
	if (is.isDummyShell()) return LARGE_FLOAT_NUMBER;

	// get the smallest exponent
	const Double* e = is.getExp();
	Double smallestAlpha = e[0];
	UInt smallAlphaPos = 0;
	for(UInt i=0; i<is.getNPrim(); i++) {
		if (e[i]<smallestAlpha) {
			smallestAlpha = e[i];
			smallAlphaPos = i;
		}
	}

	// do the job for each sub shell
	Double result = ZERO;
	UInt lmin = is.getLmin();
	UInt lmax = is.getLmax();
	const Double* c = is.getCoe();
	for(UInt iShell=lmin;iShell<=lmax;iShell++) {

		// firstly, we try to search the outmost spreaded primitive
		// as a approximation for the initial r
		// we note, that negative coeff is not important, we can revise
		// it as + sign
		UInt np = (iShell-lmin)*is.getNPrim(); 
		Double coeff = c[np+smallAlphaPos];
		Double r0 = ZERO;
		if (coeff < ZERO) {
			coeff = fabs(coeff);
		}
		Double t1 = MINUS_ONE*log(thresh/coeff);
		if (t1 > 0) {
			r0 = sqrt(t1/smallestAlpha);
		}else{
			string infor = "the coefficients is: " + boost::lexical_cast<string>(coeff);
			Excep excep("AtomShellSize","setShellRadius",EXCEPTION_SHELL_RADIUS_COEFF_TOO_SMALL,infor);
			handleExcep(excep);
		}

		// now we use Newton-Raphson method to get the real r...
		Double r;
		UInt count;
		Double logThresh = log(thresh);
		for(count=0; count<1000; count++) {

			// f and df
			Double r2  = r0*r0; 
			Double f0 = ZERO;
			Double fd = ZERO;
			for(UInt i=0; i<is.getNPrim(); i++){
				Double x  = c[np+i]*exp(MINUS_ONE*e[i]*r2);
				f0       += x;
				fd       += MINUS_TWO*e[i]*r0*x;
			}

			// if the coefficients are negative then the basis set
			// value could be negative, too. In this case, we take
			// the opposite value
			Double f,df;
			if (f0 > ZERO) {
				f  = iShell*log(r0) + log(f0) - logThresh; 
				df = iShell/r0 + fd/f0;
			}else{
				f0 = MINUS_ONE*f0;
				fd = MINUS_ONE*fd;
				f  = iShell*log(r0) + log(f0) - logThresh; 
				df = iShell/r0 + fd/f0;
			}	

			// calculate the dr
			Double dr = f/df;
			r = r0 - dr;
			if (fabs(dr) < thresh) {
				break;
			}else{
				r0 = r;
			}
		}

		// check the result
		if (count == 1000) {
			string infor = "radius calculation is failed";
			Excep excep("AtomShellSize","setShellRadius",EXCEPTION_SHELL_RADIUS_NOT_CONVERGE,infor);
			handleExcep(excep);
		}
		if (result < r) result = r;
	}
	return result;
}

void AtomShellSize::print() const 
{
	printf("Atomic number: %-5d with radius:  %-12.6f\n", (Int)atomic, atomRadius);
	for(UInt iShell=0; iShell<shellRadius.size(); iShell++) {
		Double r = getShellRadius(iShell);
		printf("Shell Radius:  %-6d  %-12.6f\n", (Int)(iShell+1), r);
	}
}

MolShellSize::MolShellSize(const MolShell& ms, const Double& thresh) 
{
	// set up basic information about atom types
	UInt atomTypes[MAX_ATOMS_TYPE];
	for(UInt i=0; i<MAX_ATOMS_TYPE; i++) atomTypes[i] = 0;

	// calculate radius and size
	UInt nAtomTypes = ms.getNAtomShellTypes();
   atomShellSizes.reserve(nAtomTypes);
	for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {
		const AtomShell& iAShell = ms.getAtomShell(iAtom);
		UInt Z = iAShell.getAtomic();
		if (isGhostAtom(Z)) continue; 
		if (atomTypes[Z] > 0) continue;
		AtomShellSize atomShellSize(iAShell,thresh);
		atomShellSizes.push_back(atomShellSize);
		atomTypes[Z] = 1;
	}
}

const AtomShellSize& MolShellSize::getAtomShellSizeInfor(const UInt& atomic) const
{
	for(UInt i=0; i<atomShellSizes.size(); i++) {
		const AtomShellSize& atomShellSize = atomShellSizes[i];
		if (atomShellSize.getAtomic() == atomic) {
			return atomShellSize;
		}
	}

	// if we approach here, then we are unable to find
	// this atom data
	string infor = "the given atomic number is: " + boost::lexical_cast<string>(atomic); 
	Excep excep("MolShellSize","getAtomShellSizeInfor",
					EXCEPTION_UNKNOWN_ATOM_TYPE_ATOM_SHELL_SIZE ,infor);
	handleExcep(excep);
	return atomShellSizes[0];
}

UInt MolShellSize::getNAtomTypes() const
{
	UInt count = 0;
	UInt possibleAtomTypes[MAX_ATOMS_TYPE];
	for(UInt i=0; i<MAX_ATOMS_TYPE; i++) possibleAtomTypes[i] = 0;
	for(UInt i=0; i<getNAtomShellData(); i++) {
		UInt atomic = atomShellSizes[i].getAtomic();
		if (possibleAtomTypes[atomic] > 0) {
			continue;
		}else{
			possibleAtomTypes[atomic] = 1;
			count++;
		}
	}
	return count;
}

void MolShellSize::getAtomTypeInfor(UIntVec& infor) const
{
	// by using this function, we assume that you reserve
	// space for the infor array
	UInt possibleAtomTypes[MAX_ATOMS_TYPE];
	for(UInt i=0; i<MAX_ATOMS_TYPE; i++) possibleAtomTypes[i] = 0;
	for(UInt i=0; i<getNAtomShellData(); i++) {
		UInt atomic = atomShellSizes[i].getAtomic();
		if (possibleAtomTypes[atomic] > 0) {
			continue;
		}else{
			infor.push_back(atomic);
		}
	}
}

void MolShellSize::print() const
{
	cout << endl;
	cout << "*******************************************" << endl;
	cout << "*        MolShellSize Object Print        *" << endl;
	cout << "*******************************************" << endl;
	for(UInt i=0; i<atomShellSizes.size(); i++) {
		const AtomShellSize& atomShellSize = atomShellSizes[i];
		atomShellSize.print();
	}
	cout << endl;
}
