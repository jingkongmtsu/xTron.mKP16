/**
 * functions and classes associated with digestutil.h
 * \author Fenglai Liu
 */
#include<cmath>
#include "shellprop.h"
#include "matrix.h"
#include "shell.h"
#include "blas.h"
#include "digestutil.h"
using namespace shellprop;
using namespace matrix;
using namespace shell;
using namespace blas;
using namespace digestutil;

void digestutil::scale4BodyInts(const UInt& LBra1Min, const UInt& LBra1Max, 
		const UInt& LBra2Min, const UInt& LBra2Max, 
		const UInt& LKet1Min, const UInt& LKet1Max, 
		const UInt& LKet2Min, const UInt& LKet2Max, 
		const Double* bra1Scale, const Double* bra2Scale, 
		const Double* ket1Scale, const Double* ket2Scale,
		Double* I) 
{
	// do we return it?
	// for L<=P, we do nothing here
	if (LBra1Max<=1 && LBra2Max<=1 && LKet1Max<=1 && LKet2Max<=1) return;

	// get dimension
	UInt nKet2Bas = getCartBas(LKet2Min,LKet2Max);
	UInt nKet1Bas = getCartBas(LKet1Min,LKet1Max);
	UInt nBra2Bas = getCartBas(LBra2Min,LBra2Max);
	UInt nBra1Bas = getCartBas(LBra1Min,LBra1Max);

	// make the dimension more clear
	// to access the integral element, we do not need the last dimension
	UInt d1 = nBra1Bas;
	UInt d2 = d1*nBra2Bas;
	UInt d3 = d2*nKet1Bas;

	// now let's scale the result
	for(UInt lBas=0; lBas<nKet2Bas; lBas++) {
		for(UInt kBas=0; kBas<nKet1Bas; kBas++) {
			for(UInt jBas=0; jBas<nBra2Bas; jBas++) {
				for(UInt iBas=0; iBas<nBra1Bas; iBas++) {
					Double N = bra1Scale[iBas]*bra2Scale[jBas]*ket1Scale[kBas]*ket2Scale[lBas];
					I[iBas+jBas*d1+kBas*d2+lBas*d3] *= N;
				}
			}
		}
	}
}

void digestutil::postProcessJKMatrix(const MolShell& s, Mtrx& K) 
{
	// firstly, we do add transpose
	K.addTranspose();

	// the diagonal shell block data is two times larger
	// therefore we scale it back
	for(UInt iAtom=0; iAtom<s.getNAtomShells(); iAtom++) {
		const AtomShell& atomShell = s.getAtomShell(iAtom);
		for(UInt iShell=0; iShell<atomShell.getNShell(); iShell++) {
			const Shell& shell = atomShell.getShell(iShell);
			UInt nCarBas = shell.getNCarBas();
			UInt offset  = shell.getBasisIndex(0,TYPE_CART);
			for(UInt j=0; j<nCarBas; j++) {
				for(UInt i=0; i<nCarBas; i++) {
					K(offset+i,offset+j) *= 0.5E0;
				}
			}
		}
	}
}

