/**
 * \file    denmtrx2.cpp
 * \brief   describe the functions SCF guess generation for density matrices
 * \author  Fenglai Liu 
 */
#include "mo.h"
#include "shell.h"
#include "excep.h"
#include "gints2d.h"
#include "oneemtrx.h"
#include "molecule.h"
#include "gintsinfor.h"
#include "atomdenmtrx.h"
#include "scfparam.h"
#include "scf.h"
#include "denmtrx.h"
using namespace mo;
using namespace shell;
using namespace excep;
using namespace gints2d;
using namespace oneemtrx;
using namespace molecule;
using namespace gintsinfor;
using namespace atomdenmtrx;
using namespace scfparam;
using namespace scf;
using namespace denmtrx;

void DenMtrx::coreGuess(const MolShell& ms, const Molecule& mol)
{
	// based on the shell information, we will 
	// form the core hamiltonian as well as the 
	// orthogonal matrix
	UInt order = 0;
	UInt nBas  = ms.getNBas();
	GIntsInfor ginfor(infor,mol);

	// set up the matrix
	Mtrx Fock(nBas,nBas);
	GInts2D ki(ms,ms,ginfor,KINETIC,order);
	ki.doMtrx(ms,ms,mol,Fock,false);
	GInts2D nai(ms,ms,ginfor,NUCLEAR_ATTRACTION,order);
	nai.doMtrx(ms,ms,mol,Fock,false);

	// now get the orthogonal matrix
	GInts2D ov(ms,ms,ginfor,TWO_BODY_OVERLAP,order);
	Mtrx twoOv(nBas,nBas);
	ov.doMtrx(ms,ms,mol,twoOv,false);
	Mtrx ortho(nBas,nBas);
	ortho.formOrthoMatrix(twoOv,ginfor.getLinearDepThresh()); 

   // now form the mo 
	MO mo(infor,ms,mol,getNSpin());
	mo.formMO(ortho,Fock);

   // fill in the density matrix
   formDenMtrx(mo);
}

void DenMtrx::coreGuess(const MolShell& ms, const Molecule& mol, const OneEMtrx& oneEMtrx)
{
	// loading the core and ortho
	// using ortho as tmp scratch matrix
	// to load in nai, too
	UInt nBas  = ms.getNBas();
	Mtrx core(nBas,nBas);
	Mtrx ortho(nBas,nBas);
	oneEMtrx.getM(KINETIC,core);
	oneEMtrx.getM(NUCLEAR_ATTRACTION,ortho);
	core.add(ortho);

	// now load in ortho
	ortho.set(ZERO);
	UInt nRow,nCol; 
	oneEMtrx.getDim(ORTHOGONAL_MATRIX,nRow,nCol);
	ortho.reset(nRow,nCol,true);
	oneEMtrx.getM(ORTHOGONAL_MATRIX,ortho);

   // now form the mo 
	MO mo(infor,ms,mol,getNSpin());
	mo.formMO(ortho,core);

   // fill in the density matrix
   formDenMtrx(mo);
}

void DenMtrx::atomicDenMtrxGuess(const MolShell& ms, const Molecule& mol) 
{
	// for the atomic density matrix guess, the density matrix
	// itself should be in it's normal type too
	// 
	// the atomic density matrix data generated also in normal type
	// so we need to check it here
	UInt type = TYPE_NORM;
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		const Mtrx& P = getMtrx(iSpin);

		// the density matrix should be square
		UInt nRow  = P.getRow();
		UInt nCol  = P.getCol();
		if (nRow != nCol) {
			string info = "in using atomic density matrix method the density matrix should be square";
			Excep excep("DenMtrx","atomicDenMtrxGuess",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
			handleExcep(excep);
		}

		// now check this matrix
		UInt type0 = -1;
		if(nRow == ms.getNBas()) {
			type0 = TYPE_NORM;
		}	

		// see whether type is defined
		if (type != type0) {
			string info = "density matrix must be in normal type(TYPE_NORM) for using the atomic density matrix guess";
			Excep excep("DenMtrx","atomicDenMtrxGuess",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
			handleExcep(excep);
		}
	}

	// firstly form the atom density matrix data
	AtomDenMtrx den(infor,mol);
	den.formAtomDenMtrxInSCF(mol,ms); 

	// now let's consider the spin state
	// if this is single spin state, it's benefit
	// that we can average alpha and beta density matrix
	if (getNSpin() == 1) {
		den.averageSpinDenMtrx();
	}

	// now let's stuff the result atom density matrix data 
	// to the diagonal block of density matrix
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		Mtrx& denMtrx = getMtrx(iSpin);
		UInt pos = 0;
		for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {

			// get the atomic number
			const AtomShell& atomShell = ms.getAtomShell(iAtom);
			UInt atomic = atomShell.getAtomic();

			// get the density matrix data
			const DenMtrx& atomDenMtrx  = den.getAtomDenMtrx(atomic);
			const Mtrx& spinAtomDenMtrx = atomDenMtrx.getMtrx(iSpin); 

			// get the position
			// because it's on diagonal block, so rowPos = colPos
			UInt rowPos = pos;
			UInt colPos = pos;
			denMtrx.addPieceInDirectSum(spinAtomDenMtrx,rowPos,colPos);
			pos = pos + atomShell.getNBas(type);
		}
	}
}

void DenMtrx::formSCFGuess(const MolShell& ms, const Molecule& mol, const SCF& scf)
{
	// get the method
	const SCFParam& par = scf.getSCFParam();
	UInt method = par.getGuess();

	// now let's see the choice
	if (method == SCF_GUESS_CORE) {
		coreGuess(ms,mol,scf.getOneEMtrx());
	}else if (method == SCF_GUESS_ATOMIC) {
		atomicDenMtrxGuess(ms,mol);
	}else if (method == SCF_GUESS_READ) {
		MO mo(infor,ms,mol,getNSpin());
		mo.recover(par.getGuessMOPath());
		formDenMtrx(mo);
	}else{
		string info = "the given guess method can not be recognized";
		Excep excep("DenMtrx","formSCFGuess",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
		handleExcep(excep);
	}
}
