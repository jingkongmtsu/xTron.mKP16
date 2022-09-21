/**
 * \file    fock.cpp
 * \brief   describe the Fock matrix, how to build it and how to handle it in SCF
 * \author  Fenglai Liu 
 */
#include<string>
#include<boost/lexical_cast.hpp>
#include "globalinfor.h"
#include "molecule.h"
#include "shell.h"
#include "matrix.h"
#include "gintsinfor.h"
#include "xcintsinfor.h"
#include "gints4d.h"
#include "xcints.h"
#include "denmtrx.h"
#include "scfparam.h"
#include "oneemtrx.h"
#include "scfintscontroller.h"
#include "fock.h"
using namespace globalinfor;
using namespace molecule;
using namespace shell;
using namespace matrix;
using namespace gintsinfor;
using namespace xcintsinfor;
using namespace gints4d;
using namespace xcints;
using namespace denmtrx;
using namespace scfparam;
using namespace oneemtrx;
using namespace scfintscontroller;
using namespace fock;
using namespace std;

Fock::Fock(const SCFParam& param, const MolShell& ms):SpinMatrix(param.getNSpin(),ms.getNBas(),ms.getNBas()),
	nRow(ms.getNBas()),nCol(ms.getNBas())
{  }

Double Fock::addCoreMtrx(const SCFParam& param, const OneEMtrx& oneEMtrx, 
		const MolShell& ms,const DenMtrx& denMtrx)
{
	// reload core matrix 
	UInt nBas = ms.getNBas();
	Mtrx core(nBas,nBas);
	oneEMtrx.getM(KINETIC,core);
	oneEMtrx.getM(NUCLEAR_ATTRACTION,core);
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		Mtrx& F = getMtrx(iSpin);
		F.add(core);
	}

	// consider the core Hamiltonian energy in normal way
	bool symmMatrix = true;
	Double ecore = ZERO;
	for(UInt iSpin=0; iSpin<denMtrx.getNSpin(); iSpin++) {
		const Mtrx& den  = denMtrx.getMtrx(iSpin);
		ecore += core.dotProduct(den,symmMatrix);
	}

	// consider the close shell case
	if (param.closeShell()) ecore *= TWO;
	return ecore;
}

Double Fock::addJKMtrx(const SCFParam& param, const SCFIntsController& intsController,
		const Molecule& mol, const MolShell& ms, DenMtrx& den)
{
	// let's see whether the user demanded to separate the JK
	if (param.doJKApart()) {

		// firstly, do the coulomb part
		UInt job    = COULOMB;
		Double eJ   = ZERO;
		GInts4D gints(intsController.getGIntsInfor(),job);
		gints.doJKMtrx(ms,den,*this,eJ,param.printSCFTimingData(),
				param.printJKMatrix(),param.saveJKOnDisk());
		printf("Coulomb energy is: %-16.10f\n", eJ);

		// now let's see whether we have exchange part
		Double eK = ZERO;
		const XCFunc& xcfunc = param.getXCFunc();
		if (xcfunc.isHybrid()) {
			job = EXCHANGE;
			GInts4D gints1(intsController.getGIntsInfor(),job);
			gints1.doJKMtrx(ms,den,*this,eK,param.printSCFTimingData(),
					param.printJKMatrix(),param.saveJKOnDisk());
			printf("Exchange energy is: %-16.10f\n", eK);
		}

		// now return the energy
		return eJ+eK;
	}

	// this is the normal way to do JK
	// let's set the JK job
	// in default, we do JK together except that
	// only J is needed
	UInt job   = COULOMB_EXCHANGE;
	const XCFunc& xcfunc = param.getXCFunc();
	if (! xcfunc.isHybrid()) {
		job = COULOMB;
	}

	// now let's do the job
	// the fock matrix itself will be passed in
	Double eJK  = ZERO;
	GInts4D gints(intsController.getGIntsInfor(),job);
	gints.doJKMtrx(ms,den,*this,eJK,param.printSCFTimingData(),
			param.printJKMatrix(),param.saveJKOnDisk());
	return eJK;
}

Double Fock::addXCMtrx(const SCFParam& param, const SCFIntsController& intsController,
		const Molecule& mol, const MolShell& ms, const DenMtrx& den)
{
	// now let's do the job
	// the fock matrix itself will be passed in
	UInt order  = 0;
	UInt job    = GROUND_STATE_DFT; 
	XCInts xcints(param.getGlobalInfor(),param.getXCFunc(),ms,
			intsController.getXCIntsInfor(),intsController.getGIntsInfor(),job,order);
	xcints.doXCMtrx(ms,den,*this,param.printSCFTimingData(),param.printXCMatrix());
	return xcints.getExc();
}


