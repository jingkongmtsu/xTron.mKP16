/**
 * \file    scfderiv.cpp
 * \brief   class mornitoring the first and second derivatives calculation for SCF
 * \author  Fenglai Liu 
 * \note
 */
#include <iostream>
#include "globalinfor.h"
#include "molecule.h"
#include "shell.h"
#include "denmtrx.h"
#include "gints4deriv.h"
#include "scfparam.h"
#include "scfderiv.h"
using namespace globalinfor;
using namespace molecule;
using namespace shell;
using namespace denmtrx;
using namespace gints4deriv;
using namespace scfparam;
using namespace scfderiv;
using namespace std;

SCFDeriv::SCFDeriv(const Molecule& mol, UInt order):jobOrder(order)
{
	// check the job order
	// only 1 or 2 is available
	if (jobOrder != 1 && jobOrder != 2) {
		string infor = "input job order is not supported for scf deriv calculation"; 
		Excep excep("SCFDeriv","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	// now let's initialize the result
	UInt nAtoms = mol.getNAtoms();
	if (jobOrder == 1) {
		derivData.init(nAtoms,3);
	}else if (jobOrder == 2) {
		derivData.init(nAtoms*3,nAtoms*3);
	}
}

void SCFDeriv::doDeriv(const SCFParam& param, const Molecule& mol, 
		const MolShell& ms, DenMtrx& den)
{
	// right now we do not have one electron deriv 
	// and the DFT part contribution, so we directly
	// go to the two electron part

	// do JK deriv
	UInt job   = COULOMB_EXCHANGE;
	const XCFunc& xcfunc = param.getXCFunc();
	if (! xcfunc.isHybrid()) {
		job = COULOMB;
	}
	const GIntsInfor& infor = param.getGIntsInfor();
	GInts4Deriv gints(infor,job,jobOrder);
	gints.doJKDeriv(ms,mol,den,derivData,true,false);
}

