/**
 * \file    scf.cpp
 * \brief   class mornitoring the whole SCF calculation
 * \author  Fenglai Liu 
 * \note
 */
#include <iostream>
#include "globalinfor.h"
#include "excep.h"
#include "molecule.h"
#include "shell.h"
#include "mo.h"
#include "vdwinfor.h"
#include "xcfunc.h"
#include "xcintsinfor.h"
#include "gints2d.h"
#include "denmtrx.h"
#include "atomdenmtrx.h"
#include "scfmacro.h"
#include "scfconv.h"
#include "fock.h"
#include "wfn.h"
#include "scfenergyconv.h"
#include "scf.h"
#include "fracspininfor.h"
#ifdef WITH_GDM
#include "gdm.h"
using namespace gdm;
#endif
using namespace globalinfor;
using namespace excep;
using namespace molecule;
using namespace shell;
using namespace mo;
using namespace vdwinfor;
using namespace xcfunc;
using namespace xcintsinfor;
using namespace gints2d;
using namespace denmtrx;
using namespace atomdenmtrx;
using namespace scfmacro;
using namespace scfconv;
using namespace scfenergyconv;
using namespace fock;
using namespace wfn;
using namespace scf;
using namespace std;
using namespace fracspininfor;


SCF::SCF(const GlobalInfor& infor, const Molecule& mol, 
		const MolShell& ms):eNulRep(mol.getNuclearRepulsion()),
	eCore(ZERO),eJK(ZERO),eXC(ZERO),eTotal(ZERO),param(infor,mol),
	intsController(param),oneEMtrx(infor,param.getSec(),param.useFile()),
	mo(infor,ms,mol,param.getNSpin()) 
{ 
	// pre-compute the one E matrix
	oneEMtrx.formMtrx(intsController.getGIntsInfor(),mol,ms,ms,
			param.printSCFTimingData(),param.printCoreMatrix()); 

	// on the other hand, we may also need the free atom density matrix
	const XCFunc& xcfunc = param.getXCFunc(); 
	if (xcfunc.useAtomDenMtrxData()) {
		AtomDenMtrx den(infor,mol);
		den.formAtomDenMtrxInSCF(mol,ms); 
		den.writeToDisk();
	}
}

SCF::SCF(const GlobalInfor& infor, const SCFParam& param0, const Molecule& mol, 
		const MolShell& ms):eNulRep(mol.getNuclearRepulsion()),
	eCore(ZERO),eJK(ZERO),eXC(ZERO),eTotal(ZERO),param(param0),
	intsController(param),oneEMtrx(infor,param.getSec(),param.useFile()),
	mo(infor,ms,mol,param.getNSpin()) 
{  

	// pre-compute the one E matrix
	oneEMtrx.formMtrx(intsController.getGIntsInfor(),mol,ms,ms,
			param.printSCFTimingData(),param.printCoreMatrix()); 

	// on the other hand, we may also need the free atom density matrix
	const XCFunc& xcfunc = param.getXCFunc(); 
	if (xcfunc.useAtomDenMtrxData()) {
		AtomDenMtrx den(infor,mol);
		den.formAtomDenMtrxInSCF(mol,ms); 
		den.writeToDisk();
	}
}

void SCF::updateMO(const MO& mos) 
{
	// the input mo and mo here should have same nSpin
	if (mos.getNSpin() != mo.getNSpin()) {
		string infor = "input mo has different nSpin with the mo data inside"; 
		Excep excep("SCF","updateMO",EXCEPTION_SCF_INPUT_INFOR_INVALID,infor);
		handleExcep(excep);
	}

	// now let's do update
	for(UInt iSpin=0; iSpin<mos.getNSpin(); iSpin++) {
		const Mtrx& oldMO = mos.getMtrx(iSpin);
		Mtrx& newMO = mo.getMtrx(iSpin);
		newMO = oldMO;
	}
}

void SCF::doSCF(const GlobalInfor& infor, const Molecule& mol, 
		const MolShell& ms, DenMtrx& den)
{
	SCFConv conv(infor,param);
	doSCFConv(infor, mol, ms, den, conv);
}

void SCF::doSCFConv(const GlobalInfor& infor, const Molecule& mol, 
		const MolShell& ms, DenMtrx& den, SCFConv& conv)
{
	// initilize the data
	//SCFConv conv(infor,param);
	Fock fock(param,ms);

	// for using GDM, we will impose the use the energy 
	// difference check for SCF convergence
	SCFEnergyConvController energyDiff(infor,mol);

#ifdef WITH_GDM

	// let's initialize the GDM
	GDM gdm(param,mol,ms);

#endif

	// do we print out the initial density matrix?
	if (param.printInitialDensityMatrix()) {
		den.print("initial density matrix before SCF");
	}

	// here loop begins
	for(UInt iter=0; iter<param.getMaxSCFCycles(); iter++) {

		// form the fock 
		// core and JK part
		if (iter>0) {
			fock.initialize(ZERO); 
		}
		eTotal  = eNulRep;
		eCore   = fock.addCoreMtrx(param,oneEMtrx,ms,den);
		eJK     = fock.addJKMtrx(param,intsController,mol,ms,den);
		eTotal += eCore;
		eTotal += eJK;

		// do we have DFT functional?
		// add in XC part
		const XCFunc& xcfunc = param.getXCFunc();
		if (! xcfunc.withoutDFTFunc()) {
			eXC = fock.addXCMtrx(param,intsController,mol,ms,den);
			eTotal += eXC;
		}else{
			eXC = ZERO;
		}

		// debug, let's print out the calculated fock matrix
		if (param.printFockMatrix()) {
			fock.print("fock matrix after adding all of components");
		}

		// now form new Fock based on scf convergence
		// algorithm
		conv.convSCF(iter,eTotal,den,oneEMtrx,fock);
		if (param.printSCFConvClass()) {
			conv.print();
		}

		// debug, let's print out the scf fock matrix
		if (param.printSCFFockMatrix()) {
			fock.print("new fock matrix derived from scfConv class for generating new MO/density matrix");
		}

		// print out energy components
		if (! param.isSCFSilent()) {

			// get the SCF functional name
			const XCFunc& xcfunc = param.getXCFunc();
			string funcName = xcfunc.getFuncName();

			// now print out the stuff
			printf("DIIS error %-16.10f\n", conv.DIISError()); 
			if (conv.scfEnds()) {
				printf("%s iter is: %d, @@@@@ scf converged energy is: %-16.10f \n", funcName.c_str(), (Int)iter, eTotal);
			}else{
				printf("%s iter is: %d, @@@@@ total energy is: %-16.10f \n", funcName.c_str(), (Int)iter, eTotal);
			}
			// print out the most lowest energy for the whole scf
			Double lowestE;
			UInt lowestEIndex;
			conv.getLowestEnergyInfor(lowestE,lowestEIndex);
			printf("The most lowest energy so far is: %-16.10f, scf index is %d\n",lowestE, (Int)lowestEIndex);
			printf("NucRepulsion %-16.10f\n", eNulRep);
			printf("core %-16.10f\n", eCore);
			printf("JK   energy is: %-16.10f\n", eJK);
			printf("XC   %-16.10f\n", eXC);
			printf("MTAG %s %3d %3d %-16.10f %-16.10f\n",funcName.c_str(), (Int)iter, (Int)lowestEIndex, eTotal, conv.DIISError()); 
		}

#ifdef WITH_GDM

		// let's do GDM
		if (useGDM(conv.getJob())) {
			gdm.doGDM(iter,mo,fock,eTotal);
			cout << "GDM is in use" << endl;

			// we also need to update the MO energy
			SpinMatrix tmp(fock);
			MO newMO(mo);
			newMO.formMO(ms,mol,oneEMtrx,tmp);
			for(UInt iSpin=0; iSpin<newMO.getNSpin(); iSpin++) {
				const DoubleVec& moE0 = newMO.getEnergyVec(iSpin);
				DoubleVec& moE1       = mo.getEnergyVec(iSpin);
				moE1 = moE0;
			}
		}
#else
		if (useGDM(conv.getJob())) {
			string infor = "you can not use gdm job here, since the GDM job is not defined. Please add gdm library and recompile";
			Excep excep("SCF","doSCF",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
			handleExcep(excep);
		}
#endif

		// form the MO
		// we will reform the mo if SCF ends
		// and update density matrix if SCF is not ended
		if (! useGDM(conv.getJob())) {
			mo.formMO(ms,mol,oneEMtrx,fock);
		}
    // let's see whether we have any fractional spin job?
		// if so we need to further scale the fractional number of mo
		FracSpinInfor fracSpinInfor(infor,mol);
		if (fracSpinInfor.doFracSpinJob()) {
		  mo.formFracSpinMO(fracSpinInfor);
		}
		//
		// form the density matrix
		den.formDenMtrx(mo);

		// is the SCF proceed to an end?
		if (conv.scfEnds()) {
			doPostSCF(infor,mol,ms,den);
			break;
		}else if (useGDM(conv.getJob()) && energyDiff.hasEnergyConv(conv)) {
			doPostSCF(infor,mol,ms,den);
			break;
		}else{
			if (param.printDensityMatrixInSCF()) {
				den.print("density matrix during scf process");
			}
		}
	}
}

void SCF::doPostSCF(const GlobalInfor& infor, const Molecule& mol,
		const MolShell& ms, const DenMtrx& den) const
{
	// firstly let's see whether we do post SCF?
	// if not we just returned here
	if (! param.doPostSCFWork()) return;

	// let's check whether we need to print out stuff
	// or save the results?
	// print or save MO?
	if (param.useFile())	mo.writeToDisk();
	if (param.printMOData()) mo.print(param.getMOPrintOption());

	// form density matrix
	if (param.printResultDensityMatrix()) {
		den.print("result density matrix");
	}
	if (param.saveResultDenMtrxOnDisk()) {
		den.writeToDisk(param.getResultDenMtrxPath());
	}

	// let's see whether we print out the wfx file
	WFN wfx(infor,mol);
	wfx.dump(eTotal,infor,mol,ms,mo);

	// do we do dipole moment calculation after SCF ends?
	if (param.doDipole()) {

		// we need to set up the momentum integral matrix
		MtrxVec momIntMtrx(3);
		UInt nBas = ms.getNBas();
		for(UInt i=0; i<3; i++) {
			Mtrx& intMtrx = momIntMtrx[i];
			intMtrx.init(nBas,nBas);
		}

		// now let's do the calculation
		GInts2D gints2d(ms,ms,param.getGIntsInfor(),MOM_P,0);
		gints2d.doMultiMtrx(ms,ms,momIntMtrx,false);

		// let's combine the density matrix to produce
		// the dipole moment
		// we multiply -1 for the negative of electron density
		Double dipole[3];
		dipole[0] = ZERO;
		dipole[1] = ZERO;
		dipole[2] = ZERO;
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			const Mtrx& denMtrx = den.getMtrx(iSpin);
			for(UInt iM=0; iM<3; iM++) {
				const Mtrx& intMtrx = momIntMtrx[iM];
				dipole[iM] += MINUS_ONE*intMtrx.dotProduct(denMtrx,true);
			}
		}

		// if this is close shell, we need to scale the
		// result by two so that to add both alpha and 
		// beta part
		if (den.getNSpin() == 1) {
			dipole[0] *= TWO;
			dipole[1] *= TWO;
			dipole[2] *= TWO;
		}

		// we have to add atomic contribution to the 
		// dipole moment
		// use the geometry of the molecule
		Double atomicDipole[3];
		atomicDipole[0] = ZERO;
		atomicDipole[1] = ZERO;
		atomicDipole[2] = ZERO;
		for(UInt iAtom=0; iAtom<mol.getNAtoms(); iAtom++) {
			const Atom& atom = mol.getAtom(iAtom);
			UInt atomic      = atom.getAtomic();
			const Double* xyz= atom.getXYZ();
			atomicDipole[0] += xyz[0]*atomic;
			atomicDipole[1] += xyz[1]*atomic;
			atomicDipole[2] += xyz[2]*atomic;
		}

		// now let's add them together
		dipole[0] += atomicDipole[0];
		dipole[1] += atomicDipole[1];
		dipole[2] += atomicDipole[2];

		// convert it into  debye unit
		const Double UNIT_TO_DEBYE = 2.541765E0;
		dipole[0] *= UNIT_TO_DEBYE;
		dipole[1] *= UNIT_TO_DEBYE;
		dipole[2] *= UNIT_TO_DEBYE;

		// calculate total dipole
		Double totalDipole = sqrt(dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2]);

		// now let's print out the dipole moment
		// alpha state
		printf("========================================\n");
		printf("Dipole Moment(Debye) in post SCF process\n");
		printf("========================================\n");
		printf("dipole moment: %-15.8f  %-15.8f  %-15.8f\n", dipole[0], dipole[1], dipole[2]);
		printf("total dipole moment: %-15.8f\n", totalDipole);
	}

	// do we do the VDW calculation for the post SCF?
	if (param.doVDWInPostSCF()) {
		VDWInfor vdwInfor(infor,mol);
		XCIntJobInfor xcJobInfor(infor,param.getGIntsInfor(),param.getXCIntsInfor(),GROUND_STATE_DFT,0);
		vdwInfor.doXDM(mol,ms,xcJobInfor,den);
		printf("VDW energy: %-15.8f\n", vdwInfor.vdwEnergy());
	}
}
