/**
 * \file    atomdenmtrx.cpp
 * \brief   classes used to generate the SCF converged density matrix etc. for free atom
 * \author  Fenglai Liu and Jing Kong
 */
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>
#include "parameterparsing.h"
#include "globalinfor.h"
#include "textread.h"
#include "molecule.h"
#include "element.h"
#include "shell.h"
#include "excep.h"
#include "xcfunc.h"
#include "mo.h"
#include "scf.h"
#include "atomdenmtrx.h"
#include "scfconv.h"

using namespace boost::filesystem;
using namespace parameterparsing;
using namespace globalinfor;
using namespace textread;
using namespace molecule;
using namespace element;
using namespace shell;
using namespace excep;
using namespace xcfunc;
using namespace mo;
using namespace scf;
using namespace atomdenmtrx;
using namespace std;
using namespace scfconv;

AtomDenMtrx::AtomDenMtrx(const GlobalInfor& infor0, const Molecule& mol):infor(infor0),
	section(mol.getSec()),scfEXMethod("HF"),scfECMethod("NONE"),maxSCFCycles(100),scfRunInSilence(true)
{
	// generate atom types information
	// and initilize the atom density matrix data
	UInt nTypes = atomTypes.size();
	atomTypes.reserve(nTypes);
	mol.getAtomTypes(atomTypes);
	atomDenMtrices.reserve(nTypes);

	// now let's read in the parameters
	string input = infor.getInputFile();
	ParameterParsing pp(input,"atomdenmtrx",section);
	if (pp.hasAnyParameters()) {

		// scf method for running atom density matrix data
		string key = "scf_method";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			LineParse l(value);
			if (l.getNPieces() == 1) {
				string tmpValue = l.findValue(0);
				scfEXMethod = tmpValue;
			}else if (l.getNPieces() == 2) {
				string tmpValue = l.findValue(0);
				scfEXMethod = tmpValue;
				tmpValue = l.findValue(1);
				scfECMethod = tmpValue;
			}else{
				string infor = "the value for scf_method is invalid, only one/two names should be given";
				Excep excep("AtomDenMtrx","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// SCF cycles
		key = "max_scf_cycles";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "max_scf_cycles can not be processed, not an integer";
				Excep excep("AtomDenMtrx","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp > 0) {
				maxSCFCycles = tmp;
			}else{
				string infor = "invalid number of maximum SCF cycles set";
				Excep excep("AtomDenMtrx","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
		}

		// do we keep en eye on SCF process?
		key = "scf_debug";
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				scfRunInSilence = false;
			}
		}
	}
}

void AtomDenMtrx::formAtomDenMtrxInSCF(const Molecule& mol, const MolShell& ms)
{
	// get the functional information
	// for the calculation
	XCFunc xcfunc(scfEXMethod,scfECMethod);

	// now loop over free atom to get the density matrix
	// for the free atom, charge is always 0
	UInt nTypes = atomTypes.size();
	for(UInt i=0; i<nTypes; i++) {

		// the atomic number for this type
		UInt atomic = atomTypes[i];
		if (isGhostAtom(atomic)) continue;

		// set up the one atom molecule
		// also we update the section number
		Molecule molAtom(atomic);
		molAtom.updateSec(section);

		// now create the one atom shell molshell
		UInt iAtom = ms.getNAtomShells()+1;
		for(UInt j=0; j<ms.getNAtomShells(); j++) {
			const AtomShell& as = ms.getAtomShell(j);
			if (as.getAtomic() == atomic) {
				iAtom = j;
				break;
			}
		}
		if (iAtom == ms.getNAtomShells()+1) {
			string info = "can not find the corresponding atom shell data for the free atom(atomic number): " 
				+ boost::lexical_cast<string>(atomic);
			Excep excep("AtomDenMtrx","formAtomDenMtrxInSCF",EXCEPTION_ATOM_DENSITY_MATRICES_ERROR,info);
			handleExcep(excep);
		}
		MolShell atomShell(ms.getAtomShell(iAtom));

		// set the scf parameters
		// we note that because to avoid overwriting the 
		// files for the exsiting SCF section, we force
		// everything here inside memory
		SCFParam param(infor,molAtom);
		param.keepInMem();

		// let's disable post SCF process to avoid possible trouble
		param.disablePostSCF();

		// we need to allocate enough SCF cycles for it for run
		UInt newSCFCycles = maxSCFCycles;
		param.changeMaxSCFCycles(newSCFCycles);

		// if scf is not in debug mode, we do it in silence way
		if (scfRunInSilence) {
			param.keepSCFSilence();
		}

		// update the xcfunc information
		param.updateXCFunc(xcfunc);

		// now let's do the SCF
		// here we always starting from core guess
		SCF scf(infor,param,molAtom,atomShell);
		DenMtrx den(infor,molAtom,atomShell,atomShell,param.getNSpin());
		den.coreGuess(atomShell,molAtom);
		//scf.doSCF(infor,molAtom,atomShell,den);
		SCFConv scfConv(infor, param, "diis");
		scf.doSCFConv(infor,molAtom,atomShell,den,scfConv);

		// now push the result density matrices into the result list
		atomDenMtrices.push_back(den);
	}
}

void AtomDenMtrx::writeToDisk() const
{
	// form the general path
	const string& scratchDir = infor.getScratchDir();
	path priPath(scratchDir.c_str());
	string sec = boost::lexical_cast<string>(section);
	path sect(sec.c_str());
	string denp = "atomdenmtrx";
	path denpath(denp.c_str());
	priPath /= sect;
	priPath /= denpath;

	// remove old data if there's same name folder exists
	// else create new one
	if (exists(priPath)) {
		remove_all(priPath);
	}
	create_directories(priPath);

	// for each type of atom, we construct the density matrix
	UInt nTypes = atomTypes.size();
	for(UInt i=0; i<nTypes; i++) {

		// the atomic number for this type
		UInt atomic = atomTypes[i];
		if (isGhostAtom(atomic)) continue;

		// form the path for each atom
		string atom = boost::lexical_cast<string>(atomic);
		path atomPath(atom.c_str());
		path p(priPath);
		p /= atomPath;

		// remove old data if there's same name folder exists
		// else create new one
		if (exists(p)) {
			remove_all(p);
		}
		create_directories(p);

		// now let's write the data
		const DenMtrx& d = getAtomDenMtrx(atomic);
		d.writeToDisk(p.string());
	}
}

void AtomDenMtrx::recover(const MolShell& ms) 
{
	// form the general path
	const string& scratchDir = infor.getScratchDir();
	path priPath(scratchDir.c_str());
	string s = boost::lexical_cast<string>(section);
	path sect(s.c_str());
	string denp = "atomdenmtrx";
	path denpath(denp.c_str());
	priPath /= sect;
	priPath /= denpath;

	// checking it's existance
	if (! exists(priPath)) {
		string info = "the corresponding atomic density matrix folder does not exist: " + priPath.string();
		Excep excep("AtomDenMtrx","constructor",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}

	// set up atom type information
	UInt nTypes = MAX_ATOMS_TYPE;
	atomTypes.reserve(nTypes);

	// here we will recover the atom information from the given file structure
	WordConvert w;
	for (directory_iterator iter(priPath); iter!=directory_iterator(); ++iter) {
		path p(iter->path());
		if (is_directory(p)) {
			// because the name uses atomic number
			// now let's retreieve it
			path name = p.filename();
			UInt atomic = -1;
			if (! w.toUInt(name.string(),atomic)) {
				string info = "failed to read in atom information from the given path: " + p.string();
				Excep excep("AtomDenMtrx","constructor",EXCEPTION_DIR_MISSING,info);
				handleExcep(excep);
			}
			atomTypes.push_back(atomic);
		}
	}

	// finally according to the atom types information
	// reserve density matrix and let's read it from hard disk
	nTypes = atomTypes.size();
	atomDenMtrices.reserve(nTypes);
	for(UInt i=0; i<nTypes; i++) {

		// the atomic number for this type
		UInt atomic = atomTypes[i];
		if (isGhostAtom(atomic)) continue;

		// set up the one atom molecule
		// also we update the section number
		Molecule molAtom(atomic);
		molAtom.updateSec(section);

		// now create the one atom shell molshell
		UInt iAtom = ms.getNAtomShells()+1;
		for(UInt j=0; j<ms.getNAtomShells(); j++) {
			const AtomShell& as = ms.getAtomShell(j);
			if (as.getAtomic() == atomic) {
				iAtom = j;
				break;
			}
		}
		if (iAtom == ms.getNAtomShells()+1) {
			string info = "can not find the corresponding atom shell data for the free atom(atomic number): " 
				+ boost::lexical_cast<string>(atomic);
			Excep excep("AtomDenMtrx","recover",EXCEPTION_ATOM_DENSITY_MATRICES_ERROR,info);
			handleExcep(excep);
		}
		MolShell atomShell(ms.getAtomShell(iAtom));

		// now build the density matrix
		DenMtrx den(infor,molAtom,atomShell,atomShell,molAtom.getNSpin());

		// form the path for each atom
		// now let's read in the data
		string atomicNum = boost::lexical_cast<string>(atomic);
		path atomPath(atomicNum.c_str());
		path p(priPath);
		p /= atomPath;
		den.recover(p.string());

		// finally let's push back the data
		atomDenMtrices.push_back(den);
	}
}

void AtomDenMtrx::averageSpinDenMtrx()
{
	UInt nAtomDenMtrx = atomTypes.size();
	for(UInt i=0; i<nAtomDenMtrx; i++) {
		DenMtrx& denMtrx = atomDenMtrices[i];
		if (denMtrx.getNSpin() == 2) {
			Mtrx& alpha = denMtrx.getMtrx(0);
			const Mtrx& beta = denMtrx.getMtrx(1);
			alpha.add(beta);
			alpha.scale(HALF);
		}
	}
}

const DenMtrx& AtomDenMtrx::getAtomDenMtrx(const UInt& atomic) const
{

	// obtain the position
	UInt pos = atomTypes.size()+1;
	for(UInt i=0; i<atomTypes.size(); i++) {
		if (atomTypes[i] == atomic) pos = i;
	}
	if (pos == atomTypes.size()+1) {
		string info = "can not find the corresponding density matrices for the free atom(atomic number): " 
			+ boost::lexical_cast<string>(atomic);
		Excep excep("AtomDenMtrx","formAtomDenMtrxInSCF",EXCEPTION_ATOM_DENSITY_MATRICES_ERROR,info);
		handleExcep(excep);
	}

	// return the result
	return atomDenMtrices[pos];
}

void AtomDenMtrx::print() const
{
	for(UInt i=0; i<atomTypes.size(); i++) {
		cout << "For Atom " << atomTypes[i] << endl;
		atomDenMtrices[i].print("density matrix");
	}
}

