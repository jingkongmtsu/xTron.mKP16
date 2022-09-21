/**
 * \file    mo.cpp
 * \brief   describing the MO generated in the SCF procedure
 * \author  Fenglai Liu and Jing Kong
 */
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
#include "filerw.h"
#include "fracspininfor.h"
#include "shell.h"
#include "gints2d.h"
#include "oneemtrx.h"
#include "gintsinfor.h"
#include "molecule.h"
#include "filerw.h"
#include "textread.h"
#include "matrix.h"
#include "blas.h"
#include "blas1.h"
#include "mo.h"
using namespace boost::filesystem;
using namespace filerw;
using namespace fracspininfor;
using namespace shell;
using namespace gints2d;
using namespace oneemtrx;
using namespace gintsinfor;
using namespace matrix;
using namespace molecule;
using namespace textread;
using namespace blas;
using namespace mo;
using namespace std;

MO::MO(const GlobalInfor& infor0, const MolShell& ms, 
		const Molecule& mol, const UInt& nSpin):SpinMatrix(nSpin,ms.getNBas(),ms.getNBas()),
	infor(infor0),section(mol.getSec()),shellCode(ms.getCode()),
	nAlpha(mol.getNEle(0)),nBeta(mol.getNEle(1)),alphaEnergy(ms.getNBas(),ZERO) 
{ 
	// do we need to initilize beta part?
	if (nSpin == 2) {
		betaEnergy.assign(ms.getNBas(),ZERO);
	}
}

void MO::solveFock(const UInt& iSpin, const Mtrx& ortho, Mtrx& Fock)
{ 
	// get the corresponding MO
	Mtrx& mo = getMtrx(iSpin);

	// additional check for MO dimension
	if(Fock.getRow() != Fock.getCol()) {
		string info = "input Fock matrix must be square, row == col";
		Excep excep("MO","solveFock",EXCEPTION_FOCK_SOLVING_ERROR,info);
		handleExcep(excep);
	}	
	if(Fock.getRow() != mo.getRow()) {
		string info = "Fock matrix dimension does not match the mo dimension";
		Excep excep("MO","solveFock",EXCEPTION_FOCK_SOLVING_ERROR,info);
		handleExcep(excep);
	}

	// mornitoring that whether the Fock matrix 
	// is with redundant basis set
	// then we reform the number of orbital - MO's column
	bool withRedundantBasis = false;
	if (ortho.getCol()!=ortho.getRow()) withRedundantBasis = true;
	UInt nMO  = ortho.getCol();
	UInt nBas = Fock.getRow();
	if (nMO > mo.getCol()) {
		string info = "if with redundant basis set, why nMO from orthogonal matrix > MO dimension?";
		Excep excep("MO","solveFock",EXCEPTION_FOCK_SOLVING_ERROR,info);
		handleExcep(excep);
	}

	// reform the data
	if (mo.getCol() != nMO) {
		bool reformData = false;
		mo.reset(nBas,nMO,reformData);
	}

	// transform the Fock matrix with orhogonal matrix
	// we note that only lower triangular part of Fock is used
	// so the Fock matrix must be symmetrical here
	// F' = S^{T}FS
	//Fock.print("original fock matrix");
	mo.symMatrixMult(Fock,ortho,'L',ONE,ZERO);
	//mo.print("FS in solving mo", 5);
	if (withRedundantBasis) {
		bool reformData = false;
		Fock.reset(nMO,nMO,reformData);
	}
	Fock.mult(ortho,mo,true,false,ONE,ZERO);
	//ortho.print("orthogonal matrix used in solving the Fock matrix",5);
	//Fock.print("final fock matrix S^{T}FS before diagonalize",5);

	// diagonalize F'
	if (iSpin == 0) {
		alphaEnergy.assign(nMO,ZERO);
		Fock.doEigen(alphaEnergy);
	}else{
		betaEnergy.assign(nMO,ZERO);
		Fock.doEigen(betaEnergy);
	}
	//Fock.print("initial mo",5);

	// form MO coefficients
	mo.mult(ortho,Fock,false,false,ONE,ZERO);
	//mo.print("result mo",5);

	// finally we need to restore the Fock's dimension
	if (withRedundantBasis) {
		bool reformData = true;
		Fock.reset(nBas,nBas,reformData);
	}
}

void MO::formMO(const Mtrx& ortho, Mtrx& Fock)
{
	// solve Fock matrix for alpha part
	solveFock(0,ortho,Fock);

	// copy beta part information
	if (getNSpin() == 2) {
		const Mtrx& alpha = getMtrx(0);
		Mtrx& beta = getMtrx(1);
		beta = alpha;
		betaEnergy = alphaEnergy;
	}
}

void MO::formMO(const MolShell& ms, const Molecule& mol, const OneEMtrx& oneEMtrx, SpinMatrix& Fock)
{
	// get the ortho matrix
	// we note that oneEMtrx stores the 
	// ortho matrix as a full matrix
	UInt nRow,nCol; 
	oneEMtrx.getDim(ORTHOGONAL_MATRIX,nRow,nCol);
	Mtrx ortho(nRow,nCol);
	oneEMtrx.getM(ORTHOGONAL_MATRIX,ortho);

	// solve Fock matrix for each part
	for(UInt iSpin=0; iSpin<getNSpin();iSpin++) {
		solveFock(iSpin,ortho,Fock.getMtrx(iSpin));
	}
}

void MO::writeToDisk() const
{
	// get the mo path 
	const string& scratchDir = infor.getScratchDir();
	path p(scratchDir.c_str());
	string sec = boost::lexical_cast<string>(section);
	path sect(sec.c_str());
	string mop = "mo";
	path mopath(mop.c_str());
	p /= sect;
	p /= mopath;

	// remove old data if there's same name folder exists
	// else create new one
	if (exists(p)) {
		remove_all(p);
	}
	create_directories(p);

	// firstly let's record the mo dimenion data so 
	// that when recovering we can double check it
	path folder(p);
	path record("record.txt");
	folder /= record;
	string file = folder.string();

	// open the output stream
	std::ofstream outf;
	const char * output_file = file.c_str();
	outf.open(output_file,std::ofstream::out);

	// print out the comments and data
	// print out the data
	string line = "# this file records dimension data for the given mo";
	outf << line << endl;
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		line = boost::lexical_cast<string>(mtrxList[iSpin].getRow()) + "  ";
		line = line + boost::lexical_cast<string>(mtrxList[iSpin].getCol());
		outf << line << endl;
	}
	outf.close();

	// now let's begin to write current MO data
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {

		// mo file name
		string name = "a.bin";
		if (iSpin == 1) {
			name = "b.bin";
		}

		// now do file writing
		const DoubleVec& vec = mtrxList[iSpin].getVec();
		FileReadWrite fileRW1(p.string(),name);
		fileRW1.write(&vec.front(),vec.size());

		// mo energy file name
		name = "ea.bin";
		if (iSpin == 1) {
			name = "eb.bin";
		}

		// now do file writing
		FileReadWrite fileRW2(p.string(),name);
		if (iSpin == 0) {
			fileRW2.write(&alphaEnergy.front(),alphaEnergy.size());
		}else{
			fileRW2.write(&betaEnergy.front(),betaEnergy.size());
		}
	}
}

void MO::recover(const string& moPath) 
{
	// mo dir location
	path p(moPath.c_str());

	// test that whether we have the mo data for the path
	if (! exists(p)) {
		string info = "the corresponding mo file does not exist: " + p.string();
		Excep excep("MO","recover",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}

	// let's read in the mo dimension data first
	// now let's see whether we have the record.txt?
	path folder(p);
	path record("record.txt");
	folder /= record;
	if (! exists(folder)) {
		string info = "record.txt does not exist for the mo, it contains MO dimension data"; 
		Excep excep("MO","recover",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}
	string file = folder.string();

	// let's recover the data, firstly open it
	std::ifstream inf;
	const char * input_file = file.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string info = "can not open the record.txt for recovering mo data";
		Excep excep("MO","recover",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}

	// set up dimension
	UInt nRow1 = 0;
	UInt nCol1 = 0;
	UInt nRow2 = 0;
	UInt nCol2 = 0;

	// recover the data, get the nSpin state
	string line;
	WordConvert w;
	inf.seekg(0,ios::beg);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isCom() || l.isEmp()) continue;
		if (l.getNPieces() == 2) {
			string v1  = l.findValue(0);
			string v2  = l.findValue(1);
			if (nRow1 == 0) {
				if (! w.toUInt(v1,nRow1)) {
					string info = "nRow reading from record.txt for the first mo is invalid";
					Excep excep("MO","recover",EXCEPTION_MO_ERROR,info);
					handleExcep(excep);
				}
				if (! w.toUInt(v2,nCol1)) {
					string info = "nCol reading from record.txt for the first mo is invalid";
					Excep excep("MO","recover",EXCEPTION_MO_ERROR,info);
					handleExcep(excep);
				}
			}else{
				if (! w.toUInt(v1,nRow2)) {
					string info = "nRow reading from record.txt for the second mo is invalid";
					Excep excep("MO","recover",EXCEPTION_MO_ERROR,info);
					handleExcep(excep);
				}
				if (! w.toUInt(v2,nCol2)) {
					string info = "nCol reading from record.txt for the second mo is invalid";
					Excep excep("MO","recover",EXCEPTION_MO_ERROR,info);
					handleExcep(excep);
				}
			}
		}
	}
	inf.close();

	// now let's possibly reset the mo dimension
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		if (iSpin == 0) {
			if (nRow1 != mtrxList[iSpin].getRow()) {
				string info = "something wrong with the MO basis set dimension, it should not changed. " 
					"We fould dimension reading from record.txt does not match the current basis dimension";
				Excep excep("MO","recover",EXCEPTION_MO_ERROR,info);
				handleExcep(excep);
			}
			if (nCol1 != mtrxList[iSpin].getCol()) {
				bool reformData = false;
				mtrxList[iSpin].reset(nRow1,nCol1,reformData);
			}
		}else{
			if (nRow2 != mtrxList[iSpin].getRow()) {
				string info = "something wrong with the MO basis set dimension, it should not changed. " 
					"We fould dimension reading from record.txt does not match the current basis dimension";
				Excep excep("MO","recover",EXCEPTION_MO_ERROR,info);
				handleExcep(excep);
			}
			if (nCol2 != mtrxList[iSpin].getCol()) {
				bool reformData = false;
				mtrxList[iSpin].reset(nRow2,nCol2,reformData);
			}
		}
	}

	// now begin to read in
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {

		// mo file name
		string name = "a.bin";
		if (iSpin == 1) {
			name = "b.bin";
		}

		// now do file writing
		// the mo file must exist
		// else we need to report an error
		DoubleVec& vec = mtrxList[iSpin].getVec();
		FileReadWrite fileRW1(p.string(),name);
		fileRW1.read(&vec.front(),vec.size());

		// mo energy file name
		name = "ea.bin";
		if (iSpin == 1) {
			name = "eb.bin";
		}

		// the mo energy data is not necessary
		// so if we have, we will use the data file
		// if we do not have it, we will initialize
		// the data to be zero
		string emoPath = moPath + "/" + name;
		path ep(emoPath.c_str());
		if (! exists(ep)) {
			if (iSpin == 0) {
				alphaEnergy.assign(alphaEnergy.size(),ZERO);
			}else{
				betaEnergy.assign(betaEnergy.size(),ZERO);
			}
		}else{
			FileReadWrite fileRW2(p.string(),name);
			if (iSpin == 0) {
				fileRW2.read(&alphaEnergy.front(),alphaEnergy.size());
			}else{
				fileRW2.read(&betaEnergy.front(),betaEnergy.size());
			}
		}
	}
}

Double MO::getLUMOHOMODiff(UInt iSpin) const 
{
	// get the homo and lumo index
	UInt iLUMO = getLUMOIndex(iSpin);
	UInt iHOMO = iLUMO-1;

	// we may have a problem that 
	// we do not have virtual orbitals
	if (iLUMO>=getNOrb(iSpin)) {
		string info = "because no virtual orbitals in result mo, we can not get HOMO LUMO difference";
		Excep excep("MO","getLUMOHOMODiff",EXCEPTION_NO_VIRTUAL_ORBITAL,info);
		handleExcep(excep);
		return ZERO;
	}

	// let's get the E
	Double EHOMO = ZERO;
	Double ELUMO = ZERO;
	if (iSpin == 0) {
		EHOMO = alphaEnergy[iHOMO];
		ELUMO = alphaEnergy[iLUMO];
	}else{
		EHOMO = betaEnergy[iHOMO];
		ELUMO = betaEnergy[iLUMO];
	}
	return (ELUMO - EHOMO);
}

void MO::formFracSpinMO(const FracSpinInfor& infor)
{
	// now let's do the scaling work for each spin state
	for(UInt iSpin=0; iSpin<2; iSpin++) {

		// now get the value
		Int nScaledMO = infor.getNScaledMO(iSpin);
		Int begin     = infor.getMOBeginIndex(iSpin);

		// whether we do scaling for this spin state?
		if (nScaledMO == 0) continue;

		// now get the mo data
		Mtrx& spinMO = getMtrx(iSpin);

		// check the fractional spin information
		if (nScaledMO>(Int)spinMO.getCol()) {
			string info0 = "for the spin state alpha, ";
			if (iSpin == 1) {
				info0 = "for the spin state beta, ";
			}
			string info  = info0 + "the number of scaled mo > all of mo number, invalid situation";
			Excep excep("MO","formFracSpinMO",EXCEPTION_INVALID_CASE_FRAC_SPIN,info);
			handleExcep(excep);
		}

		// check the fractional spin information
		if (begin>=(Int)spinMO.getCol()) {
			string info0 = "for the spin state alpha, ";
			if (iSpin == 1) {
				info0 = "for the spin state beta, ";
			}
			string info  = info0 + "the top mo index for scaling is invalid, it's larger than the possible mo index";
			Excep excep("MO","formFracSpinMO",EXCEPTION_INVALID_CASE_FRAC_SPIN,info);
			handleExcep(excep);
		}

		// check the fractional spin information
		if (begin+1<nScaledMO) {
			string info0 = "for the spin state alpha, ";
			if (iSpin == 1) {
				info0 = "for the spin state beta, ";
			}
			string info  = info0 + "the top mo index for scaling is less than the number of scaling mo, invalid situation";
			Excep excep("MO","formFracSpinMO",EXCEPTION_INVALID_CASE_FRAC_SPIN,info);
			handleExcep(excep);
		}

		// now let's print out the information
		if (iSpin == 0) {
			printf("for alpha spin mo\n");
		}else{
			printf("for beta  spin mo\n");
		}
		printf("mo begin index: %d, number of mos for scaled: %d\n", begin, nScaledMO);
		Int end = begin - nScaledMO;
		for(Int iMO=begin; iMO>end; iMO--) {
			Double val = infor.getScaleVal(iSpin,begin-iMO);
			vscal(spinMO.getPtr(0,iMO),val,spinMO.getCol());
		}   
	}

	// now let's modify the other information
	const Molecule& newMol = infor.getMol();
	nAlpha = newMol.getNAlpha();
	nBeta  = newMol.getNBeta();
	printf("number of electrons are reset, new alpha electron occupation number %d, "
			"and new beta electron occupation number %d\n", (Int)nAlpha, (Int)nBeta);
	printf("please note that the occupation number includes the fractional spin "
			"so that occulation number does not equal to the original number of electrons\n");

	// remove all of energy data
	alphaEnergy.assign(alphaEnergy.size(),ZERO);
	betaEnergy.assign(betaEnergy.size(),ZERO);

}

void MO::print(UInt level) const 
{
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {

		if (level == 1 || level == 3) {
			cout << "number of occupied mo: " << getNOcc(iSpin) << endl;
			cout << "eigen value of MO: " << endl;
			const DoubleVec& e = getEnergyVec(iSpin);

			// now print out the eigen values
			UInt nMOinBatch = 6;
			UInt nMOBatch   = e.size()/nMOinBatch;
			if (e.size() -nMOBatch*nMOinBatch != 0) nMOBatch += 1;
			for(UInt iBatch=0; iBatch<nMOBatch; iBatch++) {

				// calculate the mo begin and mo end index for this batch
				UInt moBegin = iBatch*nMOinBatch;
				UInt moEnd   = (iBatch+1)*nMOinBatch;
				if (moEnd>e.size()) moEnd = e.size();
				for(UInt iMO=moBegin; iMO<moEnd; iMO++) {
					printf("%-12d", (Int)iMO);
				}
				printf("\n");
				for(UInt iMO=moBegin; iMO<moEnd; iMO++) {
					printf("%-12.6f", e[iMO]);
				}
				printf("\n");
			}
		}

		// now print out the whole MO data
		if (level==2 || level ==3) {
			string t;
			if (iSpin == 0) {
				t = "Alpha MO Coefficients: ";
			}else{
				t = "Beta MO Coefficients: ";
			}	
			const Mtrx& M = mtrxList[iSpin];
			M.print(t);
		}
	}
}

