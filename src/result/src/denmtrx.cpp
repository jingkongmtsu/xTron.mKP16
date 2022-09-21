/**
 * \file    denmtrx.cpp
 * \brief   describe the functions of density matrices except guess generation
 * \author  Fenglai Liu 
 */
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include "mo.h"
#include "shell.h"
#include "filerw.h"
#include "excep.h"
#include "textread.h"
#include "molecule.h"
#include "denmtrx.h"
using namespace boost::filesystem;
using namespace std;
using namespace mo;
using namespace shell;
using namespace filerw;
using namespace excep;
using namespace textread;
using namespace molecule;
using namespace denmtrx;

DenMtrx::DenMtrx(const GlobalInfor& infor0, const Molecule& mol, 
					const MolShell& rs, const MolShell& cs, 
					const UInt& nSpin):SpinMatrix(nSpin,rs.getNCarBas(),cs.getNCarBas()),
	infor(infor0),section(mol.getSec()),nAlpha(mol.getNAlpha()),nBeta(mol.getNBeta()) 
{
	// for density matrix, we need to set it to Cartesian dimension
	// since we may need to transform the data into the Cart. dimension
	// however, it's normal dimension is the nBas
	UInt nRow = rs.getNBas();
	UInt nCol = cs.getNBas();
	bool reformData = false;
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		mtrxList[iSpin].reset(nRow,nCol,reformData);
	}
}

void DenMtrx::formDenMtrx(const MO& mos) 
{
	// now form density matrix for each spin state
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {

		// additional checking 
		Mtrx& den      = getMtrx(iSpin);
		if(mos.getNBas(iSpin) != den.getRow() || mos.getNBas(iSpin) != den.getCol()) {
			string info = "density matrix dimension does not match the given mo";
			Excep excep("DenMtrx","formDenMtrx",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
			handleExcep(excep);
		}	
		if(mos.getNOcc(iSpin)>mos.getNOrb(iSpin)) {
			string info = "input mo has occupied orbital larger than the number of orbital?";
			Excep excep("DenMtrx","formDenMtrx",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
			handleExcep(excep);
		}

		// real work
		bool trans1    = false;
		bool trans2    = true;
		UInt nrow      = mos.getNBas(iSpin);
		UInt nocc      = mos.getNOcc(iSpin);
		UInt rowMO1Pos = 0;
		UInt colMO1Pos = 0;
		UInt rowMO2Pos = 0;
		UInt colMO2Pos = 0;
		const Mtrx& mo = mos.getMtrx(iSpin);
		den.mult(mo,mo,rowMO1Pos,colMO1Pos,rowMO2Pos,colMO2Pos,nrow,nocc,nrow,nocc,
				trans1,trans2,ONE,ZERO);
	}
}

void DenMtrx::writeToDisk() const
{
	// get the densitym matrix path
	const string& scratchDir = infor.getScratchDir();
	path p(scratchDir.c_str());
	string sec = boost::lexical_cast<string>(section);
	path sect(sec.c_str());
	string denp = "denmtrx";
	path denpath(denp.c_str());
	p /= sect;
	p /= denpath;

	// remove old data if there's same name folder exists
	// else create new one
	if (exists(p)) {
		remove_all(p);
	}
	create_directories(p);

	// firstly let's record the dimenion data so 
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
	string line = "# this file records dimension data for the given density matrix";
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
	}
}

void DenMtrx::writeToDisk(const string& denPath) const
{
	// remove old data if there's same name folder exists
	// else create new one
	path p(denPath.c_str());
	if (exists(p)) {
		remove_all(p);
	}
	create_directories(p);

	// firstly let's record the dimenion data so 
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
	string line = "# this file records dimension data for the given density matrix";
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
	}
}

void DenMtrx::recover(const string& denPath) 
{
	// mo dir location
	path p(denPath.c_str());

	// test that whether we have the mo data for the path
	if (! exists(p)) {
		string info = "the corresponding density matrix file does not exist: " + p.string();
		Excep excep("DenMtrx","recover",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}

	// let's read in the dimension data first
	// now let's see whether we have the record.txt?
	path folder(p);
	path record("record.txt");
	folder /= record;
	if (! exists(folder)) {
		string info = "record.txt does not exist for the density matrix class, it contains dimension data"; 
		Excep excep("DenMtrx","recover",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}
	string file = folder.string();

	// let's recover the data, firstly open it
	std::ifstream inf;
	const char * input_file = file.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string info = "can not open the record.txt for recovering density matrix data";
		Excep excep("DenMtrx","recover",EXCEPTION_FILE_MISSING,info);
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
					string info = "nRow reading from record.txt for the first density matrix is invalid";
					Excep excep("DenMtrx","recover",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
					handleExcep(excep);
				}
				if (! w.toUInt(v2,nCol1)) {
					string info = "nCol reading from record.txt for the first density matrix is invalid";
					Excep excep("DenMtrx","recover",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
					handleExcep(excep);
				}
			}else{
				if (! w.toUInt(v1,nRow2)) {
					string info = "nRow reading from record.txt for the second density matrix is invalid";
					Excep excep("DenMtrx","recover",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
					handleExcep(excep);
				}
				if (! w.toUInt(v2,nCol2)) {
					string info = "nCol reading from record.txt for the second density matrix is invalid";
					Excep excep("DenMtrx","recover",EXCEPTION_DENSITY_MATRICES_FORM_ERROR,info);
					handleExcep(excep);
				}
			}
		}
	}
	inf.close();

	// now let's possibly reset the density matrix dimension
	for(UInt iSpin=0; iSpin<getNSpin(); iSpin++) {
		if (iSpin == 0) {
			if (nRow1 != mtrxList[iSpin].getRow() || nCol1 != mtrxList[iSpin].getCol()) {
				bool reformData = false;
				mtrxList[iSpin].reset(nRow1,nCol1,reformData);
			}
		}else{
			if (nRow2 != mtrxList[iSpin].getRow() || nCol2 != mtrxList[iSpin].getCol()) {
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
		// the binary file must exist
		// else we need to report an error in fileRW
		//
		// we note that here we do not use the vector size
		// because the density matrix is initialized in Cartesian dimension,
		// see the constructor above. Therefore we use the row and column
		// instead
		DoubleVec& vec = mtrxList[iSpin].getVec();
		UInt len = mtrxList[iSpin].getRow()*mtrxList[iSpin].getCol();
		FileReadWrite fileRW1(p.string(),name);
		fileRW1.read(&vec.front(),len);
	}
}

