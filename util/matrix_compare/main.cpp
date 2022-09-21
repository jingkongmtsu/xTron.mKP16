//
// testing the xc matrix 
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string> 
#include<boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "blas1.h"
#include "excep.h"
#include "filerw.h"
#include "textread.h"
#include "globalinfor.h"
#include "molecule.h"
#include "shell.h"
#include "filerw.h"
#include "spinmatrix.h"
using namespace std;
using namespace boost::filesystem;
using namespace blas;
using namespace excep;
using namespace filerw;
using namespace textread;
using namespace globalinfor;
using namespace molecule;
using namespace shell;
using namespace filerw;
using namespace spinmatrix;

int main(int argc, char* argv[])
{
	// read in the parameter for input file and output file 
	string inf  = argv[1];
	string outf = argv[2];
	freopen (outf.c_str(),"w",stdout);

	// also read in other parameters
	string p1   = argv[3];
	string p2   = argv[4];

	// do we have density matrix path also defined?
	bool doDigestion = false;
	string p3;
	if (argc>5) {
		doDigestion = true;
		p3 = argv[5];
	}

	// set the matrix
	// result 1 is the standard data for comparison
	// result 2 is the testing matrix data
	GlobalInfor glob1(inf);
	Molecule mol(inf,1);
	MolShell s(inf,mol);
	printf("number of basis sets %d\n", s.getNBas());
	SpinMatrix result1(mol.getNSpin(),s.getNBas(),s.getNBas());
	SpinMatrix result2(mol.getNSpin(),s.getNBas(),s.getNBas());

	// get the path and read in the data
	// you need to choose which one is the folder you saved 
	string name = "a.bin";	
	UInt found = p1.find("gints");
	bool standard_data_is_gints = false;
	if (found != std::string::npos) {

		// let's see which one exist
		bool gotIt = false;
		path d(p1.c_str());
		for(UInt i=0; i<3; i++) {

			// set the name
			// only one file can exist
			name = "COULOMB_EXCHANGE_INTEGRALS_a.bin";
			if (i == 1) name = "EXCHANGE_INTEGRALS_a.bin";
			if (i == 2) name = "COULOMB_INTEGRALS_a.bin";
			path d1(name.c_str());
			path d2(d);
			d2 /= d1;

			// test that whether we have the mo data for the path
			if (! exists(d2)) {
				continue;
			}else{
				gotIt = true;
				break;
			}
		}

		// whether we issue an error
		if (! gotIt) {
			string info = "the corresponding gints data file does not exist: " + d.string();
			Excep excep("util","matrix_compare",EXCEPTION_FILE_MISSING,info);
			handleExcep(excep);
		}
		standard_data_is_gints = true;
	}
	string path1= p1 + "/" + name;
	FileReadWrite fileRW1(path1);
	Mtrx& M1 = result1.getMtrx(0);
	fileRW1.read(M1.getPtr(0,0),s.getNBas()*s.getNBas());

	// get the path and read in the data
	name = "a.bin";	
	string path2= p2 + "/" + name;
	FileReadWrite fileRW2(path2);
	Mtrx& M2 = result2.getMtrx(0);
	fileRW2.read(M2.getPtr(0,0),s.getNBas()*s.getNBas());

	// read in the density matrix data
	// if you want to read in density matrix, the data should be
	// in same dimension of the results
	if (doDigestion) {

		// set the density matrix data
		SpinMatrix denMtrx(mol.getNSpin(),s.getNBas(),s.getNBas());
		name = "a.bin";	
		string path3 = p3 + "/" + name;
		FileReadWrite fileRW3(path3);
		Mtrx& M3 = denMtrx.getMtrx(0);
		fileRW3.read(M3.getPtr(0,0),s.getNBas()*s.getNBas());

		// now let's calculate the energy 1
		bool symmMatrix = true;
		const Mtrx& alphaDenMtrx = denMtrx.getMtrx(0);
		Double eJK1 = M1.dotProduct(alphaDenMtrx,symmMatrix);
		Double eJK2 = M2.dotProduct(alphaDenMtrx,symmMatrix);
		if (standard_data_is_gints) {
			printf("the gints4d energy is %-16.10f\n", eJK1);
		}else{
			printf("the standard xc energy is %-16.10f\n", eJK1);
		}
		printf("the xc energy is %-16.10f\n", eJK2);
	}

	// compare the difference
	Double thresh = 1.0E-6;
	bool failed = M1.diffComp(M2,thresh,true,true);

}
