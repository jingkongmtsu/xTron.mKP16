//
// note:
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string> 
#include<vector>
#include<cmath>
#include<boost/lexical_cast.hpp>
#include "tbb/tbb.h"
#include "libgen.h"
#include "filerw.h"
#include "excep.h"
#include "textread.h"
#include "globalinfor.h"
#include "gintsinfor.h"
#include "molecule.h"
#include "matrix.h"
#include "denmtrx.h"
#include "shell.h"
#include "gints4d.h"
using namespace excep;
using namespace textread;
using namespace molecule;
using namespace filerw;
using namespace globalinfor;
using namespace gintsinfor;
using namespace shell;
using namespace matrix;
using namespace denmtrx;
using namespace gints4d;
using namespace tbb;
using namespace std;

void gints4d_test(bool doLargeJob)
{
	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "             GINTS 4D TESTING                 " << endl;
	if (doLargeJob) {
		cout << "       testing the timing performace          " << endl;
	}else{
		cout << "   testing the accuracy of 4D integrals etc.  " << endl;
	}
	cout << "----------------------------------------------" << endl;
	string inf = "gints4d/input.in";
	if (doLargeJob) inf = "gints4d/large.in";
	bool inDebug = true;
	if (doLargeJob) inDebug = false;

	// testing that how many jobs we have
	UInt nMol = 0;

	// open the file
	ifstream in;
	const char * input_file = inf.c_str();
	in.open(input_file,ios::in);
	if (!in) {
		string infor = "failed to open the input file";
		Excep excep("gints4d_test","large job test",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// counting that how many sections we have
	string key = "molecule";
	string line;
	WordConvert w;
	in.seekg(0,ios::beg);
	while(getline(in,line)) {
		LineParse l(line);
		if (l.isSec() && w.compare(l.getSecName(), key)) {
			nMol++;
		}
	}
	in.close();

	// now this is the real test section
	for(UInt iMol=1; iMol<=nMol; iMol++) {

		// form the molecule information etc.
		GlobalInfor globalInfor(inf);
		Molecule molecule(inf,iMol);
		MolShell s(inf,molecule);
		GIntsInfor infor(globalInfor,molecule);

		// set up density matrix
		UInt nBas = s.getNBas();
		DenMtrx den(globalInfor,molecule,s,s,molecule.getNSpin());

		// set density matrix
		// for large job, we use core guess
		// else we make all of density matrix to be one
		// that's what is in the standard data file
		if (doLargeJob) {
			//den.coreGuess(s,molecule);
		}else{
			Mtrx& D = den.getMtrx(0);
			string fileName = "gints4d/d" + boost::lexical_cast<string>(iMol) + ".txt";
			FileReadWrite fileRW(fileName);
			fileRW.readMatrixFromTextFile(D.getPtr(),D.getRow(),D.getCol());
		}

		// for large job, print out nBas
		if (doLargeJob) {
			cout << "nBas " << s.getNBas() << " nCarBas " << s.getNCarBas() << endl;
		}

		// set the job list
		UInt jobList[3];
		UInt nJobs = 3;
		jobList[0] = COULOMB;
		jobList[1] = EXCHANGE;
		jobList[2] = COULOMB_EXCHANGE;
		for(UInt iJob=0; iJob<nJobs; iJob++) {

			// set up results
			SpinMatrix M(molecule.getNSpin(),nBas,nBas);

			// work code
			// we do not test energy here
			Double energy = ZERO;
			UInt job = jobList[iJob];
			GInts4D gints(infor,job);
			gints.doJKMtrx(s,den,M,energy,true,false,false);

			// get the alpha matrix - actually we only
			// conpare the alpha part, since all of examples
			// here are for close shell
			Mtrx& M0 = M.getMtrx(0);
			M0.copyLowToUpper();

			// get the standard data file
			// and compare the results
			if (inDebug) {

				// get the name of the standard data file
				Mtrx N(nBas,nBas);
				if (job == COULOMB_EXCHANGE) {
					for(UInt x=0; x<2; x++) {
						string name = "j";
						if (x==1) name = "k";
						name = name + boost::lexical_cast<string>(iMol) + ".txt";
						string fileName = "gints4d/" + name;
						Mtrx tmp(nBas,nBas);
						FileReadWrite fileRW(fileName);
						fileRW.readMatrixFromTextFile(tmp.getPtr(),nBas,nBas);
						N.add(tmp);
					}
				}else{
					string name = "k";
					if (job == COULOMB) name = "j";
					name = name + boost::lexical_cast<string>(iMol) + ".txt";
					string fileName = "gints4d/" + name;
					FileReadWrite fileRW(fileName);
					fileRW.readMatrixFromTextFile(N.getPtr(),nBas,nBas);
				}

				// debug test
				/*
				M0.print("data matrix");
				N.print("reference matrix");
				*/

				// get the job name
				string jobName = "COULOMB";
				if (job == EXCHANGE) {
					jobName = "EXCHANGE";
				}else if (job == COULOMB_EXCHANGE) {
					jobName = "COULOMB_EXCHANGE";
				}

				// do the comparison work
				Double thresh = 1.0E-5;
				bool failed = M0.diffComp(N,thresh,false);
				if (failed) {
					cout << "the calculation of " << jobName << " for job " << iMol << " is failed with accuracy of " << thresh << endl;
				}else{
					cout << "the calculation of " << jobName << " for job " << iMol << " is passed with accuracy of " << thresh << endl;
				}
			}
			cout << endl;
		}
	}
}
