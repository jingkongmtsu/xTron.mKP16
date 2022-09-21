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
#include "gints4deriv.h"
using namespace excep;
using namespace textread;
using namespace molecule;
using namespace filerw;
using namespace globalinfor;
using namespace gintsinfor;
using namespace shell;
using namespace matrix;
using namespace denmtrx;
using namespace gints4deriv;
using namespace tbb;
using namespace std;

void gints4deriv_test(bool doLargeJob)
{
	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "           GINTS 4D DERIV TESTING             " << endl;
	if (doLargeJob) {
		cout << "       testing the timing performace          " << endl;
	}else{
		cout << "   testing the accuracy of 4D integrals deriv " << endl;
	}
	cout << "----------------------------------------------" << endl;
	string folder = "gints4deriv/large/";
	if (! doLargeJob) folder = "gints4deriv/";
	string inf = folder + "input.in";

	// testing that how many jobs we have
	UInt nMol = 0;
	ifstream in;
	const char * input_file = inf.c_str();
	in.open(input_file,ios::in);
	if (!in) {
		string infor = "failed to open the input file";
		Excep excep("gints4deriv_test","input molecules count",EXCEPTION_FILE_MISSING,infor);
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

		// for large job, print out nBas
		if (doLargeJob) {
			cout << "nBas " << s.getNBas() << " nCarBas " << s.getNCarBas() << endl;
		}

		// set density matrix
		// for large job, just make up a guess
		// for accurary test, we read in the guess
		if (doLargeJob) {
			den.coreGuess(s,molecule);
		}else{
			for(UInt iSpin=0; iSpin<molecule.getNSpin(); iSpin++) {
				Mtrx& D = den.getMtrx(iSpin);
				string fileName = folder + "d" + boost::lexical_cast<string>(iMol) + ".txt";
				if (molecule.getNSpin() == 2) {
					fileName = folder + "d" + boost::lexical_cast<string>(iMol) + "_a.txt";
					if (iSpin == 1) {
						fileName = folder + "d" + boost::lexical_cast<string>(iMol) + "_b.txt";
					}
				}
				FileReadWrite fileRW(fileName);
				fileRW.readMatrixFromTextFile(D.getPtr(),D.getRow(),D.getCol());
			}
		}

		// set the job list
		UInt jobList[3];
		UInt nJobs = 3;
		jobList[0] = COULOMB;
		jobList[1] = EXCHANGE;
		jobList[2] = COULOMB_EXCHANGE;
		for(UInt order=1; order<=2; order++) {

			// set up results
			Mtrx M;
			if (order == 1) {
				M.init(molecule.getNAtoms(),3);
			}else if (order == 2) {
				M.init(molecule.getNAtoms()*3, molecule.getNAtoms()*3);
			}

			for(UInt iJob=0; iJob<nJobs; iJob++) {

				// work code
				// we do not test energy here
				UInt job = jobList[iJob];
				GInts4Deriv gints(infor,job,order);
				M.set(ZERO);
				gints.doJKDeriv(s,molecule,den,M,true,false);
				if (order == 2) {
					M.copyLowToUpper();
				}
				M.print("result derivatives");

				// set the standard data
				Mtrx N;
				UInt nAtoms = molecule.getNAtoms();
				if (order == 1) {
					N.init(nAtoms,3);
				}else if (order == 2) {
					N.init(3*nAtoms,3*nAtoms);
				}

				// get the standard data file
				// and compare the results
				// get the name of the standard data file
				UInt nFiles = 1;
				if (job == COULOMB_EXCHANGE) nFiles = 2;
				for(UInt iFile=0; iFile<nFiles; iFile++) {

					// set the file name
					string name = "k";
					if (job == COULOMB) name = "j";
					if (iFile==1) name = "j"; 
					if (order == 1) {
						name = name + boost::lexical_cast<string>(iMol) + "g.txt";
					}else if (order == 2) {
						name = name + boost::lexical_cast<string>(iMol) + "h.txt";
					}
					string fileName = folder + name;

					// set temp matrix
					Mtrx N0;
					if (order == 1) {
						if (! doLargeJob) {
							N0.init(3,nAtoms);
							FileReadWrite fileRW(fileName);
							fileRW.readMatrixFromTextFile(N0.getPtr(),3,nAtoms);
							N0.transpose(true);
						}else{
							N0.init(nAtoms,3);
							FileReadWrite fileRW(fileName);
							fileRW.readMatrixFromTextFile(N0.getPtr(),nAtoms,3);
						}
					}else if (order == 2) {
						N0.init(3*nAtoms,3*nAtoms);
						FileReadWrite fileRW(fileName);
						fileRW.readMatrixFromTextFile(N0.getPtr(),3*nAtoms,3*nAtoms);
						N0.copyLowToUpper();
					}

					// add it to N
					N.add(N0);
				}

				// set job name
				string name = "COULOMB";
				if (job == EXCHANGE) {
					name = "EXCHANGE";
				}else if (job == COULOMB_EXCHANGE) {
					name = "COULOMB_EXCHANGE";
				}

				// do the comparison work
				Double thresh = 1.0E-6;
				if (doLargeJob) thresh = 1.0E-5;
				bool failed = M.diffComp(N,thresh,false);
				if (failed) {
					cout << "the calculation of " << name << " is failed with accuracy of " << thresh << endl;
				}else{
					cout << "the calculation of " << name << " is passed with accuracy of " << thresh << endl;
				}
			}
			cout << endl;
		}
	}
}
