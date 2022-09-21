//
// testing the xc matrix 
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string> 
#include<vector>
#include<cmath>
#include<boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include "excep.h"
#include "filerw.h"
#include "textread.h"
#include "globalinfor.h"
#include "molecule.h"
#include "shell.h"
#include "xcintsinfor.h"
#include "gintsinfor.h"
#include "integraljobs.h"
#include "filerw.h"
#include "xcfunc.h"
#include "denmtrx.h"
#include "xcints.h"
using namespace std;
using namespace excep;
using namespace filerw;
using namespace textread;
using namespace globalinfor;
using namespace molecule;
using namespace shell;
using namespace xcintsinfor;
using namespace gintsinfor;
using namespace integraljobs;
using namespace filerw;
using namespace xcfunc;
using namespace denmtrx;
using namespace xcints;

void xcmtrx_test(bool testTiming, bool withMIC)
{
	// get the input file
	string inf = "xcmtrx/input.in";
	if (testTiming) {
		inf = "xcmtrx/timing.in";
		if (withMIC) {
			inf = "xcmtrx/timing_mic.in";
		}
	}else{
		if (withMIC) {
			inf = "xcmtrx/input_mic.in";
		}
	}

	// let's see how many molecules in the input file
	// also we get in the standard energy for compare
	UInt nMol = 0;
	bool doJob = true;
	if (doJob) {

		// open the file
		ifstream in;
		const char * input_file = inf.c_str();
		in.open(input_file,ios::in);
		if (!in) {
			string infor = "failed to open the input file";
			Excep excep("xcints_test","xcmtrx",EXCEPTION_FILE_MISSING,infor);
			handleExcep(excep);
		}

		// counting that how many sections we have
		string key = "molecule";
		nMol = 0;
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
	}

	// print out testing heading
	cout << "----------------------------------------------" << endl;
	cout << "          XCINTS TESTING FOR XCMTRX           " << endl;
	cout << "----------------------------------------------" << endl;

	// set the jobs;
	vector<UInt> jobList;
	jobList.push_back(GROUND_STATE_DFT);
	jobList.push_back(EXCHANGE);
	UInt order = 0;

	// now loop over jobs
	GlobalInfor globalInfor(inf);
	for(UInt iJob=0; iJob<jobList.size(); iJob++) {
		UInt job = jobList[iJob];

		// print heading
		cout << endl << endl;
		if (job == GROUND_STATE_DFT) {
			cout << "----------------------------------------------" << endl;
			cout << " TESTING FOR GROUND_STATE_DFT JOB             " << endl;            
			cout << "----------------------------------------------" << endl;
		}else if (job == EXCHANGE) {
			cout << "----------------------------------------------" << endl;
			cout << " TESTING FOR EXCHANGE JOB                     " << endl;            
			cout << "----------------------------------------------" << endl;
		}else if (job == COULOMB) {
			cout << "----------------------------------------------" << endl;
			cout << " TESTING FOR COULOMB JOB                      " << endl;            
			cout << "----------------------------------------------" << endl;
		}
		cout << endl << endl;

		for(UInt iMol=1; iMol<=nMol; iMol++) {

			// heading
			cout << endl << endl;
			cout << "----------------------------------------------" << endl;
			cout << "   for molecule " << iMol                       << endl;            
			cout << "----------------------------------------------" << endl;

			// form the molecule information etc.
			Molecule molecule(inf,iMol);
			MolShell s(inf,molecule);
			printf("number of basis sets %d\n", s.getNBas());
			XCIntsInfor infor(globalInfor,molecule);
			GIntsInfor ginfor(globalInfor,molecule);
			MolShellSize molShellSize(s,infor.getThresh());
			//molShellSize.print();
			XCFunc xcfunc(inf,molecule.getSec());
			XCInts xcints(globalInfor,xcfunc,s,infor,ginfor,job,order);

			// set up density matrix
			DenMtrx denMtrx(globalInfor,molecule,s,s,molecule.getNSpin()); 
			//if (testTiming) {
			//	Double val = ONE/s.getNBas();
			//	denMtrx.initialize(val);
			//}else{
				denMtrx.coreGuess(s,molecule);
			//}

			// set up result
			SpinMatrix result(molecule.getNSpin(),s.getNBas(),s.getNBas());

			// now do calculation
			boost::timer::auto_cpu_timer t;
			if (isDFTJob(job)) {
				xcints.doXCMtrx(s,denMtrx,result,true,false);
			}else{
				xcints.doJKMtrx(s,denMtrx,result,true,false);
			}

			// after the calculation, we need to do the test
			if (! testTiming) {
				SpinMatrix result0(molecule.getNSpin(),s.getNBas(),s.getNBas());
				for(UInt i=0; i<molecule.getNSpin(); i++) {
					Mtrx&       M0 = result0.getMtrx(i);
					const Mtrx& M1 = result.getMtrx(i);

					// set the name tag
					string tag = "xcmtrx";
					if (job == EXCHANGE) tag = "k";
					if (job == COULOMB)  tag = "j";
					if (job == COULOMB_EXCHANGE) tag = "jk";

					// read in data for M0
					string name = "xcmtrx/" + tag + boost::lexical_cast<string>(iMol) +".txt";
					FileReadWrite fileRW(name);
					fileRW.readMatrixFromTextFile(M0.getPtr(), M0.getRow(), M0.getCol());

					// now let's do compare
					Double thresh = 1.0E-5;
					bool failed = M1.diffComp(M0,thresh,true,true);
					if (failed) {
						cout << "the xcmtrx calculation for geometry " << iMol << " is failed with accuracy of " << thresh << endl;
					}else{
						cout << "the xcmtrx calculation for geometry " << iMol << " is successful with accuracy of " << thresh << endl;
					}
				}
			}
		}
	}
}
