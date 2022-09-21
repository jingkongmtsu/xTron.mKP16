//
// set up testing file for exrho
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
#include "integraljobs.h"
#include "xcfunc.h"
#include "espints.h"
#include "denmtrx.h"
#include "sigshellpairinfor.h"
#include "xcints.h"
using namespace std;
using namespace excep;
using namespace filerw;
using namespace textread;
using namespace globalinfor;
using namespace molecule;
using namespace shell;
using namespace xcintsinfor;
using namespace integraljobs;
using namespace xcfunc;
using namespace espints;
using namespace denmtrx;
using namespace sigshellpairinfor;
using namespace xcints;

void exrho_test(bool doLargeJob, int jobTest, bool withMIC)
{
	// firstly let's see whether we do large job on exrho
	if (jobTest == 1) {

		// set input file
		string inf = "exrho/exrho_input.in";
		if (withMIC) {
			inf = "exrho/exrho_input_mic.in";
		}

		UInt nGrids = 48;
		DoubleVec grid(3*nGrids);
		for(UInt i=0; i<nGrids; i++) {
			Double t = (100.0E0/nGrids)*i;
			grid[3*i+0] = 0.1E0+t;
			grid[3*i+1] = 0.2E0+t;
			grid[3*i+2] = 0.3E0+t;
		}

		// set the job
		vector<UInt> jobList;
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

			UInt nMol  = 1;
			for(UInt iMol=1; iMol<=nMol; iMol++) {

				// set up basic information
				Molecule molecule(inf,iMol);
				MolShell s(inf,molecule);
				GIntsInfor infor(globalInfor,molecule);
				XCIntsInfor xcInfor(globalInfor,molecule);
				UInt nCarBas = s.getNCarBas();
				UInt nBas = s.getNBas();

				// set up denPhi
				// the data here is faked
				SpinMatrix denPhi(molecule.getNSpin(),nCarBas,nGrids);
				denPhi.initialize(ONE);

				// set up pMax infor
				// all of value is one
				//Mtrx pMaxInfor(s.getNShell(),nGrids);
				//pMaxInfor.set(ONE);

				// set up result matrix
				SpinMatrix result(molecule.getNSpin(),nGrids,nBas);
				cout << "geometry: "<< iMol << " with total basis sets " << nBas << endl;
				cout << "Cartesian number of basis sets: " << s.getNCarBas() << endl;

				// let's see what is the length of cart den etc.
				GIntJobInfor jobInfor(infor,job,order);
				SigMolShellPairInfor sp(s,s,jobInfor); 
				UInt len = sp.getSigSPVecLen(s,s,TYPE_CART);

				// set up Coulomb stuff
				DoubleVec cartDen(len,ONE);
				DoubleVec wRho(nGrids,ONE);
				DoubleVec halfJRhoVec(nGrids,ZERO);
				DoubleVec halfJRhoVecMtrx(len,ZERO);

				// initialize the MIC env
				XCFunc xcfunc(inf,molecule.getSec());
				XCInts xcints(globalInfor,xcfunc,s,xcInfor,infor,job,order);
				XCIntJobInfor xcJobInfor(globalInfor,infor,xcInfor,job,order);

				// working code
				// we only do exchange energy density right now
				{
					boost::timer::auto_cpu_timer t;
					ESPInts exrho(infor,job);
					exrho.doMtrx(xcJobInfor,s,sp,grid,denPhi,cartDen,wRho,result,halfJRhoVec,halfJRhoVecMtrx);
				}
			}

			// finally let's return
			return;
		}
	}

	// get the input file
	string inf = "exrho/input.in";
	if (doLargeJob) {
		inf = "exrho/large_input.in";
	}

	// let's see how many molecules in the input file
	// also we get in the standard energy for compare
	UInt nMol = 0;
	if (true) {

		// open the file
		ifstream in;
		const char * input_file = inf.c_str();
		in.open(input_file,ios::in);
		if (!in) {
			string infor = "failed to open the input file";
			Excep excep("xcints_test","exho_test",EXCEPTION_FILE_MISSING,infor);
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

	// get the ex energy 
	DoubleVec ex;
	ex.reserve(2*nMol);
	if (! doLargeJob) {

		// open the file
		ifstream in;
		const char * input_file = inf.c_str();
		in.open(input_file,ios::in);
		if (!in) {
			string infor = "failed to open the input file";
			Excep excep("xcints_test","main",EXCEPTION_FILE_MISSING,infor);
			handleExcep(excep);
		}

		// counting that how many sections we have
		WordConvert w;
		string line;
		while(getline(in,line)) {
			LineParse l(line);
			if (l.isCom() && line.find("energy") != string::npos) {
				for(UInt i=0; i<l.getNPieces(); i++) {
					string s= l.findValue(i);
					Double x = ZERO;
					if (w.toDouble(s,x)) {
						ex.push_back(x);
					}
				}
			}
		}
		in.close();
	}

	// print out testing heading
	cout << "----------------------------------------------" << endl;
	cout << "          XCINTS TESTING FOR EXRHO            " << endl;
	if (doLargeJob) {
		cout << "          FOR TIMING TEST ONLY                " << endl;
	}
	cout << "----------------------------------------------" << endl;
	UInt order = 0;
	UInt job   = GROUND_STATE_DFT;
	for(UInt iMol=1; iMol<=nMol; iMol++) {

		// heading
		cout << endl << endl;
		cout << "----------------------------------------------" << endl;
		cout << "   for molecule " << iMol                       << endl;            
		cout << "----------------------------------------------" << endl;
		
		// form the molecule information etc.
		GlobalInfor globalInfor(inf);
		Molecule molecule(inf,iMol);
		MolShell s(inf,molecule);
		XCIntsInfor infor(globalInfor,molecule);
		GIntsInfor  ginfor(globalInfor,molecule);
		MolShellSize molShellSize(s,infor.getThresh());
		//molShellSize.print();
		XCFunc xcfunc(inf,molecule.getSec());
		XCInts xcints(globalInfor,xcfunc,s,infor,ginfor,job,order);

		// set up density matrix
		// all of density matrix value are set to one
		DenMtrx denMtrx(globalInfor,molecule,s,s,molecule.getNSpin()); 
		denMtrx.initialize(ONE);

		// set up result - this is actually not used
		// we only care about the exchange energy calculated
		SpinMatrix result(molecule.getNSpin(),s.getNBas(),s.getNBas());

		// now do calculation
		// the result exchange energy will be 
		if (true) {
			boost::timer::auto_cpu_timer t;
			xcints.doXCMtrx(s,denMtrx,result,false,false);
		}
		if (! doLargeJob) {
			Double ae = ex[2*(iMol-1)];
			Double be = ex[2*(iMol-1)+1];
			printf("standard alpha exchange energy: %-20.14f\n", ae);
			printf("standard beta  exchange energy: %-20.14f\n", be);
		}
	}
}
