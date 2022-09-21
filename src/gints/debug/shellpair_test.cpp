#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include<cmath>
#include "tbb/tbb.h"
#include "libgen.h"
#include "globalinfor.h"
#include "gintsinfor.h"
#include "molecule.h"
#include "shell.h"
#include "denmtrx.h"
#include "denmtrxinfor.h"
#include "sigshellpairinfor.h"
#include "shellpair.h"
using namespace molecule;
using namespace globalinfor;
using namespace gintsinfor;
using namespace shell;
using namespace denmtrx;
using namespace denmtrxinfor;
using namespace sigshellpairinfor;
using namespace shellpair;
using namespace tbb;
using namespace std;

void shellpair_test(bool doLargeJob)
{
	cout << "----------------------------------------------" << endl;
	cout << "SHELL PAIR INFOR DATA TESTING                 " << endl;
	if (doLargeJob) {
		cout << "testing the timing performace                 " << endl;
	}else{
		cout << "testing the accuracy of shell pairs etc.      " << endl;
	}
	cout << "----------------------------------------------" << endl;
	string inf = "shellpair/input.in";
	Double thresh = 1.0E-14;
	if (doLargeJob) {

		// form the sig infor
		GlobalInfor globalInfor(inf);
		Molecule molecule(inf,2);
		MolShell s(inf,molecule);
		GIntsInfor ginfor(globalInfor,molecule);
		GIntJobInfor jobInfor(ginfor,COULOMB,0);
		DenMtrx den(globalInfor, molecule, s, s, molecule.getNSpin());
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			Mtrx& M = den.getMtrx(iSpin);
			M.set(0.1E0);
		}
		DenMtrxInfor denInfor(globalInfor,s,s,den);
		cout << "total number of basis sets: " << s.getNBas() << endl;
		tick_count t0 = tick_count::now();
		SigMolShellPairInfor sp(s, s, denInfor, jobInfor); 
		tick_count t1 = tick_count::now();
		Double t = (t1-t0).seconds();
		printf("%s  %-12.6f\n", "SigMolShellPairInfor forming time with seconds ", t);
		UInt nAtoms = s.getNAtomShells();
		UInt nTolAtomPairs = nAtoms*(nAtoms+1)/2;
		UInt nSigAtomPairs = sp.getNSigAtomShellPairs();
		cout << "total atom pairs: " << nTolAtomPairs << endl;
		cout << "total sig atom shell pairs: " << nSigAtomPairs << endl;
		cout << "difference is " << nTolAtomPairs - nSigAtomPairs << endl;
		cout << "percentage is " << nSigAtomPairs/nTolAtomPairs << endl;
		sp.print(1);
		sp.statPrint(thresh,s,s);

	}else {

		// form the sig shell pair information
		GlobalInfor globalInfor(inf);
		Molecule molecule(inf,1);
		MolShell s(inf,molecule);
		GIntsInfor ginfor(globalInfor,molecule);
		GIntJobInfor jobInfor(ginfor,EXCHANGE,0);
		DenMtrx den(globalInfor, molecule, s, s, molecule.getNSpin());
		for(UInt iSpin=0; iSpin<den.getNSpin(); iSpin++) {
			Mtrx& M = den.getMtrx(iSpin);
			M.set(0.1E0);
		}
		DenMtrxInfor denInfor(globalInfor,s,s,den);
		SigMolShellPairInfor sp(s, s, denInfor, jobInfor); 
		UInt nAtoms = s.getNAtomShells();
		UInt nTolAtomPairs = nAtoms*(nAtoms+1)/2;
		UInt nSigAtomPairs = sp.getNSigAtomShellPairs();
		cout << "total atom pairs: " << nTolAtomPairs << endl;
		cout << "total sig atom shell pairs: " << nSigAtomPairs << endl;
		cout << "difference is " << nTolAtomPairs - nSigAtomPairs << endl;
		molecule.print();
		s.print(3);
		sp.statPrint(thresh,s,s);
		sp.print(3);
	}
}

