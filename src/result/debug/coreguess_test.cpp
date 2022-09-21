//
// note:
// set up core guess test 
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include<cmath>
#include "libgen.h"
#include "filerw.h"
#include "globalinfor.h"
#include "molecule.h"
#include "matrix.h"
#include "denmtrx.h"
#include "shell.h"
using namespace molecule;
using namespace filerw;
using namespace globalinfor;
using namespace shell;
using namespace matrix;
using namespace denmtrx;
using namespace std;

void coreguess_test()
{
	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "             CORE GUESS TESTING               " << endl;
	cout << "----------------------------------------------" << endl;
	string inf = "core/input.in";

	UInt nMol = 2;
	for (UInt iMol=0; iMol<nMol; iMol++) {

		// form the molecule information etc.
		GlobalInfor globalInfor(inf);
		Molecule molecule(inf,iMol+1);
		MolShell s(inf,molecule);

		// now form the density matrix
		UInt nSpin = molecule.getNSpin();
		UInt nBas  = s.getNBas();
		DenMtrx den(globalInfor,molecule,s,s,nSpin); 
		den.coreGuess(s,molecule);

		// now let's compare
		for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

			// main name
			string name = "denmtrx1";
			if (iMol == 1) name = "denmtrx2";

			// do we have other suffix?
			if (nSpin == 2) {
				if (iSpin == 0) {
					name = name + "_a";
				}else{
					name = name + "_b";
				}
			}
			name = name + ".txt";
			string fileName = "core/" + name;
			Mtrx N(nBas,nBas);
			FileReadWrite fileRW(fileName);
			fileRW.readMatrixFromTextFile(N.getPtr(),nBas,nBas);

			// now compare
			Double thresh = 1.0E-7;
			const Mtrx& M = den.getMtrx(iSpin);
			bool failed = M.diffComp(N,thresh,false);
			if (failed) {
				cout << "the calculation of " << name << " is failed with accuracy of " << thresh << endl;
			}else{
				cout << "the calculation of " << name << " is passed with accuracy of " << thresh << endl;
			}
		}
	}
}
