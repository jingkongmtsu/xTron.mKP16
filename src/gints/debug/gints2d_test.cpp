//
// note:
// set up gints2d test for KI, OV and NAI (Nov. 2014)
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include<cmath>
#include "tbb/tbb.h"
#include "libgen.h"
#include "filerw.h"
#include "globalinfor.h"
#include "gintsinfor.h"
#include "molecule.h"
#include "matrix.h"
#include "shell.h"
#include "gints2d.h"
using namespace molecule;
using namespace filerw;
using namespace globalinfor;
using namespace gintsinfor;
using namespace shell;
using namespace matrix;
using namespace gints2d;
using namespace tbb;
using namespace std;

void gints2d_test(bool doLargeJob)
{
	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "             GINTS 2D TESTING                 " << endl;
	if (doLargeJob) {
		cout << "       testing the timing performace          " << endl;
	}else{
		cout << "   testing the accuracy of 2D integrals etc.  " << endl;
	}
	cout << "----------------------------------------------" << endl;
	string inf = "gints2d/input.in";
	if (doLargeJob) inf = "gints2d/dna.in";

	// now let's organize the job list
	bool doEnergy = true;

	// testing the integral matrix - as a part of Fock matrix
	if(doEnergy) {

		// now do each job
		UInt nMol = 2;
		if (doLargeJob) nMol = 3;
		for(UInt iMol=1; iMol<=nMol; iMol++) {

			// form the molecule information etc.
			GlobalInfor globalInfor(inf);
			Molecule molecule(inf,iMol);
			MolShell s(inf,molecule);
			GIntsInfor infor(globalInfor,molecule);

			// form the job list
			UInt order = 0;
			UIntVec jobList;
			jobList.push_back(KINETIC);
			jobList.push_back(NUCLEAR_ATTRACTION);
			jobList.push_back(TWO_BODY_OVERLAP);

			// set up result matrix
			UInt nBas = s.getNBas();
			Mtrx M(nBas,nBas);
			cout << "geometry: "<< iMol << " with total basis sets " << nBas << endl;
			cout << "Cartesian number of basis sets: " << s.getNCarBas() << endl;

			// now loop over the jobss
			for(UInt iJob=0; iJob<jobList.size(); iJob++) {

				// print out job
				if (jobList[iJob]==KINETIC)            cout << "job is: kinetic         " << endl;
				if (jobList[iJob]==NUCLEAR_ATTRACTION) cout << "job is: nai             " << endl;
				if (jobList[iJob]==TWO_BODY_OVERLAP)   cout << "job is: two body overlap" << endl;
				if (jobList[iJob]==MOM_P) cout << "job is: momentum integral for dipole calculation" << endl;

				// set the multi-Matrix
				UInt nM = 0;
				if (jobList[iJob] == MOM_P) nM = 3;
				MtrxVec Mlist(nM);
				if (nM>0) {
					for(UInt i=0; i<nM; i++) {
						Mlist[i].init(nBas,nBas);
					}
				}

				// working code
				M.set(ZERO);
				GInts2D gints(s,s,infor,jobList[iJob],order);
				if (jobList[iJob] == MOM_P) {
					gints.doMultiMtrx(s,s,Mlist,true);
				}else{
					gints.doMtrx(s,s,molecule,M,true);
				}
				cout << endl;

				// now do debug work
				if (! doLargeJob && jobList[iJob] != MOM_P) {

					// get the standard data file
					string name = "ki";
					if (jobList[iJob]==NUCLEAR_ATTRACTION) name = "nai";
					if (jobList[iJob]==TWO_BODY_OVERLAP)   name = "ov";
					if (iMol == 2) {
						name = name + "2";
					}
					name = name + ".txt";
					string fileName = "gints2d/" + name;
					Mtrx N(nBas,nBas);
					FileReadWrite fileRW(fileName);
					fileRW.readMatrixFromTextFile(N.getPtr(),nBas,nBas);

					// now compare
					Double thresh = 1.0E-7;
					M.copyLowToUpper();
					bool failed = M.diffComp(N,thresh,false);
					if (failed) {
						cout << "the calculation of " << name << " is failed with accuracy of " << thresh << endl;
					}else{
						cout << "the calculation of " << name << " is passed with accuracy of " << thresh << endl;
					}
				}
			}
		}
	}
}
