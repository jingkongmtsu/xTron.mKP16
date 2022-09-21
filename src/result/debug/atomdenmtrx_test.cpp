//
// note:
// set up test for atom density matrices
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string> 
#include<vector>
#include<cmath>
#include<boost/lexical_cast.hpp>
#include "libgen.h"
#include "globalinfor.h"
#include "excep.h"
#include "textread.h"
#include "molecule.h"
#include "shell.h"
#include "atomdenmtrx.h"
using namespace molecule;
using namespace excep;
using namespace textread;
using namespace globalinfor;
using namespace shell;
using namespace atomdenmtrx;
using namespace std;

void atomDenMtrx_test()
{
	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "         ATOM DENSITY MATRICES TESTING        " << endl;
	cout << "----------------------------------------------" << endl;
	string inf = "atomdenmtrx/input.in";
	GlobalInfor globalInfor(inf);
	globalInfor.print();

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
			Excep excep("atomdenmtrx_test","atomdenmtrx_test",EXCEPTION_FILE_MISSING,infor);
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

	// now let's do the atom density matrix test
	for (UInt iMol=1; iMol<=nMol; iMol++) {

		// read in standard energy
		DoubleVec standard;
		vector<string> atomNames;
		cout << endl;
		cout << "*****************************************************" << endl;
		cout << "*  Molecule: " << boost::lexical_cast<string>(iMol)    << endl;
		cout << "*****************************************************" << endl;

		// open the file
		ifstream in;
		const char * input_file = inf.c_str();
		in.open(input_file,ios::in);
		if (!in) {
			string infor = "failed to open the input file";
			Excep excep("atomdenmtrx_test","atomdenmtrx_test",EXCEPTION_FILE_MISSING,infor);
			handleExcep(excep);
		}

		// counting that how many sections we have
		WordConvert w;
		string line;
		string line2;
		UInt count = 0;
		bool getit = false;
		while(getline(in,line)) {
			LineParse l(line);
			if (l.isCom() && line.find("atom energy result") != string::npos) {
				count = count + 1;
				if (count == iMol) {
					getit = true;
					while(getline(in,line2)) {
						LineParse l2(line2);
						if (l2.isCom() && line2.find("end") != string::npos) break;
						atomNames.push_back(l2.findValue(0));
						string s= l2.findValue(1);
						Double x = ZERO;
						if (w.toDouble(s,x)) {
							standard.push_back(x);
						}else{
							string infor = "failed to get the standard energy data";
							Excep excep("atomdenmtrx_test","atomdenmtrx_test",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
							handleExcep(excep);
						}
					}
					break;
				}
			}
		}
		in.close();

		// real calculation
		Molecule molecule(inf,iMol);
		MolShell s(inf,molecule);
		AtomDenMtrx den(globalInfor,molecule);
		den.formAtomDenMtrxInSCF(molecule,s); 
		for(UInt i=0; i<atomNames.size(); i++) {
			printf("standard energy for atom  %s is: %-14.7f\n", atomNames[i].c_str(), standard[i]);
		}
		den.print();

		// write it into files
		den.writeToDisk();

		// now let's try to recover it
		AtomDenMtrx den2(globalInfor,molecule.getSec());
		den2.recover(s);
		cout << "the atom density matrix after recovered " << endl;
		den2.print();
	}

}
