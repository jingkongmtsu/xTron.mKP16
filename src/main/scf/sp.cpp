//
// This is the template to set up a single point SCF calculation
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
#include "libgen.h"
#include "globalinfor.h"
#include "scfparam.h"
#include "excep.h"
#include "textread.h"
#include "molecule.h"
#include "denmtrx.h"
#include "shell.h"
#include "scf.h"
using namespace molecule;
using namespace excep;
using namespace textread;
using namespace globalinfor;
using namespace scfparam;
using namespace shell;
using namespace denmtrx;
using namespace scf;
using namespace std;

void doSP(const string& inf, const string& outf)
{
	// set the inf and outf
	freopen (outf.c_str(),"w",stdout);

	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "          SINGLE POINT SCF CALCULATION        " << endl;
	cout << "----------------------------------------------" << endl;

	// let's see how many molecules in the input file
	// right now we only support two cases:
	// 1  only one molecule;
	// 2  two molecules, one is for generating guess 
	//    and the second is for real calculation
	UInt nMol = 0;
	bool doJob = true;
	if (doJob) {

		// open the file
		ifstream in;
		const char * input_file = inf.c_str();
		in.open(input_file,ios::in);
		if (!in) {
			string infor = "failed to open the input file";
			Excep excep("scf_sp","scf_sp",EXCEPTION_FILE_MISSING,infor);
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

		// if nMols is not correct, report an error
		if (nMol != 1 && nMol != 2) {
			string infor = "we only support nMol=1 or nMol=2 cases";
			Excep excep("scf_sp","scf_sp",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
			handleExcep(excep);
		}
	}

	// this is real job
	GlobalInfor globalInfor(inf);
	globalInfor.print();
	Molecule molecule(inf,1);
	MolShell s(inf,molecule);
	cout << "number of basis sets " << s.getNBas() << endl;

	// set up scf
	SCF scf(globalInfor,molecule,s);

	// now form the input density matrix
	UInt nSpin = molecule.getNSpin();
	DenMtrx den(globalInfor,molecule,s,s,nSpin); 
	den.formSCFGuess(s,molecule,scf);

	// do we save the initial density matrix?
	SCFParam par(globalInfor,molecule);
	if (par.saveInitDenMtrxOnDisk()) {
		const string& p = par.getInitDenMtrxPath(); 
		den.writeToDisk(p);
	}

	// do scf
	if (true) {
		boost::timer::auto_cpu_timer t;
		scf.doSCF(globalInfor,molecule,s,den);
	}

	if (nMol == 2) {
		// heading for B05
		cout << "-----------Second SCF begins-----------" << endl;

		// now prepare the real calculation
		Molecule molecule2(inf,2);
		MolShell s2(inf,molecule2);

		// now form the density matrix
		nSpin = molecule2.getNSpin();
		DenMtrx den2(globalInfor,molecule2,s2,s2,nSpin); 
		den2.formDenMtrx(scf.getMO());

		// do scf
		SCF scf2(globalInfor,molecule2,s2);
#ifdef WITH_GDM
		scf2.updateMO(scf.getMO());
#endif
		if (true) {
			boost::timer::auto_cpu_timer t;
			scf2.doSCF(globalInfor,molecule2,s2,den2);
		}
	}
}
