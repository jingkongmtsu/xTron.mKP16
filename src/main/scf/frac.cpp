//
// This is the template to set up fractional spin calculation
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
#include "excep.h"
#include "textread.h"
#include "molecule.h"
#include "denmtrx.h"
#include "fracspininfor.h"
#include "shell.h"
#include "scf.h"
#include "parameterparsing.h"
using namespace molecule;
using namespace excep;
using namespace textread;
using namespace globalinfor;
using namespace shell;
using namespace denmtrx;
using namespace fracspininfor; 
using namespace scf;
using namespace std;
using namespace parameterparsing;

void fracSpin(const string& inf, const string& outf)
{
	// set the inf and outf
	freopen (outf.c_str(),"w",stdout);

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
			Excep excep("frac","frac",EXCEPTION_FILE_MISSING,infor);
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
		if (nMol != 4 && nMol != 3) {
			string infor = "we must have three/four sections set up for fractional spin calculation";
			Excep excep("frac","frac",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
			handleExcep(excep);
		}
	}

	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "          INITIAL GUESS CALCULATION           " << endl;
	cout << "----------------------------------------------" << endl;

	// this is real job
	GlobalInfor globalInfor(inf);
	Molecule molecule0(inf,1);
	MolShell s0(inf,molecule0);
	cout << "number of basis sets " << s0.getNBas() << endl;

	// set up scf
	SCF scf0(globalInfor,molecule0,s0);

	// now form the input density matrix
	UInt nSpin = molecule0.getNSpin();
	DenMtrx den0(globalInfor,molecule0,s0,s0,nSpin); 
	den0.formSCFGuess(s0,molecule0,scf0);

	// do scf
	if (true) {
		boost::timer::auto_cpu_timer t;
		scf0.doSCF(globalInfor,molecule0,s0,den0);
	}

	// heading for second SCF run
	cout << "----------------------------------------------" << endl;
	cout << "             SECOND SCF CALCULATION           " << endl;
	cout << "----------------------------------------------" << endl;

	// now prepare the real calculation
	Molecule molecule(inf,2);
	MolShell s(inf,molecule);

	// now set up density matrix
	DenMtrx den(globalInfor,molecule,s,s,nSpin); 
	den.formDenMtrx(scf0.getMO());

	// now it's SCF
	SCF scf(globalInfor,molecule,s);
#ifdef WITH_GDM
	scf.updateMO(scf0.getMO());
#endif
	if (true) {
		boost::timer::auto_cpu_timer t;
		scf.doSCF(globalInfor,molecule,s,den);
	}

	// now let's see whether we have three sections or four sections
	if (nMol == 3) {

		// heading for third SCF run
		cout << "----------------------------------------------" << endl;
		cout << "             THIRD SCF CALCULATION            " << endl;
		cout << "----------------------------------------------" << endl;

		// now prepare the real calculation
		Molecule molecule3(inf,3);
		MolShell s3(inf,molecule3);

		// let's see whether we have any fractional spin job?
		// if no, then the newMol is same with molecule2
		// else we need to use the modified molecule
		FracSpinInfor fracSpinInfor(globalInfor,molecule3);
		const Molecule& newMol = fracSpinInfor.getMol();
		newMol.print();

		// now form the density matrix
		// we use the new molecule information so that to get
		// the scalied electron occupation number
		nSpin = newMol.getNSpin();
		DenMtrx den3(globalInfor,newMol,s3,s3,nSpin); 

		// form density matrix
		// we copy the alpha mo from previous geometry
		MO mo3(globalInfor,s3,molecule3,nSpin);
		const MO& oldMOData = scf.getMO();
		for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
			const Mtrx& oldMO = oldMOData.getMtrx(iSpin);
			Mtrx& newMO = mo3.getMtrx(iSpin);
			newMO = oldMO;
		}
		if (fracSpinInfor.doFracSpinJob()) {
			mo3.formFracSpinMO(fracSpinInfor);
			den3.formDenMtrx(mo3);
		}else{
			string infor = "for the third SCF run you must define fractional spin job information";
			Excep excep("frac","frac",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
			handleExcep(excep);
		}

		// do scf
		SCF scf3(globalInfor,newMol,s3);
		if (true) {
			boost::timer::auto_cpu_timer t;
			scf3.doSCF(globalInfor,newMol,s3,den3);
		}
	}else{

		// heading for second SCF run
		cout << "----------------------------------------------" << endl;
		cout << "              THRID SCF CALCULATION           " << endl;
		cout << "----------------------------------------------" << endl;

		// now prepare the real calculation
		Molecule molecule2(inf,3);
		MolShell s2(inf,molecule2);

		// set up the mo
		// the beta part we copy the alpha part
		// so that to pretend to be ROHF calculation
		nSpin = molecule2.getNSpin();
		MO mo2(globalInfor,s2,molecule2,nSpin);
		const MO& oldMO = scf.getMO();
		const Mtrx& alphaOldMO = oldMO.getMtrx(0);
		Mtrx& alphaNewMO = mo2.getMtrx(0);
		alphaNewMO.copyMatrix(alphaOldMO);
		if (oldMO.getNSpin() == 2) {
			cout << "copy beta mo to alpha mo" << endl;
			Mtrx& betaNewMO = mo2.getMtrx(1);
			betaNewMO.copyMatrix(alphaOldMO);
		}
		mo2.print(3);

		// now set up density matrix
		DenMtrx den2(globalInfor,molecule2,s2,s2,nSpin); 
		den2.formDenMtrx(mo2);

		// now it's SCF
		SCF scf2(globalInfor,molecule2,s2);
		if (true) {
			boost::timer::auto_cpu_timer t;
			scf2.doSCF(globalInfor,molecule2,s2,den2);
		}

		// heading for third SCF run
		cout << "----------------------------------------------" << endl;
		cout << "             FOURTH SCF CALCULATION           " << endl;
		cout << "----------------------------------------------" << endl;

		// now prepare the real calculation
		Molecule molecule3(inf,4);
		MolShell s3(inf,molecule3);

		// let's see whether we have any fractional spin job?
		// if no, then the newMol is same with molecule2
		// else we need to use the modified molecule
		FracSpinInfor fracSpinInfor(globalInfor,molecule3);
		const Molecule& newMol = fracSpinInfor.getMol();
		newMol.print();

		// now form the density matrix
		// we use the new molecule information so that to get
		// the scalied electron occupation number
		nSpin = newMol.getNSpin();
		DenMtrx den3(globalInfor,newMol,s3,s3,nSpin); 

		// form density matrix
		// we copy the alpha mo from previous geometry
		// also we note that the beta mo is same with alpha in this case
		MO mo3(globalInfor,s3,molecule3,nSpin);
		Mtrx& newAlphaMO = mo3.getMtrx(0);
		newAlphaMO.copyMatrix(alphaNewMO);
		Mtrx& newBetaMO  = mo3.getMtrx(1);
		newBetaMO.copyMatrix(alphaNewMO);
		if (fracSpinInfor.doFracSpinJob()) {
			mo3.formFracSpinMO(fracSpinInfor);
			mo3.print(3);
			den3.formDenMtrx(mo3);
		}else{
			string infor = "for the third SCF run you must define fractional spin job information";
			Excep excep("frac","frac",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
			handleExcep(excep);
		}

		// do scf
		SCF scf3(globalInfor,newMol,s3);
		if (true) {
			boost::timer::auto_cpu_timer t;
			scf3.doSCF(globalInfor,newMol,s3,den3);
		}
	}
}
