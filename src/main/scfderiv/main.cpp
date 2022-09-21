//
// This is the template to do SCF deriv calculation
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
#include "mo.h"
#include "scfparam.h"
#include "globalinfor.h"
#include "excep.h"
#include "textread.h"
#include "molecule.h"
#include "denmtrx.h"
#include "shell.h"
#include "scfderiv.h"
using namespace mo;
using namespace molecule;
using namespace scfparam;
using namespace excep;
using namespace textread;
using namespace globalinfor;
using namespace shell;
using namespace denmtrx;
using namespace scfderiv;
using namespace std;

int main(int argc, char* argv[])
{
	// set the inf and outf
	string inf   = argv[1];
	string outf  = argv[2];
	freopen (outf.c_str(),"w",stdout);

	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "           SCF DERIVATIVES CALCULATION        " << endl;
	cout << "----------------------------------------------" << endl;

	// this is real job
	GlobalInfor globalInfor(inf);
	globalInfor.print();
	Molecule molecule(inf,1);
	MolShell s(inf,molecule);
	cout << "number of basis sets " << s.getNBas() << endl;

	// form the scf param
	SCFParam par(globalInfor,molecule);

	// now reading the mo
	UInt nSpin = molecule.getNSpin();
	MO mo(globalInfor,s,molecule,nSpin);
	mo.recover(par.getGuessMOPath());

	// now let's form density matrix
	DenMtrx den(globalInfor,molecule,s,s,nSpin); 
	den.formDenMtrx(mo); 

	// do scf derivatives
	UInt order = 1;
	if (true) {
		boost::timer::auto_cpu_timer t;
		SCFDeriv deriv(molecule,order);
		deriv.doDeriv(par,molecule,s,den);
	}
}
