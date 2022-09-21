#include<cstdio>
#include<cstdlib>
#include "molecule.h"
#include "shell.h"
#include "shellsize.h"
using namespace molecule;
using namespace shell;
using namespace shellsize;

//
// testing C2P/P2C and basis set normalization scale
//
extern void cptrans_test(bool withPara);
extern void cptrans_para_test();

Int main(int argc, char* argv[])
{
	// setting the testing jobs
	bool testShell  = false;
	bool testCPTrans= false;
	bool testCPTrans_para= false;
	bool testShellSize   = false;

	// parse the input parameter
	for(Int i=1; i<argc; i++) {
		string com = argv[i];
		if (com == "shell" ) testShell = true;
		if (com == "c2p" )   testCPTrans = true;
		if (com == "c2p_para") testCPTrans_para = true;
		if (com == "shellsize") testShellSize = true;
	}

	// test shell 
	if (testShell) {
		string inf = "shell_test/input.in";
		for(UInt i=1; i<=3; i++) {
			Molecule molecule(inf,i);
			MolShell ms(inf,molecule,0,true);
			ms.print(3);
		}
	}

	// test shell size
	if (testShellSize) {
		string inf = "shell_test/input.in";
		Molecule molecule(inf,1);
		MolShell ms(inf,molecule,0);
		MolShellSize size(ms,1.0E-8);
		size.print();
	}

	// test cptrans etc.
	if (testCPTrans) {
		bool withPara = true;
		cptrans_test(withPara);
	}

	// test cptrans with parallel 
	if (testCPTrans_para) {
		cptrans_para_test();
	}

	return 0;
}
