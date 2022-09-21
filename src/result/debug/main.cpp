#include<string> 
using namespace std;

///
/// testing core guess 
///
extern void coreguess_test();

///
/// testing free atom density matrix etc.
///
extern void atomDenMtrx_test();

int main(int argc, char* argv[])
{

	// job definition
	bool testingCoreGuess  = false;
	bool testingAtomDenMtrx  = false;

	// parse the input parameter
	for(int i=1; i<argc; i++) {
		string com = argv[i];
		if (com == "coreguess") testingCoreGuess = true;
		if (com == "atomdenmtrx") testingAtomDenMtrx = true;
	}

	/////////////////////////////////
	//    test the core guess      //
	/////////////////////////////////
	if (testingCoreGuess) {
		coreguess_test();
	}

	//////////////////////////////////
	//    test the atom density etc.//
	//////////////////////////////////
	if (testingAtomDenMtrx) {
		atomDenMtrx_test();
	}

	return 0;
}
