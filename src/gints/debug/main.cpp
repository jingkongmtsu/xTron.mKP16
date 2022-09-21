#include<string> 
using namespace std;

///
/// shell pair test function
///
extern void shellpair_test(bool largeJob);

///
/// gints 2D test
///
extern void gints2d_test(bool doLargeJob);

///
/// integral infor test
///
extern void integralinfor_test();

///
/// gints 4D test
///
extern void gints4d_test(bool doLargeJob);

///
/// gints 4D derivatives test
///
extern void gints4deriv_test(bool doLargeJob);

int main(int argc, char* argv[])
{

	// job definition
	bool testingShellPair       = false;
	bool testingLargeShellPair  = false;
	bool testingGInts2D         = false;
	bool testingGInts4D         = false;
	bool testingGInts4Deriv     = false;
	bool testingGInts4DPara     = false;
	bool testingIntegralInfor   = false;
	bool testingGInts4Deriv_large = false;

	// parse the input parameter
	for(int i=1; i<argc; i++) {
		string com = argv[i];
		if (com == "shellpair")        testingShellPair          = true;
		if (com == "largeshellpair")   testingLargeShellPair     = true;
		if (com == "gints2d")          testingGInts2D            = true;
		if (com == "gints4d")          testingGInts4D            = true;
		if (com == "gints4deriv")      testingGInts4Deriv        = true;
		if (com == "gints4d_para")     testingGInts4DPara        = true;
		if (com == "integral_infor")   testingIntegralInfor      = true;
		if (com == "gints4deriv_large")  testingGInts4Deriv_large = true;
	}

	/////////////////////////////////
	// test the shell pair class   //
	/////////////////////////////////
	if (testingShellPair || testingLargeShellPair) {
		bool doLargeJob = false;
		if (testingLargeShellPair) doLargeJob = true;
		shellpair_test(doLargeJob);
	}

	/////////////////////////////////
	// test the gints2d class      //
	/////////////////////////////////
	if (testingGInts2D) {
		bool doLargeJob = false;
		gints2d_test(doLargeJob);
	}

	/////////////////////////////////
	// test the integral infor     //
	// reading                     //
	/////////////////////////////////
	if (testingIntegralInfor) {
		integralinfor_test();
	}

	/////////////////////////////////
	// test the gints4d class      //
	/////////////////////////////////
	if (testingGInts4D) {
		bool doLargeJob = false;
		gints4d_test(doLargeJob);
	}else if (testingGInts4DPara) {
		bool doLargeJob = true;
		gints4d_test(doLargeJob);
	}

	/////////////////////////////////
	// test the gints4deriv class  //
	/////////////////////////////////
	if (testingGInts4Deriv || testingGInts4Deriv_large) {
		bool doLargeJob = false;
		if (testingGInts4Deriv_large) doLargeJob = true;
		gints4deriv_test(doLargeJob);
	}

	return 0;
}
