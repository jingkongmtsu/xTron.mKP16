#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string> 
#include<vector>
#include<cmath>
using namespace std;

//
// testing the exchange energy density forming
//
extern void exrho_test(bool doLargeJob, int job, bool withMIC);

//
// testing the rsp integrals for MIC
//
extern void espints_test();

//
// testing the forming of xc mtrx
//
extern void xcmtrx_test(bool testTiming, bool withMIC);

int main(int argc, char* argv[])
{
	// job definition
	bool testingXCMtrx      = false;
	bool withMIC            = false;
	bool testingXCMtrxPara  = false;
	bool testingEXRho       = false;
	bool testingEXRhoPara   = false;
	bool testingEXRhoTiming = false;
	bool testingESPInts     = false;

	// parse the input parameter
	for(int i=1; i<argc; i++) {
		string com = argv[i];
		if (com == "xcmtrx")        testingXCMtrx      = true;
		if (com == "xcmtrx_para")   testingXCMtrxPara  = true;
		if (com == "with_mic")      withMIC            = true;
		if (com == "exrho")         testingEXRho       = true;
		if (com == "espints")       testingESPInts     = true;
		if (com == "exrho_para")    testingEXRhoPara   = true;
		if (com == "exrho_timing")  testingEXRhoTiming = true;
	}

	if (testingXCMtrx || testingXCMtrxPara) {
		bool testTiming = false;
		if (testingXCMtrxPara) testTiming = true;
		xcmtrx_test(testTiming,withMIC);
	}

	if (testingEXRho || testingEXRhoPara || testingEXRhoTiming) {
		int job = 0;
		if (testingEXRhoTiming) job = 1;
		bool doLargeJob = false;
		if (testingEXRhoPara) doLargeJob = true;
		exrho_test(doLargeJob,job,withMIC);
	}

	if (testingESPInts) {
		espints_test();
	}

	return 0;
}
