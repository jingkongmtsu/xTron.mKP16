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
#include "libgen.h"
#include "excep.h"
#include "textread.h"
using namespace excep;
using namespace textread;
using namespace std;

// external functions
extern void doSP(const string& inf, const string& outf);
extern void fracSpin(const string& inf, const string& outf);

int main(int argc, char* argv[])
{
	// set the inf and outf
	string inf   = argv[1];
	string outf  = argv[2];

	// let's see whether this is a fractional spin case
	bool doJob = true;
	bool isFracSpinJob = false;
	if (doJob) {

		// open the file
		ifstream in;
		const char * input_file = inf.c_str();
		in.open(input_file,ios::in);
		if (!in) {
			string infor = "failed to open the input file";
			Excep excep("main","main",EXCEPTION_FILE_MISSING,infor);
			handleExcep(excep);
		}

		// counting that how many sections we have
		string key = "frac_spin";
		string line;
		WordConvert w;
		in.seekg(0,ios::beg);
		while(getline(in,line)) {
			LineParse l(line);
			if (l.isSec() && w.compare(l.getSecName(), key)) {
				isFracSpinJob = true;
			}
		}
		in.close();
	}

	// whether this is fractional spin calculation?
	if (isFracSpinJob) {
		fracSpin(inf,outf);
	}else{
		doSP(inf,outf);
	}

	// everything good, return
	return 0;
}
