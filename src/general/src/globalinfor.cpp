/**
 * cpp file for recording global information 
 * \author fenglai liu 
 */
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "excep.h"
#include "textread.h"
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include "parameterparsing.h"
#include "globalinfor.h"
using namespace std;
using namespace excep;
using namespace textread;
using namespace boost::filesystem;
using namespace parameterparsing;
using namespace globalinfor;

void GlobalInfor::formScratchDir() 
{
	// create an unique folder name for current job
	path uniquePath = unique_path();

	// get the scracth dir
	char* scratch = getenv("EMUL_SCRATCH_DIR");
	if (! scratch) {
		Excep excep("GlobalInfor","formScratchDir",EXCEPTION_DIR_MISSING,
				"EMUL_SCRATCH_DIR is not defined");
		handleExcep(excep);
	}

	// form the real scratch path
	path scratchPath(scratch);
	scratchPath /= uniquePath;

	// let's see whether the path exist
	// if so we issue an error
	scrDir = scratchPath.string();
	if (exists(scratchPath)) {
		printf("something wrong with the scratch dir %s\n", scrDir.c_str());
		Excep excep("GlobalInfor","formScratchDir",EXCEPTION_DIR_ALREADY_EXIST,
				"the scratch dir already existed, please remove the dir first");
		handleExcep(excep);
	}
}

void GlobalInfor::inforParse()
{
	// read in keywords
	UInt section = 0;
	bool isClusterSection = false;
	ParameterParsing pp(input,"global_infor",section,isClusterSection);
	UInt nKey = pp.getKeyNumber();
	WordConvert w;
	for(UInt iKey=0; iKey<nKey; iKey++) {

		// get the key name
		string key = pp.getKey(iKey);

		// whether it do multi-threads
		if (w.compare("MULTI_THREADS",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  multi_threads
			/// group    global_infor
			/// values   true or false (T or F also works)
			/// 
			/// this keyword defines whether you want to use multi-threads mode
			/// for the real calculation. If it's false, the program will be 
			/// running in serial mode. the Default value is true.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				enableMultiThreads = true;
			}else if (value == "FALSE" || value == "F") {
				enableMultiThreads = false;
				nCPUThreads=1;
			}else{
				string infor = "only true or false allowed to process multi_threads";
				Excep excep("GlobalInfor","inforParse",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("CPU_THREADS_NUMBER",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  cpu_threads_number
			/// group    global_infor
			/// values   set the number of threads you really want to use
			/// 
			/// the value you set can be less than the number of CPU core,
			/// it must be >= 1. If you set 1 then the program running 
			/// in serial mode. We use the following function 
			/// boost::thread::hardware_concurrency() to detect the possible
			/// number of threads, and if the threads number larger than
			/// this value it will be set to the value that the above function
			/// generates.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (w.toUInt(value,tmp)) {
				if (tmp>0) {
					nCPUThreads = tmp;
				}else{
					string infor = "the cpu_threads_number can not be 0";
					Excep excep("GlobalInfor","inforParse",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
					handleExcep(excep);
				}
			}else{
				string infor = "only integer allowed to process cpu_threads_number";
				Excep excep("GlobalInfor","inforParse",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else {
			cout << key << endl;
			Excep excep("Global_infor","inforParse",EXCEPTION_INPUT_PARAMETER_INVALID,
					"the keyword value in the input section can not be parsed");
			handleExcep(excep);
		}
	}

	// check nCPUThreads with multi-threads mode
	// since only multi-threads are enabled we can read in the information
	if (! enableMultiThreads && nCPUThreads>1) {
		nCPUThreads=1;

		// issue a warning message
		string infor = "reset nCPUThreads to be 1 since the multi-threads mode is turning off";
		Excep excep("GlobalInfor","inforParse",EXCEPTION_IMPROPER_THREAD_NUM_SETTING,infor);
		handleExcep(excep);
	}

	// we need to re-check the number of threads that user set
	// it should not exceed the nProcs 
	// here is the situation is multi-threads are not set, so in
	// default it's true
	UInt nProcs = boost::thread::hardware_concurrency();
	if (nCPUThreads == 0) {
		nCPUThreads=nProcs;
	} else if(nCPUThreads>nProcs) {
		nCPUThreads=nProcs;

		// issue a warning message
		string infor1 = "you can not set nCPUThreads more than the core number,";
		string infor2 = " we reset it using boost::thread::hardware_concurrency()";
		string infor  = infor1+infor2;
		Excep excep("GlobalInfor","inforParse",EXCEPTION_IMPROPER_THREAD_NUM_SETTING,infor);
		handleExcep(excep);
	}
}

GlobalInfor::GlobalInfor(const string& in):input(in),enableMultiThreads(true),nCPUThreads(0)        
{
	formScratchDir();
	inforParse();

	// finally, let's print out the whole input file
	cout << endl;
	cout << "*******************************************" << endl;
	cout << "*     The Original Input File Content     *" << endl;
	cout << "*******************************************" << endl;

	// open the input file
	std::ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "failed to open the input file";
		Excep excep("GlobalInfor","constructor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// now let's print out each line of input
	string line;
	while(getline(inf,line)) {
		cout << line << endl;
	}

	// now close the input
	inf.close();
}

void GlobalInfor::print() const 
{
	cout << endl;
	cout << "*******************************************" << endl;
	cout << "*            GlobalInfor Print            *" << endl;
	cout << "*******************************************" << endl;
	printf("input   file is        : %s\n", input.c_str());
	printf("scratch file is        : %s\n", scrDir.c_str());
	printf("number of CPU threads  : %d\n", (Int)nCPUThreads);
	if (enableMultiThreads) {
		printf("we use multi-threads \n");
	}else{
		printf("we use single-threads\n");
	}
}
