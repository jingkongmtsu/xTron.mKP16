//
// testing log:
//
// Feb. 2014: 
// all of files in general are well tested
// textread;
// parameterparsing;
// globalinfor;
// xcfunc;
// filerw
//
// May 2014:
// test localMemScr
//
// June 2014:
// test all of components after we revised the scalablevec.h
// now we only stick to the TBB allocator
//
// Dec. 2014:
// add in histdataman test, and passed
// however, we did not test the situation to recover record.txt 
// from contructor
//
//
#include<cstdio>
#include<string>
#include<iostream>
#include<fstream>
#include<cmath>
#include "libgen.h"
#include "textread.h"
#include "parameterparsing.h"
#include "globalinfor.h"
#include "xcfunc.h"
#include "filerw.h"
#include "matrix.h"
#include "localmemscr.h"
#include "histdataman.h"
using namespace std;
using namespace textread;
using namespace xcfunc;
using namespace parameterparsing;
using namespace globalinfor;
using namespace filerw;
using namespace matrix;
using namespace localmemscr;
using namespace histdataman;

Int main(int argc, char* argv[])
{
	//
	// set input and output
	//
	string input = "input.in";
	string outf  = "test.out";
	freopen (outf.c_str(),"w",stdout);

	//
	// set test parameter
	//
	bool testTextRead  = false;
	bool testPP        = false;
	bool testGInfor    = false;
	bool testxcfunc    = false;
	bool testRW        = false;
	bool testLocalMemScr = false;
	bool testHistDataMan = false;
	for(Int i=1; i<argc; i++) {
		string com = argv[i];
		if (com == "textread"  ) testTextRead    = true;
		if (com == "parameterparsing"  ) testPP  = true;
		if (com == "globalinfor"  ) testGInfor   = true;
		if (com == "xcfunc"  )      testxcfunc   = true;
		if (com == "filerw"  )      testRW       = true;
		if (com == "localmemscr"  ) testLocalMemScr = true;
		if (com == "histdataman"  ) testHistDataMan = true;
	}

	// global infor
	if (testGInfor) {
		GlobalInfor gInfor(input);
		gInfor.print();
	}

	// textread 
	if (testTextRead) {
		ifstream inf;
		const char * input_file = input.c_str();
		inf.open(input_file,ios::in);
		string line;
		while(getline(inf,line)) {
			LineParse l(line);
			cout << "line is " << line << endl;
			if (l.isEmp()) {
				cout << "empty line" << endl;
			} else if (l.isCom()) {
				cout << "comment line" << endl;
				for(UInt i=0; i<l.getNPieces(); i++) {
					cout << "  piece in the line: " << l.findValue(i) << endl;
				}
			} else if (l.isSec()) {
				cout << "section line " << l.getSecName() << endl;
			} else if (l.isShlSign()) {
				cout << " shell sign " << endl;
			}else{
				for(UInt i=0; i<l.getNPieces(); i++) {
					cout << "  piece in the line: " << l.findValue(i) << endl;
				}
			}
		}
		inf.close();
	}

	// parameter parsing
	if (testPP) {
		ParameterParsing pp(input,"xcfunc",1);
		pp.print();
	}

	// xcfunc
	if (testxcfunc) {
		XCFunc xcfunc1(input,1);
		xcfunc1.print();
	}

	// fileRW
	// this testing is only on linux!
	if (testRW) {
		UInt len = 1000000;
		DoubleVec testing(len);
		for(UInt i=0; i<len; i++) testing[i] = i+1;
		FileReadWrite fileRW("/tmp", "testing");
		fileRW.write(&testing[0],len);
		DoubleVec testing2(len);
		fileRW.read(&testing2[0],len);
		Double diff = ZERO;
		for(UInt i=0; i<len; i++) diff += testing2[i] - testing[i];
		if (fabs(diff)<THRESHOLD_MATH) {
			cout << "file read and write success" << endl;
		}else{
			cout << "file read and write failed" << endl;
		}
	}

	// localmemscr
	if(testLocalMemScr) {
		LocalMemScr scr(10000);
		Double* head = NULL;
		for(UInt i=0; i<130; i++) {
			Double* t = scr.getNewMemPos(i+1);
			if (i>0) {
				cout << "memory position " << t-head << endl;
			}else{
				head = t;
			}		
		}
	}

	// testing histdataman
	if(testHistDataMan) {

		// test with both file on disk or on the memory
		GlobalInfor infor(input);
		UInt sec   = 1;
		UInt nSpin = 2;
		for(UInt times=0; times<2; times++) {

			// set up hist object
			bool inFile = true;
			if (times == 1) inFile = false;
			HistDataMan hist(infor,"test",sec,nSpin,inFile);
			if (inFile) {
				cout << "testing histdataman with file saved on disk" << endl;
			}else{
				cout << "testing histdataman with file saved on memory" << endl;
			}

			// set up data matrix
			UInt n = 100;
			Mtrx a1(n,n);
			Mtrx b1(n,n);
			for(UInt i=0; i<n; i++) {
				for(UInt j=0; j<n; j++) {
					a1(j,i) = i+j+1;
					b1(j,i) = a1(j,i) + 1;
				}
			}
			n = 200;
			Mtrx a2(n,n);
			Mtrx b2(n,n);
			for(UInt i=0; i<n; i++) {
				for(UInt j=0; j<n; j++) {
					a2(j,i) = i+j+2;
					b2(j,i) = a2(j,i) + 1;
				}
			}
			n = 300;
			Mtrx a3(n,n);
			Mtrx b3(n,n);
			for(UInt i=0; i<n; i++) {
				for(UInt j=0; j<n; j++) {
					a3(j,i) = i+j+3;
					b3(j,i) = a3(j,i) + 1;
				}
			}

			// save/read data test
			// 0  matrix - symmetrical 
			// 1  matrix - asymmetrical 
			// 2  vector
			for(UInt i=0; i<3; i++) {

				// set the symm status
				bool inSymm = true;
				if (i>0) inSymm = false;
				if (i == 0) {
					cout << "histdataman test on matrix symmetrical save/read" << endl;
				}else if (i == 1) {
					cout << "histdataman test on matrix asymmetrical save/read" << endl;
				}else{
					cout << "histdataman test on vector save/read" << endl;
				}

				// store data
				if (i==2) {
					hist.storeData(a1.getVec(),0);
					hist.storeData(b1.getVec(),1);
					hist.storeData(a2.getVec(),0);
					hist.storeData(b2.getVec(),1);
					hist.storeData(a3.getVec(),0);
					hist.storeData(b3.getVec(),1);
				}else{
					hist.storeData(a1,inSymm,0);
					hist.storeData(b1,inSymm,1);
					hist.storeData(a2,inSymm,0);
					hist.storeData(b2,inSymm,1);
					hist.storeData(a3,inSymm,0);
					hist.storeData(b3,inSymm,1);
				}

				// now get data
				UInt n   = 100;
				UInt sec = 0+i*3;
				bool thresh = 1.0E-8;
				Mtrx t(n,n);
				if (i==2) {
					hist.retrieveData(t.getVec(),sec,0);
				}else{
					hist.retrieveData(t,sec,0);
				}
				if(t.diffComp(a1,thresh,inSymm)) {
					cout << "histdataman a1 matrix failed" << endl;
				}
				t.set(ZERO);
				if (i==2) {
					hist.retrieveData(t.getVec(),sec,1);
				}else{
					hist.retrieveData(t,sec,1);
				}
				if(t.diffComp(b1,thresh,inSymm)) {
					cout << "histdataman b1 matrix failed" << endl;
				}

				n = 200;
				sec = 1+i*3;
				Mtrx t1(n,n);
				if (i==2) {
					hist.retrieveData(t1.getVec(),sec,0);
				}else{
					hist.retrieveData(t1,sec,0);
				}
				if(t1.diffComp(a2,thresh,inSymm)) {
					cout << "histdataman a2 matrix failed" << endl;
				}
				t1.set(ZERO);
				if (i==2) {
					hist.retrieveData(t1.getVec(),sec,1);
				}else{
					hist.retrieveData(t1,sec,1);
				}
				if(t1.diffComp(b2,thresh,inSymm)) {
					cout << "histdataman b2 matrix failed" << endl;
				}

				n = 300;
				sec = 2+i*3;
				Mtrx t2(n,n);
				if (i==2) {
					hist.retrieveData(t2.getVec(),sec,0);
				}else{
					hist.retrieveData(t2,sec,0);
				}
				if(t2.diffComp(a3,thresh,inSymm)) {
					cout << "histdataman a3 matrix failed" << endl;
				}
				t2.set(ZERO);
				if (i==2) {
					hist.retrieveData(t2.getVec(),sec,1);
				}else{
					hist.retrieveData(t2,sec,1);
				}
				if(t2.diffComp(b3,thresh,inSymm)) {
					cout << "histdataman b3 matrix failed" << endl;
				}
			}
		}
	}

	return 0;

}
