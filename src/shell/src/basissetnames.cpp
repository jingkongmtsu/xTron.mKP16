/**
 * CPP files corresponding to the basissetnames.h
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include<fstream>
#include<cstdlib>
#include "textread.h"
#include "excep.h"
#include "basissetnames.h"
using namespace std;
using namespace textread;
using namespace excep;
using namespace basissetnames;

string basissetnames::getBasisSetFileNameFromInput(const string& inputBasisSetName)
{
	// get the basis file env
	char* path = getenv("EMUL_BASIS_FILE_DIR");
	if (path==NULL) {
		Excep excep("basissetnames","getBasisSetFileNameFromInput",EXCEPTION_DIR_MISSING, 
				"EMUL_BASIS_FILE_DIR is not defined");
		handleExcep(excep);
	}

	// set up the setting file path
	// where we can find the xcfunc definition
	string settingPath(path);
	string setting_file = settingPath + "/" + "basis.conf";
	string input_file   = setting_file;

	// now open the resource
	const char* source = input_file.c_str();
	ifstream inf;
	inf.open(source,ios::in);
	if (!inf) {
		Excep excep("basissetnames","getBasisSetFileNameFromInput",EXCEPTION_FILE_MISSING,input_file);
		handleExcep(excep);
	}

	// now let's go to search this file
	string line;
	WordConvert w;
	bool find = false;
	string fileName = "NOT_VALID_BASIS_SET_FILE_NAME";
	while(getline(inf,line)) {
		LineParse l(line);

		// for basis.conf, pure comment line or empty line, we just skip them
		if(l.isCom() || l.isEmp()) continue;

		// pre-check 
		UInt n = l.getNPieces();
		if(n != 2) {
			cout << line << endl;
			Excep excep("basissetnames","getBasisSetFileNameFromInput",EXCEPTION_INPUT_PARAMETER_INVALID,
					"This line of content is invalid to parse");
			handleExcep(excep);
		}

		// if this is the name we want to get, then we get it
		string name = l.findValue(0);
		if (w.compare(name,inputBasisSetName)) {
			fileName = l.findValue(1);
			find  = true;
			break;
		}	
	}

	// close the file
	inf.close();

	// let's see whether we get the basis set file name 
	if (! find) {
		Excep excep("basissetnames","getBasisSetFileNameFromInput",EXCEPTION_FAIL_LOCATE_BASIS_SET_NAME,
				inputBasisSetName);
		handleExcep(excep);
	}

	// finally return the file name
	return fileName;
}

bool basissetnames::isValidBasisSetName(const string& name)
{
	// get the basis file env
	char* path = getenv("EMUL_BASIS_FILE_DIR");
	if (path==NULL) {
		Excep excep("basissetnames","isValidBasisSetName",EXCEPTION_DIR_MISSING,"EMUL_BASIS_FILE_DIR is not defined");
		handleExcep(excep);
	}

	// set up the setting file path
	// where we can find the xcfunc definition
	string settingPath(path);
	string setting_file = settingPath + "/" + "basis.conf";
	string input_file   = setting_file;

	// now open the resource
	const char* source = input_file.c_str();
	ifstream inf;
	inf.open(source,ios::in);
	if (!inf) {
		Excep excep("basissetnames","isValidBasisSetName",EXCEPTION_FILE_MISSING,input_file);
		handleExcep(excep);
	}

	// search this file
	string line;
	bool find = false;
	WordConvert w;
	while(getline(inf,line)) {
		LineParse l(line);

		// for basis.conf, pure comment line or empty line, we just skip them
		if(l.isCom() || l.isEmp()) continue;

		// pre-check 
		UInt n = l.getNPieces();
		if(n != 2) {
			cout << line << endl;
			Excep excep("basissetnames","isValidBasisSetName",EXCEPTION_INPUT_PARAMETER_INVALID,
					"This line of content is invalid to parse");
			handleExcep(excep);
		}

		// if this is the name we want to get, then we get it
		string name0 = l.findValue(0);
		if (w.compare(name,name0)) {
			find  = true;
			break;
		}	
	}

	// close the file
	inf.close();

	// finally return the find status
	return find;
}

