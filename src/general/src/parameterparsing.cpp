/**
 * \file   parameterparsing.cpp
 * \author Fenglai Liu 
 */
#include "excep.h"
#include "textread.h"
#include "parameterparsing.h"
#include<iostream>
#include<fstream>
#include <boost/lexical_cast.hpp>
using namespace excep;
using namespace textread;
using namespace parameterparsing;

ParameterParsing::ParameterParsing(const string& input, 
		const string& className, const UInt& sec, bool isCluster)
{

	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		Excep excep("ParameterParsing","constructor",EXCEPTION_FILE_MISSING,input);
		handleExcep(excep);
	}

	// get the section key word - so that we know which section
	// we are after
	string secKey = "molecule";
	if (sec == 0) {
		secKey = "cluster";
		if (! isCluster) secKey = "global_infor";
	}

	// find that section
	UInt n = 0;
	bool find = false;
	string line;
	WordConvert w;
	inf.seekg(0,ios::beg);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isSec() && w.compare(l.getSecName(), secKey)) {
			if (sec == 0) {
				find = true;
				break;
			}else{
				n++;
				if (sec == n) {
					find = true;
					break;
				}
			}
		}
	}

	// checking the finding status
	// for global section, user may want to use it's default value
	// therefore if it's not found this is OK
	if (isCluster || sec > 0) {
		if (! find) {
			string infor = "section number is ";
			infor = infor + boost::lexical_cast<string>(sec);
			if (secKey == "global_infor") infor = "section is global_infor";
			Excep excep("ParameterParsing","constructor",EXCEPTION_SECTION_SEARCH_FAILED,infor);
			handleExcep(excep);
		}
	}

	// now from here we begin to search the class definition with 
	// the given geometry section. This is only applied to non-global
	// sections. Make sure the section we get is only for this geometry
	// if we did not find something, that means user did not
	// provide any custom data
	if (sec >= 0 && isCluster) {
		find  = false;
		while(getline(inf,line)) {
			LineParse l(line);
			if (l.isSec() && w.compare(l.getSecName(), className)) {
				find = true;
				break;
			}else if (l.isSec() && w.compare(l.getSecName(), "molecule")) {
				break;
			}else if (l.isSec() && w.compare(l.getSecName(), "cluster")) {
				break;
			}
		}
		if (! find) return;
	}

	// now we enter into the user definition for this class
	// here we have some special cases
	// where the value corresponding to the keyword must retain
	// it's case. 
	// we name it one by one here
	while(getline(inf,line)) {
		LineParse l(line);

		// where we stop reading?
		if (l.isSec()) break;
		if (l.isCom() || l.isEmp()) continue;
		if (l.getNPieces() < 2) {
			cout << line << endl;
			Excep excep("ParameterParsing","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
					"The input parameter keyword should have its value defined");
			handleExcep(excep);
		}

		// debug choice
		// always has the keyword of "DEBUG"
		if (w.compare(l.findValue(0), "debug")) {

			// in this case, there could be two situations:
			// 1  debug is only for the given class
			//    the value attached with debug is the debug level indicator
			// 2  debug is only for one class
			//    the value attached with debug is the class name with debug level 0
			// 
			// or only one class need to be debugged
			if (l.getNPieces() == 2) {
				string value   = l.findValue(1);
				if (w.isInt(value)) {
					string keyword = className;
					w.capitalize(keyword);
					parameters.insert(pair<string,string>(keyword,value));
				}else{
					string keyword = value;
					value = "0";
					w.capitalize(keyword);
					parameters.insert(pair<string,string>(keyword,value));
				}
			}else{

				// then in this case, we have more classes need to debug printing
				// one is the class name, the other is the debug level
				// if the debug level is not explicitly given, it's always zero
				for(UInt i=1; i<l.getNPieces(); i++) {
					string keyword = l.findValue(i);
					string value   = "0";
					if (w.isInt(keyword)) continue; // this is the debug level indicator (number)
					if (i<l.getNPieces() - 1) {
						string tmpValue = l.findValue(i+1);
						if (w.isInt(tmpValue)) value = tmpValue;
					}
					w.capitalize(keyword);
					parameters.insert(pair<string,string>(keyword,value));
				}
			}

		}else if (l.getNPieces() == 2) {

			// this is the first key-value section: one key to one value
			string keyword = l.findValue(0);
			string value   = l.findValue(1);

			// special case
			w.capitalize(keyword);
			if (keyword != "GUESS_MO_PATH" && keyword != "WFN_PATH" 
					&& keyword != "RESULT_DENSITY_MATRIX_PATH"
					&& keyword != "INITIAL_DENSITY_MATRIX_PATH"
					&& keyword != "GUESS_DENSITY_MATRIX_PATH") {
				w.capitalize(value);
			}
			parameters.insert(pair<string,string>(keyword,value));

		}else {

			// this is another key-value section: one key and multiple value
			// here we will wrap the values together and waiting for further
			// processing
			string keyword = l.findValue(0);
			string value   = l.findValue(1);
			for(UInt i=2; i<l.getNPieces(); i++) {
				value = value + " " + l.findValue(i);
			}
			w.capitalize(keyword);
			w.capitalize(value);
			parameters.insert(pair<string,string>(keyword,value));
		}
	}

	// close the input file
	inf.close();

}

bool ParameterParsing::hasKeyDefined(string oriKey) const {
	WordConvert w;
	w.capitalize(oriKey);
	map<string,string>::const_iterator it;
	it = parameters.find(oriKey);
	if (it == parameters.end()) return false;
	return true;
}

string ParameterParsing::getValue(string oriKey) const {
	WordConvert w;
	w.capitalize(oriKey);
	map<string,string>::const_iterator it;
	it = parameters.find(oriKey);
	if (it == parameters.end()) {
		return "NONE";
	}else{
		return it->second;
	}
}

string ParameterParsing::getKey(const UInt& index) const {

	// firstly we need to check the index number
	// it should be inside the range
	if (index>=getKeyNumber()) {
		Excep excep("ParameterParsing","getKey",EXCEPTION_INPUT_PARAMETER_INVALID, 
				"The input index is out of the range");
		handleExcep(excep);
	}

	// now let's find out the key matches the input index
	UInt n = 0;
	map<string,string>::const_iterator it;
	for (it = parameters.begin(); it != parameters.end(); ++it) {
		if (n == index) {
			return it->first;
		}
		n++;
	}

	// we should not return here
	// in case just return "NONE"
	return "NONE";
}

void ParameterParsing::print() const {
	cout << "List of parameters captured:" << endl;
	map<string,string>::const_iterator it;
	for (it = parameters.begin(); it != parameters.end(); ++it) {
		cout << "Key " << it->first << " Value " << it->second << endl;
	}
}
