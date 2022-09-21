/**
 * cpp file for text process
 * \author fenglai liu and jing kong
 */
#include "textread.h"
#include "excep.h"
#include <boost/algorithm/string.hpp>   // string handling
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <cctype>
using namespace textread;
using namespace excep;
using namespace boost;

LineParse::LineParse(const string& input):isEmpty(false),isComment(false),
	isSection(false),isShellSign(false),number(0),line(input) 
{
	trim(line);
	if(line.empty()){ 
		isEmpty = true;
	} else if(line.at(0) == '#' || line.at(0) == '!') {
		isComment = true;

		// some times people may have multiple comment sign
		// let's count it
	   UInt n = 0;
		for(UInt i=0; i<line.size(); i++) {
			if(line.at(i) == '#' || line.at(i) == '!') n++;
		}

		// sometimes we put useful data in comment section
		// therefore we also analyze the commment sections
		string s=line.substr(n);
		trim(s);
		split(pieces, s, is_any_of(" ") || is_any_of(";") || is_any_of(",") 
				|| is_any_of("="), token_compress_on);
		number = pieces.size();
	} else if(line.compare(0,4,"****") == 0) {
		isShellSign = true;
	} else if(line.at(0) == '%') {
		isSection = true;	
		string s = line.substr(1);
		vector<string>sstr;
		split(sstr, s, is_any_of("#") || is_any_of("!"));

		// get the section name and other meaningful stuff
		string s1 = sstr[0];
		trim(s1);
		split(pieces, s1, is_any_of(" ") || is_any_of(";") || is_any_of(","), token_compress_on);
		number = pieces.size();
		sectionName=pieces[0];
		trim(sectionName);
	} else { 
		// two steps: 
		// first one is to check whether it contains some comment 
		// and get rid of the comment part
		// second, split the substr without a comment
		vector<string>sstr;
		split(sstr, line, is_any_of("#") || is_any_of("!"));
		string s=sstr[0];
		trim(s);
		split(pieces, s, is_any_of(" ") || is_any_of(";") || is_any_of(",") 
				|| is_any_of("="), token_compress_on);
		number = pieces.size();
	}
}

string LineParse::findValue(UInt i) const 
{
	if (i>=number) {
		Excep excep("LineParse","findValue",EXCEPTION_VECTOR_OVERFLOW,
				"Pieces array is overflowed");
		handleExcep(excep);
	}
	return pieces.at(i);
}

bool WordConvert::toInt(const string& s, Int& x) {
	try {
		x = lexical_cast<Int>(s);
	} catch(boost::bad_lexical_cast& error) {
		return false;
	}
	return true;
}

bool WordConvert::toUInt(const string& s, UInt& x) {
	bool issize_t = isUInt(s);
	if (! issize_t) {
		x = -1;
		return false;
	}
	try {
		x = lexical_cast<UInt>(s);
	} catch(boost::bad_lexical_cast& error) {
		return false;
	}
	return true;
}

bool WordConvert::toDouble(string s, Double& x) {
	try {
		// for some basis set files, they may contain the fortan
		// type of exponents like 1.0D-2. We have to replace the 
		// D into E so that to make it understood by program
		for (UInt i=0; i<s.length(); i++) {
			if (s.at(i) == 'D' || s.at(i) == 'd') {
				s.at(i) = 'E';
			}
		}
		x = lexical_cast<Double>(s);
	} catch(boost::bad_lexical_cast& error) {
		return false;
	}
	return true;
}

bool WordConvert::compare(string s1, string s2){
	to_upper(s1);
	to_upper(s2);
	trim(s1);
	trim(s2);
	if (s1 == s2) {
		return true;
	} else {
		return false;
	}
}

void WordConvert::capitalize(string& s){
	to_upper(s);
}

bool WordConvert::isInt(const string& s) {
	// do not use it anymore...
	// since the compiler always give warning message...
	/*
		Int x;
		try {
		x = lexical_cast<Int>(s);
		} catch(boost::bad_lexical_cast& error) {
		return false;
		}
		return true;
		*/
	if (s[0] == '+' || s[0] == '-') {
		for(UInt i=1;i<s.size();i++)
			if (!isdigit(s[i])) return false;
	}else{
		for(UInt i=0;i<s.size();i++)
			if (!isdigit(s[i])) return false;
	}
	return true;
}

bool WordConvert::isUInt(const string& s) {
	if (s[0] == '+') {
		for(UInt i=1;i<s.size();i++)
			if (!isdigit(s[i])) return false;
	}else{
		for(UInt i=0;i<s.size();i++)
			if (!isdigit(s[i])) return false;
	}
	return true;
}
