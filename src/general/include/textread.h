/**
 * \file  textread.h
 * \brief classes for reading and parsing the text file
 *
 * Basically, the idea is like this:
 * 
 * All of text file are treated as some ensemble of lines, and each line is
 * some ensemble of words. Hence we have two level of classes, one is used 
 * to describe the "line", the other is used to describe the "word".
 * 
 * What is the word?
 * 
 * Word is the components inside the line, they are separated with each other
 * by comma, space etc. 
 *
 * Algorithm:
 * 
 * Here we consider the "#" as symbols of comments. Hence all of words
 * following the "#", they are considered as comments. ("!" is the comment symbol
 * too)
 * 
 * For spliting the line, we consider the ",", space, "=", ";" as the delimiters
 * to get the words
 * 
 * There's another situation, that we may want to search for a section of "data"
 * or keywords, so we use "%" sign or "****" sign to show that the data section
 * begins (% used in the input, and **** used in the shell data part). More 
 * information please see the sample file.
 *
 * \author Fenglai Liu and Jing Kong
 */

#ifndef TEXTREAD_H
#define TEXTREAD_H
#include "libgen.h"
#include <string>
#include <vector>

namespace textread {

	using namespace std;

	/**
	 * \class LineParse
	 * \brief LineParse is used to parsing a line of information from the given text file.
	 */
	class LineParse {

		private:

			bool           isEmpty;             ///< whether this is empty line
			bool           isComment;           ///< whether this is comment line
			bool           isSection;           ///< beginning of a section
			bool           isShellSign;         ///< whether shell data section begins or ends
			UInt           number;              ///< number of pieces getting from the line
			string         line;                ///< the line read in from text file
			string         sectionName;         ///< the section name
			vector<string> pieces;              ///< pieces getting from parsing the line

		public:

			///
			/// constructor
			///
			LineParse(const string& input); 

			///
			/// destructor
			///
			~LineParse() { };

			///
			/// find the word for the given value
			///
			string findValue(UInt i)  const; 

			///
			/// whether this is a section sign?
			///
			bool isSec() const {
				return isSection;
			};

			///
			/// whether this is a shell sign?
			///
			bool isShlSign() const {
				return isShellSign;
			};

			///
			/// get section name
			///
			string getSecName() const {
				return sectionName;
			};

			///
			/// get number of pieces for this line
			///
			UInt getNPieces () const {
				return number;
			};

			///
			/// whether this is a comment line?
			///
			bool isCom() const {
				return isComment;
			};

			///
			/// whether it's an empty line?
			///
			bool isEmp() const {
				return isEmpty;
			};
	};

	/**
	 * \class  WordConvert
	 *
	 * WordConvert is used to convert the word into different forms. For example,
	 * to convert a string into double or int. This is a utility class
	 */
	class WordConvert {

		public:
			WordConvert(){ };
			~WordConvert(){ };
			bool toUInt(const string& s, UInt& x);
			bool toInt(const string& s, Int& x);
			bool toDouble(string s, Double& x);
			bool compare(string s1, string s2);
			void capitalize(string& s);
			bool isInt(const string& s);
			bool isUInt(const string& s);
	};

}


#endif
