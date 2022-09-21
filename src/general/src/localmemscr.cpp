/**
 * cpp file for manipulating local memory for single thread
 * \author fenglai liu and jing kong
 */
#include <cstdio>
#include <boost/lexical_cast.hpp>
#include "excep.h"
#include "localmemscr.h"
using namespace excep;
using namespace localmemscr;

Double* LocalMemScr::getNewMemPos(const UInt& len) 
{
	// firstly check that whether there's no memory usr
	// if in this case, things is very simple
	if (lastMemLen == 0) {

		// check that whether the length is less than the 
		// whole memory length
		if (len>memHolder.size()) {
			string infor = "required memory is larger than what we have in total";
			Excep excep("LocalMeMScr","getNewMemPos",EXCEPTION_VECTOR_OVERFLOW,infor);
			handleExcep(excep);
		}

		// now it's safe to do the work
		lastMemLen = len;
		return &memHolder[0];
	}

	// secondly, let's go to see that wether we can
	// append a new memory position to the memory user
	UInt lastPos    = lastMemUser;
	UInt lastLen    = lastMemLen;
	UInt newPos     = lastPos+lastLen;
	if (newPos>=memHolder.size() || newPos+len>=memHolder.size()) {
		string infor = "The memory poition causing collapse is " + 
			boost::lexical_cast<string>(newPos) + " it's length is " + boost::lexical_cast<string>(len);
		Excep excep("LocalMeMScr","getNewMemPos",EXCEPTION_VECTOR_OVERFLOW,infor);
		handleExcep(excep);
	}
	lastMemUser = newPos;
	lastMemLen  = len;
	return &memHolder[newPos];
}

