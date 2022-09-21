/**
 * CPP files corresponding to the block matrix list class
 * \author fenglai liu
 */
#include<string>
#include<algorithm>
#include<iostream>
#include "boost/lexical_cast.hpp"
#include "excep.h"
#include "blockmatrixlist.h"
using namespace std;
using namespace excep;
using namespace blockmatrixlist;

BlockMtrxList::BlockMtrxList(const UInt& lengthLimit, const UInt& incLengthFactor, 
					const UInt& nRowLimit, const UInt& nColLimit):lenLimit(lengthLimit),curPos(0)
{
	// check the length limit, it should be > 0
	if (lenLimit == 0) {
		string infor = "input length limit is 0 for non-default constructor";
		Excep excep("BlockMtrxList","constructor",EXCEPTION_BLOCK_MATRIX_LIST_CONSTRUCTOR_ERROR,infor);
		handleExcep(excep);
	}

	// now let's initialize the whole stuff
	UInt len = lengthLimit + incLengthFactor;
	BlockMtrx A(nRowLimit,nColLimit);
	for(UInt iBlock=0; iBlock<len; iBlock++) {
		blockList.push_back(A);
	}
}

void BlockMtrxList::initMem(const UInt& lengthLimit, const UInt& incLengthFactor, 
					const UInt& nRowLimit, const UInt& nColLimit)
{
	// set the value
	lenLimit = lengthLimit;
	curPos   = 0;

	// check the length limit, it should be > 0
	if (lenLimit == 0) {
		string infor = "input length limit is 0 for non-default constructor";
		Excep excep("BlockMtrxList","initMem",EXCEPTION_BLOCK_MATRIX_LIST_CONSTRUCTOR_ERROR,infor);
		handleExcep(excep);
	}

	// now let's initialize the whole stuff
	UInt len = lengthLimit + incLengthFactor;
	BlockMtrx A(nRowLimit,nColLimit);
	for(UInt iBlock=0; iBlock<len; iBlock++) {
		blockList.push_back(A);
	}
}

void BlockMtrxList::clear()
{
	// for the simplest case, just delete whatever inside
	if (lenLimit == 0) {
		blockList.clear();
		return;
	}

	// now clear the content
	for(UInt iBlock=0; iBlock<len(); iBlock++) {
		BlockMtrx& B = getBlockMatrix(iBlock);
		B.resetBlockData();
	}

	// finally reset the curpos
	curPos = 0;
}

void BlockMtrxList::update(const BlockMtrx& A)
{

	// this is for the simplest way to update list
	if (lenLimit == 0) {

		// whether the list is empty?
		if (blockList.empty()) {
			blockList.push_back(A);
			return;
		}

		// now we need to do a search
		BlockMatrices::iterator it=find(blockList.begin(),blockList.end(),A);

		// append/update it
		if (it == blockList.end()) {
			blockList.push_back(A);
		}else{
			it->add(A);
		}

		// now return
		return;
	}

	// this is for list has been already been initialized
	for(UInt iBlock=0; iBlock<len(); iBlock++) {
		BlockMtrx& B = getBlockMatrix(iBlock);
		if (B == A) {
			B.add(A);
			return;
		}
	}

	// now let's append the matrix to the end
	// we need to consider the case that memory is not 
	// enough, in this case we still need to push back 
	// the data block
	if (curPos < blockList.size()) {
		BlockMtrx& B = getBlockMatrix(curPos);
		B = A;
	}else{
		blockList.push_back(A);
	}
	curPos++;
}

void BlockMtrxList::merge(const BlockMtrxList& A)
{
	for(UInt iBlock=0; iBlock<A.len(); iBlock++) {
		const BlockMtrx& B = A.getBlockMatrix(iBlock);
		update(B);
	}
}

void BlockMtrxList::updateMatrix(Mtrx& A) const
{
	bool inTranspose = false;
	for(UInt iBlock=0; iBlock<len(); iBlock++) {
		const BlockMtrx& B = getBlockMatrix(iBlock);
		B.updateData(A,inTranspose);
	}
}

bool BlockMtrxList::hasIt(const BlockMtrx& A) const
{
	// whether the list is empty?
	if (blockList.empty()) {
		return false;
	}

	// now we need to do a search if we do not have length limit
	if (lenLimit == 0) {
		BlockMatrices::const_iterator it=find(blockList.begin(),blockList.end(),A);
		if (it == blockList.end()) {
			return false;
		}
		return true;
	}

	// search for the case with length limit
	for(UInt iBlock=0; iBlock<len(); iBlock++) {
		const BlockMtrx& B = getBlockMatrix(iBlock);
		if (B == A) return true;
	}
	return false;
}

void BlockMtrxList::print(const string& title) const
{
	cout << title << endl;
	for(UInt iBlock=0; iBlock<len(); iBlock++) {
		const BlockMtrx& B = getBlockMatrix(iBlock);
		string n = boost::lexical_cast<string>(iBlock);
		string t = "block matrix: " + n;
		B.blockPrint(t);
	}
}

