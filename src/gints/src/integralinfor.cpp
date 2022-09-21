/**
 * cpp file associated with integralinfor.h
 * \author Fenglai Liu 
 */
#include <boost/lexical_cast.hpp>
#include<iostream>
#include<fstream>
#include<string>
#include<cstdio>
#include "excep.h"
#include "textread.h"
#include "shellprop.h"
#include "integraljobs.h"
#include "integralinfor.h"
using namespace std;
using namespace excep;
using namespace textread;
using namespace shellprop;
using namespace integraljobs;
using namespace integralinfor;

void ElemDerivInfor::debugPrint() const 
{
	// first pos
	string pos1 = "NULL";
	if(firstDerivPos == BRA1) {
		pos1 = "BRA1";
	}else if(firstDerivPos == BRA2) {
		pos1 = "BRA2";
	}else if (firstDerivPos == KET1) {	
		pos1 = "KET1";
	}else if (firstDerivPos == KET2) {	
		pos1 = "KET2";
	}else if (firstDerivPos == OPERATOR) {	
		pos1 = "OPER";
	}

	// second pos
	string pos2 = "NULL";
	if(secondDerivPos == BRA1) {
		pos2 = "BRA1";
	}else if(secondDerivPos == BRA2) {
		pos2 = "BRA2";
	}else if (secondDerivPos == KET1) {	
		pos2 = "KET1";
	}else if (secondDerivPos == KET2) {	
		pos2 = "KET2";
	}else if (secondDerivPos == OPERATOR) {	
		pos2 = "OPER";
	}

	// first dir
	string dir1 = "NULL";
	if(firstDerivDir == GINTS_DERIV_X) {
		dir1 = "X";
	}else if(firstDerivDir == GINTS_DERIV_Y) {
		dir1 = "Y";
	}else if (firstDerivDir == GINTS_DERIV_Z) {	
		dir1 = "Z";
	}

	// second dir
	string dir2 = "NULL";
	if(secondDerivDir == GINTS_DERIV_X) {
		dir2 = "X";
	}else if(secondDerivDir == GINTS_DERIV_Y) {
		dir2 = "Y";
	}else if (secondDerivDir == GINTS_DERIV_Z) {	
		dir2 = "Z";
	}

	// now print 
	printf("1st deriv pos is: %-4s, dir is: %-4s, 2ed deriv pos is: %-4s, dir is: %-4s\n", 
			pos1.c_str(), dir1.c_str(),pos2.c_str(), dir2.c_str());
}

RedundantDerivInfor::RedundantDerivInfor(const UInt& oper, 
		const UInt& redundantPos0,const ElemDerivInfor& deriv):ElemDerivInfor(deriv),
	redundantPos(redundantPos0),rhsDerivLength(0) 
{ 
	//
	// for energy calculation, or one body integral case
	// we do not have redundant derivatives information defined
	//
	UInt nBody = getOperOrder(oper);
	if (getOrder() == 0 || nBody == 1) {
		return;
	}

	// let's check that whether the 
	// redundant position is defined in base class
	// it's required that the first deriv position
	// should be equal to the redundantPos
	if (firstDerivPos != redundantPos) {
		string infor = "it's required that the first derivative pos is redundant";
		Excep excep("RedundantDerivInfor","constructor",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
		handleExcep(excep);
	}		

	//
	// initialize the array
	//
	ElemDerivInfor deriv0(*this);
	for(UInt i=0; i<MAX_NUM_RHS_DERIV; i++) {
		rhsDeriv[i]    = deriv0;
		rhsDerivPos[i] = -1;
		rhsCoefs[i]    = MINUS_ONE;
	}

	//
	// one thing to remind: for operator of NAI,
	// there will be derivatives on the operator 1/|r-C|
	// therefore the redundant one is derivatives on C 
	// the redundant position is OPERATOR
	//

	// we do the order 1
	if (getOrder() >= 1) {

		// the translation invariance shift the nuclear
		// the dirivative direction of x, y z are same
		// therefore the RHS will be a linear combination
		// of all possible positions except the redundant one
		rhsDerivLength = 0; 
		for(UInt i=0; i<nBody; i++) {

			// set the name of position
			UInt pos = BRA1;
			if (i==1) {
				pos = BRA2;
			}else if (i==2) {
				pos = KET1;
			}else if (i==3) {
				pos = KET2;
			}

			// whether the position is same with redundant pos?
			if (pos == redundantPos) continue;

			// now update the pos
			rhsDeriv[rhsDerivLength].update1stPos(pos);
			rhsDerivLength++; 
		}
	}

	// now let's deal with second order deriv
	// 
	// if second deriv pos is not redundant, then
	// here we do not need to change anything on RHS
	// derived from the first order derivatives
	//
	// however, if it's same with the redundant
	// position, the RHS need to be re-evaluated
	if (getOrder() >= 2 && secondDerivPos == redundantPos) {

		// reset the deriv infor
		// however, here in these cases all of coefs
		// are positive
		ElemDerivInfor deriv0(*this);
		for(UInt i=0; i<MAX_NUM_RHS_DERIV; i++) {
			rhsDeriv[i] = deriv0;
			rhsCoefs[i] = ONE;
		}

		// now generate the derivatives	
		// remember, the derivatives direction of x, y z etc.
		// does not alter
		UInt index = 0;
		for(UInt i=0; i<nBody; i++) {

			// set the name of position
			UInt pos1 = BRA1;
			if (i==1) {
				pos1 = BRA2;
			}else if (i==2) {
				pos1 = KET1;
			}else if (i==3) {
				pos1 = KET2;
			}

			// whether the position is same with redundant pos?
			if (pos1 == redundantPos) continue;

			// over loop the second position
			for(UInt j=0; j<nBody; j++) {

				// set the name of position
				UInt pos2 = BRA1;
				if (j==1) {
					pos2 = BRA2;
				}else if (j==2) {
					pos2 = KET1;
				}else if (j==3) {
					pos2 = KET2;
				}

				// whether the position is same with redundant pos?
				if (pos2 == redundantPos) continue;

				// now let's form the derivInfor
				// for translational invariance the direction does 
				// not alter
				ElemDerivInfor rhs(deriv0);
				rhs.update1stPos(pos1);
				rhs.update2edPos(pos2);

				// is there same term already in the rhs array?
				bool hasThisRHS = false;
				UInt rhsIndex = -1;
				for(UInt k=0; k<index; k++) {
					if (rhsDeriv[k] == rhs) {
						hasThisRHS = true;
						rhsIndex = k;
						break;
					}
				}

				// now let's update
				if (hasThisRHS) {
					rhsCoefs[rhsIndex] += ONE;
				}else{

					// double check the index range
					if (index >= MAX_NUM_RHS_DERIV) {
						string infor = "index is out of boundary for 2ed deriv, firstDerivPos == secondDerivPos";
						Excep excep("RedundantDerivInfor","constructor",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
						handleExcep(excep);
					}

					// now update the term
					rhsDeriv[index] = rhs;
					index++;
				}
			}
		}

		// the final rhs length 
		rhsDerivLength = index;
	}
}

void RedundantDerivInfor::updatePos(const ElemDerivInforVec& derivVec)
{
	for(UInt i=0; i<rhsDerivLength; i++) {

		// search the position
		bool gotIt = false;
		for(UInt j=0; j<derivVec.size(); j++) {
			const ElemDerivInfor& deriv = derivVec[j];
			if (deriv == rhsDeriv[i]) {
				rhsDerivPos[i] = j;
				gotIt = true;
				break;
			}
		}

		// double check
		if (! gotIt) {
			cout << "rhs derivative term " << endl;
			rhsDeriv[i].debugPrint();
			string infor = "fail to get the identical deriv term in the input array";
			Excep excep("RedundantDerivInfor","updatePos",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
			handleExcep(excep);
		}
	}
}

void RedundantDerivInfor::print() const 
{
	// now print out the redundant pos
	string pos1 = "NULL";
	if(redundantPos == BRA1) {
		pos1 = "BRA1";
	}else if(redundantPos == BRA2) {
		pos1 = "BRA2";
	}else if (redundantPos == KET1) {	
		pos1 = "KET1";
	}else if (redundantPos == KET2) {	
		pos1 = "KET2";
	}else if (redundantPos == OPERATOR) {	
		pos1 = "OPER";
	}
	cout << "redundant derivative position is: " << pos1 << endl;

	// now it's RHS
	cout << "the LHS term information: " << endl;
	debugPrint();
	cout << "the RHS term information: " << endl;
	for(UInt i=0; i<rhsDerivLength; i++) {
		printf("term: %-1d, coef: %-2.1f, pos: %-d\n", (Int)(i+1), rhsCoefs[i], (Int)rhsDerivPos[i]);
		rhsDeriv[i].debugPrint();
	}
	cout << endl;
}

////////////////////////////////////////////////////////////////////
//                 #### SingleIntegralInfor ####                  // 
////////////////////////////////////////////////////////////////////
void SingleIntegralInfor::initMemInfor(const string& fileName)
{
	// let me open the file to locate the code section
	ifstream inf;
	const char * input_file = fileName.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "can not open the data file: " + fileName;
		Excep excep("SingleIntegralInfor","initMemInfor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// locate the code
	string line;
	inf.seekg(0,ios::beg);
	bool getIt = false;
	string scode = boost::lexical_cast<string>(LCode);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.getNPieces() >= 2) {
			string val = l.findValue(0);
			if (val == scode) {
				getIt = true;

				// read in the memory size
				// should be following the LCode
				string v = l.findValue(1);
				UInt memSize = boost::lexical_cast<UInt>(v);
				memAllocLen = memSize;

				// now stop reading
				break;
			}
		}
	}
	if (! getIt) {
		string infor = "failed to get the memory data section for the LCode: " + scode;
		Excep excep("SingleIntegralInfor","initMemInfor",EXCEPTION_DATA_SECTION_NOT_FOUND,infor);
		handleExcep(excep);
	}

	// finally close the file
	inf.close();
}

void SingleIntegralInfor::initDerivInfor(const string& fileName)
{
	// let me open the file to locate the code section
	ifstream inf;
	const char * input_file = fileName.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "can not open the data file: " + fileName;
		Excep excep("SingleIntegralInfor","initDerivInfor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// locate the code
	string line;
	inf.seekg(0,ios::beg);
	bool getIt = false;
	string scode = boost::lexical_cast<string>(LCode);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.getNPieces() > 0) {
			string val = l.findValue(0);
			if (val == scode) {
				getIt = true;
				break;
			}
		}
	}
	if (! getIt) {
		string infor = "failed to get the deriv data section for the LCode: " + scode;
		Excep excep("SingleIntegralInfor","initDerivInfor",EXCEPTION_DATA_SECTION_NOT_FOUND,infor);
		handleExcep(excep);
	}

	// reserve space for the array
	// only nuclear attraction does not have 
	// redundant position, because NAI's redundant 
	// position is the nuclei center
	UInt nBody = getOperOrder(oper);
	if (oper != NUCLEAR_ATTRACTION) nBody -= 1; 
	UInt maxNum = 0;
	if (order == 1) {
		maxNum = nBody*3;      // each position, like BRA1 on derivatives x, y and z
	}else if (order == 2) {
		UInt n = nBody*3;
		maxNum = ((1+n)*n)/2;  
	}else {
		string infor = "derivatives order is invalid, right now only support order=1 and order=2";
		Excep excep("SingleIntegralInfor","initDerivInfor", EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
	}
	derivInforArray.reserve(maxNum);

	// now let's read in the data
	WordConvert w;
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isEmp()) break;

		// whether this is line contains the position information?
		string p1 = l.findValue(0);
		w.capitalize(p1);
		if (p1 == "BRA1" || p1 == "BRA2" || p1 == "KET1" || p1 == "KET2") {

			// get the pos1
			UInt pos1 = NULL_POS;
			if (p1 == "BRA1") {
				pos1 = BRA1;
			}else if (p1 == "BRA2") {
				pos1 = BRA2;
			}else if (p1 == "KET1") {
				pos1 = KET1;
			}else if (p1 == "KET2") {
				pos1 = KET2;
			}else {
				cout << "the read in position is: " << p1 << endl;
				string infor = "this position is invalid";
				Excep excep("SingleIntegralInfor","initDerivInfor",
						EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
				handleExcep(excep);
			}

			// let's see whether we have pos2
			UInt pos2 = NULL_POS;
			if (l.getNPieces() >= 2) {
				string p2 = l.findValue(1);
				w.capitalize(p2);
				if (p2 == "BRA1") {
					pos2 = BRA1;
				}else if (p2 == "BRA2") {
					pos2 = BRA2;
				}else if (p2 == "KET1") {
					pos2 = KET1;
				}else if (p2 == "KET2") {
					pos2 = KET2;
				}else {
					cout << "the read in position is: " << p2 << endl;
					string infor = "this position is invalid";
					Excep excep("SingleIntegralInfor","initDerivInfor",
							EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
					handleExcep(excep);
				}
			}

			// now let's read in dir
			string line2;
			while(getline(inf,line2)) {

				// we break when we read in the end sign
				if (line2.find("##") != std::string::npos) {
					break;
				}

				// now prepare for read in
				LineParse l(line2);

				// read in first dir
				UInt dir1 = NO_DERIV;
				string p1 = l.findValue(0);
				if (p1 == "X" || p1 == "x") {
					dir1 = GINTS_DERIV_X;
				}else if (p1 == "Y" || p1 == "y") {
					dir1 = GINTS_DERIV_Y;
				}else if (p1 == "Z" || p1 == "z") {
					dir1 = GINTS_DERIV_Z;
				}else {
					cout << "the read in direction is: " << p1 << endl;
					string infor = "this direction is invalid";
					Excep excep("SingleIntegralInfor","initDerivInfor",
							EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
					handleExcep(excep);
				}

				// read in second dir
				UInt dir2 = NO_DERIV;
				if (l.getNPieces() >= 2) {
					string p2 = l.findValue(1);
					if (p2 == "X" || p2 == "x") {
						dir2 = GINTS_DERIV_X;
					}else if (p2 == "Y" || p2 == "y") {
						dir2 = GINTS_DERIV_Y;
					}else if (p2 == "Z" || p2 == "z") {
						dir2 = GINTS_DERIV_Z;
					}else {
						cout << "the read in direction is: " << p2 << endl;
						string infor = "this direction is invalid";
						Excep excep("SingleIntegralInfor","initDerivInfor",
								EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
						handleExcep(excep);
					}
				}

				// finally update data
				ElemDerivInfor deriv;
				deriv.update1stPos(pos1);
				deriv.update2edPos(pos2);
				deriv.update1stDir(dir1);
				deriv.update2edDir(dir2);
				derivInforArray.push_back(deriv);
			}
		}
	}
	
	// close the file
	inf.close();

	// initialize the redundant derivatives based on
	// what we have read
	// and it finishes all of work
	formRedundantDerivInfor();
}

void SingleIntegralInfor::formRedundantDerivInfor()
{
	// firstly check whether deriv array has something
	if (derivInforArray.size() == 0) {
		string infor = "no unique derivatives information, so it's impossible to form redundant ones";
		Excep excep("SingleIntegralInfor","formRedundantDerivInfor",
				EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
		handleExcep(excep);
	}

	// look for the redundant position
	// here for NAI the redundant position is just operator
	// else we need to check
	UInt redundantPos = NULL_POS;
	if (oper == NUCLEAR_ATTRACTION) {
		redundantPos = OPERATOR;
	}else{

		// firstly we set up vectors to preserve 
		// the position information
		UInt posList[4];
		for(UInt i=0; i<4; i++) posList[i] = 0;
		for(UInt iDeriv=0; iDeriv<derivInforArray.size(); iDeriv++) {
			const ElemDerivInfor& infor = derivInforArray[iDeriv];

			// we only check the first derivative position
			// because it always exist
			// also, because the first pos in general identical
			// to the second one, so it's impossible that 
			// the first derivatives go through some position;
			// but the higher derivatives does not
			// therefore we only check the firstderivPos
			UInt pos = infor.getDerivPos();
			if (pos == BRA1) {
				posList[0] = 1;
			}else if (pos == BRA2) {
				posList[1] = 1;
			}else if (pos == KET1) {
				posList[2] = 1;
			}else if (pos == KET2) {
				posList[3] = 1;
			}
		}

		// now check the result
		// there's only one redundant position
		// which is the position does not show
		// in the unique derivatives
		UInt count = 0;
		UInt nBody = getOperOrder(oper);
		for(UInt i=0; i<nBody; i++) {
			if (posList[i] == 0) count++;
		}
		if (count != 1) {
			if (count == 0) {
				string infor = "can not find the redundant derivative position based on unique ones";
				Excep excep("SingleIntegralInfor","formRedundantDerivInfor",
						EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
			}else{
				string infor = "why there are multiple redundant derivative position exist?";
				Excep excep("SingleIntegralInfor","formRedundantDerivInfor",
						EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
			}
		}

		// now let's identify which is the redundant position
		for(UInt i=0; i<nBody; i++) {

			// get the pos
			UInt pos = BRA1;
			if (i==1) {
				pos = BRA2;
			}else if (i==2) {
				pos = KET1;
			}else if (i==3) {
				pos = KET2;
			}

			// now assign it
			if (posList[i] == 0) {
				redundantPos = pos;
				break;
			}
		}
	}

	// let's compute how many redundant derivatives we have
	// reserve space to form redundant derivatives information
	UInt nRedDeriv = 0;
	UInt nBody = getOperOrder(oper);
	if (order == 1) {
		nRedDeriv = 3;  // derivatives on X, Y and Z
	}else if (order == 2) {
		UInt n = nBody;
		if (oper == NUCLEAR_ATTRACTION) n = nBody+1; // plus the operator
		nRedDeriv = 3*(n*3); // x, y and z for redudant pos times all of other derivatives
	}else {
		string infor = "derivatives order is invalid, right now only support order=1 and order=2";
		Excep excep("SingleIntegralInfor","formRedundantDerivInfor",
				EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
	}
	redDerivInforArray.reserve(nRedDeriv);

	// form the derivatives
	// set up variables first
	UInt pos1 = NULL_POS;
	UInt pos2 = NULL_POS;
	UInt dir1 = NO_DERIV;
	UInt dir2 = NO_DERIV;

	// order 1
	if (order == 1) {
		pos1 = redundantPos;
		for(UInt i=0; i<3; i++) {

			// generate the basic infor
			dir1 = GINTS_DERIV_X;
			if (i == 1) dir1 = GINTS_DERIV_Y;
			if (i == 2) dir1 = GINTS_DERIV_Z;
			ElemDerivInfor deriv;
			deriv.update1stPos(pos1);
			deriv.update1stDir(dir1);

			// now form the redundant one
			RedundantDerivInfor redDeriv(oper,redundantPos,deriv);
			redDeriv.updatePos(derivInforArray);
			redDerivInforArray.push_back(redDeriv);
		}
	}else if (order == 2) {

		// firstly, deal with (redundantPos, *)
		// * is not same with redundant position
		pos1 = redundantPos;
		for(UInt jPos=0; jPos<nBody; jPos++) {
			if (jPos == 0) {
				pos2 = BRA1;
			}else if (jPos == 1) {
				pos2 = BRA2;
			}else if (jPos == 2) {
				pos2 = KET1;
			}else if (jPos == 3) {
				pos2 = KET2;
			}
			if (pos2 == pos1) continue;
			for(UInt i=0; i<3; i++) {
				dir1 = GINTS_DERIV_X;
				if (i == 1) dir1 = GINTS_DERIV_Y;
				if (i == 2) dir1 = GINTS_DERIV_Z;
				for(UInt j=0; j<3; j++) {
					dir2 = GINTS_DERIV_X;
					if (j == 1) dir2 = GINTS_DERIV_Y;
					if (j == 2) dir2 = GINTS_DERIV_Z;

					// now form the basic infor
					ElemDerivInfor deriv;
					deriv.update1stPos(pos1);
					deriv.update1stDir(dir1);
					deriv.update2edPos(pos2);
					deriv.update2edDir(dir2);

					// now form the redundant one
					RedundantDerivInfor redDeriv(oper,redundantPos,deriv);
					redDeriv.updatePos(derivInforArray);
					redDerivInforArray.push_back(redDeriv);
				}
			}
		}

		// now consider the (redundantPos, redundantPos)
		pos1 = redundantPos;
		pos2 = redundantPos;
		for(UInt i=0; i<3; i++) {
			dir1 = GINTS_DERIV_X;
			if (i == 1) dir1 = GINTS_DERIV_Y;
			if (i == 2) dir1 = GINTS_DERIV_Z;
			for(UInt j=i; j<3; j++) {
				dir2 = GINTS_DERIV_X;
				if (j == 1) dir2 = GINTS_DERIV_Y;
				if (j == 2) dir2 = GINTS_DERIV_Z;

				// now form the basic infor
				ElemDerivInfor deriv;
				deriv.update1stPos(pos1);
				deriv.update1stDir(dir1);
				deriv.update2edPos(pos2);
				deriv.update2edDir(dir2);

				// now form the redundant one
				RedundantDerivInfor redDeriv(oper,redundantPos,deriv);
				redDeriv.updatePos(derivInforArray);
				redDerivInforArray.push_back(redDeriv);
			}
		}
	}

	// final check
	// we will generate all of possible derivatives
	// and see whether unique+redundant is able to
	// match the list
	if (order == 1) {
		check1stDeriv();
	}else if (order == 2) {
		check2edDeriv();
	}
}

void SingleIntegralInfor::check1stDeriv() const
{
	// set variable
	UInt pos1 = NULL_POS;
	UInt pos2 = NULL_POS;
	UInt dir1 = NO_DERIV;
	UInt dir2 = NO_DERIV;

	// now generate possible derivatives
	UInt nBody = getOperOrder(oper);
	UInt len   = 0;
	if (oper == NUCLEAR_ATTRACTION) nBody += 1; // plus the operator
	for(UInt iPos=0; iPos<nBody; iPos++) {

		// set pos1
		if (iPos == 0) {
			pos1 = BRA1;
		}else if (iPos == 1) {
			pos1 = BRA2;
		}else if (iPos == 2) {
			pos1 = KET1;
		}else if (iPos == 3) {
			pos1 = KET2;
		}

		// for NAI, iPos is 2 is the operator
		if (oper == NUCLEAR_ATTRACTION && iPos == 2) {
			pos1 = OPERATOR;
		}

		// now generate the dir
		for(UInt i=0; i<3; i++) {
			dir1 = GINTS_DERIV_X;
			if (i == 1) dir1 = GINTS_DERIV_Y;
			if (i == 2) dir1 = GINTS_DERIV_Z;

			// now form the basic infor
			ElemDerivInfor deriv;
			deriv.update1stPos(pos1);
			deriv.update1stDir(dir1);
			deriv.update2edPos(pos2);
			deriv.update2edDir(dir2);
			len++;

			// can we find it in the deriv vector?
			bool findIt = false;
			for(UInt iDeriv=0; iDeriv<derivInforArray.size(); iDeriv++) {
				const ElemDerivInfor& infor = derivInforArray[iDeriv];
				if (deriv == infor) {
					findIt = true;
					break;
				}
			}

			// if find it let's continue
			if (findIt) continue;

			// reset and try to find it in redundant ones
			findIt = false;
			for(UInt iDeriv=0; iDeriv<redDerivInforArray.size(); iDeriv++) {
				const RedundantDerivInfor& infor = redDerivInforArray[iDeriv];
				if (deriv == infor) {
					findIt = true;
					break;
				}
			}

			// if find it let's continue
			if (findIt) continue;

			// now we need to issue an error
			cout << "missing derivatives information:" << endl;
			deriv.debugPrint();
			string infor = "the first order derivatives check failed";
			Excep excep("SingleIntegralInfor","check1stDeriv", 
					EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
		}
	}

	// finally check the length
	UInt nDeriv = redDerivInforArray.size() + derivInforArray.size();
	if (nDeriv != len) {
		string infor = "the number of 2ed derivatives does not equal to the calculated ones";
		Excep excep("SingleIntegralInfor","check2edDeriv",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
	}
}

void SingleIntegralInfor::check2edDeriv() const
{
	// set variable
	UInt pos1 = NULL_POS;
	UInt pos2 = NULL_POS;
	UInt dir1 = NO_DERIV;
	UInt dir2 = NO_DERIV;

	// now generate possible derivatives
	UInt nBody = getOperOrder(oper);
	UInt len   = 0;
	for(UInt iPos=0; iPos<nBody; iPos++) {

		// set pos1
		if (iPos == 0) {
			pos1 = BRA1;
		}else if (iPos == 1) {
			pos1 = BRA2;
		}else if (iPos == 2) {
			pos1 = KET1;
		}else if (iPos == 3) {
			pos1 = KET2;
		}

		for(UInt jPos=iPos; jPos<nBody; jPos++) {

			// set pos2
			if (jPos == 0) {
				pos2 = BRA1;
			}else if (jPos == 1) {
				pos2 = BRA2;
			}else if (jPos == 2) {
				pos2 = KET1;
			}else if (jPos == 3) {
				pos2 = KET2;
			}

			// now generate the dir
			for(UInt i=0; i<3; i++) {
				dir1 = GINTS_DERIV_X;
				if (i == 1) dir1 = GINTS_DERIV_Y;
				if (i == 2) dir1 = GINTS_DERIV_Z;

				// for the second pos, if it's same
				// with the first one, we need to
				// take care of symmetry
				UInt begin = 0;
				if (pos1 == pos2) begin = i;
				for(UInt j=begin; j<3; j++) {
					dir2 = GINTS_DERIV_X;
					if (j == 1) dir2 = GINTS_DERIV_Y;
					if (j == 2) dir2 = GINTS_DERIV_Z;

					// now form the basic infor
					ElemDerivInfor deriv;
					deriv.update1stPos(pos1);
					deriv.update1stDir(dir1);
					deriv.update2edPos(pos2);
					deriv.update2edDir(dir2);
					len++;

					// can we find it in the deriv vector?
					bool findIt = false;
					for(UInt iDeriv=0; iDeriv<derivInforArray.size(); iDeriv++) {
						const ElemDerivInfor& infor = derivInforArray[iDeriv];
						if (deriv == infor) {
							findIt = true;
							break;
						}
					}

					// if find it let's continue
					if (findIt) continue;

					// reset and try to find it in redundant ones
					findIt = false;
					for(UInt iDeriv=0; iDeriv<redDerivInforArray.size(); iDeriv++) {
						const RedundantDerivInfor& infor = redDerivInforArray[iDeriv];
						if (deriv == infor) {
							findIt = true;
							break;
						}
					}

					// if find it let's continue
					if (findIt) continue;

					// now we need to issue an error
					cout << "missing derivatives information:" << endl;
					deriv.debugPrint();
					string infor = "the second order derivatives check failed";
					Excep excep("SingleIntegralInfor","check2edDeriv", 
							EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
				}
			}
		}
	}

	// now we may need to do additional check
	if (oper == NUCLEAR_ATTRACTION) {
		pos1 = OPERATOR;
		UInt nBody = 3;
		for(UInt jPos=0; jPos<nBody; jPos++) {

			// set pos2
			if (jPos == 0) {
				pos2 = BRA1;
			}else if (jPos == 1) {
				pos2 = BRA2;
			}else if (jPos == 2) {
				pos2 = OPERATOR;
			}

			// now generate the dir
			for(UInt i=0; i<3; i++) {
				dir1 = GINTS_DERIV_X;
				if (i == 1) dir1 = GINTS_DERIV_Y;
				if (i == 2) dir1 = GINTS_DERIV_Z;

				// for the second pos, if it's same
				// with the first one, we need to
				// take care of symmetry
				UInt begin = 0;
				if (pos1 == pos2) begin = i;
				for(UInt j=begin; j<3; j++) {
					dir2 = GINTS_DERIV_X;
					if (j == 1) dir2 = GINTS_DERIV_Y;
					if (j == 2) dir2 = GINTS_DERIV_Z;

					// now form the basic infor
					ElemDerivInfor deriv;
					deriv.update1stPos(pos1);
					deriv.update1stDir(dir1);
					deriv.update2edPos(pos2);
					deriv.update2edDir(dir2);
					len++;

					// can we find it in the deriv vector?
					bool findIt = false;
					for(UInt iDeriv=0; iDeriv<derivInforArray.size(); iDeriv++) {
						const ElemDerivInfor& infor = derivInforArray[iDeriv];
						if (deriv == infor) {
							findIt = true;
							break;
						}
					}

					// if find it let's continue
					if (findIt) continue;

					// reset and try to find it in redundant ones
					findIt = false;
					for(UInt iDeriv=0; iDeriv<redDerivInforArray.size(); iDeriv++) {
						const RedundantDerivInfor& infor = redDerivInforArray[iDeriv];
						if (deriv == infor) {
							findIt = true;
							break;
						}
					}

					// if find it let's continue
					if (findIt) continue;

					// now we need to issue an error
					cout << "missing derivatives for NAI:" << endl;
					deriv.debugPrint();
					string infor = "the second order derivatives check failed";
					Excep excep("SingleIntegralInfor","check2edDeriv", 
							EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
				}
			}
		}
	}

	// finally check the length
	UInt nDeriv = redDerivInforArray.size() + derivInforArray.size();
	if (nDeriv != len) {
		string infor = "the number of 2ed derivatives does not equal to the calculated ones";
		Excep excep("SingleIntegralInfor","check2edDeriv",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
	}
}

void SingleIntegralInfor::print() const
{
	cout << "*********************************************************" << endl;
	cout << "integral LCode is: " << LCode << endl;
	cout << "integral job is: " << oper << endl;
	cout << "derivative order is: " << order << endl;
	cout << "memory length allocated by the integral engine: " << memAllocLen << endl;
	cout << endl;

	if (order>0) {
		cout << "now it's unique derivatives information" << endl;
		for(UInt iDeriv=0; iDeriv<derivInforArray.size(); iDeriv++) {
			const ElemDerivInfor& infor = derivInforArray[iDeriv];
			infor.debugPrint();
		}
		cout << endl;
		cout << "now it's redundant derivatives information" << endl;
		for(UInt iDeriv=0; iDeriv<redDerivInforArray.size(); iDeriv++) {
			const RedundantDerivInfor& infor = redDerivInforArray[iDeriv];
			infor.print();
		}
		cout << endl;
	}
}

////////////////////////////////////////////////////////////////////
//                       #### IntegralInfor ####                  // 
////////////////////////////////////////////////////////////////////
IntegralInfor::IntegralInfor(const UInt& oper0, const UInt& order0):oper(oper0),
	order(order0),maxL(0),maxMemLenForSP(0),maxMemLenForD(0),
	maxMemLenForF(0),maxMemLenForG(0),maxMemLenForH(0),maxMemLenForI(0),
	nElemDerivs(0),nRedDerivs(0)
{
	// get the operator name for the setting file
	string operName = getOperNameInSettingFile(oper); 

	// we need to get the setting folder dir
	char* path = getenv("EMUL_SETTING_FILE");
	if (path==NULL) {
		Excep excep("IntegralInfor","constructor",EXCEPTION_DIR_MISSING,
				"EMUL_SETTING_FILE is not defined");
		handleExcep(excep);
	}
	string p(path);

	// now set the setting file name
	// here we need to be careful
	// for the 2ERI and 3ERI job setting, the name of memfile
	// and the derivFile is different
	// you need to take care of this in the future
	string memFile   = p + "/mem_infor_" + operName + "_" + boost::lexical_cast<string>(order) + ".txt";
	string derivFile = p + "/deriv_infor_" + operName + "_" + boost::lexical_cast<string>(order) + ".txt";
	
	// let's count how many integral we have
	// we do it by searching memory file
	ifstream inf;
	const char * input_file = memFile.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		cout << "please check the file in the setting file folder you defined in .bashrc etc." << endl;
		string infor = "can not open the data file: " + memFile;
		Excep excep("IntegralInfor","constructor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// let's count
	inf.seekg(0,ios::beg);
	UInt nCodes = 0;
	string line;
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.getNPieces() > 0) nCodes += 1;
	}
	inf.clear();

	// now let's form the integral infor 
	// for each LCode
	// we also check the maxL over there
	integralsInfor.reserve(nCodes); 
	inf.seekg(0,ios::beg);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.getNPieces() >= 2) {
			string val = l.findValue(0);
			LInt code = boost::lexical_cast<LInt>(val);
			SingleIntegralInfor singleIntInfor(code,oper,order);
			integralsInfor.push_back(singleIntInfor);

			// read in the second value, which is the memory length
			val = l.findValue(1);
			UInt memLength = boost::lexical_cast<UInt>(val);

			// double check the length of the line
			UInt len = l.getNPieces();
			if (len<=2) {
				cout << line << endl;
				string infor = "this line in the memory setting file is not valid, "
					"we do not have angular momentum infor here";
				Excep excep("IntegralInfor","constructor",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
				inf.close();
				handleExcep(excep);
			}

			// now check the other information on the line, 
			// which is related to the L
			// we start from the 3rd element
			UInt maxAng = 0;
			for(UInt i=2; i<len; i++) {
				string val = l.findValue(i);
				UInt L = boost::lexical_cast<UInt>(val);
				UInt lmin = -1;
				UInt lmax = -1;
				decodeL(L,lmin,lmax);
				if (lmax>maxAng) maxAng = lmax;
			}

			// now let's see the maxL
			if (maxAng>maxL) maxL = maxAng;

			// also set up the memory useage for each shell
			if (maxAng == 0 || maxAng == 1) {
				if (maxMemLenForSP<memLength) maxMemLenForSP = memLength;  
			}else if (maxAng == 2) {
				if (maxMemLenForD<memLength) maxMemLenForD = memLength;  
			}else if (maxAng == 3) {
				if (maxMemLenForF<memLength) maxMemLenForF = memLength;  
			}else if (maxAng == 4) {
				if (maxMemLenForG<memLength) maxMemLenForG = memLength;  
			}else if (maxAng == 5) {
				if (maxMemLenForH<memLength) maxMemLenForH = memLength;  
			}else if (maxAng == 6) {
				if (maxMemLenForI<memLength) maxMemLenForI = memLength;  
			}else{
				cout << line << endl;
				string infor = "this line in the memory setting file is not valid, "
					"the angular momentum value is higher than what we can record(maxL = 6)";
				Excep excep("IntegralInfor","constructor",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
				inf.close();
				handleExcep(excep);
			}
		}
	}

	// finally close the file
	inf.close();

	// now let's form the mem infor and deriv infor
	for(UInt i=0; i<integralsInfor.size(); i++) {

		// update the memory information
		integralsInfor[i].initMemInfor(memFile);

		// update derivatives information
		if (order>0) {
			integralsInfor[i].initDerivInfor(derivFile);
		}
	}

	// generate the number of elementary derivatives and redundant derivatives
	// we note, that such number should be same for all of integral codes
	// we will check it here, too
	if (order>0) {

		// get the number
		const SingleIntegralInfor& intInfor0 = integralsInfor[0];
		nElemDerivs = intInfor0.getElemDerivInforLen();
		nRedDerivs  = intInfor0.getRedDerivInforLen();

		// now let's check all of other integral infor
		for(UInt i=1; i<integralsInfor.size(); i++) {
			const SingleIntegralInfor& intInfor = integralsInfor[i];
			UInt nElemDerivs1 = intInfor.getElemDerivInforLen();
			UInt nRedDerivs1  = intInfor.getRedDerivInforLen();
			if (nElemDerivs != nElemDerivs1 || nRedDerivs != nRedDerivs1) {
				string infor = "the number of elementary derivatives and redundant derivatives "
					"are different from one integral code to another, which is wrong";
				cout << "the original integral code information, where we get nElemDerivs etc. " << endl;
				intInfor0.print();
				cout << "the integral code information, which get different nElemDerivs etc. " << endl;
				intInfor.print();
				Excep excep("IntegralInfor","constructor",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
			}
		}
	}
}

const SingleIntegralInfor& IntegralInfor::getIntegralInfor(const LInt& code) const
{
	// now search the integral infor
	for(UInt i=0; i<integralsInfor.size(); i++) {
		LInt LCode = integralsInfor[i].getLCode();
		if (LCode == code) {
			return integralsInfor[i];
		}
	}

	// we should not get here
	// so we issue an error
	cout << "the LCode for searching is: " << code << endl;
	string infor = "can not find the SingleIntegralInfor for the given LCode";
	Excep excep("IntegralInfor","getIntegralInfor",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
	handleExcep(excep);
	return integralsInfor[0];
}

UInt IntegralInfor::getMaxMemLen(UInt L) const 
{
	// we check whether L is larger than the maxL
	// in this case we do not have any record data
	// so report error
	if (L>maxL) {
		string infor = "the input L is larger than the maximum angular momentum of integral engine";
		Excep excep("IntegralInfor","getMaxMemLen",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	// now let's return the memory usage for each shell type
	if (L == 0 || L == 1) {
		return maxMemLenForSP;  
	}else if (L == 2) {
		return maxMemLenForD;  
	}else if (L == 3) {
		return maxMemLenForF;  
	}else if (L == 4) {
		return maxMemLenForG;  
	}else if (L == 5) {
		return maxMemLenForH;  
	}else if (L == 6) {
		return maxMemLenForI;  
	}

	// we should not arrive here
	// in case we get here, return wrong value
	// and report error
	string infor = "we did not identify the memory length correctly";
	Excep excep("IntegralInfor","getMaxMemLen",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
	handleExcep(excep);
	return -1;
}

UInt IntegralInfor::getMemUsage(const LInt& code) const
{
	// now loop over to the the integral infor
	for(UInt i=0; i<integralsInfor.size(); i++) {
		LInt LCode = integralsInfor[i].getLCode();
		if (LCode == code) {
			return integralsInfor[i].getMemInfor();
		}
	}

	// we should not arrive here
	// fail to get the L Code matching the input code
	string infor = "fail to get the integral infor to match the given LCode";
	Excep excep("IntegralInfor","getMemUsage",EXCEPTION_DERIV_INFORMATION_VIOLATION,infor);
	handleExcep(excep);
	return -1;
}

void IntegralInfor::print() const 
{
	cout << "integral job is: " << oper << endl;
	cout << "derivative order is: " << order << endl;
	cout << "maximum L supported by the integral engine: " << maxL << endl;
	if (maxL>=1) {
		cout << "maximum memory length allocated by the integral engine for maxL = S/P: " << maxMemLenForSP << endl;
	}
	if (maxL>=2) {
		cout << "maximum memory length allocated by the integral engine for maxL = D  : " << maxMemLenForD << endl;
	}
	if (maxL>=3) {
		cout << "maximum memory length allocated by the integral engine for maxL = F  : " << maxMemLenForF << endl;
	}
	if (maxL>=4) {
		cout << "maximum memory length allocated by the integral engine for maxL = G  : " << maxMemLenForG << endl;
	}
	if (maxL>=5) {
		cout << "maximum memory length allocated by the integral engine for maxL = H  : " << maxMemLenForH << endl;
	}
	if (maxL>=6) {
		cout << "maximum memory length allocated by the integral engine for maxL = I  : " << maxMemLenForI << endl;
	}
	cout << "number of integrals is: " << integralsInfor.size() << endl;
	cout << endl;
	for(UInt i=0; i<integralsInfor.size(); i++) {
		integralsInfor[i].print();
	}
}
