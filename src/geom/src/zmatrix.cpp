/**
 * \file zmatrix.cpp
 * \brief  CPP file corresponding to the zmatrix.h
 * \author Fenglai Liu 
 */
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cmath>
#include <boost/lexical_cast.hpp>
#include "element.h"
#include "excep.h"
#include "textread.h"
#include "molecule.h"
#include "geomutil.h"
#include "zmatrix.h"
using namespace element;
using namespace excep;
using namespace textread;
using namespace molecule;
using namespace geomutil;
using namespace zmatrix;

void ZMatrixItem::print() const {
	if (isNormalAtom()) {
		printf("%-6d %-6d %-12.6f %-6d %-12.6f %-6d %-12.6f\n",(Int)atomic,(Int)atom2+1,bond,(Int)atom3+1,
				angle,(Int)atom4+1,dihedral);
	}else if (is1stAtom()) {
		printf("%-6d\n",(Int)atomic);
	}else if (is2edAtom()) {
		printf("%-6d %-6d %-12.6f\n",(Int)atomic,(Int)atom2+1,bond);
	}else if (is3rdAtom()) {
		printf("%-6d %-6d %-12.6f %-6d %-12.6f\n",(Int)atomic,(Int)atom2+1,bond,(Int)atom3+1,angle);
	}
}

UInt ZMatrix::getSecLoc(ifstream& inf, UInt section) const {

	// get the key word
	string key = "molecule";
	if (section == 0) key = "cluster";

	// additional check on section
	if (section>=MAX_MOLECULE_SECTION_NUMBER) {
		string infor = "given section is: " + boost::lexical_cast<string>(section); 
		Excep excep("ZMatrix","getSecLoc",EXCEPTION_ILLEGAL_GEOM_SECTION_INDEX,infor);
		handleExcep(excep);
	}

	// find the keyword
	UInt n = 0;
	bool find = false;
	string line;
	WordConvert w;
	inf.seekg(0,ios::beg);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isSec() && w.compare(l.getSecName(), key)) {
			if (section == 0) {
				find = true;
				break;
			}else{
				n++;
				if (section == n) {
					find = true;
					break;
				}
			}
		}
	}

	// checking the finding status
	if (! find) {
		string infor = "section number is ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("ZMatrix","getSecLoc",EXCEPTION_SECTION_SEARCH_FAILED,infor);
		handleExcep(excep);
	}

	// remember where we are
	UInt loc = inf.tellg(); 
	return loc;
}

void ZMatrix::readInVarValues(ifstream& inf, UInt loc, 
		vector<string>& vars, DoubleVec& vals) const 
{

	// go to the given location
	inf.seekg(loc,ios::beg);

	// the following line is the charge and multiplicity
	// let's by pass it
	string line;
	getline(inf,line);

	// let's try to bypass the geometry section 
	bool hasVar = true;
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isEmp() || l.isCom()) break;
		if (l.isSec()) {
			hasVar = false;
			break;
		}
	}

	// if we do not have variables, just return
	if (! hasVar) return;

	// now let's read in the variable and it's value
	WordConvert w;
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isSec()) break;

		// for the variable definition, only two data fields
		if(l.getNPieces() == 2) {

			// the variable name
			string var = l.findValue(0);
			w.capitalize(var);
			vars.push_back(var);

			// read in double value
			Double val = ZERO;
			if (!w.toDouble(l.findValue(1), val)) {
				string infor = line;
				cout << "the input data field is not a double value" << endl;
				Excep excep("ZMatrix","readInVarValues",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			vals.push_back(val);
		}else{
			string infor = line;
			cout << "the given line of data can not be processed" << endl;
			Excep excep("ZMatrix","readInVarValues",EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR,infor);
			handleExcep(excep);
		}
	}
}

void ZMatrix::readInGeom(ifstream& inf, UInt loc,
		const vector<string>& vars, const DoubleVec& vals) 
{
	// go to the given location
	inf.seekg(loc,ios::beg);

	// the following line is the charge and multiplicity
	// let's by pass it
	string line;
	getline(inf,line);

	// now let's read in the geometry data
	// please be aware that the values in the geometry
	// set could be variables, which means it's "string"
	UInt atom1Index = 0;
	WordConvert w;
	while(getline(inf,line)) {

		// begin the real parsing
		LineParse l(line);

		// where we stop reading?
		if (l.isEmp() || l.isCom() || l.isSec()) break;

		// read in the atomic number, here the read in string may be
		// the atomic number or the atom symbol
		string s = l.findValue(0);
		UInt atomic = -1;
		if (isAtomicSymbol(s)) {
			atomic = getAtomicNumber(s);
		}else{
			if (w.isUInt(s)) {
				w.toUInt(s,atomic);
				if (atomic>MAX_ATOMS_TYPE) {
					string infor = line;
					Excep excep("ZMatrix","readInGeom",INVALID_ATOMIC_NUMBER,infor);
					handleExcep(excep);
				}
			} else {
				string infor = line;
				Excep excep("ZMatrix","readInGeom",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
		} 

		// read in the data field
		UInt atom2Index = -1;
		UInt atom3Index = -1;
		UInt atom4Index = -1;
		Double bondVal  = ZERO;
		Double angleVal = ZERO;
		Double dihedralVal = ZERO;
		for(UInt iField=1; iField<l.getNPieces(); iField++) {
			
			// read in the atom index
			if (iField == 1 || iField == 3 || iField == 5) {

				// read in data
				UInt atomIndex = -1;
				string s = l.findValue(iField);
				if (w.isUInt(s)) {
					w.toUInt(s,atomIndex);
					if (atomIndex == 0) {
						string infor = "the atom index should be at least 1, see the line: "+ line;
						Excep excep("ZMatrix","readInGeom",EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR,infor);
						handleExcep(excep);
					}else if (atomIndex>=atom1Index+1) {
						string infor = "the atom index should < the first atom index, see the line: "+ line;
						Excep excep("ZMatrix","readInGeom",EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR,infor);
						handleExcep(excep);
					}
				} else {
					string infor = line;
					cout << "the input data field should be atom index" << endl;
					Excep excep("ZMatrix","readInGeom",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
					handleExcep(excep);
				}
				atomIndex = atomIndex-1;

				// now let's assign the data
				if (iField == 1) {
					atom2Index = atomIndex;
				}else if (iField == 3) {
					atom3Index = atomIndex;
				}else {
					atom4Index = atomIndex;
				}
			}else{

				// let's see whether it data is a string
				string s = l.findValue(iField);
				Double val = ZERO;
				if (! w.toDouble(s, val)) {
					
					// firstly change it to upper case
					w.capitalize(s);
					
					// now let's search the input vars
					UInt pos = vars.size() + 1;
					for(UInt iVar=0; iVar<vars.size(); iVar++){
						if (vars[iVar] == s) {
							pos = iVar;
							break;
						}
					}

					// double check whether we got the variable defined?
					if(pos == vars.size() + 1) {
						string infor = "in this line the variable " + s + " is not defined below: "+ line;
						Excep excep("ZMatrix","readInGeom",EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR,infor);
						handleExcep(excep);
					}

					// now let's reset the val
					val = vals[pos];
				}

				// now let's double check
				if (iField == 2) {
					if (val<ZERO) {
						string infor = "the given bond distance is less than 0";
						Excep excep("ZMatrix","readInGeom",EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR,infor);
						handleExcep(excep);
					}else{
						bondVal = val;
					}
				} else if (iField == 4) {
					if(val<ZERO || val>180.0E0) {
						string infor = "the given angle should be between 0 and 180";
						Excep excep("ZMatrix","readInGeom",EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR,infor);
						handleExcep(excep);
					}else{
						angleVal = val;
					}
				} else{
					dihedralVal = val;
				}
			}
		}

		// now let's make the z matrix item
		ZMatrixItem zmatrix(atomic,atom2Index,atom3Index,atom4Index,
				bondVal,angleVal,dihedralVal);
		zmatrixList.push_back(zmatrix);
		atom1Index++;
	}
}

ZMatrix::ZMatrix(const string& input, UInt section)
{
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "section number is ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("ZMatrix","Constructor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// now read in molecule data
	// get the section location
	UInt loc = getSecLoc(inf,section);

	// we may have variables defined in geometry section
	// so let's read it first
	vector<string> vars;
	DoubleVec vals;
	readInVarValues(inf,loc,vars,vals); 

	// now read in the z matrix data
	readInGeom(inf,loc,vars,vals); 

	// close the input file
	inf.close();
}

void ZMatrix::loadInToXYZ(DoubleVec& coord) const
{
	// determine how much atoms we have
	UInt n = getNAtoms();

	// let's check that whether coord has enough memory
	// to hold the result
	if (coord.size()/3 < n) {
		string infor = "the input xyz vector's length is less than nAtoms";
		Excep excep("ZMatrix","loadInToXYZ",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	// the first atom is in origin
	coord[0] = ZERO;
	coord[1] = ZERO;
	coord[2] = ZERO;

	// the second atom, placed on X axis, and always > 0
	if (n >= 2) {
		const ZMatrixItem& item = getZMatrix(1);
		Double bond = item.getBondDistance();
		coord[3] = bond;
		coord[4] = ZERO;
		coord[5] = ZERO;
	} else {
		return;
	}

	// the third atom,placed on the xy plane, y>0
	if (n >= 3) {
		const ZMatrixItem& item = getZMatrix(2);
		UInt atom2   = item.getAtom2Index(); 
		Double bond  = item.getBondDistance();
		Double angle = item.getAngleVal();
		Double x = ZERO;
		Double y = ZERO;
		Double z = ZERO;
		if (atom2 == 0) {
			// the bond formed is between this atom and the first atom
			x = bond*cos(angle*DEGREE_TO_RADIAN);
			y = bond*sin(angle*DEGREE_TO_RADIAN);
		}else{
			// the bond formed is between this atom and the second atom
			const ZMatrixItem& item1 = getZMatrix(1);
			Double bond1 = item1.getBondDistance();
			x = bond1-bond*cos(angle*DEGREE_TO_RADIAN);
			y = bond*sin(angle*DEGREE_TO_RADIAN);
		}
		coord[6] = x;
		coord[7] = y;
		coord[8] = z;
	} else {
		return;
	}

	// let me set the precision for transformation
	// we think that 1.0E-10 should be good enough
	const Double MEASURE_PREC = 1.0E-10;

	// deal with the other atoms
	if (n > 3) {
		for (UInt iAtom=3; iAtom<n; iAtom++) {

			//
			// here we have such sequence of atom:
			// atom   1  (bondlength)  2  (angle)  3  (dihedral)
			//the atom's coordinate is what we want to seek (x,y,z)
			//

			//
			// firstly, let's derive all of these cordinates
			//
			const ZMatrixItem& item = getZMatrix(iAtom);
			UInt atom2 = item.getAtom2Index(); // atom connecting with bond
			UInt atom3 = item.getAtom3Index(); // atom connecting with angle
			UInt atom4 = item.getAtom4Index(); // atom connecting with dihedral
			Double x1  = coord[3*atom2+0];
			Double y1  = coord[3*atom2+1];
			Double z1  = coord[3*atom2+2];
			Double x2  = coord[3*atom3+0];
			Double y2  = coord[3*atom3+1];
			Double z2  = coord[3*atom3+2];
			Double x3  = coord[3*atom4+0];
			Double y3  = coord[3*atom4+1];
			Double z3  = coord[3*atom4+2];

			// set variables
			Double angleX  = ZERO;
			Double cosalfa = ZERO;
			bool labelX = false;
			Double angleY  = ZERO;
			Double cosbeta = ZERO;
			bool labelY = false;
			Double angleZ  = ZERO;
			Double cosgama = ZERO;
			bool labelZ = false;		// used in rotation around axis

			// move the second atom to the origin, so as the first atom and the third atom
			x1 -= x2;
			y1 -= y2;
			z1 -= z2;
			x3 -= x2;
			y3 -= y2;
			z3 -= z2;
			Double distance = sqrt(x1*x1+y1*y1+z1*z1);

			// move the line between 1 and 2 to the xz plane, z>0
			Double dis = sqrt(x1*x1+y1*y1);
			if(fabs(dis)>MEASURE_PREC) {
				cosgama = x1/dis;
				if (y1<ZERO) {
					//counter-clockwise
					angleZ =  acos(cosgama);
				}else {
					//clockwise
					angleZ = -acos(cosgama);
				}
				rotationByZ(x1, y1, angleZ);
				rotationByZ(x3, y3, angleZ);
			} else {
				labelZ = true;
			}

			// move the line between 1 and 2 to the x axis, x>0
			dis = sqrt(x1*x1+z1*z1);
			if(fabs(dis) > MEASURE_PREC) {
				cosbeta = x1/dis;
				if (z1>ZERO) {
					//counter-clockwise
					angleY =  acos(cosbeta);
				}else {
					//clockwise
					angleY = -acos(cosbeta);
				}
				rotationByY (x1, z1, angleY);
				rotationByY (x3, z3, angleY);
			} else {
				labelY = true;
			}

			// move the line between 2 and 3 to the xy plane, y>0
			dis = sqrt(y3*y3+z3*z3);
			if(fabs(dis) > MEASURE_PREC) {
				cosalfa = y3/dis;
				if (z3>ZERO) {
					// clockwise
					angleX = -acos(cosalfa);
				}else {
					// counter-clockwise
					angleX = acos(cosalfa);
				}
				rotationByX (y3, z3, angleX);
			} else {
				labelX = true;
			}

			// restrict the dihedral within -180 to 180
			Double dihedral = item.getDihedralVal();
			if (dihedral <-180.0E0) {
				while (true) {
					if (dihedral > -180.0E0 || fabs(dihedral+180.0E0)<MEASURE_PREC) {
						break;
					} else {
						dihedral = dihedral + 360.0E0;
					}
				}
			} else if (dihedral > 180.0E0) {
				while (true) {
					if (dihedral < 180.0E0 || fabs(dihedral-180.0E0)<MEASURE_PREC) {
						break;
					} else {
						dihedral = dihedral - 360.0E0;
					}
				}
			}

			// calculate the x, y, z value 
			Double bond  = item.getBondDistance();
			Double angle = item.getAngleVal();
			Double z = bond*sin(angle*DEGREE_TO_RADIAN)*sin(dihedral*DEGREE_TO_RADIAN);
			Double y = bond*sin(angle*DEGREE_TO_RADIAN)*cos(dihedral*DEGREE_TO_RADIAN);
			Double x = distance-bond*cos(angle*DEGREE_TO_RADIAN);

			// move the x,y,z value back to the original coordinate system
			if (! labelX) {
				angleX = -angleX;
				rotationByX (y, z, angleX);
			}
			if (! labelY) {
				angleY = -angleY;
				rotationByY (x, z, angleY);
			}
			if (! labelZ) {
				angleZ = -angleZ;
				rotationByZ (x, y, angleZ);
			}
			x += x2;
			y += y2;
			z += z2;

			// write the result back
			coord[3*iAtom+0] = x;
			coord[3*iAtom+1] = y;
			coord[3*iAtom+2] = z;
		}
	}
}

void ZMatrix::print() const {
	cout << "zmatrix data " << endl;
	for (UInt i=0; i<getNAtoms(); i++) {
		zmatrixList[i].print();
	}
	cout << endl;
}


