/**
 * \file molecule.cpp
 * \brief  CPP files corresponding to the molecule.h
 * \author Fenglai Liu and Jing Kong
 */
#include<iostream>
#include<fstream>
#include<cstdio>
#include <boost/lexical_cast.hpp>
#include "element.h"
#include "molecule.h"
#include "excep.h"
#include "textread.h"
#include "zmatrix.h"
using namespace element;
using namespace excep;
using namespace molecule;
using namespace zmatrix;
using namespace textread;

///////////////////////////////////////////////////
// ####  private function for molecule class
///////////////////////////////////////////////////

UInt Molecule::getSecLoc(ifstream& inf) {

	// get the key word
	string key = "molecule";
	if (section == 0) key = "cluster";

	// additional check on section
	if (section>=MAX_MOLECULE_SECTION_NUMBER) {
		string infor = "given section is: " + boost::lexical_cast<string>(section); 
		Excep excep("Molecule","getSecLoc",EXCEPTION_ILLEGAL_GEOM_SECTION_INDEX,infor);
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
				// let's check the data unit
				for(UInt iPiece=0; iPiece<l.getNPieces(); iPiece++) {
					string choice = l.findValue(iPiece);
					if (w.compare(choice,"bohr")) {
						unitInAngstrom = false;
					}else if (w.compare(choice,"angstrom"))  {
						unitInAngstrom = true;
					}
				}

				find = true;
				break;
			}else{
				n++;
				if (section == n) {
					// let's check the data unit
					for(UInt iPiece=0; iPiece<l.getNPieces(); iPiece++) {
						string choice = l.findValue(iPiece);
						if (w.compare(choice,"bohr")) {
							unitInAngstrom = false;
						}else if (w.compare(choice,"angstrom"))  {
							unitInAngstrom = true;
						}
					}

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
		Excep excep("Molecule","getSecLoc",EXCEPTION_SECTION_SEARCH_FAILED,infor);
		handleExcep(excep);
	}

	// remember where we are
	UInt loc = inf.tellg(); 
	return loc;
}

void Molecule::readChargeMult(ifstream& inf, const UInt& loc) 
{
	// go to the given location
	inf.seekg(loc,ios::beg);

	// directly read in data, it should be just under
	// the title of molecule if it's defined
	string line;
	getline(inf,line);
	LineParse l(line);
	WordConvert w;
	UInt number = l.getNPieces();

	// double check and read in data
	if(number != 2) {
		string infor = "Something wrong with the line of charge and spin multilicity";
		Excep excep("Molecule","readChargeMult",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}
	
	// charge and multiplicity
	bool convert = w.toInt(l.findValue(0), charge);
	if (! convert) {
		string infor = "Charge can not be read";
		Excep excep("Molecule","readChargeMult",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
		handleExcep(excep);
	}
	convert = w.toUInt(l.findValue(1), multiplicity);
	if (! convert) {
		string infor = "multiplicity can not be read";
		Excep excep("Molecule","readChargeMult",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
		handleExcep(excep);
	}
}

void Molecule::readAtomData(const string& input, const UInt& loc)
{
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "section number is ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("Molecule","readAtomData",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// go to the given location
	// we need to bypass the first line of data
	// since it's spin multiplicity and charge
	inf.seekg(loc,ios::beg);
	string line;
	getline(inf,line);

	// let's see whether they are the xyz format data?
	// we check the first line of data
	bool isXYZ = true;
	getline(inf,line);
	LineParse l(line);
	if(l.getNPieces() < 4) {
		// for xyz format, it's at least 
		// "atom x y z" 
		// so must be at least 4 data fields
		isXYZ = false;
	}

	// read in atom data
	if (isXYZ) {
		readAtomDataInXYZ(inf,loc);
	}

	// now close the file
	inf.close();

	// if this is not xyz, this is zmatrix format
	if (! isXYZ) {
		readAtomDataInZMatrix(input);
	}
}

void Molecule::readAtomDataInZMatrix(const string& input) 
{
	// form the zmatrix data 
	ZMatrix zmatrix(input,section);

	// now let's transform the data into the xyz
	DoubleVec xyz(3*zmatrix.getNAtoms());
	zmatrix.loadInToXYZ(xyz);

	// now let's push the geometry data into the data set
	Double coord[3];
	for(UInt i=0; i<zmatrix.getNAtoms(); i++) {
		const ZMatrixItem& item = zmatrix.getZMatrix(i);
		
		// get the information
		UInt atomic = item.getAtomic();
		Double x = xyz[3*i+0];
		Double y = xyz[3*i+1];
		Double z = xyz[3*i+2];

		// scale the result to bohr
		if (isInputGeomDataInAng()) {
			coord[0] = x*ANGSTROM_TO_BOHR;
			coord[1] = y*ANGSTROM_TO_BOHR;
			coord[2] = z*ANGSTROM_TO_BOHR;
		}

		// add in new data
		Atom newAtom(coord,atomic);
		atoms.push_back(newAtom);
	}
}

void Molecule::readAtomDataInXYZ(ifstream& inf, const UInt& loc) 
{
	// go to the given location
	// we need to bypass the first line of data
	inf.seekg(loc,ios::beg);
	string line;
	getline(inf,line);

	// read the data
	WordConvert w;
	Double coord[3];
	UInt atomic = 0;
	while(getline(inf,line)) {

		// begin the real parsing
		LineParse l(line);

		// where we stop reading?
		if (l.isEmp() || l.isCom() || l.isSec()) break;

		// checking the data field range, should be at least 4
		if(l.getNPieces() < 4) {
			string infor = line;
			Excep excep("Molecule","readAtomDataInXYZ",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}

		// read in the atomic number, here the read in string may be
		// the atomic number or the atom symbol
		string s = l.findValue(0);
		if (isAtomicSymbol(s)) {
			atomic = getAtomicNumber(s);
		}else{
			if (w.isUInt(s)) {
				w.toUInt(s,atomic);
				if (atomic>MAX_ATOMS_TYPE) {
					string infor = line;
					Excep excep("Molecule","readAtomDataInXYZ",INVALID_ATOMIC_NUMBER,infor);
					handleExcep(excep);
				}
			} else {
				string infor = line;
				Excep excep("Molecule","readAtomDataInXYZ",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
		} 

		// read in the coordinates
		if (!w.toDouble(l.findValue(1), coord[0])) {
			string infor = line;
			Excep excep("Molecule","readAtomDataInXYZ",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
			handleExcep(excep);
		}
		if (!w.toDouble(l.findValue(2), coord[1])) {
			string infor = line;
			Excep excep("Molecule","readAtomDataInXYZ",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
			handleExcep(excep);
		}
		if (!w.toDouble(l.findValue(3), coord[2])) {
			string infor = line;
			Excep excep("Molecule","readAtomDataInXYZ",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
			handleExcep(excep);
		}

		// scale the result to bohr
		if (isInputGeomDataInAng()) {
			coord[0] *= ANGSTROM_TO_BOHR;
			coord[1] *= ANGSTROM_TO_BOHR;
			coord[2] *= ANGSTROM_TO_BOHR;
		}

		// add in new data
		Atom newAtom(coord,atomic);
		atoms.push_back(newAtom);
	}
}

void Molecule::setNAlphaBeta() 
{
	// get the total electrons
	Int tote = 0;
	for (UInt i=0; i<getNAtoms(); i++) tote += getAtomic(i); 
	tote -= charge;

	// here we need to check the multiplicity provided
	Int diff = multiplicity - 1; // difference between nalpha - nbeta
	if ((tote + diff)%2 != 0 || diff < 0) {
		string infor = "Multiplicity conflicts with the given molecular system";
		Excep excep("Molecule","setNAlphaBeta",EXCEPTION_INVALID_NALPHA_NBETA,infor);
		handleExcep(excep);
	} 
	if ((tote - diff)%2 != 0 || tote - diff < 0) {
		string infor = "Multiplicity conflicts with the given molecular system";
		Excep excep("Molecule","setNAlphaBeta",EXCEPTION_INVALID_NALPHA_NBETA,infor);
		handleExcep(excep);
	}
	
	// set the alpha and beta
	nalpha = 0.5*(tote+diff);
	nbeta  = 0.5*(tote-diff);
}

void Molecule::updateAtomIndexInfor()
{
	UInt index = 0;
	for(UInt i=0; i<getNAtoms(); i++) {
		atoms[i].updateAtomIndex(index);
		index++;
	}
}

bool Molecule::checkDistance() const
{
	// set an constant for checking the distance
	// 0.00001 should be good enough
	Double smallR = 0.00001E0;
	for(UInt i=0; i<getNAtoms(); i++) {
		for(UInt j=i; j<getNAtoms(); j++) {
			if (j == i) continue;
			Double R = getDistance(i,j);
			if (R<smallR) return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////
// ####  constructor or constructor like functions
///////////////////////////////////////////////////
Molecule::Molecule(const string& input, UInt whichGeom):unrestricted(false),unitInAngstrom(true),
	section(whichGeom),multiplicity(1),nalpha(0),nbeta(0),charge(0)
{
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "section number is ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("Molecule","Constructor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// now read in molecule data
	// get the section location
	UInt location = getSecLoc(inf);

	// read the charge and spin multiplicity data
	readChargeMult(inf,location); 

	// close the input file
	inf.close();

	// get the atom data
	// we will open the file inside again
	readAtomData(input,location);

	// now determine the nalpha and nbeta 
	setNAlphaBeta();

	// reform the atom index
	updateAtomIndexInfor();

	// finally check the distance between atoms
	bool everythingGood = checkDistance();
	if (! everythingGood) {
		string infor = "section number is ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("Molecule","Constructor",EXCEPTION_ATOMS_POS_OVERLAP,infor);
		handleExcep(excep);
	}
}

Molecule::Molecule(const UInt& atomic):unrestricted(false),unitInAngstrom(true),
	section(-1),multiplicity(1),nalpha(0),nbeta(0),charge(0)
{
	multiplicity = getGroundStateAtomMult(atomic);
	Double xyz[3];
	xyz[0] = ZERO;
	xyz[1] = ZERO;
	xyz[2] = ZERO;
	Atom atom(xyz,atomic);
	atom.updateAtomIndex(0);
	atoms.reserve(1);
	atoms.push_back(atom);
	setNAlphaBeta();
}

///////////////////////////////////////////////////
// ####  public functions
///////////////////////////////////////////////////
void Molecule::getAtomTypes(UIntVec& atomTypes) const {

	// obtain the types information
	UInt possibleAtomTypes[MAX_ATOMS_TYPE];
	for(UInt i=0; i<MAX_ATOMS_TYPE; i++) possibleAtomTypes[i] = 0;
	for(UInt i=0; i<getNAtoms(); i++) {
		UInt atomic = getAtomic(i);
		if (possibleAtomTypes[atomic] > 0) {
			continue;
		}else{
			possibleAtomTypes[atomic] = 1;
		}
	}

	// count the number of types
	UInt count = 0;
	for(UInt i=0; i<MAX_ATOMS_TYPE; i++) {
		if (possibleAtomTypes[i] > 0) count++;
	}
	atomTypes.assign(count,0);

	// analyze the result
	count = 0;
	for(UInt i=0; i<MAX_ATOMS_TYPE; i++) {
		if (possibleAtomTypes[i] > 0) {
			atomTypes[count] = i; 
			count++;
		}
	}
}

UInt Molecule::getNAtomTypes() const {

	// obtain the types information
	UInt count = 0;
	UInt possibleAtomTypes[MAX_ATOMS_TYPE];
	for(UInt i=0; i<MAX_ATOMS_TYPE; i++) possibleAtomTypes[i] = 0;
	for(UInt i=0; i<getNAtoms(); i++) {
		UInt atomic = getAtomic(i);
		if (possibleAtomTypes[atomic] > 0) {
			continue;
		}else{
			possibleAtomTypes[atomic] = 1;
			count++;
		}
	}
	return count;
}

const Double* Molecule::getXYZ(UInt atom) const {
	const Atom& a = getAtom(atom);
	return a.getXYZ();
}

UInt Molecule::getAtomic(UInt atom) const {
	const Atom& a = getAtom(atom);
	return a.getAtomic();
}

string Molecule::getSymbol(UInt atom) const {
	UInt atomicNum = getAtomic(atom);
	return getAtomicSymbol(atomicNum);
}

void Molecule::getAtomDistanceMatrix(Double* dmatrix, const UInt& order) const{

	if (order == 0) {
		UInt natoms = getNAtoms();
		for(UInt i=0; i<natoms; i++) {
			for(UInt j=i; j<natoms; j++) {
				const Double* xyz1 = getXYZ(i); 
				const Double* xyz2 = getXYZ(j); 
				Double x2 = (xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]); 
				Double y2 = (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]); 
				Double z2 = (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]); 
				Double r = sqrt(x2+y2+z2);
				dmatrix[j+i*natoms] = r;
				if (i != j) dmatrix[i+j*natoms] = r;
			}
		}
	}else {
		string infor = "This order is not available yet";
		Excep excep("Molecule","getAtomDistanceMatrix",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
		handleExcep(excep);
	}
}

Double Molecule::getDistance(const UInt& atom1, const UInt& atom2) const {
	const Double* xyz1 = getXYZ(atom1); 
	const Double* xyz2 = getXYZ(atom2); 
	Double x2 = (xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]); 
	Double y2 = (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]); 
	Double z2 = (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]); 
	Double r = sqrt(x2+y2+z2);
	return r;
}

Double Molecule::getNuclearRepulsion() const {
	Double E = ZERO;
	for(UInt i=0; i<getNAtoms(); i++) {
		for(UInt j=i; j<getNAtoms(); j++) {
			if (j == i) continue;
			UInt ZA = getAtomic(i);
			UInt ZB = getAtomic(j);
			Double R = getDistance(i,j);
			E += ZA*ZB/R;
		}
	}
	return E;
}

void Molecule::forceUNRestricted()
{
	// set the unrestricted to be true
	unrestricted = true;

	// we need to check whether we can really do it
	// for the case that nalpha != nbeta, it's natrually
	// unrestricted. However, for the single electron
	// case, you can not do it
	if (isSingEleSystem()) {
		string infor = "for single electron system you can not force the unrestricted state"; 
		Excep excep("Molecule","forceUNRestricted",EXCEPTION_ILLEGAL_SET_UNRESTRICTED,infor);
		handleExcep(excep);
	}

}

void Molecule::print() const {

	cout << "Molecule information for section " << section << endl;
	cout << "charge         : "        << charge  << endl;	
	cout << "multiplicity   : "        << multiplicity  << endl;	
	cout << "nelectron      : "        << getTotalEle()  << endl;	
	cout << "nalpha         : "        << nalpha << endl;	
	cout << "nbeta          : "        << nbeta  << endl;	
	if (unrestricted) {
		cout << "the molecule is forced to the unrestricted " << endl;	
	}
	if (unitInAngstrom) {
		cout << "the input geometry data is in unit of angstrom " << endl;	
	}else{
		cout << "the input geometry data is in unit of Bohr " << endl;	
	}
	for (UInt i=0; i<getNAtoms(); i++) {
		const Atom& a = getAtom(i);
		Int atomic    = a.getAtomic();
		Int index     = a.getAtomIndex();
		const Double* xyz = a.getXYZ();
		string sym  = getAtomicSymbol(atomic); 
		printf("%-3s  %-6d  %-3d  %-10.5f  %-10.5f  %-10.5f\n", sym.c_str(), index, atomic,
				xyz[0]*BOHR_TO_ANGSTROM,xyz[1]*BOHR_TO_ANGSTROM,xyz[2]*BOHR_TO_ANGSTROM);
	}
}


