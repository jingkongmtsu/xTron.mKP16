/**
 * CPP files corresponding to the parseshell.h
 * \author  Fenglai liu and Jing kong
 */
#include "excep.h"
#include "element.h"
#include "molecule.h"
#include "textread.h"
#include "basissetnames.h"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include "parseshell.h"
using namespace std;
using namespace boost;
using namespace excep;
using namespace element;
using namespace molecule;
using namespace textread;
using namespace basissetnames;
using namespace parseshell;

/***********************************************************************************
 *               utility function for the working classes here                     *
 ***********************************************************************************/
UInt parseshell::findBasisSetStyle(const string& input, const UInt& dataLoc)
{
	// set up the default basis set style
	UInt basisStyle = INVALID_BASIS_INPUT;             

	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "can not open the data file for parsing the basis set style";
		Excep excep("parseshell","findBasisSetStyle",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// firstly, move the pointer to the proper beginning
	// if the data location is > 0
	if (dataLoc > 0) {
		inf.seekg(dataLoc,ios::beg);
	}

	// this variable is used to define the hit of keywords
	// for detecting the nwchem style
	// finding basis, +1 
	// finding end,   +1 
	UInt nwchemStyleHitKeys = 0;

	// check the basis set style
	WordConvert w;
	string line;
	while(getline(inf,line)) {

		// we omit the empty lines and comment lines 
		LineParse l(line);
		if (l.isEmp() || l.isCom()) continue;

		// for the Gaussian style of input, the head of the data
		// section starts with atom symbol with an "0"
		// we will check it here
		if (l.getNPieces() == 2) {
			string k0 = l.findValue(0);
			string k1 = l.findValue(1);
			if (k1 == "0" && isAtomicSymbol(k0)) {
				basisStyle = GAUSSIAN_BASIS_INPUT; 
				break;
			}
		}

		// for the NWChem style of input, the head of the data
		// starts with keyword "basis", and it will close with 
		// the end symbol
		// here we will check whether we see the "basis" keyword 
		if (l.getNPieces() > 1 && ! l.isSec()) {
			string k0 = l.findValue(0);
			if (w.compare(k0,"basis")) {
				nwchemStyleHitKeys += 1;
			}
		}

		// for the NWChem style of input, the head of the data
		// starts with keyword "basis", and it will close with 
		// the end symbol
		// here we will check whether we see the "end" keyword 
		if (l.getNPieces() == 1 && ! l.isSec()) {
			string k0 = l.findValue(0);
			if (w.compare(k0,"end")) {
				nwchemStyleHitKeys += 1;
			}
		}

		// if this is NWChem style data, now we identified it
		// then exit the loop
		if (nwchemStyleHitKeys == 2) {
			basisStyle = NWCHEM_BASIS_INPUT; 
			break;
		}

		// if this is in the user input file, we should stop when
		// the basis set section ends
		if (dataLoc > 0 && l.isSec()) break;
	}

	// close the file
	inf.close();

	// now return the style
	return basisStyle; 
}

/***********************************************************************************
 *               functions related to the RawShell class                           *
 ***********************************************************************************/
void RawShell::splitSPShellData(const RawShell& spShell)
{
	// check that whether it's a SP shell
	if (spShell.shellType != "SP") {
		Excep excep("RawShell","splitSPShellData",EXCEPTION_INPUT_PARAMETER_INVALID,"the input rawshell is not a SP shell");
		handleExcep(excep);
	}

	// also check the dimension data etc.
	if (spShell.ngau != ngau) {
		Excep excep("RawShell","splitSPShellData",EXCEPTION_INPUT_PARAMETER_INVALID,
				"the input rawshell number of Gaussians does not equal the new rawshell data here");
		handleExcep(excep);
	}

	// finally check the shell type here
	if (shellType != "S" && shellType != "P") {
		Excep excep("RawShell","splitSPShellData",EXCEPTION_INPUT_PARAMETER_INVALID,
				"the new rawshell's shell type is not S or P shell, something wrong here?");
		handleExcep(excep);
	}

	// we only need to reset the exp and coefficient array
	// firstly set the exp array data
	for(UInt i=0;i<ngau;i++) {
		exp[i] = spShell.exp[i];
	}

	// coefficients we need to see whether it's S or P
	if (shellType == "S") {
		for(UInt i=0;i<ngau;i++) {
			coe[i] = spShell.coe[i];
		}
	}else{
		for(UInt i=0;i<ngau;i++) {
			coe[i] = spShell.coe[i+ngau];
		}
	}
}

/***********************************************************************************
 *               functions related to the RawAtomShell class                       *
 ***********************************************************************************/
void RawAtomShell::readRawAtomShellData(const string& input, const UInt& dataLoc) 
{
	// firstly, let's see whether this is the user input file
	// or the basis set data file
	string inputFile(input);
	if (dataLoc == 0) {

		// get the basis file env
		char* path = getenv("EMUL_BASIS_FILE_DIR");
		if (path==NULL) {
			Excep excep("RawAtomShell","readRawAtomShellData",EXCEPTION_DIR_MISSING, "EMUL_BASIS_FILE_DIR is not defined");
			handleExcep(excep);
		}

		// set up the setting file path
		string settingPath(path);
		string basisFile = settingPath + "/" + inputFile + ".bas";
		inputFile = basisFile;
	}

	// check the possible error
	UInt basisStyle = findBasisSetStyle(inputFile,dataLoc);
	if (basisStyle == INVALID_BASIS_INPUT) {
		string infor = "we can not identify the shell data style, basis set format are not supported.";
		Excep excep("RawAtomShell","readRawAtomShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
		handleExcep(excep);
	}

	// finally, depending on the basis set data style;
	// we will read in the raw atom shell data
	// we note that we have already checked the basis set style;
	// so we do not check it again here
	if (basisStyle == NWCHEM_BASIS_INPUT) {
		readNWChemShellData(inputFile,dataLoc);
	}else if (basisStyle == GAUSSIAN_BASIS_INPUT) {
		readGaussianShellData(inputFile,dataLoc);
	}
}

void RawAtomShell::readGaussianShellData(const string& input, const UInt& dataLoc)
{
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "can not open file for raw atom shell reading";
		Excep excep("RawAtomShell","readGaussianShellData",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// firstly, move the pointer to the proper beginning
	// if the data location is > 0
	// the location is related to the section part in 
	// the input file, for reading data from basis set
	// files the dataLoc is 0
	if (dataLoc > 0) {
		inf.seekg(dataLoc,ios::beg);
	}

	// atomic symbol
	string symbol= getAtomicSymbol(atomic);

	// proceed to the given location
	// now let's go to find the given atom basis set part
	WordConvert w;
	string line;
	bool find = false;
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.getNPieces() != 2) {
			continue;
		}else if (l.isEmp() || l.isCom()){
			continue;
		}else if (dataLoc > 0 && l.isSec()) {
			break;
		}else{
			string k0 = l.findValue(0);
			string k1 = l.findValue(1);
			if (k1 == "0" && w.compare(k0,symbol)) {
				find = true;
				break;
			}
		}
	}
	if (! find) {
		string infor = "For atom : " + symbol;
		Excep excep("RawAtomShell","readGaussianShellData",EXCEPTION_RAW_ATOM_SHELL_DATA_NOT_FOUND,infor);
		handleExcep(excep);
	}

	// process each shell data for the given atom
	while(getline(inf,line)) {

		// each atom's data should be saperated with "****" sign
		// but for the basis set files, it could end with empty line
		// for the general basis set input, it could end with "end" 
		LineParse l(line);
		if (l.isShlSign() || l.isEmp() || l.isSec()){
			break;
		}

		// firstly, let's check the shell type information
		// we list all possible shell types that we can find
		string shelltype = l.findValue(0);
		w.capitalize(shelltype);
		if (shelltype != "S" && shelltype != "SP" &&
				shelltype != "SPD" && shelltype != "P" &&	
				shelltype != "D" && shelltype != "F" &&	
				shelltype != "G" && shelltype != "H" &&	
				shelltype != "I") {
			string infor = "For atom : " + symbol + "shell type is: " + shelltype + ", can not identify the shell type";
			Excep excep("RawAtomShell","readGaussianShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
			handleExcep(excep);
		}

		// now read in number of gaussians
		UInt ngau = 0;
		if (!w.toUInt(l.findValue(1),ngau)) {
			string infor = "For atom : " + symbol; 
			Excep excep("RawAtomShell","readGaussianShellData",
					EXCEPTION_RAW_ATOM_SHELL_NGAU_INVALID,infor);
			handleExcep(excep);
		}

		// now read in scale factor
		Double scalef = ZERO;
		if (!w.toDouble(l.findValue(2),scalef)) {
			string infor = "For atom : " + symbol; 
			Excep excep("RawAtomShell","readGaussianShellData",
					EXCEPTION_RAW_ATOM_SHELL_SCALE_FAC_INVALID,infor);
			handleExcep(excep);
		}

		// let's form the shell type of raw shell
		// also we figure out the number of types in the shell data
		string type(shelltype);
		UInt nTypeGau = 1;
		if (shelltype == "SPD") type = "SP";
		if (shelltype == "SP" || shelltype == "SPD") nTypeGau = 2;

		// find the pure/Cartesian format
		// we need a saperate flag to signal the shell state,
		// spherical(in other words, pure) or Cartesian
		bool pureState = true;
		if(l.getNPieces() >= 4) {
			string x = l.findValue(3);
			if (w.compare(x,"C")) {
				pureState = false;
			}
		}

		// now put all information into raw shell
		RawShell rs(pureState,scalef,nTypeGau,ngau,type);

		// to read in each shell's coe and exp
		DoubleVec coe2(ngau,ZERO); // for storing the P coe in the SP shell data
		DoubleVec coe3(ngau,ZERO); // for storing the D coe in the SPD shell data
		for (UInt i=0; i<ngau; i++) {

			// now read in each gaussian data
			getline(inf,line);
			LineParse l(line);

			// let's check the data length
			if(l.getNPieces() < 2) {
				string infor = "the data fields for the shell is not valid, at least two data fields, coefficients and exp";
				Excep excep("RawAtomShell","readGaussianShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
				handleExcep(excep);
			}else if(w.compare(shelltype,"SP") && l.getNPieces() != 3) {
				string infor = "for SP shell the number of data fields is not valid, should be exp. coeff.(S) coeff.(P)";
				Excep excep("RawAtomShell","readGaussianShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
				handleExcep(excep);
			}else if(w.compare(shelltype,"SPD") && l.getNPieces() != 4) {
				string infor = "for SPD shell the data is not valid, should be exp. coeff.(S) coeff.(P) coeff.(D)";
				Excep excep("RawAtomShell","readGaussianShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
				handleExcep(excep);
			}

			// exponential factor
			Double e = ZERO;
			if (!w.toDouble(l.findValue(0),e)) {
				Excep excep("RawAtomShell","readGaussianShellData",
						EXCEPTION_RAW_ATOM_SHELL_EXP_DATA_INVALID,line);
				handleExcep(excep);
			}
			rs.exp[i] = e;

			// coefficient
			Double c = ZERO;
			if (!w.toDouble(l.findValue(1),c)) {
				Excep excep("RawAtomShell","readGaussianShellData",
						EXCEPTION_RAW_ATOM_SHELL_COEF_DATA_INVALID,line);
				handleExcep(excep);
			}
			rs.coe[i] = c;

			// read in additional coefficient for SP/SPD
			if (w.compare(shelltype,"SP") || w.compare(shelltype,"SPD")) {
				Double c = ZERO;
				if (!w.toDouble(l.findValue(2),c)) {
					Excep excep("RawAtomShell","readGaussianShellData",
							EXCEPTION_RAW_ATOM_SHELL_COEF_DATA_INVALID,line);
					handleExcep(excep);
				}
				coe2[i] = c;
			}

			// SPD
			if (w.compare(shelltype,"SPD")) {
				Double c = ZERO;
				if (!w.toDouble(l.findValue(3),c)) {
					Excep excep("RawAtomShell","readGaussianShellData",
							EXCEPTION_RAW_ATOM_SHELL_COEF_DATA_INVALID,line);
					handleExcep(excep);
				}
				coe3[i] = c;
			}
		}

		// additional treatment for the SP/SPD shell, attach the P shell data
		if (w.compare(shelltype,"SP") || w.compare(shelltype,"SPD")) {
			for(UInt i=0;i<ngau;i++) {
				rs.coe[ngau+i] = coe2[i];
			}
		}

		// finally appending this shell
		shells.push_back(rs);

		// here we do additional treatment for the SPD shell
		// for SPD shell, we just split it into the SP shell and another D shell
		if (shelltype == "SPD") {
			nTypeGau = 1;
			RawShell rawShellD(pureState,scalef,nTypeGau,ngau,type);
			for(UInt i=0;i<ngau;i++) {
				rawShellD.exp[i] = rs.exp[i];
				rawShellD.coe[i] = coe3[i];
			}
			shells.push_back(rawShellD);
		}
	}

	// close the file
	inf.close();
}

void RawAtomShell::readNWChemShellData(const string& input, const UInt& dataLoc)
{
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "can not open file for raw atom shell reading";
		Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// firstly, move the pointer to the proper beginning
	// if the data location is > 0
	// the location is related to the section part in 
	// the input file, for reading data from basis set
	// files the dataLoc is 0
	if (dataLoc > 0) {
		inf.seekg(dataLoc,ios::beg);
	}

	// atomic symbol
	string symbol= getAtomicSymbol(atomic);
	//printf("this is for atom: %s\n", symbol.c_str());

	// here we need to set up a location variable to store the 
	// location of the file
	UInt loc = 0;

	// for nwchem format shell data, the pure or Cartesian style of 
	// basis set is defined in the head of each atom shell section
	bool pureState = true;

	// proceed to the given location
	// now let's go to find the given atom basis set part
	WordConvert w;
	string line;
	bool find = false;
	while(getline(inf,line)) {

		// check whether it's comment or the empty lines
		LineParse l(line);
		if (l.isEmp() || l.isCom()){
			continue;
		}

		// if it's in the user input file, when section ends
		// then we stop too
		if (dataLoc>0 && l.isSec()) break;

		// whether it's the basis set section begins
		// store that possible location
		string k0 = l.findValue(0);
		if (w.compare(k0,"BASIS")) {

			// detect the possible shell status?
			// it's in the third data field
			// the default stat is Cartesian, follow the rules in the NWChem program
			string stat = "CARTESIAN";
			if (l.getNPieces() >= 3) {
				stat = l.findValue(2);
				w.capitalize(stat);
			}

			// save a copy of location
			loc = inf.tellg();

			// now let's see the next line
			string line2;
			getline(inf,line2);
			LineParse l2(line2);
			string atom = l2.findValue(0);
			if (w.compare(atom,symbol)) {

				// now for the found atom, write the pure status back
				if (stat == "SPHERICAL") {
					pureState = true;
				}else if (stat == "CARTESIAN") {
					pureState = false;
				}else{
					string infor = "the pure/Cartesian status here in the shell format is invalid: " + stat;
					Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
					handleExcep(excep);
				}

				// now let's get out
				find = true;
				break;
			}
		}
	}
	if (! find) {
		string infor = "For atom : " + symbol;
		Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_RAW_ATOM_SHELL_DATA_NOT_FOUND,infor);
		handleExcep(excep);
	}

	// now let's move the data to the atom shell data section
	// it's position is recorded by the loc
	inf.seekg(loc,ios::beg);

	// process each shell data for the given atom
	while(getline(inf,line)) {

		// the end of shell data could be end with empty line, 
		// section name; or the end symbol
		LineParse l(line);
		if (l.isEmp() || l.isSec()){
			break;
		}

		// check whether it's end with "end" symbol
		if (l.getNPieces() == 1) {
			string stat = l.findValue(0);
			if (w.compare(stat,"END")) {
				break;
			}
		}

		// firstly, let's check the shell type information
		// we list all possible shell types that we can find
		// it seems that nwchem does not support the SPD type of shell
		// therefore we will issue an error if we have SPD shell defined here
		string shelltype = l.findValue(1);
		w.capitalize(shelltype);
		if (shelltype != "S" && shelltype != "SP" && shelltype != "P" &&	
				shelltype != "D" && shelltype != "F" &&	
				shelltype != "G" && shelltype != "H" && shelltype != "I") {
			string infor = "For atom : " + symbol + "shell type is: " + shelltype + ", can not identify this shell type";
			Excep excep("RawAtomShell","readNWChemShellData", EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
			handleExcep(excep);
		}

		// debug code
		//if (symbol == "C") {
		//	printf("the data file lines %s\n", line.c_str());
		//}

		// we need to save a location
		UInt location = inf.tellg();

		// now read in number of gaussians
		// also because NWChem support the general format shell data
		// therefore it may contains multiple shells
		//  
		// for SP shells, it does not support general format of shell data
		// so for SP shells, it's only for one shell section
		UInt ngau    = 0;
		UInt nShells = 1;
		string line2;
		while(getline(inf,line2)) {
			LineParse l2(line2);

			// this is debug code
			//if (symbol == "C") {
			//	printf("the data file lines for counting ngau: %s\n", line2.c_str());
			//}

			// we need to test that whether this data is double
			// if not, we will stop
			string data = l2.findValue(0);
			Double x = ZERO;
			if (! w.toDouble(data,x)) {
				break;
			}

			// determine the number of shell data
			ngau++;
			if (shelltype != "SP") {
				UInt ns = l2.getNPieces() - 1;
				if (ns != nShells) nShells = ns;
			}
		}

		// we should have at least one Gaussian function inside
		if (ngau == 0) {
			string infor = "why the number of Gaussian primitive functions is 0?";
			Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
			handleExcep(excep);
		}

		// set up the number of Gaussians type
		UInt nTypeGau = 1;
		if (shelltype == "SP") nTypeGau = 2;

		// now let's move back to the data beginning part
		// to read in shell's coe and exp
		inf.seekg(location,ios::beg);
		DoubleVec expFac(ngau,ZERO); 
		DoubleVec coe(nTypeGau*ngau*nShells,ZERO); 
		for (UInt i=0; i<ngau; i++) {

			// now read in each gaussian data
			getline(inf,line2);
			LineParse l2(line2);

			// this is debug code
			//if (shelltype == "SP") {
			//	printf("the data file lines for reading the data: %s\n", line2.c_str());
			//}

			// let's check the data length
			if(l2.getNPieces() -1 !=  nShells && shelltype != "SP") {
				cout << line2 << endl;
				cout << " the number of data fields " << l2.getNPieces() << endl;
				cout << " number of shells          " << nShells << endl;
				Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,
						"the number of data fileds does not equal to nShells");
				handleExcep(excep);
			}
			if(w.compare(shelltype,"SP") && l2.getNPieces() != 3) {
				cout << line2 << endl;
				Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_RAW_ATOM_SHELL_INVALID,
						"the number of data fields for the SP shell must be 3 (two coefficients + one exp)");
				handleExcep(excep);
			}

			// exponential factor
			Double e = ZERO;
			if (!w.toDouble(l2.findValue(0),e)) {
				Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_RAW_ATOM_SHELL_EXP_DATA_INVALID,line2);
				handleExcep(excep);
			}
			expFac[i] = e;

			// coefficient array
			for (UInt j=1; j<l2.getNPieces(); j++) {
				Double c = ZERO;
				if (!w.toDouble(l2.findValue(j),c)) {
					Excep excep("RawAtomShell","readNWChemShellData",EXCEPTION_RAW_ATOM_SHELL_COEF_DATA_INVALID,line2);
					handleExcep(excep);
				}
				coe[i+(j-1)*ngau] = c;
			}
		}

		// now put them into raw shell
		for (UInt j=0; j<nShells; j++) {

			// set up the shell data
			RawShell rs(pureState,ONE,nTypeGau,ngau,shelltype);
			for(UInt ip=0;ip<ngau;ip++) {
				rs.exp[ip] = expFac[ip];
				rs.coe[ip] = coe[ip+j*ngau];
			}

			// for SP shell, we need to add another coeff array
			// for this case, the array of coe must be 2*ngau
			if (shelltype == "SP") {
				for(UInt ip=0;ip<ngau;ip++) {
					rs.coe[ip+ngau] = coe[ip+ngau];
				}
			}

			// finally push the data inside
			shells.push_back(rs);
		}
	}

	// close the file
	inf.close();
}

/***********************************************************************************
 *               functions related to the RawMolShell class                        *
 ***********************************************************************************/
UInt RawMolShell::parsingLoc(const string& input, UInt section, UInt part) const
{
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "section number is ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("RawMolShell","parsingLoc",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// get the key word
	string key = "molecule";
	if (section == 0) key = "cluster";

	// find the keyword of molecule section
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

	// checking the finding
	if (! find) {
		string infor = "section number is ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("RawMolShell","parsingLoc",EXCEPTION_SECTION_SEARCH_FAILED,infor);
		handleExcep(excep);
	}

	// to see whether this is aux basis or common basis?
	key = "basis";
	if (part > 0) {
		key = "aux_basis";
	}

	// furthermore, we search the basis set part appending to the geom definition 
	// that is to say, basis set definition should be like this:
	// molecule(cluster)...
	// basis set definition
	// another molecule...
	n = 0;
	find = false;
	UInt pos = -1;
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isSec() && w.compare(l.getSecName(), key)) {
			if (part == 0) {
				find = true;
				pos = inf.tellg();
				break;
			}else{
				n++;
				if (part == n) {
					find = true;
					pos = inf.tellg();
					break;
				}
			}
		}

		// in case we searched another molecule/cluster section
		bool isMolSec = l.isSec() && w.compare(l.getSecName(), "molecule");
		bool isCulSec = l.isSec() && w.compare(l.getSecName(), "cluster");
		if (isMolSec || isCulSec) {
			string infor = "searching molecule section " + boost::lexical_cast<string>(section); 
			Excep excep("RawMolShell","parsingLoc",EXCEPTION_SHELL_DATA_FAIL_TO_BE_FOUND,infor);
			handleExcep(excep);
		}
	}

	// further checking
	if (! find) {
		if (part>0) {
			string infor = "We did not fine the aux basis set section: " +  
				boost::lexical_cast<string>(section) + " part: " + 
				boost::lexical_cast<string>(part); 
			Excep excep("RawMolShell","parsingLoc",EXCEPTION_SHELL_DATA_FAIL_TO_BE_FOUND,infor);
			handleExcep(excep);
		}else{
			string infor = "We did not fine the normal basis set section: " +  
				boost::lexical_cast<string>(section); 
			Excep excep("RawMolShell","parsingLoc",EXCEPTION_SHELL_DATA_FAIL_TO_BE_FOUND,infor);
			handleExcep(excep);
		}
	}

	// close input file
	inf.close();

	// return the value
	return pos;
}

void RawMolShell::parsingBasisSetInfor(const string& input, const UInt& loc, 
		vector<string>& basisSetNames, UIntVec& elementNames) const
{
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string infor = "can not open the user input file";
		Excep excep("RawMolShell","parsingBasisSetInfor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// let's read the basis set data information
	inf.seekg(loc,ios::beg);
	string line;
	WordConvert w;
	while(getline(inf,line)) {

		// now let's read in the basis set string
		// we omit the empty lines and comment lines 
		// and we break when the section ends
		LineParse l(line);
		if (l.isEmp() || l.isCom()) continue;
		if (l.isSec()) break;

		//
		// let's try to read in the possible basis set data
		// two formats supported, the simplest one is that only
		// basis set data provided for all of atoms, like:
		// %basis
		// 6-31G**
		// %end
		//
		// or like:
		// %basis
		// H  6-31G**
		// O  cc-pvdz
		// %end
		//
		if (l.getNPieces() == 2) {
			string s  = l.findValue(0);
			string basisSetName = l.findValue(1);

			// firstly let's check whether it's atom symbol
			// or atomic number
			UInt atomic = 0;
			if (isAtomicSymbol(s)) {
				atomic = getAtomicNumber(s);
			}else if (w.isUInt(s)) {
				UInt v = 0;
				w.toUInt(s,v);
				if (v > 0 && v<MAX_ATOMS_TYPE) {
					atomic = v;
				}
			}

			// if the basis set name or the atomic number is not valid,
			// we need to issue an error
			if (atomic == 0 || ! isValidBasisSetName(basisSetName)) {
				cout << line << endl;
				string infor = "we can not process the basis set data name or the atomic symbol/number "
					"you provided in the above line"; 
				Excep excep("RawMolShell","parsingBasisSetInfor",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
				handleExcep(excep);
			}

			// now let's push the data inside
			elementNames.push_back(atomic);
			basisSetNames.push_back(basisSetName);

		}else if (l.getNPieces() == 1) {

			// for this situation, we should only have one line of information
			string s  = l.findValue(0);
			if (isValidBasisSetName(s)) {
				basisSetNames.push_back(s);
			}else{
				cout << line << endl;
				string infor = "we can not process the basis set data name you provided in the above line"; 
				Excep excep("RawMolShell","parsingBasisSetInfor",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
				handleExcep(excep);
			}
			break;
		}
	}

	// close input file
	inf.close();
}

RawMolShell::RawMolShell(const string& input, const Molecule& mol, 
		UInt part, bool debug_print)
{
	// get the data section
	UInt location = parsingLoc(input,mol.getSec(),part);

	// let's prepare the atom types data for reading the raw shell data
	UInt nAtomTypes = mol.getNAtomTypes();
	UIntVec elements(nAtomTypes);
	mol.getAtomTypes(elements);

	// let's try to identify whether the basis set is contained in
	// the input file
	// if the basis set data is defined in the standard basis set file,
	// here the basis style will be "INVALID_BASIS_INPUT"
	UInt basisStyle = findBasisSetStyle(input,location);
	if (basisStyle != INVALID_BASIS_INPUT) {
		for(UInt i=0; i<elements.size(); i++) {
			UInt atomic  = elements[i];
			if (isGhostAtom(atomic)) continue;
			RawAtomShell atomShell(input,location,atomic);
			molShells.push_back(atomShell);
		}
	}else{

		// now read possible basis set names if it's not user defined input file
		// we note that the sequence of elements in the elementNames may not
		// same with the array of elements, so that's the reason we have two arrays
		vector<string> basisSetNames;
		UIntVec elementNames;
		basisSetNames.reserve(nAtomTypes);
		elementNames.reserve(nAtomTypes);
		parsingBasisSetInfor(input,location,basisSetNames,elementNames);

		// firstly, if the basis set data is defined per each element, let's check
		// whether the element names cover all of possible elements in the molecule
		if (elementNames.size() > 0) {
			for(UInt i=0; i<elementNames.size(); i++) {
				UInt ele  = elementNames[i];
				bool find = false;
				for(UInt j=0; j<elements.size(); j++) {
					if (ele == elements[j]) {
						find = true;
						break;
					}
				}
				if (!find) {
					string infor = "did you forget the define the basis set data for this element: " 
						+ boost::lexical_cast<string>(ele);
					Excep excep("RawMolShell","constructor",EXCEPTION_RAW_ATOM_SHELL_INVALID,infor);
					handleExcep(excep);
				}
			}
		}

		// now let's read in the basis set data
		molShells.reserve(elements.size());
		for(UInt i=0; i<elements.size(); i++) {

			// atom information
			UInt atomic  = elements[i];
			UInt loc     = 0;  // since we read basis set data file, the location is set to 0
			if (isGhostAtom(atomic)) continue;

			// get the basis set data name
			// if the elementNames is larger than 0; we need to get the corresponding
			// element's basis set data name
			string basisName = basisSetNames[0];
			UInt pos = 0;
			if (elementNames.size() > 0) {
				for(UInt j=0; j<elementNames.size(); j++) {
					if (elementNames[j] == atomic) {
						pos = j;
						break;
					}
				}
			}
			basisName = basisSetNames[pos];

			// get the basis set file name
			string basisFileName = getBasisSetFileNameFromInput(basisName);

			// finally form the atom shell data
			RawAtomShell atomShell(basisFileName,loc,atomic);
			molShells.push_back(atomShell);
		}
	}

	// debug printing
#ifdef DEBUG
	if (debug_print) print();
#endif
}

const RawAtomShell& RawMolShell::getAtomShell(UInt atomic) const
{
	// search for the atom type information
	for(UInt i=0; i<getNAtomShell(); i++) {
		const RawAtomShell& ras = molShells[i];
		if (ras.atomic == atomic) return ras;
	}

	// now we did not find the matched atomic number
	// throw an error
	string infor = "for atom whose atomic number is: " + boost::lexical_cast<string>(atomic); 
	Excep excep("RawMolShell","getAtomShell", 
			EXCEPTION_RAW_ATOM_SHELL_DATA_NOT_FOUND,infor);
	handleExcep(excep);
	return molShells[0];
}

bool RawMolShell::hasHighAngShellData() const
{
	// loop over the shells data to see whether 
	// we have shells whose angular momentum > F
	for(UInt i=0; i<getNAtomShell(); i++) {
		const RawAtomShell& ras = molShells[i];
		for(UInt iShell=0; iShell<ras.getNShell(); iShell++) {
			const RawShell& s = ras.getShell(iShell);
			if (s.shellType == "G" || s.shellType == "H" || s.shellType == "I") {
				return true;
			}
		}
	}

	// now let's return false
	return false;
}

bool RawMolShell::hasSPShellData() const
{
	// loop over the shells data to see whether 
	// we have SP shells 
	for(UInt i=0; i<getNAtomShell(); i++) {
		const RawAtomShell& ras = molShells[i];
		for(UInt iShell=0; iShell<ras.getNShell(); iShell++) {
			const RawShell& s = ras.getShell(iShell);
			if (s.shellType == "SP") {
				return true;
			}
		}
	}

	// now let's return false
	return false;
}

void RawMolShell::print() const 
{
	cout << "***************************************" << endl;
	cout << "*           RawMolShell               *" << endl;
	cout << "***************************************" << endl;
	UInt length = molShells.size();
	for (UInt i=0;i<length;i++) {
		const RawAtomShell& atomShell = molShells[i];
		string sym = getAtomicSymbol(atomShell.atomic);
		printf ("%3s  \n", sym.c_str());
		UInt len = atomShell.shells.size();
		for (UInt j=0;j<len;j++) {
			const RawShell& sh = atomShell.shells[j];
			const char* usePure = "PURE";
			if (!sh.isPure) usePure = "CART";
			printf ("%3s  %4d  %6.4f  %4s  \n", sh.shellType.c_str(), 
					(Int)sh.ngau, sh.scale, usePure);
			for (UInt k=0;k<sh.ngau;k++) {
				if (sh.shellType == "SP"){
					printf ("%14.7f  %14.7f  %14.7f  \n", sh.exp[k], sh.coe[k],sh.coe[sh.ngau+k]);
				}else{
					printf ("%14.7f  %14.7f  \n", sh.exp[k], sh.coe[k]);
				}
			}
		}
		printf ("%4s  \n", "****");
	}	
}


