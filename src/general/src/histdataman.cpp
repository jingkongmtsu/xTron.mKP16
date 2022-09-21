/**
 * cpp file for the historical data manipulation class
 * \author fenglai liu 
 */
#include<iostream>
#include<fstream>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>  
#include "excep.h"
#include "filerw.h"
#include "matrix.h"
#include "textread.h"
#include "globalinfor.h"
#include "histdataman.h"
using namespace excep;
using namespace filerw;
using namespace matrix;
using namespace textread;
using namespace boost;
using namespace boost::filesystem;
using namespace globalinfor;
using namespace histdataman;

string HistDataMan::setLocation(const GlobalInfor& infor, const string& name, 
		const UInt& section) const
{
	// form the data folder name
	string dir(infor.getScratchDir());
	string sec = boost::lexical_cast<string>(section);
	string moduleName(name);
	to_lower(moduleName);

	// now use path to connect all of string together
	path p0(dir.c_str());
	path p1(sec.c_str());
	path p2(moduleName.c_str());
	p0 /= p1;
	p0 /= p2;
	//cout << "location name " << p0.string() << endl;
	return p0.string();
}

HistDataMan::HistDataMan(const GlobalInfor& infor, const string& name, const UInt& sec, 
		const UInt& nSpins, bool inFile):nSpin(nSpins),location("NONE")
{
	// initilize the data when using file/memory
	if (inFile) {

		// take a look that whether there exists the same name 
		// file or folder. If has, destroy it; else we create a new one
		// all of these done in silence
		location = setLocation(infor,name,sec);
		path data(location.c_str());

		// does it exist?
		if (exists(data)) {
			remove_all(data);
		}

		// now create a new one
		// we note that if the scratch dir and geometry folder
		// is not created, the dir creation will create them first
		create_directories(data);
	}

	// roughly initilize the status vector
	// this is used only for reserving some space
	UInt n = 100;
	status.reserve(N_FIELD*n);
}

HistDataMan::HistDataMan(const GlobalInfor& infor, const string& name, 
		const UInt& sec):nSpin(0),location("NONE")
{
	// initilize the data folder location
	location = setLocation(infor,name,sec);
	path folder(location.c_str());

	// does it exist?
	// we may need to report an error here
	if (! exists(folder)) {
		string infor = "the given location does not exist: " + location; 
		Excep excep("HistDataMan","constructor",EXCEPTION_DIR_MISSING,infor);
		handleExcep(excep);
	}

	// now let's see whether we have the record.txt?
	path record("record.txt");
	folder /= record;
	if (! exists(folder)) {
		string infor = "record.txt in the given location does not exist: " + location; 
		Excep excep("HistDataMan","constructor",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}
	string file = folder.string();

	// let's recover the data, firstly open it
	std::ifstream inf;
	const char * input_file = file.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		string info = "can not open the record.txt for location: " + location;
		Excep excep("HistDataMan","constructor",EXCEPTION_FILE_MISSING,info);
		handleExcep(excep);
	}

	// recover the data, get the nSpin state
	string line;
	WordConvert w;
	inf.seekg(0,ios::beg);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isCom() || l.isEmp()) continue;
		if (l.getNPieces() == 2) {
			string key  = l.findValue(0);
			if (w.compare(key, "nspin")) {
				string value = l.findValue(1);
				if (value == "1") {
					nSpin = 1;
				}else if (value == "2") {
					nSpin = 2;
				}else{
					string info = "nSpin is invalid in the record.txt for location: " + location;
					Excep excep("HistDataMan","constructor",EXCEPTION_HISTDATAMAN_ERROR,info);
					handleExcep(excep);
				}
			}
		}
	}

	// roughly initilize the status vector
	// this is used only for reserving some space
	UInt n = 100;
	status.reserve(N_FIELD*n);

	// now rewind the file let's get the data
	inf.seekg(0,ios::beg);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isCom() || l.isEmp()) continue;
		string key  = l.findValue(0);
		if (w.compare(key, "data_section")) {

			// check the number of fields
			if (l.getNPieces() != N_FIELD+1) {
				string info = "number of data is invalid in the record.txt for location: " + location;
				Excep excep("HistDataMan","constructor",EXCEPTION_HISTDATAMAN_ERROR,info);
				handleExcep(excep);
			}

			// now let's read one by one
			// it's number of row 
			string value = l.findValue(1);
			UInt nRow = 0;
			if (w.toUInt(value,nRow)) {
				status.push_back(nRow);
			}else{
				string info = "nRow in data section is invalid in the record.txt for location: " + location;
				Excep excep("HistDataMan","constructor",EXCEPTION_HISTDATAMAN_ERROR,info);
				handleExcep(excep);
			}

			// it's number of col
			value = l.findValue(2);
			UInt nCol = 0;
			if (w.toUInt(value,nCol)) {
				status.push_back(nCol);
			}else{
				string info = "ncol in data section is invalid in the record.txt for location: " + location;
				Excep excep("HistDataMan","constructor",EXCEPTION_HISTDATAMAN_ERROR,info);
				handleExcep(excep);
			}

			// whether the data is symmetrized
			value = l.findValue(3);
			UInt symm = 0;
			if (w.toUInt(value,symm)) {
				if (symm == SYMM_DATA || symm == ASYMM_DATA) {
					status.push_back(symm);
				}else{
					string info = "symm status is not in allowed value in reading from record.txt";
					Excep excep("HistDataMan","constructor",EXCEPTION_HISTDATAMAN_ERROR,info);
					handleExcep(excep);
				}
			}else{
				string info = "symm status is invalid in the record.txt for location: " + location;
				Excep excep("HistDataMan","constructor",EXCEPTION_HISTDATAMAN_ERROR,info);
				handleExcep(excep);
			}

			// finally, for nSpin == 2 case; we need to record another copy
			// this is for the beta part
			if (nSpin == 2) {
				status.push_back(nRow);
				status.push_back(nCol);
				status.push_back(symm);
			}
		}
	}

	// now close the file
	inf.close();
}

void HistDataMan::updateRecord(const UInt& nRow, const UInt& nCol, const UInt& symmStatus) const
{

	// do we have the record.txt yet?
	path folder(location.c_str());
	path record("record.txt");
	folder /= record;
	string file = folder.string();
	if (! exists(folder)) {

		// open the output stream
		std::ofstream outf;
		const char * output_file = file.c_str();
		outf.open(output_file,std::ofstream::out);

		// print out the comments
		string line = "# this file records data for location: " + location;
		outf << line << endl;
		line = "nspin = " + boost::lexical_cast<string>(nSpin);
		outf << line << endl;
		outf.close();
	}

	// let's open the file, perhaps again
	std::ofstream outf;
	const char * output_file = file.c_str();
	outf.open(output_file,std::ofstream::app);

	// print out the data
	string line = "data_section  ";
	line = line + boost::lexical_cast<string>(nRow) + "  ";
	line = line + boost::lexical_cast<string>(nCol) + "  ";
	line = line + boost::lexical_cast<string>(symmStatus) + "  ";
	outf << line << endl;
	outf.close();
}

string HistDataMan::getNewDataFileName(UInt iSpin) const
{
	// get the name
	// number of data starts from 0,1,2,3,4,5....
	// data name starts from 0,1,2,3,4....
	UInt n = nData();
	UInt newDataSection = 0;
	if (n > 0) newDataSection = n;
	string name = boost::lexical_cast<string>(newDataSection);
	if (nSpin == 2) {
		string suffix = "a";
		if (iSpin == 1) suffix = "b";
		name += "_";
		name += suffix;
	}
	name += ".bin";
	return name;
}

string HistDataMan::getDataFileName(const UInt& dataSectionNumber, UInt iSpin) const
{
	// check the given data section number
	UInt n = nData();
	if(dataSectionNumber>=n) {
		string info = "invalid data section number given for location: " + location;
		Excep excep("HistDataMan","getDataFileName",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}

	// get the name
	string name = boost::lexical_cast<string>(dataSectionNumber);
	if (nSpin == 2) {
		string suffix = "a";
		if (iSpin == 1) suffix = "b";
		name += "_";
		name += suffix;
	}
	name += ".bin";
	return name;
}

void HistDataMan::storeData(const DoubleVec& oriData, UInt iSpin)
{
	// store data
	if (useFile()) {
		string name = getNewDataFileName(iSpin);
		FileReadWrite fileRW(location,name);
		fileRW.write(&oriData.front(),oriData.size());
	}else{
		dataList.push_back(oriData);
	}

	// save data status
	status.push_back(oriData.size());
	status.push_back(0);
	status.push_back(ASYMM_DATA);

	// finally, update the record.txt
	if (useFile()) {
		updateRecord(oriData.size(),0,ASYMM_DATA);
	}
}	

void HistDataMan::storeData(const Mtrx& oriData, bool isSymm, UInt iSpin)
{
	// for symm matrix, we only store its lower triangular part
	if (isSymm) {

		// check whether row == col
		if(oriData.getRow()!=oriData.getCol()) {
			string info = "input symmetry status conflicts with input matrix(row!=col) for given location: " + location;
			Excep excep("HistDataMan","storeData",EXCEPTION_HISTDATAMAN_ERROR,info);
			handleExcep(excep);
		}	

		// form symmetrical form of vector
		UInt n = oriData.getRow();
		UInt len = n*(n+1)/2;
		DoubleVec vecData(len,ZERO);
		oriData.formSymmVec(vecData);

		// save it to file
		if (useFile()) {
			string name = getNewDataFileName(iSpin);
			FileReadWrite fileRW(location,name);
			fileRW.write(&vecData.front(),len);
		}else{
			dataList.push_back(vecData);
		}

		// save data length
		status.push_back(n);
		status.push_back(n);
		status.push_back(SYMM_DATA);

		// finally, update the record.txt
		if (useFile()) {
			updateRecord(n,n,SYMM_DATA);
		}
	} else{	

		// save the whole matrix
		const DoubleVec& vec = oriData.getVec();
		if (useFile()) {
			string name = getNewDataFileName(iSpin);
			FileReadWrite fileRW(location,name);
			fileRW.write(&vec.front(),vec.size());
		}else{
			dataList.push_back(vec);
		}

		// save data length
		status.push_back(oriData.getRow());
		status.push_back(oriData.getCol());
		status.push_back(ASYMM_DATA);

		// finally, update the record.txt
		if (useFile()) {
			updateRecord(oriData.getRow(),oriData.getCol(),ASYMM_DATA);
		}
	}
}	

void HistDataMan::retrieveData(DoubleVec& data, const UInt& dataSectionNumber, 
		UInt iSpin) const
{
	// additional check 
	if(dataSectionNumber>=nData()) {
		string info = "invalid data section number given for location: " + location;
		Excep excep("HistDataMan","retrieveData, vector",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}
	if(data.size()!=vecDim(dataSectionNumber,iSpin)) { 
		string info = "input vector size conflicts with data stored for location: " + location;
		Excep excep("HistDataMan","retrieveData, vector",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}

	// set up tmp vector
	UInt len = vecDim(dataSectionNumber,iSpin);
	DoubleVec vec(len);

	// retrieve data
	if (useFile()) {
		string name = getDataFileName(dataSectionNumber,iSpin);
		FileReadWrite fileRW(location,name);
		fileRW.read(&vec.front(),len);
	}else{
		retrieveDataFromMem(vec,dataSectionNumber,iSpin);
	}

	// load data
	for(UInt i=0; i<len; i++) {
		data[i] += vec[i];
	}
}	

void HistDataMan::retrieveData(Mtrx& M, const UInt& dataSectionNumber, UInt iSpin) const
{
	// additional check 
	if(dataSectionNumber>=nData()) {
		string info = "invalid data section number given for location: " + location;
		Excep excep("HistDataMan","retrieveData, matrix",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}
	if(M.getRow()!=rowDim(dataSectionNumber,iSpin)) {
		string info = "input matrix row conflicts with data stored for location: " + location;
		Excep excep("HistDataMan","retrieveData, matrix",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}	
	if(M.getCol()!=colDim(dataSectionNumber,iSpin)) {
		string info = "input matrix col conflicts with data stored for location: " + location;
		Excep excep("HistDataMan","retrieveData, matrix",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}	

	// whether the data is stored in symmetrized way
	bool isSymm = false;
	if (isSymmData(dataSectionNumber)) isSymm = true;

	// for symm matrix, we only store its lower triangular part
	if (isSymm) {

		// check 
		if(M.getRow()!=M.getCol()) {
			string info = "input symmetry status conflicts with input matrix(row!=col) for given location: " + location;
			Excep excep("HistDataMan","retrieveData",EXCEPTION_HISTDATAMAN_ERROR,info);
			handleExcep(excep);
		}	

		// form symmetrical form of vector
		UInt n = M.getRow();
		UInt len = n*(n+1)/2;
		DoubleVec vecData(len,ZERO);

		// loading data from file/memory
		if (useFile()) {
			string name = getDataFileName(dataSectionNumber,iSpin);
			FileReadWrite fileRW(location,name);
			fileRW.read(&vecData.front(),vecData.size());
		}else{
			retrieveDataFromMem(vecData,dataSectionNumber,iSpin);
		}

		// finally load in data into matrix
		M.loadSymmVec(vecData);

	}else{	

		// set up tmp vector and load in data
		DoubleVec vec(rowDim(dataSectionNumber,iSpin)*colDim(dataSectionNumber,iSpin));
		if (useFile()) {
			string name = getDataFileName(dataSectionNumber,iSpin);
			FileReadWrite fileRW(location,name);
			fileRW.read(&vec.front(),vec.size());
		}else{
			retrieveDataFromMem(vec,dataSectionNumber,iSpin);
		}

		// add it into the matrix
		DoubleVec& oriVec = M.getVec();
		for(UInt i=0; i<oriVec.size(); i++) {
			oriVec[i] += vec[i];
		}
	}
}	

void HistDataMan::retrieveDataFromMem(DoubleVec& vec, const UInt& dataSectionNumber,
		UInt iSpin) const
{
	// additional check
	if (dataSectionNumber>=nData()) {
		string info = "invalid data section number given for location: " + location;
		Excep excep("HistDataMan","retrieveDataFromMem",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}

	// get the dimension data and check
	UInt rDim = rowDim(dataSectionNumber,iSpin);
	UInt cDim = colDim(dataSectionNumber,iSpin);
	UInt dim  = rDim;

	// dim
	// for vector we just use rDim
	// for matrix we need to know whether it's symmetrical or not
	if (cDim>0) {
		if (isSymmData(dataSectionNumber,iSpin)) {
			dim = rDim*(rDim+1)/2;
		}else{
			dim *= cDim;
		}
	}

	// check
	if(vec.size()!=dim) {
		string info = "input vector length conflicts with data stored for location: " + location;
		Excep excep("HistDataMan","retrieveDataFromMem",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}

	// consider the section of data in terms of spin states
	// for spin polarized situation,
	// for alpha, the data section number is like 0 2 4 6 8...
	// for beta,  the data section number is like 1 3 5 7 9...
	UInt section = sect(dataSectionNumber,iSpin);

	// now loop over to get the data
	UInt count = 0;
	bool getData = false;
	for(HistDataList::const_iterator it= dataList.begin(); it != dataList.end(); ++it) {
		if (count == section) {
			vec = *it;
			getData = true;
			break;
		}
		count++;
	}

	// check whether we have got the data
	if(! getData) {
		string info = "fail to get data from memory in location: " + location;
		Excep excep("HistDataMan","retrieveDataFromMem",EXCEPTION_HISTDATAMAN_ERROR,info);
		handleExcep(excep);
	}
}

