/**
 * \file    filerw.cpp
 * \brief   cpp file for file reading and writing
 * \author  Fenglai Liu and Jing Kong
 */
#include <boost/filesystem.hpp>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "textread.h"
#include "excep.h"
#include "filerw.h"
using namespace excep;
using namespace textread;
using namespace boost::filesystem;
using namespace filerw;

FileReadWrite::FileReadWrite(const string& dirName, const string& dataFileName)
{
	// composing the path for file name
	path dir(dirName.c_str());
	path file(dataFileName.c_str());
	dir /= file;
	fileName = dir.string();
}

void FileReadWrite::read(Double* data, const UInt& dataLength) const
{
	UInt dataTypeLength = sizeof(Double);
	FILE *fp = fopen(fileName.c_str(),"rb");  
	if (fp!=NULL) {  
		UInt status = fread(&data[0],dataTypeLength,dataLength,fp);  
		if (status != dataLength) {
			Excep excep("FileReadWrite","read",EXCEPTION_BIN_FILE_READ_ERROR,fileName);
			handleExcep(excep);
		}
	}else{
		Excep excep("FileReadWrite","read",EXCEPTION_FILE_MISSING,fileName);
		handleExcep(excep);
	}	
	fclose(fp);
}

void FileReadWrite::write(const Double* data, const UInt& dataLength) const
{
	UInt dataTypeLength = sizeof(Double);
	FILE *fp = fopen(fileName.c_str(),"wb");  
	if (fp!=NULL) {  
		UInt status = fwrite(&data[0],dataTypeLength,dataLength,fp);  
		if (status != dataLength) {
			Excep excep("FileReadWrite","write",EXCEPTION_BIN_FILE_WRITE_ERROR,fileName);
			handleExcep(excep);
		}
	}else{
		Excep excep("FileReadWrite","write",EXCEPTION_FILE_WRITE_FAIL,fileName);
		handleExcep(excep);
	}	
	fclose(fp);
}

void FileReadWrite::readMatrixFromTextFile(Double* array, UInt row, UInt col) const
{
	// open the input file
	std::ifstream inf;
	const char * input_file = fileName.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		Excep excep("FileReadWrite","readMatrixFromTextFile",EXCEPTION_FILE_MISSING,fileName);
		handleExcep(excep);
	}

	// read in data
	string line;
	UInt colCount = 0;
	UInt rowCount = 0;
	WordConvert w;
	while(getline(inf,line)) {

		// basic input
		LineParse l(line);
		if (l.isCom()) continue;
		if (l.isEmp()) break;

		// test that whether the given line is index line
		// the index line has all data to be int since they are index
		// on the other hand, if only one data in the line,
		// it must be the index line since we always has 
		// the index before the data line
		UInt ncolperline = l.getNPieces();
		if (ncolperline == 1) continue;
		if (ncolperline > 1) {
			string t1 = l.findValue(0);
			string t2 = l.findValue(1);
			if (w.isInt(t1) && w.isInt(t2)) continue;
		}

		// read in the concencutive data block
		// we start from the second data since the first
		// one is the row index
		for(UInt i=1; i<ncolperline; i++) {
			Double tmp = ZERO;
			w.toDouble(l.findValue(i),tmp);
			array[rowCount+(colCount+i-1)*row] = tmp;
		}
		rowCount++;

		// reset the row, also incremental the col count
		if (rowCount == row) {
			rowCount = 0;
			colCount += ncolperline - 1;
		}
	}
	inf.close();
}

