/**
 * \file  filerw.h
 * \author Fenglai Liu and Jing Kong
 */

#ifndef FILERW_H
#define FILERW_H
#include "libgen.h"
#include<string>

namespace filerw {

	using namespace std;

	/**
	 * \class FileReadWrite
	 * \brief class for reading and writing a data file
	 *
	 * Since we use the size_t (UInt), therefore it should not have 
	 * the 2G file limit for the file functions here (for 64 bit system)
	 * for 32 bit system, the large file limit is still there
	 *
	 * Currently we support the double precision vector reading/writing
	 * in the raw type, Double* 
	 */
	class FileReadWrite {

		private:

			string fileName;   ///< file name for reading or writing

		public:

			/**
			 * normal file reading/writing constructor
			 * normal means that the data is stored in binary way
			 */
			FileReadWrite(const string& dirName, const string& dataFileName); 

			/**
			 * only for debug use
			 * in this way, the data file is provided externally
			 */
			FileReadWrite(const string& input):fileName(input) { };

			/**
			 * destructor
			 */
			~FileReadWrite() { };

			/**
			 * read in data from the file into the data
			 */
			void read(Double* data,const UInt& dataLength) const;

			/**
			 * write data from the vector into the file
			 */
			void write(const Double* data, const UInt& dataLength) const;

			/**
			 * only for debug purpose
			 * we read matrix data from the text file
			 */
			void readMatrixFromTextFile(Double* array, UInt row, UInt col) const;
	};
}


#endif
