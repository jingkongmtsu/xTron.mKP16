/**
 * \file   histdataman.h
 * \brief  class for manipulating large data 
 * \author Fenglai Liu 
 */

#ifndef HISTDATAMAN_H
#define HISTDATAMAN_H
#include<vector>
#include<list>
#include<string>
#include "libgen.h"
#include "scalablevec.h"

namespace globalinfor {
	class GlobalInfor;
}

namespace matrix {
	class Mtrx;
}

namespace histdataman {

	using namespace globalinfor;
	using namespace matrix;
	using namespace std;

	///
	/// number of data field in the status vector
	/// namely, for each data field, it's row -> col -> whether it's symmetrical
	///
	const UInt N_FIELD = 3;

	///
	/// define the data list 
	///
	typedef std::list<DoubleVec,tbb::scalable_allocator<DoubleVec> > HistDataList;

	///
	/// define the symmetry data status
	/// symm_data: for matrix only half of data is stored, length is (n+1)*n/2
	/// asymm_data: the whole matrix is stored
	/// for vector, it's also defined as asymmetrical
	///
	const UInt SYMM_DATA   = 0;
	const UInt ASYMM_DATA  = 1;

	/**
	 * \class Historical data manipulation
	 *
	 * Some important rule set up here:
	 * if the data is spin-polarized, 
	 * then the data stored must obey in this order:
	 * alpha, beta, alpha, beta.....
	 * We rely on this order for retrieving data correctly!
	 *
	 * All of spin states are default to be alpha (0)
	 */
	class HistDataMan {

		protected:

			UInt    nSpin;       ///< number of spin states for the data
			string location;     ///< where the data folder resides, if none; then data
		                        ///< is stored in memory	

			///
			/// status for each data section
			/// each data status contains three kind of field:
			/// row of data, col of data and whether the data is symmetry
			/// row and col is used for storing matrix data dimension
			/// for vector, col is set to 0
			/// if the given matrix is symmetrical, then the final data field is SYMM_DATA 
			/// if it's vector or asymmetrized matrix, then final data field is ASYMM_DATA
			///
			/// dimension for status vector is like this:
			/// (N_FIELD,nSpin,nData)
			///
			/// if everything is stored on files, then in the folder we will have a txt file
			/// named as "record.txt" to repeat the information contained in status.
			///
			UIntVec status;  

			///
			/// if we keep all of data in memory, then we will use
			/// the following data
			///
			HistDataList dataList;

			///
			/// get the location name
			/// \param    infor:     where is scratch dir
			/// \param    name :     data folder name
			/// \param    section:   molecular section
			///
			string setLocation(const GlobalInfor& infor, const string& name, 
					const UInt& section) const;

			///
			/// create the new data file name based on nData and iSpin
			///
			string getNewDataFileName(UInt iSpin = 0) const;

			///
			/// used when you want to retrieve data for given section
			/// and spin state. This is only for getting files from 
			/// disk
			///
			string getDataFileName(const UInt& dataSectionNumber, UInt iSpin = 0) const;

			///
			/// update the record file inside the folder
			///
			void updateRecord(const UInt& nRow, const UInt& nCol, const UInt& symmStatus) const;

			///
			/// function used to retrieve data for the data list for the 
			/// given section and given spin state
			///
			void retrieveDataFromMem(DoubleVec& vec, const UInt& dataSectionNumber,
					UInt iSpin = 0) const;

		public:

			/**
			 * constructor used to set up the object and save data to files etc.
			 * according to the given name
			 * \param  infor   used to retrieve the sratch folder
			 * \param  name    the data folder name
			 * \param  sec     the molecular section number 
			 * \param  nSpins  data is spin-polarized? then nSpins = 2; else nSpins = 1
			 * \param  inFile  everything is saved on File?
			 */ 
			HistDataMan(const GlobalInfor& infor, const string& name, const UInt& sec, 
					const UInt& nSpins, bool inFile);

			/**
			 * constructor used to set up the constructor and recovering data from file
			 * according to the given name. In this case, the data must be already
			 * stored in files on disk
			 * \param  infor   used to retrieve the sratch folder
			 * \param  name    the data folder name
			 * \param  sec     the molecular section number 
			 */ 
			HistDataMan(const GlobalInfor& infor, const string& name, const UInt& sec);

			/**
			 * destructor
			 */
			~HistDataMan() { };

			/////////////////////////////////////////////////////
			//        functions to retrieve/store data         //
			/////////////////////////////////////////////////////

			/**
			 * store the given data in the given spin states
			 * data in vector form
			 */
			void storeData(const DoubleVec& oriData, UInt iSpin = 0);

			/**
			 * store the given data in the given spin states
			 * data in matrix form
			 * isSymm indicates that whether matrix is symmetrical
			 */
			void storeData(const Mtrx& oriData, bool isSymm, UInt iSpin = 0);

			/**
			 * retrieve Data in the given spin state and section
			 */
			void retrieveData(DoubleVec& data, const UInt& dataSectionNumber, 
					UInt iSpin = 0) const;

			/**
			 * retrieve Data for matrix in the given spin state and section
			 */
			void retrieveData(Mtrx& M, const UInt& dataSectionNumber, UInt iSpin = 0) const;

			/////////////////////////////////////////////////////
			//  functions to retrieve data status information  //
			/////////////////////////////////////////////////////

			/**
			 * get the current total number of data 
			 * the number is spin state indepedent
			 */
			UInt nData() const {
				return status.size()/(N_FIELD*nSpin);
			};

			/**
			 * get the index of data section when it's stored in memory
			 * this is used to locate the data when it's saved in memory
			 */
			UInt sect(const UInt& section, const UInt& iSpin) const {
				return nSpin*section+iSpin;
			};

			/**
			 * get the row dimension for matrix type of data
			 */
			UInt rowDim(const UInt& section, UInt iSpin = 0) const {
				UInt pos = N_FIELD*(section*nSpin+iSpin);
				return status[pos];
			};

			/**
			 * get the col dimension for matrix type of data
			 */
			UInt colDim(const UInt& section, UInt iSpin = 0) const {
				UInt pos = N_FIELD*(section*nSpin+iSpin)+1;
				return status[pos];
			};

			/**
			 * get the dimension for vector type of data
			 */
			UInt vecDim(const UInt& section, UInt iSpin = 0) const {
				UInt pos = N_FIELD*(section*nSpin+iSpin);
				return status[pos];
			};

			/**
			 * whether the given data is stored in symmetrical way
			 * that is to say, for symmetrical way only half of it is stored
			 */
			bool isSymmData(const UInt& section, UInt iSpin = 0) const {
				UInt pos = N_FIELD*(section*nSpin+iSpin)+2;
				if (status[pos] == ASYMM_DATA) return false;
				return true;
			};

			/**
			 * do we use file to store the data?
			 */
			bool useFile() const {
				if (location == "NONE") return false;
				return true;
			};

	};
}


#endif
