/**
 * \file    globalinfor.h
 * \author  Fenglai Liu 
 */
#ifndef GLOBALINFOR_H
#define GLOBALINFOR_H
#include "libgen.h"
#include <string>

namespace globalinfor {

	using namespace std;

	/**
	 * \class   GlobalInfor
	 * \brief   gathering global setting information for running job
	 */
	class GlobalInfor {

		private:

			//
			// job related infor
			//
			string input;      ///< input file where we can located
			string scrDir;     ///< scratch space for calculation

			//
			// control related to the parallel
			//
			bool enableMultiThreads;  ///< do we open multi-threads mode(if it's disabled we do not use MIC and GPGPU)
			UInt nCPUThreads;         ///< number of CPU threads that assigned by user or use default one

			///////////////////////////////////////
			//          member functions         //
			///////////////////////////////////////

			/**
			 * forming the scratch dir folder name for current job
			 */
			void formScratchDir();

			/**
			 * parsing the information defined by user input file
			 */
			void inforParse();

		public:

			/**
			 * by giving the input file, we are able to know everything 
			 */
			GlobalInfor(const string& in);

			/**
			 * destructor
			 */
			~GlobalInfor() { };

			/**
			 * return the input file
			 */
			const string& getInputFile() const {
				return input;
			};

			/**
			 * return the scratch dir
			 */
			const string& getScratchDir() const {
				return scrDir;
			};

			/**
			 * do we use multi-threads?
			 */
			bool useMultiThreads() const {
				return enableMultiThreads;
			};

			/**
			 * to enable the multi-threads as program requires
			 */
			void enableMultiThreadsMode() {
				enableMultiThreads = true;
			};

			/**
			 * to disable the multi-threads as program requires
			 */
			void disableMultiThreadsMode() {
				enableMultiThreads = false;
			};

			/**
			 * get the number of threads
			 */
			UInt getNCPUThreads() const {
				return nCPUThreads;
			};

			/**
			 * whether we use GPGPU?
			 */
			bool useGPGPU() const { return false; };

			/**
			 * debug print
			 */
			void print() const;
	};

}

#endif
