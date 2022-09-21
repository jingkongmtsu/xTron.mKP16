/**
 * \file    localmemscr.h
 * \author  Fenglai Liu 
 */
#ifndef LOCALMEMSCR_H
#define LOCALMEMSCR_H
#include "libgen.h"
#include "scalablevec.h"

namespace localmemscr {

	using namespace std;

	/**
	 * \class   LocalMemScr
	 * \brief   manipulating local memory for single thread
	 *
	 * Currently the local memory scratch is only for double type of data
	 *
	 * The implementation of this class is very simple. However, it's not 
	 * flexible; for example; if you have a lot of places which need 
	 * memory, and these places are totally independent with each other;
	 * such class will cause memory waste since it will always hold the 
	 * memory until reset is called.
	 *
	 * Use it with your caution!
	 */
	class LocalMemScr {

		private:

			UInt lastMemUser;      ///< the last memory user position
			UInt lastMemLen;       ///< for the last memory user, it's required mem length
			DoubleVec memHolder;   ///< holding memory for the scratch

		public:

			/**
			 * initilize the local memory sratch
			 */
			LocalMemScr(const UInt& len):lastMemUser(0),lastMemLen(0),
			memHolder(len,ZERO) { };

			/**
			 * destructor
			 */
			~LocalMemScr() { };

			/**
			 * return a new memory position according to the length
			 * it required
			 */
			Double* getNewMemPos(const UInt& len); 

			/**
			 * reset the local memory scratch
			 * to erase all of memory users information
			 */
			void reset() {
				lastMemUser = 0;
				lastMemLen  = 0;
				memHolder.assign(memHolder.size(),ZERO);
			};
	};

}

#endif
