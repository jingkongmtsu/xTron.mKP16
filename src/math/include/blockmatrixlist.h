/**
 * \file   blockmatrixlist.h
 * \author Fenglai Liu 
 */
#ifndef BLOCKMATRIXLIST_H
#define BLOCKMATRIXLIST_H
#include "libgen.h"
#include<vector>
#include "blockmatrix.h"
using namespace std;
using namespace blockmatrix;

namespace blockmatrixlist {

	// define the block matrix list
	typedef vector<BlockMtrx,tbb::scalable_allocator<BlockMtrx> > BlockMatrices;

	///
	/// class for block matrix list 
	///
	/// mathematically, this could be used for sparse matrix 
	/// this class is actually a wrapper class for BlockMatrices data
	///
	/// to use this class, we have two ways.
	///
	/// the first way is very simple, you just contruct an empty list using 
	/// the default constructor, then adding in data when you have it.
	///
	/// however, the first method will cause the inside block matrix data
	/// memory rellocation when you call the clear function. This is not good
	/// when you want to frequently use the block matrix list to store data,
	/// then move the data, clear it then use it again.
	///
	/// if you want to frequently use the block matrix list, it's better that
	/// you know the matrix length limit first, together with the block matrix
	/// size limit (just like we use the block matrix list to store the tmp result
	/// in gints4d.cpp, see the usage over there). In this way we will initialize
	/// the block matrix list to hold the memory first, then you can add in data,
	/// move the data to other places; clear data and use the list again.
	///
	/// in the second way you must have a length limit value, and we will use 
	/// the curPos to track the last block matrix position. If the length limit 
	/// is reached, then it means the list is full and you need to move the data out.
	///
	/// the length limit is only a "soft" limit. It means you can be over this limit.
	/// It's the user's responsibility to watch that whether the limit is reached.
	///
	/// a good example to use the block matrix list is in the gints4d.cpp. Please 
	/// see the code over there.
	///
	class BlockMtrxList {

		protected:

			///
			/// sometimes we impose a block matrix length limit, so that
			/// to ensure that the block matrix number in the list is less
			/// than the given limit. In default such limit is just set
			/// to 0, means no limit 
			///
			UInt lenLimit;            

			///
			/// this is the current block matrix position for filling in data
			///
			UInt curPos;

			BlockMatrices blockList;  ///< a list of block matrices

		public:

			///
			/// constructor - set up an empty list for default use
			///
			BlockMtrxList():lenLimit(0),curPos(0) { };

			///
			/// constructor - set up the list for hold the memory first
			///
			/// this constructor will form the block matrix list fully
			/// it will constructor the list and block matrix to hold the memory first
			/// later when you do add or clear it does not affects the memory hold inside
			///
			/// \param length limit: this is the length limit set for the list
			/// \param incLengthFactor: when do initialization, the list is not strickly on the length limit.
			///                         we will have additional increment on length for allocation, this is what
			///                         this value about
			/// \param nRowLimit: block matrix row limit
			/// \param nColLimit: block matrix col limit
			///
			BlockMtrxList(const UInt& lengthLimit, const UInt& incLengthFactor, 
					const UInt& nRowLimit, const UInt& nColLimit);

			///
			/// we can also initialize the memory later with this function
			/// but actually it's same with the above constructor
			///
			void initMem(const UInt& lengthLimit, const UInt& incLengthFactor,
					const UInt& nRowLimit, const UInt& nColLimit);

			///
			/// destructor
			///
			~BlockMtrxList() { };

			///
			/// copy constructor
			///
			BlockMtrxList(const BlockMtrxList& A):lenLimit(A.lenLimit),curPos(A.curPos),
			blockList(A.blockList) { }; 

			///
			/// operator =
			///
			BlockMtrxList& operator=(const BlockMtrxList& A) {
				if (this == &A) {
					return *this;
				}else{
					lenLimit  = A.lenLimit;
					curPos    = A.curPos;
					blockList = A.blockList;
					return *this;
				}
			};

			///
			/// get the actually length of list
			/// if the block list is pre-initialized, then
			/// the actually length is the curPos
			///
			UInt len() const { 
				if (lenLimit == 0) {
					return blockList.size(); 
				}
				return curPos;
			};

			///
			/// get the length limit
			///
			UInt getLenLimit() const { return lenLimit; };

			///
			/// whether the block list reaches it's length limit?
			///
			bool reachLenLimit() const {
				if (lenLimit == 0) return false;
				return (curPos >= lenLimit);
			};	

			///
			/// clear the content of the matrix list
			///
			void clear(); 

			///
			/// get the block matrix - const type
			///
			const BlockMtrx& getBlockMatrix(const UInt& i) const { return blockList[i]; }; 

			///
			/// get the block matrix - non-const type
			///
			BlockMtrx& getBlockMatrix(const UInt& i) { return blockList[i]; }; 

			///
			/// operation of update
			/// update performs two kinds of operations,
			/// 1 if the given block matrix is new to the list, we append it
			///   to the end of the list;
			/// 2 if the given block matrix already exists in the list, we 
			///   then update the old one with the new one's value
			///
			void update(const BlockMtrx& A);

			///
			/// operation of merge
			/// merge the given list into this one
			///
			void merge(const BlockMtrxList& A);

			///
			/// update this list into the given matrix
			///
			/// here the function loop over each block and add it into
			/// the given matrix
			///
			/// we note that all of the block matrices are just added
			/// into the A in "normal" order, which means no transpose
			/// status is pre-defined
			///
			void updateMatrix(Mtrx& A) const;

			///
			/// whether the list has the same Block matrix in it?
			///
			bool hasIt(const BlockMtrx& A) const;

			///
			/// print out content of block matrix
			///
			void print(const string& title) const;
	};

}

#endif

