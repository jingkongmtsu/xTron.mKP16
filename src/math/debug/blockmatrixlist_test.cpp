//
//  note:
//  testing the merge, update, updateMatrix functions, all pass
//
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<vector>
#include "omp_setting.h"
#include "scalablevec.h"
#include "blockmatrixlist.h"
using namespace blockmatrixlist;
using namespace std;
using namespace omp_setting;

void blockmatrixlist_test()
{
	// lapcck testing
	cout << "----------------------------------------------" << endl;
	cout << "block matrix list testing                     " << endl;
	cout << "test update function                          " << endl;
	cout << "test merge function                           " << endl;
	cout << "test updateMatrix function                    " << endl;
	cout << "----------------------------------------------" << endl;
	cout << endl;

	// heading
	cout << "----------------------------------------------" << endl;
	cout << "test block matrix list in default constructor " << endl;
	cout << "----------------------------------------------" << endl;
	{
		// set up empty block matrix list and add in something
		cout << "testing update function" << endl;
		BlockMtrxList list;
		for(UInt n=0; n<2; n++) {

			// initial the block matrix
			UInt row = 3;
			UInt col = 3;
			BlockMtrx x(row,col,n*3,n*3);
			for(UInt i=0; i<col; i++) {
				for(UInt j=0; j<row; j++) {
					x(j,i) = i+2*j+1;
				}
			}

			// print it
			x.blockPrint("block matrix");

			// now update list
			list.update(x);
		}
		list.print("block matrix list");
		cout << endl << endl;

		// set up empty block matrix list and add in something
		cout << "testing merge function between list2 and list" << endl;
		BlockMtrxList list2;
		for(UInt n=0; n<2; n++) {

			// initial the block matrix
			UInt row = 3;
			UInt col = 3;
			BlockMtrx x(row,col,n*3,0);
			for(UInt i=0; i<col; i++) {
				for(UInt j=0; j<row; j++) {
					x(j,i) = i+2*j+1;
				}
			}
			// now update list
			list2.update(x);
		}
		list2.print("another block matrix list which is going to add in the list");
		list.merge(list2);
		list.print("block matrix list after merging");
		cout << endl << endl;

		// construct mother matrix
		UInt rowa = 6;
		UInt cola = 6;
		Mtrx a(rowa,cola);
		cout << "initial matrix is all zero " << endl;
		list.updateMatrix(a);
		a.print("matrix after updating");
	}

	// heading
	cout << "----------------------------------------------" << endl;
	cout << "test block matrix list non-default constructor" << endl;
	cout << "----------------------------------------------" << endl;
	{
		// set up the initial block list
		BlockMtrxList list(3,2,4,4);
		cout << "the initial block list size " << list.len() << endl;
		cout << "block list length limit " << list.getLenLimit() << endl;

		// set up empty block matrix list and add in something
		cout << "testing update function" << endl;
		for(UInt n=0; n<2; n++) {

			// initial the block matrix
			UInt row = 3;
			UInt col = 3;
			BlockMtrx x(row,col,n*3,n*3);
			for(UInt i=0; i<col; i++) {
				for(UInt j=0; j<row; j++) {
					x(j,i) = i+2*j+1;
				}
			}

			// print it
			x.blockPrint("block matrix");

			// now update list
			list.update(x);
		}
		list.print("block matrix list");
		cout << endl << endl;

		// set up empty block matrix list and add in something
		cout << "testing merge function between list2 and list" << endl;
		BlockMtrxList list2;
		for(UInt n=0; n<2; n++) {

			// initial the block matrix
			UInt row = 3;
			UInt col = 3;
			BlockMtrx x(row,col,n*3,0);
			for(UInt i=0; i<col; i++) {
				for(UInt j=0; j<row; j++) {
					x(j,i) = i+2*j+1;
				}
			}
			// now update list
			list2.update(x);
		}
		list2.print("another block matrix list which is going to add in the list");
		list.merge(list2);
		list.print("block matrix list after merging");
		if (list.reachLenLimit()) {
			cout << "reaching block matrix list length limit" << endl;
		}
		cout << endl << endl;

		// construct mother matrix
		UInt rowa = 6;
		UInt cola = 6;
		Mtrx a(rowa,cola);
		cout << "initial matrix is all zero " << endl;
		list.updateMatrix(a);
		a.print("matrix after updating");
	}
}

