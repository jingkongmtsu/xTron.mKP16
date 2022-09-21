//
// this is a small function to do transpose work
// for vector form of matrix, used in testing
//
#include "libgen.h"
#include<vector>
using namespace std;

void transpose(UInt row, UInt col, vector<Double>& A)
{
	vector<Double>vTmp(col*row);
	for(UInt i=0; i<col; i++) {
		for(UInt j=0; j<row; j++) {
			vTmp[i+j*col]=A[j+i*row];
		}
	}
	for(UInt i=0; i<row*col; i++) {
		A[i] = vTmp[i];
	}
}

void lowToUpper(UInt n, vector<Double>& A)
{
	for(UInt j=0; j<n; j++) {
		for(UInt i=j; i<n; i++) {
			A[j+i*n]=A[i+j*n];
		}
	}
}
