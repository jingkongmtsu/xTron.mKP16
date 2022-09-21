#include<string> 
using namespace std;

///
/// BLAS test function
///
extern void blas_test(bool testblas1, bool testblas2, bool testblas3);

///
/// BLAS test for parallel
///
extern void para_blas_test(bool testblas1, bool testblas2, bool testblas3);

///
/// LAPACK test function
///
extern void lapack_test();

///
/// parallel LAPACK test function
///
extern void para_lapack_test();

///
/// matrix test function
///
extern void matrix_test();

///
/// blockmatrix test function
///
extern void blockmatrix_test();

///
/// blockmatrixlist test function
///
extern void blockmatrixlist_test();

///
/// matrix test for parallel
///
extern void para_matrix_test();

int main(int argc, char* argv[])
{
	// setting the testing jobs
	bool testingMatrix          = false;
	bool testingBLAS1           = false;
	bool testingBLAS2           = false;
	bool testingBLAS3           = false;
	bool testingLAPACK          = false;
	bool testingBlockMatrix     = false;
	bool testingBlockMatrixList = false;

	// parallel testing
	bool testingParaBLAS1  = false;
	bool testingParaBLAS2  = false;
	bool testingParaBLAS3  = false;
	bool testingParaLAPACK = false;
	bool testingParaMatrix = false;

	// parse the input parameter
	for(int i=1; i<argc; i++) {
		string com = argv[i];

		// serial testing
		if (com == "blas1" ) testingBLAS1 = true;
		if (com == "blas2" ) testingBLAS2 = true;
		if (com == "blas3" ) testingBLAS3 = true;
		if (com == "lapack" ) testingLAPACK = true;
		if (com == "matrix" ) testingMatrix = true;
		if (com == "blockmatrix" ) testingBlockMatrix = true;
		if (com == "blockmatrixlist" ) testingBlockMatrixList = true;

		// parallel testing
		if (com == "pblas1" ) testingParaBLAS1  = true;
		if (com == "pblas2" ) testingParaBLAS2  = true;
		if (com == "pblas3" ) testingParaBLAS3  = true;
		if (com == "plapack") testingParaLAPACK = true;
		if (com == "pmatrix") testingParaMatrix = true;
	}

	// blas
	if (testingBLAS1 || testingBLAS2 || testingBLAS3) {
		blas_test(testingBLAS1,testingBLAS2,testingBLAS3);
	}

	// parallel blas test
	// all 3 level of BLAS will be done together
	if (testingParaBLAS1 || testingParaBLAS2 || testingParaBLAS3) {
		para_blas_test(testingParaBLAS1,testingParaBLAS2,testingParaBLAS3);
	}

	// lapack
	if (testingLAPACK) {
		lapack_test();
	}

	// parallel lapack
	if (testingParaLAPACK) {
		para_lapack_test();
	}

	//matrix 
	if (testingMatrix) {
		matrix_test();
	}

	//blockmatrix 
	if (testingBlockMatrix) {
		blockmatrix_test();
	}

	//blockmatrixlist 
	if (testingBlockMatrixList) {
		blockmatrixlist_test();
	}

	//parallel matrix 
	if (testingParaMatrix) {
		para_matrix_test();
	}

	return 0;

}
