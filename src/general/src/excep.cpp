/**
 * cpp file for exception handling
 * \author fenglai liu and jing kong
 *
 * the code should be able to run on both cpu and co-processor
 */
#include <iostream>
#include <cstdlib>
#include "excep.h"
using namespace excep;

// begin the offload attribute region
//#pragma offload_attribute (push, target (mic))

string Excep::codeInfor() const 
{
	string status;
	switch(code) 
	{
		// group 1
		case EXCEPTION_VECTOR_OVERFLOW: 
			status = "vector overflow";
			break;	
		case EXCEPTION_VECTOR_DIMENSION_CONFLICT:
			status = "the dimension for two working vectors are not matching each other";
			break;	
		case EXCEPTION_DATA_IS_NOT_INITIALIZED:
			status = "the data inside vector/matrix is not initialized";
			break;	
		case EXCEPTION_MEMORY_ALLOCATION_FAILED:
			status = "to allocate memory for the data is failed";
			break;	

		// group 2
		case EXCEPTION_FILE_MISSING: 
			status = "file is missing";
			break;	
		case EXCEPTION_DIR_MISSING: 
			status = "directory is missing";
			break;	
		case EXCEPTION_BIN_FILE_READ_ERROR: 
			status = "binary file reading is failed";
			break;	
		case EXCEPTION_BIN_FILE_WRITE_ERROR: 
			status = "binary file writing is failed";
			break;	
		case EXCEPTION_FILE_WRITE_FAIL: 
			status = "fail to open a file for writing";
			break;	
		case EXCEPTION_DATA_SECTION_NOT_FOUND:
			status = "fail to find the request data section in the given file";
			break;	
		case EXCEPTION_USE_DEFAULT_FILE_PATH:
			status = "file path is not defined, so we use the default one";
			break;	
		case EXCEPTION_DIR_ALREADY_EXIST:              
			status = "the directory can not be established, it's already there";
			break;	
		case EXCEPTION_IMPROPER_FILE_NAME:
			status = "the file name is invalid, please check the file name you gave";
			break;	

		// group 3
		case EXCEPTION_SECTION_SEARCH_FAILED: 
			status = "section search is failed";
			break;	
		case EXCEPTION_INPUT_PARAMETER_INVALID: 
			status = "parameters defined in the input file is invalid";
			break;	
		case INVALID_FUNCTIONAL_NAME: 
			status = "functional name is invalid";
			break;	
		case EXCEPTION_MODULE_SEARCH_FAILED: 
			status = "the given class/module name searching in input file is failed";
			break;	
		case INVALID_FUNCTIONAL_COEFFICIENTS: 
			status = "the functional coefficients array conflicts with other definitions";
			break;	
		case EXCEPTION_CONVERT_TO_DOUBLE_FAILED:
			status = "fail to convert an value to Double type";
			break;	
		case EXCEPTION_CONVERT_TO_INT_FAILED:
			status = "fail to convert an value to Int(UInt) type";
			break;	
		case EXCEPTION_FUNCTION_NOT_AVAILABLE:
			status = "the required calculation is not available yet";
			break;	
		case EXCEPTION_INTEGER_OVERFLOW:
			status = "the integer is overflowing";
			break;	
		case EXCEPTION_INTEGER_LENGTH_UNDEFINED:
			status = "the sizeof(Int) returns undefined integer length";
			break;	
		case EXCEPTION_HISTDATAMAN_ERROR:   
			status = "error arises in HistDataMan class";
			break;	
		case EXCEPTION_IMPROPER_THREAD_NUM_SETTING: 
			status = "the threads number is improperly set, a reset is needed";
			break;	
		case INVALID_FUNCTIONAL_METHOD_SWITCH:
			status = "the method switch definition in the functional information is invalid (see xcfunc.conf)";
			break;	
		case EXCEPTION_HARDWARE_CONFLIT:
			status = "the use of hardware conficts with each other";
			break;	
		case EXCEPTION_RESULT_IS_MEANINGLESS: 
			status = "the function result is meaningless, please check the implementation and the input parameters";
			break;	

		// group 4
		case EXCEPTION_POLARIZABILITIES:
			status = "the given atom's POLARIZABILITY data is not available";
			break;
		case INVALID_ATOMIC_NUMBER:
			status = "the given atom's atomic number is invalid";
			break;
		case EXCEPTION_INVALID_NALPHA_NBETA:
			status = "the result number of alpha/beta electrons are invalid";
			break;
		case EXCEPTION_ATOMS_POS_OVERLAP:
			status = "two atoms superpose with each other";
			break;
		case EXCEPTION_ILLEGAL_GEOM_SECTION_INDEX:
			status = "section number for this geometry is invalid";
			break;	
		case EXCEPTION_ILLEGAL_SET_UNRESTRICTED: 
			status = "it's valid to set the unrestricted status in the molecule section";
			break;	
		case EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR:
			status = "the given zmatrix data from input file has an error";
			break;	

		// group 5
		case EXCEPTION_DIVIDE_BY_ZERO:
			status = "NaN is caused because dividing by zero";
			break;
		case MACRO_IN_MATH_NOT_PROPERLY_DEFINED:
			status = "macro in math is not properly defined, failed to call vendor's function";
			break;
		case EXCEPTION_CHAR_VARIABLE_INVALID:
			status = "to call vendor's function, the input char parameter is invalid";
			break;
		case EXCEPTION_ARRAY_DIMENSION_NOT_MATCH:
			status = "the input array dimension data is invalid to perform computation";
			break;
		case LAPACK_FUNCTION_FAILED:
			status = "the lapack function returns error message";
			break;
		case EXCEPTION_MATRIX_DIMENSION_CONFLICT:
			status = "dimensions for matrix are in conflict so that operation can not proceed";
			break;
		case EXCEPTION_EIGEN_VALUE_LESS_THAN_ZERO:
			status = "the calculated eigen value is less than 0";
			break;
		case EXCEPTION_LINEAR_DEPENDENCY_IN_MATRIX:
			status = "linear dependency observed in the matrix";
			break;
		case EXCEPTION_TOO_MANY_TO_PRINT_MATRIX:
			status = "too many data to be printed out in matrix, so no printing available";
			break;
		case EXCEPTION_INVALID_COLS_PRINT_MATRIX:
			status = "unmatch column number provided in printing matrix data";
			break;
		case EXCEPTION_MATRIX_RESET_WARNING:
			status = "something abnormal in matrix reset, rellocation of memory is needed";
			break;
		case EXCEPTION_MATRIX_DIMENSION_OVERFLOW:
			status = "the position to access matrix element is overflow";
			break;
		case EXCEPTION_BLOCK_MATRIX_ROW_COL_POS_CONFLICT:
			status = "row and column positions between block matrices do not match";
			break;
		case EXCEPTION_MATRIX_TRANSPOSE_ERROR:
			status = "something wrong happens in transposing the given matrix";
			break;
		case EXCEPTION_BLOCK_MATRIX_LIST_CONSTRUCTOR_ERROR:
			status = "something wrong happens in constructing block matrix list";
			break;
		case EXCEPTION_BLOCK_MATRIX_INVALID_MATRIX_PASS_IN:
			status = "invalid type of matrix passed in the function";
			break;
		case EXCEPTION_DIVIDE_BY_SMALL_NUMBER_BYPASS:
			status = "the value is divided by a small number, by we bypass the calculation";
			break;
		case EXCEPTION_LINEAR_DEPENDENCY_ERROR_IN_MATRIX:
			status = "the matrix is linear depedent so that the given operation on matrix can not perform";
			break;

		// group 6:
		case EXCEPTION_SHELL_RADIUS_COEFF_TOO_SMALL:
			status = "the coeff. is too small compare with input threshold value";
			break;
		case EXCEPTION_SHELL_RADIUS_NOT_CONVERGE: 
			status = "shell radius calculation is not converging";
			break;
		case EXCEPTION_SHELL_BASIS_INDEX_TYPE_INVALID:
			status = "input basis set index type is invalid";
			break;
		case EXCEPTION_RAW_ATOM_SHELL_DATA_NOT_FOUND:
			status = "raw atom shell data is not found";
			break;
		case EXCEPTION_RAW_ATOM_SHELL_INVALID:
			status = "the shell data read in raw atom shell data is invalid";
			break;
		case EXCEPTION_RAW_ATOM_SHELL_NGAU_INVALID:
			status = "the ngau read in raw atom shell data is invalid";
			break;
		case EXCEPTION_RAW_ATOM_SHELL_SCALE_FAC_INVALID:
			status = "the scale factor read in raw atom shell data is invalid";
			break;
		case EXCEPTION_RAW_ATOM_SHELL_EXP_DATA_INVALID:
			status = "the exponential factor read in raw atom shell data is invalid";
			break;
		case EXCEPTION_RAW_ATOM_SHELL_COEF_DATA_INVALID:
			status = "the coefficients read in raw atom shell data is invalid";
			break;
		case EXCEPTION_SHELL_DATA_FAIL_TO_BE_FOUND:
			status = "failed to find the shell data section in input file";
			break;
		case EXCEPTION_UNKNOWN_COMP_SHELL_TYPE:
			status = "we do not support this kind of composite shell";
			break;
		case EXCEPTION_BASIS_SET_INDEX_OUTOF_RANGE:
			status = "the given basis set index is out of range";
			break;
		case EXCEPTION_UNKNOWN_SHELL_SYMBOL:
			status = "the shell symbol is unknown";
			break;
		case EXCEPTION_ILLEGAL_AUX_SHELL_INDEX:
			status = "the aux shell index is illegal";
			break;
		case EXCEPTION_SHELL_NOT_EMPTY:
			status = "shell data is not empty";
			break;
		case EXCEPTION_INVALID_SHELL_DATA:
			status = "it seems the shell data is internally conflicted with each other";
			break;
		case EXCEPTION_SHELL_DIMENSION_CONFLICT:
			status = "the shell dimension is in conflict with input data";
			break;
		case EXCEPTION_RAW_SHELL_DATA_CONFLICT_MOL_DATA:
			status = "the raw shell data does not match the molecule data";
			break;
		case EXCEPTION_UNKNOWN_ATOM_TYPE_ATOM_SHELL_SIZE: 
			status = "the given atom type is missing in the atom shell size class";
			break;
		case EXCEPTION_UNSUPPORTED_ANGULAR_MOMENTUM:
			status = "the angular momentum is not supported in this program";
			break;
		case EXCEPTION_CPTRANS_DATA_NOT_AVAILABLE:
			status = "in processing the C2P/P2C transformation etc. given L is not available";
			break;
		case EXCEPTION_BLOCK_MATRIX_CPTRANS_ERROR:
			status = "in processing the C2P/P2C transformation error occurs";
			break;
		case EXCEPTION_INVALID_SHELL_DATA_FOR_COMPRESSING:
			status = "the given shell data can not be compressed";
			break;
		case EXCEPTION_FAIL_LOCATE_BASIS_SET_NAME:
			status = "we failed to locate the basis set name/file_name from the configuration file basis.conf";
			break;

		// group 7:
		case EXCEPTION_GINTS_ILLEGAL_K_DIGESTION:
			status = "fatal error happens in the exchange digestion process";
			break;
		case EXCEPTION_GINTS_SHELL_PAIR_ERROR:
			status = "fatal error happens in forming the shell pair data";
			break;
		case EXCEPTION_GINTS_LCODE_MISSING_OR_INVALID:
			status = "the L code is invalid or it can not be processed here";
			break;
		case EXCEPTION_GINTS_INVALID_INTEGRAL_JOB:
			status = "invalid given integral job in gints module";
			break;
		case EXCEPTION_GINTS_FAIL_TO_FIND_MATRIX: 
			status = "fail to find the given matrix data from storage in form of memory/disk";
			break;
		case EXCEPTION_SIG_SHELL_PAIR_INT_BOUNDARY_ERROR:
			status = "integral boundary is not correctly calculated in sigfinificant shell pair infor class";
			break;
		case EXCEPTION_DERIV_INFORMATION_VIOLATION:
			status = "integral derivatives information is in violation status, please double check it";
			break;
		case EXCEPTION_GINTS_ILLEGAL_DERIV_DIGESTION:
			status = "something wrong happens in digesting the integral derivatives";
			break;

		// group 8:
		case EXCEPTION_XCINTS_INVALID_SIGBASIS_INDEX:
			status = "the derived significant basis-set/shell index is invalid";
			break;
		case EXCEPTION_XCINTS_INVALID_GRID_CHOICE:
			status = "the grid choice in xcints module is invalid";
			break;
		case EXCEPTION_XCINTS_FAIL_TO_FIND_ATOMGRIDS:
			status = "fail to find the atom grid data according to the input atomic number";
			break;
		case EXCEPTION_XCINTS_FAIL_TO_LOCATE_BATCH_ATOM_IN_SIGATOMLIST:
			status = "fail to find the mother atom in the given sig atom list";
			break;
		case EXCEPTION_XCINTS_INVALID_XCVAR:
			status = "the DFT variable passing in or generated is invalid";
			break;
		case EXCEPTION_XCINTS_INVALID_VAR_DERIV_ORDER:
			status = "the DFT variable derivative order is invalid";
			break;
		case EXCEPTION_XCINTS_INVALID_FUNC_DERIV_ORDER:
			status = "the DFT functional derivative order is invalid";
			break;
		case EXCEPTION_XCINTS_DFTMATRIX_ERROR:           
			status = "the dft matrix gets error";
			break;
		case EXCEPTION_INVALID_XCINTS_CALCULATION: 
			status = "error arises in xcints(top module) calculation";
			break;
		case EXCEPTION_XCINTS_BATCHBASIS_ERROR:
			status = "error arises in batch basis calculation";
			break;
		case EXCEPTION_XCINTS_DEBUG_PRINT_DISABLED:
			status = "the debug print is disabled in xcints module";
			break;
		case EXCEPTION_XCINTS_TOO_MANY_BATCH_POINTS:
			status = "there are too many batch points in the given xcints calculation";
			break;

		// group 9:
		case EXCEPTION_DENSITY_MATRICES_FORM_ERROR:
			status = "error arises in forming the density matrix";
			break;
		case EXCEPTION_FOCK_SOLVING_ERROR:
			status = "error arises in solving Fock matrix to derive mo";
			break;
		case EXCEPTION_ATOM_DENSITY_MATRICES_ERROR:   
			status = "error arises in forming the SCF density matrices/MO etc. data for free atom";
			break;
		case EXCEPTION_NO_VIRTUAL_ORBITAL:
			status = "no virtual orbitals information, this causes error for the given function";
			break;
		case EXCEPTION_PMAXINFOR_ERROR: 
			status = "error arises in forming density matrix infor data";
			break;
		case EXCEPTION_INVALID_CASE_FRAC_SPIN:
			status = "there's something wrong in the fractional spin calculation";
			break;
		case EXCEPTION_MO_ERROR: 
			status = "error arises in mo class";
			break;

		// group 10:
		case EXCEPTION_HISTORY_FOCK_DIMENSION_ERROR:
			status = "fock matrix derived from history has conflict in dimension with current fock matrix";
			break;
		case EXCEPTION_SCFDIIS_ERROR:
			status = "error arises in DIIS of scf procedure";
			break;
		case EXCEPTION_SCFEADIIS_ERROR:
			status = "error arises in EDIIS/ADIIS of scf procedure";
			break;
		case EXCEPTION_SCFEADIIS_WARNING:
			status = "warning message arises in EDIIS/ADIIS of scf procedure";
			break;
		case EXCEPTION_SCFCONV_WARNING:
			status = "warning message arises during scf convergence progress";
			break;
		case EXCEPTION_SCFCONV_ERROR:    
			status = "error arises in SCFConv class during scf";
			break;
		case EXCEPTION_SCF_INPUT_INFOR_INVALID:
			status = "in the scf process the input information is invalid";
			break;
		case EXCEPTION_DIIS_CONTROLLER_ERROR:
			status = "error arises in DIIS/ADIIS/EDIIS etc. DIIS like algorithm application";
			break;
		 	
		// default
		default:
			status = "invalid exception code passed in";
			break;
	}
	return status;
}

bool Excep::doTerminate() const
{
	// warning sign
	if (code == EXCEPTION_POLARIZABILITIES || 
			code == EXCEPTION_IMPROPER_THREAD_NUM_SETTING ||
			code == EXCEPTION_LINEAR_DEPENDENCY_IN_MATRIX || 
			code == EXCEPTION_TOO_MANY_TO_PRINT_MATRIX    ||
			code == EXCEPTION_INVALID_COLS_PRINT_MATRIX   ||
			code == EXCEPTION_XCINTS_DEBUG_PRINT_DISABLED ||
			code == EXCEPTION_DIVIDE_BY_SMALL_NUMBER_BYPASS || 
			code == EXCEPTION_SCFEADIIS_WARNING ||
			code == EXCEPTION_SCFCONV_WARNING ||
			code == EXCEPTION_MATRIX_RESET_WARNING || 
			code == EXCEPTION_USE_DEFAULT_FILE_PATH) {
		return false;
	}
	return true;
}

void Excep::print() const
{
	string status = codeInfor();
	cout << "class/namespace   : " << moduleName << endl;
	cout << "function          : " << methodName << endl;
	cout << "exception status  : " << status     << endl;
	cout << "detailed infor    : " << infor      << endl;
	cout << endl;
}

void excep::handleExcep(const Excep& excep)
{
	// is it a serious error?
	bool serious = excep.doTerminate();
	if (serious) {
		cout << "************************" << endl;
		cout << "fatal error occurs: " << endl;
		cout << "************************" << endl;
	}else{
		cout << "*******************************" << endl;
		cout << "this is just a warning sign: " << endl;
		cout << "*******************************" << endl;
	}

	// now print out the infor
	// and may terminate program in abnormal way
	excep.what();
	if (serious) exit(1);
}

// end the offload attribute region
//#pragma offload_attribute (pop)

