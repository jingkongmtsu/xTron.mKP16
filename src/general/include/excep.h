/**
 * \file   excep.h
 * \brief  exception handling
 * \author Fenglai Liu and Jing Kong
 *
 * the exceptions definition will be made available on both host and co-processor
 *
 */
#ifndef EXCEP_H
#define EXCEP_H

// begin the offload attribute region
//#pragma offload_attribute (push, target (mic))

#include "libgen.h"
#include <string>

//
// standardize the exceptions, and group them together
// each group has 200 exceptions
//

// 
// group 1: exception related to STL etc.
// number 1-200
//
#define EXCEPTION_VECTOR_OVERFLOW                1
#define EXCEPTION_VECTOR_DIMENSION_CONFLICT      2
#define EXCEPTION_DATA_IS_NOT_INITIALIZED        3
#define EXCEPTION_MEMORY_ALLOCATION_FAILED       4

// 
// group 2: exception related to file handling
// number 201-400
//
#define EXCEPTION_FILE_MISSING                   201
#define EXCEPTION_DIR_MISSING                    202
#define EXCEPTION_BIN_FILE_READ_ERROR            203
#define EXCEPTION_BIN_FILE_WRITE_ERROR           204
#define EXCEPTION_FILE_WRITE_FAIL                205
#define EXCEPTION_DATA_SECTION_NOT_FOUND         206
#define EXCEPTION_USE_DEFAULT_FILE_PATH          207
#define EXCEPTION_DIR_ALREADY_EXIST              208
#define EXCEPTION_IMPROPER_FILE_NAME             209

// 
// group 3: exception related to general folder
// number 401-600
//
#define EXCEPTION_SECTION_SEARCH_FAILED          401
#define EXCEPTION_INPUT_PARAMETER_INVALID        402
#define INVALID_FUNCTIONAL_NAME                  403
#define EXCEPTION_MODULE_SEARCH_FAILED           404
#define INVALID_FUNCTIONAL_COEFFICIENTS          405
#define EXCEPTION_CONVERT_TO_DOUBLE_FAILED       406
#define EXCEPTION_CONVERT_TO_INT_FAILED          407
#define EXCEPTION_FUNCTION_NOT_AVAILABLE         408
#define EXCEPTION_INTEGER_OVERFLOW               409
#define EXCEPTION_INTEGER_LENGTH_UNDEFINED       410
#define EXCEPTION_HISTDATAMAN_ERROR              411
#define EXCEPTION_IMPROPER_THREAD_NUM_SETTING    412
#define INVALID_FUNCTIONAL_METHOD_SWITCH         413
#define EXCEPTION_HARDWARE_CONFLIT               414
#define EXCEPTION_RESULT_IS_MEANINGLESS          415

//
// group 4: exception related to geom
// number 601-800
//
#define EXCEPTION_POLARIZABILITIES               601
#define INVALID_ATOMIC_NUMBER                    602
#define EXCEPTION_INVALID_NALPHA_NBETA           603
#define EXCEPTION_ATOMS_POS_OVERLAP              604
#define EXCEPTION_ILLEGAL_GEOM_SECTION_INDEX     605
#define EXCEPTION_ILLEGAL_SET_UNRESTRICTED       606
#define EXCEPTION_ZMATRIX_READ_IN_DATA_ERROR     607

//
// group 5: exception related to MATH
// number 801-1000
//
#define EXCEPTION_DIVIDE_BY_ZERO                      801
#define MACRO_IN_MATH_NOT_PROPERLY_DEFINED            802
#define EXCEPTION_CHAR_VARIABLE_INVALID               803
#define EXCEPTION_ARRAY_DIMENSION_NOT_MATCH           804
#define LAPACK_FUNCTION_FAILED                        805
#define EXCEPTION_MATRIX_DIMENSION_CONFLICT           806
#define EXCEPTION_EIGEN_VALUE_LESS_THAN_ZERO          807
#define EXCEPTION_LINEAR_DEPENDENCY_IN_MATRIX         808
#define EXCEPTION_TOO_MANY_TO_PRINT_MATRIX            809
#define EXCEPTION_INVALID_COLS_PRINT_MATRIX           810
#define EXCEPTION_MATRIX_RESET_WARNING                811
#define EXCEPTION_MATRIX_DIMENSION_OVERFLOW           812
#define EXCEPTION_BLOCK_MATRIX_ROW_COL_POS_CONFLICT   813
#define EXCEPTION_MATRIX_TRANSPOSE_ERROR              814
#define EXCEPTION_BLOCK_MATRIX_LIST_CONSTRUCTOR_ERROR 815
#define EXCEPTION_BLOCK_MATRIX_INVALID_MATRIX_PASS_IN 816
#define EXCEPTION_DIVIDE_BY_SMALL_NUMBER_BYPASS       817
#define EXCEPTION_LINEAR_DEPENDENCY_ERROR_IN_MATRIX   818

//
// group 6: exception related to shell
// number 1001-1200
//
#define EXCEPTION_SHELL_RADIUS_COEFF_TOO_SMALL        1001
#define EXCEPTION_SHELL_RADIUS_NOT_CONVERGE           1002
#define EXCEPTION_SHELL_BASIS_INDEX_TYPE_INVALID      1003
#define EXCEPTION_RAW_ATOM_SHELL_DATA_NOT_FOUND       1004
#define EXCEPTION_RAW_ATOM_SHELL_INVALID              1005
#define EXCEPTION_RAW_ATOM_SHELL_NGAU_INVALID         1006
#define EXCEPTION_RAW_ATOM_SHELL_SCALE_FAC_INVALID    1007
#define EXCEPTION_RAW_ATOM_SHELL_EXP_DATA_INVALID     1008
#define EXCEPTION_RAW_ATOM_SHELL_COEF_DATA_INVALID    1009
#define EXCEPTION_SHELL_DATA_FAIL_TO_BE_FOUND         1010
#define EXCEPTION_UNKNOWN_COMP_SHELL_TYPE             1011
#define EXCEPTION_BASIS_SET_INDEX_OUTOF_RANGE         1012
#define EXCEPTION_UNKNOWN_SHELL_SYMBOL                1013
#define EXCEPTION_ILLEGAL_AUX_SHELL_INDEX             1014
#define EXCEPTION_SHELL_NOT_EMPTY                     1015
#define EXCEPTION_INVALID_SHELL_DATA                  1016
#define EXCEPTION_SHELL_DIMENSION_CONFLICT            1017
#define EXCEPTION_RAW_SHELL_DATA_CONFLICT_MOL_DATA    1018
#define EXCEPTION_UNKNOWN_ATOM_TYPE_ATOM_SHELL_SIZE   1019
#define EXCEPTION_UNSUPPORTED_ANGULAR_MOMENTUM        1020
#define EXCEPTION_CPTRANS_DATA_NOT_AVAILABLE          1021
#define EXCEPTION_BLOCK_MATRIX_CPTRANS_ERROR          1022
#define EXCEPTION_INVALID_SHELL_DATA_FOR_COMPRESSING  1023
#define EXCEPTION_FAIL_LOCATE_BASIS_SET_NAME          1024

//
// group 7: exception related to gints
// number 1201-1400
//
#define EXCEPTION_GINTS_ILLEGAL_K_DIGESTION         1201
#define EXCEPTION_GINTS_SHELL_PAIR_ERROR            1202
#define EXCEPTION_GINTS_LCODE_MISSING_OR_INVALID    1203
#define EXCEPTION_GINTS_INVALID_INTEGRAL_JOB        1204
#define EXCEPTION_GINTS_FAIL_TO_FIND_MATRIX         1205
#define EXCEPTION_SIG_SHELL_PAIR_INT_BOUNDARY_ERROR 1206
#define EXCEPTION_DERIV_INFORMATION_VIOLATION       1207
#define EXCEPTION_GINTS_ILLEGAL_DERIV_DIGESTION     1208

//
// group 8: exception related to xcints
// number 1401-1600
//
#define EXCEPTION_XCINTS_INVALID_SIGBASIS_INDEX                   1401
#define EXCEPTION_XCINTS_INVALID_GRID_CHOICE                      1402
#define EXCEPTION_XCINTS_FAIL_TO_FIND_ATOMGRIDS                   1403
#define EXCEPTION_XCINTS_FAIL_TO_LOCATE_BATCH_ATOM_IN_SIGATOMLIST 1404
#define EXCEPTION_XCINTS_INVALID_XCVAR                            1405
#define EXCEPTION_XCINTS_INVALID_VAR_DERIV_ORDER                  1406
#define EXCEPTION_XCINTS_INVALID_FUNC_DERIV_ORDER                 1407
#define EXCEPTION_XCINTS_DFTMATRIX_ERROR                          1408
#define EXCEPTION_INVALID_XCINTS_CALCULATION                      1409
#define EXCEPTION_XCINTS_BATCHBASIS_ERROR                         1410
#define EXCEPTION_XCINTS_DEBUG_PRINT_DISABLED                     1411
#define EXCEPTION_XCINTS_TOO_MANY_BATCH_POINTS                    1412

//
// group 9: exception related to the result module
// number 1601-1800
//
#define EXCEPTION_DENSITY_MATRICES_FORM_ERROR       1601
#define EXCEPTION_FOCK_SOLVING_ERROR                1602
#define EXCEPTION_ATOM_DENSITY_MATRICES_ERROR       1603
#define EXCEPTION_NO_VIRTUAL_ORBITAL                1604
#define EXCEPTION_PMAXINFOR_ERROR                   1605
#define EXCEPTION_INVALID_CASE_FRAC_SPIN            1606
#define EXCEPTION_MO_ERROR                          1607

//
// group 10: exception related to the scf module
// number 1801-2000
//
#define EXCEPTION_HISTORY_FOCK_DIMENSION_ERROR      1801
#define EXCEPTION_SCFDIIS_ERROR                     1802
#define EXCEPTION_SCFEADIIS_ERROR                   1803
#define EXCEPTION_SCFEADIIS_WARNING                 1804
#define EXCEPTION_SCFCONV_WARNING                   1805
#define EXCEPTION_SCFCONV_ERROR                     1806
#define EXCEPTION_SCF_INPUT_INFOR_INVALID           1807
#define EXCEPTION_DIIS_CONTROLLER_ERROR             1808

namespace excep {

	using namespace std;

	/**
	 * \class exception for base class
	 */
	class Excep {

		protected:

			UInt    code;         ///< exception code
			string  moduleName;   ///< the class/namespace name where trigger the exception
			string  methodName;   ///< the function name where trigger the exception
			string  infor;        ///< more information for the exception

		public:

			///
			/// constructor
			///
			Excep(const string& name1, const string& name2, const UInt& code0,
					const string& infor0):code(code0),moduleName(name1),
			methodName(name2),infor(infor0) { };

			///
			/// destructor
			///
			~Excep() { };

			///
			/// get the information related to the code
			///
			string codeInfor() const;

			///
			/// Shall we terminate the program, or proceed to give an warning sign?
			///
			bool doTerminate() const;

			///
			/// print out the information for the base type of exception
			///
			void print() const;

			///
			/// what is the content for exception?
			///
			void what() const { print(); };
	};


	/////////////////////////////////////////////////////////////
	/// this is the function to handle the exceptions defined  //
	/////////////////////////////////////////////////////////////
	
	///
	/// handle the base type of exception
	///
	void handleExcep(const Excep& excep);

}

// end the offload attribute region
//#pragma offload_attribute (pop)

#endif

