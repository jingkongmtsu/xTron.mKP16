///
/// cpp file corresponding to xcintsinfor.h
///
#include<iostream>
#include "textread.h"
#include "excep.h"
#include "parameterparsing.h"
#include "globalinfor.h"
#include "molecule.h"
#include "gridinfor.h"
#include "xcintsinfor.h"
using namespace textread;
using namespace excep;
using namespace parameterparsing;
using namespace molecule;
using namespace globalinfor;
using namespace gridinfor;
using namespace xcintsinfor;
using namespace std;

XCIntsInfor::XCIntsInfor(const GlobalInfor& infor, const Molecule& mol):doPartitionWeights(true),
	insigGridsCheck(true),doOddElectron(false),doXCEnergyProfing(false),
	pruneMethod(NON_PRUNE_GRID),stepFunctionType(BECKE_ORI_STEP_FUNCTION),gridChoice(STANDARD_GRID),
	nRad(50),nAng(194),nPtPerCPUThread(500),sigGridCheckThresh(1.0E-12),
	threshold(1.0E-10),funcTol(1.0E-12),nSpin(mol.getNSpin()),
	isCloseShell(mol.isCloseShell()),oddElecPar1(0.5E0),oddElecPar2(0.5E0),
	xcvar_debug(0),molegrids_debug(0),sigatombasis_debug(0),batchgrid_debug(0), 
	batchbasis_debug(0),batchvar_debug(0),batchfunc_debug(0),batchxcmtrx_debug(0),
	espints_debug(0),exrho_debug(0),printTimingForESP(false)
{
	nDen = mol.getNBeta() > 0 ? 2 : nSpin;

	// deal with input key words
	UInt section = mol.getSec();
	string input = infor.getInputFile();
	ParameterParsing pp(input,"xcints",section);
	if (! pp.hasAnyParameters()) return;

	// now let's read in the key of xcints
	UInt nKey = pp.getKeyNumber();
	WordConvert w;
	for(UInt iKey=0; iKey<nKey; iKey++) {

		// get the key name
		string key = pp.getKey(iKey);

		// now begin the read of the keywords
		if (w.compare("THRESHOLD",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  threshold
			/// group    xcints
			/// values   double precision value
			/// 
			/// this keyword defines the boundary value to determine the 
			/// significant basis functions for the numerical calculations.
			/// For example, if the threshold value is 1.0E-12, then if the
			/// basis function's contribution to the batch grid points
			/// is less than this boundary value, it will be neglected in 
			/// the following calculations.
			///
			/// the default value for this keyword is 1.0E-10.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value. not double.";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			threshold = tmp;

		}else if (w.compare("ODD_ELECTRON_PARAMETER",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  odd_electron_parameter 
			/// group    xcints
			/// values   two double precision value
			/// 
			/// this keyword defines the parameter values for performing the 
			/// "odd electron population", please see the paper below:
			///
			/// "Analyzing effects of strong electron correlation within 
			///  Kohn-Sham density-functional theory"
			///  Emil Proynov, Fenglai Liu, and Jing Kong
			///  Phys. Rev. A 88, 032510
			///
			///  for the input, you need to assign two double precision numbers,
			///  the first one is the a1 parameter, the second one is the a2 
			///  parameter; more detail please refer to the paper and the functional
			///  code becke05.f in xcfunc directory, see the function Becke05_odd_electron.
			///
			///  the default value for this keyword is (0.5, 0.5).
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			LineParse l(value);
			if (l.getNPieces() == 2) {
				for(UInt i=0; i<2; i++) {
					Double tmp = 0.0;
					string tmpValue = l.findValue(i);
					if (!w.toDouble(tmpValue,tmp)) {
						string infor = "the odd electron population parameter is invalid, not a double value";
						Excep excep("XCIntsInfor","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
						handleExcep(excep);
					}
					if (i == 0) {
						oddElecPar1 = tmp;
					}else{
						oddElecPar2 = tmp;
					}
				}
			}else{
				string infor = "the value for odd_electron_parameter is invalid, you need to set two parameters";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("FUNCTIONAL_THRESHOLD",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  functional_threshold
			/// group    xcints
			/// values   double precision value
			/// 
			/// this keyword defines the boundary value to determine the 
			/// significance of the electron density, gradient density etc.
			/// during the functional and functional derivatives calculations.
			/// For example, if the threshold value is 1.0E-12, then if the 
			/// electron density value on this batch grid all less than this 
			/// boundary value, it will be neglected in the functional 
			/// calculations.
			///
			/// the default value for this keyword is 1.0E-12.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process functional_threshold value. not double.";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			funcTol = tmp;

		}else if (w.compare("INSIGNIFICANT_GRID_CHECK_THRESHOLD",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  insignificant_grid_check_threshold
			/// group    xcints
			/// values   double precision value
			/// 
			/// this keyword defines the boundary value to determine 
			/// whether the grid point is significant for the numerical
			/// calculations. 
			///
			/// in the function AtomGrids::formQuadraturePtsWts, when we
			/// form the quadrature grid using the Eular-Maclaurine formula;
			/// we will check the weights of the radial points to see whether
			/// it is less than the threshold value. If so, we will consider
			/// this grid point to be insignificant and neglect it in the 
			/// following calculations.
			///
			/// the default value for this keyword is 1.0E-12.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value for insignificant grid check. not double.";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			sigGridCheckThresh = tmp;

		}else if (w.compare("INSIGNIFICANT_GRIDS_CHECK",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  insignificant_grids_check
			/// group    xcints
			/// values   true or false (T or F also works)
			/// 
			/// this keyword determines that whether we performs the 
			/// insignificant grid check inside the function
			/// function AtomGrids::formQuadraturePtsWts. For more 
			/// details please see the explanation of the keyword
			/// insignificant_grid_check_threshold.
			///
			/// disable this keyword may raise up a bit of calculation
			/// cost, especially for the large grid set. Usually
			/// to perform such insignificant grid check should not 
			/// affect precision, but you can disable it in terms of 
			/// debugging purpose.
			///
			/// the default value for this keyword is true.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				insigGridsCheck = true;
			}else if (value == "FALSE" || value == "F") {
				insigGridsCheck = false;
			}else{
				string infor = "only true or false allowed to process insignificant_grids_check";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("DO_ODD_ELECTRON",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  do_odd_electron
			/// group    xcints
			/// values   true or false (T or F also works)
			/// 
			/// this keyword determines that whether you do
			/// Odd electron Population. Please see the paper
			/// for more details.
			/// 
			/// "Analyzing effects of strong electron correlation within 
			///  Kohn-Sham density-functional theory"
			///  Emil Proynov, Fenglai Liu, and Jing Kong
			///  Phys. Rev. A 88, 032510
			///
			/// If you choose to do the odd electron population, it 
			/// actually will be performed in the post-SCF section.
			/// Because the odd electron population requires the DFT
			/// variables like Laplacian of electron density and 
			/// exchange energy density etc. and such variables may not
			/// exist in your DFT application (for example, B3LYP SCF
			/// calculation), therefore it's not performed in the SCF
			/// but after the SCF is finished.
			///
			/// In this sense the odd electron population is like one
			/// cycle XC functional calcualation.
			///
			/// the default value for this keyword is false.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				doOddElectron = true;
			}else if (value == "FALSE" || value == "F") {
				doOddElectron = false;
			}else{
				string infor = "only true or false allowed to process do_odd_electron";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("DO_XC_ENERGY_PROFILE",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  do_xc_energy_profile
			/// group    xcints
			/// values   true or false (T or F also works)
			/// 
			/// this keyword determines that whether you do
			/// XC energy decomposition, so that to partition
			/// XC energy into components. This skill is useful
			/// for XC functional developing.
			///
			/// currently only KP14 functional supports this 
			/// function. In the future we may add all of 
			/// functional so that the xc energy decomposition
			/// can be performed to any xc functional.
			/// 
			/// the default value for this keyword is false.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				doXCEnergyProfing	= true;
			}else if (value == "FALSE" || value == "F") {
				doXCEnergyProfing	= false;
			}else{
				string infor = "only true or false allowed to process do_xc_energy_profile";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("DO_PARTITION_WEIGHTS",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  do_partition_weights
			/// group    xcints
			/// values   true or false (T or F also works)
			/// 
			/// this keyword determines that whether we have the
			/// partition weights added to the grid system.
			/// 
			/// the partition weights is used to partition the 
			/// grid system into atom-centered grid system. 
			/// For example, the Becke weights, see the paper:
			///
			/// "A multicenter numerical integration scheme for 
			/// polyatomic molecules",
			/// The Journal of Chemical Physics, 1988, 88, 2547.
			///
			/// the default value for this keyword is true.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				doPartitionWeights = true;
			}else if (value == "FALSE" || value == "F") {
				doPartitionWeights = false;
			}else{
				string infor = "only true or false allowed to process do_partition_weights";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("GRID_POINTS",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  grid_points
			/// group    xcints
			/// values   pruned grid name, or two numbers,
			///          the first one is the radial grid points,
			///          the second one is the angular grid points
			/// 
			/// this keyword determines that grid set chosen for
			/// the XC part numerical calculations. There are two
			/// choices for selecting the grid set.
			///
			/// the first choice is to use the pruned grid set,
			/// for example; the SG1 grid etc. Currently the allowed
			/// pruned grid set are:
			/// SG1 : SG1 grid set
			/// BAKER: Baker grid set provided by PQS program
			/// 
			/// the second choice is to use the un-pruned grids.
			/// Here you need to specify the number of radial points
			/// and the number of angular points. Radial points are
			/// generated through the Eular-Maclaurine formula,
			/// angular points are generated through Lebedev formula.
			/// Therefore, you must assign two numbers for using the 
			/// un-pruned grids. For example, 
			///
			/// grid_points  =  50  194
			///
			/// it means you use 50 radial points, 194 angular points.
			///
			/// For the Lebedev points, you can only select the values
			/// below: 
			/// 6,    14,   26,   38,   50,   74,   86,   110,  
			/// 146,  170,  194,  230,  266,  302,  350,  434, 
			/// 590,  770,  974,  1202, 1454, 1730, 2030, 2354, 
			/// 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
			/// other values are not allowed.
			///
			/// the default value for this keyword  (50, 194).
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			LineParse l(value);
			if (l.getNPieces() == 1) {
				string tmpValue = l.findValue(0);
				w.capitalize(tmpValue);
				if (tmpValue == "SG1") {
					pruneMethod = SG1_GRID;
				}else if (tmpValue == "BAKER") {
					pruneMethod = BAKER_GRID;
				}else{
					string infor = "the provided prune grid method is invalid";
					Excep excep("XCIntsInfor","constructor",EXCEPTION_XCINTS_INVALID_GRID_CHOICE,infor);
					handleExcep(excep);
				}
			}else if (l.getNPieces() == 2) {
				for(UInt i=0; i<2; i++) {
					UInt tmp = 0;
					string tmpValue = l.findValue(i);
					if (!w.toUInt(tmpValue,tmp)) {
						string infor = "the number of nAng/nRad for un-pruned grid is invalid, not an integer";
						Excep excep("XCIntsInfor","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
						handleExcep(excep);
					}
					if (i == 0) {
						nRad = tmp;
					}else{
						nAng = tmp;
					}
				}
			}else{
				string infor = "the value for grid_points is invalid";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_XCINTS_INVALID_GRID_CHOICE,infor);
				handleExcep(excep);
			}

		}else if (w.compare("npt_cpu_thread",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  npt_cpu_thread
			/// group    xcints
			/// values   integer value, larger than 0
			/// 
			/// for multi-threading program performing on CPU, this keyword
			/// defines the number of grid points the user wish to apply on
			/// each CPU thread. 
			///
			/// In the XC module, each CPU thread has its own work task and
			/// finally we gather all of result together to form the final
			/// result. The value sets here is closely related to the "Batch
			/// Size" set in the atomgrids.cpp. For the working on CPU, because
			/// each CPU thread work on a single batch hence the final batch
			/// size determined dynamically over there is close to the 
			/// npt_cpu_thread.
			///
			/// the default value for this keyword is 500.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "the npt_cpu_thread setting is invalid, not an integer";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			nPtPerCPUThread = tmp;

		}else if (w.compare("GRID_QUALITY",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  GRID_QUALITY
			/// group    xcints
			/// values   three choices, STANDARD, COARSE, or FINE.
			/// 
			/// If you choose the Baker pruned grid, you have another 
			/// choice to determine what's the level of the Baker grid
			/// you want to use. "FINE" is the highest level, which uses
			/// the most number of grid points; default is "STANDARD";
			/// also you can use the "COARSE" grid.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "STANDARD") {
				gridChoice = STANDARD_GRID;
			}else if (value == "COARSE") {
				gridChoice = COARSE_GRID;
			}else if (value == "FINE") {
				gridChoice = FINE_GRID;
			}else{
				string infor = "the grid quality setting is invalid, you can use only STANDARD, COARSE or FINE";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("STEP_FUNCTION_TYPE",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  step_function_type
			/// group    xcints
			/// values   two choices, BECKE or STRATMANN
			/// 
			/// here the keywords select the step function used in 
			/// forming the partition grids. Default is "BECKE".
			///
			/// the Becke type of step function can see the paper:
			///
			/// Becke A D. 
			/// "A multicenter numerical integration scheme for polyatomic molecules" 
			/// The Journal of Chemical Physics, 1988, 88, 2547.
			///
			/// the STRATMANN type of step function see the paper:
			///
			/// R. Eric Stratmann and Gustavo E. Scuseria and Michael J. Frisch
			/// "Achieving linear scaling in exchange-correlation density functional 
			/// quadratures"
			/// Chemical Physics Letters, 257, 213 - 223, 1996
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "BECKE") {
				stepFunctionType = BECKE_ORI_STEP_FUNCTION;
			}else if (value == "STRATMANN") {
				stepFunctionType = STRATMANN_STEP_FUNCTION;
			}else{
				string infor = "the step function setting is invalid, only BECKE or STRATMANN are allowed";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else if (w.compare("xcvar",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  xcvar
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			xcvar_debug = 1;

		}else if (w.compare("batchxcmtrx",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  batchxcmtrx
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			batchxcmtrx_debug = 1;

		}else if (w.compare("batchvar",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  batchvar
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			batchvar_debug = 1;

		}else if (w.compare("espints",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  espints
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			espints_debug = 1;

		}else if (w.compare("exrho",key)) {
			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  exrho
			/// group    xcints
			/// values   in the debug keyword add exrho will make
			///          the XC code to perform debug function 
			///          to calculate exrho. This is for debugging purpose,
			///          and if you choose to debug exrho then no XC matrix etc.
			///          will be calculated. Default is not to do this debug function.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			exrho_debug = 1;

		}else if (w.compare("batchfunc",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  batchfunc
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			batchfunc_debug = 1;

		}else if (w.compare("batchgrid",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  batchgrid
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			batchgrid_debug = 1;

		}else if (w.compare("batchbasis",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  batchbasis
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			batchbasis_debug = 1;

		}else if (w.compare("molegrids",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  molegrids
			/// group    xcints
			/// values   for debugging the molegrids class you 
			/// must assign a debugging level. the level can 
			/// be 3 integers, 0, 1 or 2. 2 is the highest level
			/// and will print out everything we got.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			UInt tmp = 0;
			if (!w.toUInt(value,tmp)) {
				string infor = "molegrids_debug given is invalid, not an integer";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			if (tmp > 2) {
				string infor = "molegrids_debug level option could be only 0, 1 or 2. "
					"2 is the highest level and give you all details";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}
			molegrids_debug = tmp;

		}else if (w.compare("sigatombasis",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  sigatombasis
			/// group    xcints
			/// values   in the debug keyword add this class name
			///          and it will print out debug information 
			///          for the class. Default is no debug information
			///          print out.
			///
			/// it's worthy to note that that only under the serial
			/// mode (no multi-threading) you can print out the debug
			/// information.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			sigatombasis_debug = 1;

		}else if (w.compare("print_timing_espints",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  print_timing_espints
			/// group    xcints
			/// values   true or false (T or F also works)
			/// 
			/// this keyword is a debug option, whether you
			/// want to print out the timing data for the 
			/// espints part calculation.
			///
			/// the default value for this keyword is false.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "TRUE" || value == "T") {
				printTimingForESP = true;
			}else if (value == "FALSE" || value == "F") {
				printTimingForESP = false;
			}else{
				string infor = "only true or false allowed to process print_timing_espints";
				Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
				handleExcep(excep);
			}

		}else {
			cout << "the keyword for the xcints part is " << key << endl;
			string infor = "this keyword can not be identified";
			Excep excep("XCIntsInfor","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
			handleExcep(excep);
		}
	}
}

void XCIntsInfor::debugPrint() const
{
	cout << "*******************************************" << endl;
	cout << "*    DEBUG PRINT FOR XCIntsInfor CLASS     " << endl;
	cout << "*******************************************" << endl;
	cout << "nSpin is                                " << nSpin << endl;
	cout << "is close shell?                         " << isCloseShell << endl;
	cout << "do Partition weights?                   " << doPartitionWeights << endl;
	cout << "do insignificant grid check?            " << insigGridsCheck << endl;
	cout << "do ODD electron calculation?            " << doOddElectron << endl;
	cout << "threshold for insignificant grid check  " << sigGridCheckThresh << endl;
	cout << "functional threshold value              " << funcTol << endl;
	cout << "xcints threshold value                  " << threshold << endl;
	cout << "prune grid method                       " << pruneMethod << endl;
	cout << "grid choice                             " << gridChoice << endl;
	cout << "default nRad is                         " << nRad << endl;
	cout << "default nAng is                         " << nAng << endl;
	cout << "step function in partition weights type " << stepFunctionType << endl;
	cout << "number of points per each CPU thread:   " << nPtPerCPUThread << endl;
	if (doOddElectron) {
		cout << "ODD electron calculation parameter 1:   " << oddElecPar1 << endl;
		cout << "ODD electron calculation parameter 2:   " << oddElecPar2 << endl;
	}
}

////////////////////////////////////////////////////////////////////////////
//  ####         functions for class XCIntJobInfor                        //
////////////////////////////////////////////////////////////////////////////
XCIntJobInfor::XCIntJobInfor(const GlobalInfor& globInfor0, const GIntsInfor& ginfor0, 
		const XCIntsInfor& xcinfor, const UInt& job, const UInt& jobOrder):XCIntsInfor(xcinfor),
	globInfor(globInfor0),ginfor(ginfor0),jobName(job),order(jobOrder),
	useMultThreads(globInfor.useMultiThreads()),nCPUThreads(globInfor.getNCPUThreads()),
	initBatchSize(nPtPerCPUThread) 
{ }

bool XCIntJobInfor::doDebug(const string& name) const
{
	// for the debugging module, the exrho is an exception
	// because this choice is not for printing something,
	// but just debug exrho calculation
	// so no matter whether we use multi-threads we will
	// do debug for this choice
	if (name == "exrho") {
		if (exrho_debug > 0) {
			return true;
		}else{
			return false;
		}
	}

	// for multi-threads mode, we always return false
	if (useMultThreads) return false;

	// now let's see each key word matching
	if (name == "xcvar") {
		if (xcvar_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "molegrids") {
		if (molegrids_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "sigatombasis") {
		if (sigatombasis_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "batchbasis") {
		if (batchbasis_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "batchgrid") {
		if (batchgrid_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "batchvar") {
		if (batchvar_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "batchfunc") {
		if (batchfunc_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "batchxcmtrx") {
		if (batchxcmtrx_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else if (name == "espints") {
		if (espints_debug > 0) {
			return true;
		}else{
			return false;
		}
	} else {
		cout << "the class name you given is " << name << endl;
		cout << "this name can not be processed in doDebug of XCIntJobInfor" << endl;
		cout << "so we just return false" << endl;
	}
	return false;
}

void XCIntJobInfor::print() const
{
	cout << "*******************************************" << endl;
	cout << "*   DEBUG PRINT FOR XCIntJobInfor CLASS    " << endl;
	cout << "*******************************************" << endl;
	cout << "XCInts job name:                        " << jobName << endl;
	cout << "XCInts job order:                       " << order << endl;
	cout << "using multi-threads?                    " << useMultThreads << endl;
	cout << "number of CPU threads                   " << nCPUThreads << endl;
	cout << "initial batch size                      " << initBatchSize << endl;
	cout << endl << endl;
	debugPrint();
}
