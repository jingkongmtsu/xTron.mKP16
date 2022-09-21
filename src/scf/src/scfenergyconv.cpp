/**
 * cpp file related to scfenergyconv.h
 * \author fenglai liu 
 */
#include "excep.h"
#include "scfconv.h"
#include "molecule.h"
#include "textread.h"
#include "globalinfor.h"
#include "scfenergyconv.h"
#include "parameterparsing.h"
using namespace excep;
using namespace scfconv;
using namespace molecule;
using namespace textread;
using namespace globalinfor;
using namespace scfenergyconv;
using namespace parameterparsing;

SCFEnergyConvController::SCFEnergyConvController(const GlobalInfor& infor, 
		const Molecule& mol):nIterConvTest(5),convThresh(1.0E-10)
{
	// deal with input key words
	UInt section = mol.getSec();
	string input = infor.getInputFile();
	ParameterParsing pp(input,"scfconv",section);
	if (! pp.hasAnyParameters()) return;

	// now let's read in the key of xcints
	// here we can not check the invalid keywords, since 
	// there are other keywords also defined for scfconv section
	UInt nKey = pp.getKeyNumber();
	WordConvert w;
	for(UInt iKey=0; iKey<nKey; iKey++) {

		// get the key name
		string key = pp.getKey(iKey);

		// now begin the read of the keywords
		if (w.compare("THRESHOLD_ENERGY_CONV",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  THRESHOLD_ENERGY_CONV 
			/// group    scfconv
			/// values   double precision value
			/// 
			/// this keyword defines the threshold value to determine the 
			/// SCF convergence in terms of the energy difference. For 
			/// using the SCFEnergyConvController class, if there are
			/// n consecutive SCF total energy their energy difference
			/// abs(E(this iteration) - E(previous iteration)) is less
			/// than this threshold value, then we consider the SCF
			/// is just converged.
			///
			/// the default value for this keyword is 1.0E-10.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			Double tmp = ZERO;
			if (!w.toDouble(value,tmp)) {
				string infor = "we can not process threshold value. not double.";
				Excep excep("SCFEnergyConvController","constructor",EXCEPTION_CONVERT_TO_DOUBLE_FAILED,infor);
				handleExcep(excep);
			}
			convThresh = tmp;

		}else if (w.compare("MAX_ITERATIONS_ENERGY_CONV",key)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  MAX_ITERATIONS_ENERGY_CONV
			/// group    scfconv
			/// values   unsigned integer number
			///          
			/// this keyword defines the interations to determine the 
			/// SCF convergence in terms of the energy difference. For 
			/// using the SCFEnergyConvController class, if there are
			/// n consecutive SCF total energy(here n is defined here) 
			/// their energy difference
			/// abs(E(this iteration) - E(previous iteration)) is less
			/// than the defined threshold value, then we consider the SCF
			/// is just converged.
			///
			/// the default value for this keyword 5.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string value = pp.getValue(key);
			LineParse l(value);
			UInt tmp = 0;
			string tmpValue = l.findValue(0);
			if (!w.toUInt(tmpValue,tmp)) {
				string infor = "the number of iterations defined for MAX_ITERATIONS_ENERGY_CONV is invalid, not an integer";
				Excep excep("SCFEnergyConvController","constructor",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
				handleExcep(excep);
			}
			nIterConvTest = tmp;
		}
	}
}

bool SCFEnergyConvController::hasEnergyConv(const SCFConv& conv) const
{
	// if this is SCF cycles less than what we can
	// judged, return false
	UInt scfIndex = conv.getCurrentSCFIndex();
	if (scfIndex<=nIterConvTest) return false;

	// so let's see the energy difference
	const DoubleVec& energyList = conv.getHistEnergyData(); 
	for(UInt i=0; i<nIterConvTest; i++) {
		UInt thisIter = scfIndex-i; 
		Double diff = fabs(energyList[thisIter]-energyList[thisIter-1]);
		if (diff>convThresh) return false;
	}

	// now let's return true
	return true;
}
