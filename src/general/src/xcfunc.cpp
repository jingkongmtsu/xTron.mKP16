/**

 * CPP files corresponding to the xcfunc.h
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include<fstream>
#include <boost/lexical_cast.hpp>
#include "excep.h"
#include "textread.h"
#include "parameterparsing.h"
#include "globalinfor.h"
#include "xcfunc.h"
using namespace excep;
using namespace textread;
using namespace parameterparsing;
using namespace globalinfor;
using namespace xcfunc;

////////////////////////////////////////////////
//   memeber functions for XCFunc             //
////////////////////////////////////////////////
void XCFunc::setupXCFunc()
{
	// get the setting file env
	char* path = getenv("EMUL_SETTING_FILE");
	if (path==NULL) {
		Excep excep("XCFunc","setupXCFunc",EXCEPTION_DIR_MISSING,
				"EMUL_SETTING_FILE is not defined");
		handleExcep(excep);
	}

	// set up the setting file path
	// where we can find the xcfunc definition
	string settingPath(path);
	string setting_file = settingPath + "/" + "xcfunc.conf";
	string input_file   = setting_file;

	// now open the resource
	const char* source = input_file.c_str();
	ifstream inf;
	inf.open(source,ios::in);
	if (!inf) {
		Excep excep("XCFunc","setupXCFunc",EXCEPTION_FILE_MISSING,input_file);
		handleExcep(excep);
	}

	// check the number of functionals
	// if the correlation functional name is none, then we only perform once
	// if this is user-defined functional, w perform once two
	UInt nFunc = 2;
	if (ecName == "NONE") nFunc = 1;
	bool defineExCoe = false;
	bool defineEcCoe = false;
	for(UInt iFunc=0; iFunc<nFunc; iFunc++) {

		// get the functional name
		string functionalName = exName;
		if (iFunc == 1) functionalName = ecName;
		inf.seekg(0,ios::beg);
		bool find = false;
		string line;
		WordConvert w;
		while(getline(inf,line)) {
			LineParse l(line);
			if (l.isCom() || l.isEmp()) {
				continue;
			} else if (l.getNPieces() == 2) {
				string keyword = l.findValue(0);
				if (w.compare(keyword,"functional")) {
					string funcName = l.findValue(1);
					if (w.compare(functionalName,funcName)) {
						find = true;
						break;
					}
				}
			}
		}
		if (! find) {
			string infor = "in xcfunc.conf functional name search is failed: ";
			infor = infor + functionalName;
			Excep excep("XCFunc","setupXCFunc",EXCEPTION_MODULE_SEARCH_FAILED,infor);
			handleExcep(excep);
		}

		// now we begin to process the information given
		while(getline(inf,line)) {
			LineParse l(line);

			// shall we stop here?
			// for xcfunc.conf, pure comment line or empty line, this is the end of definition
			if(l.isCom() || l.isEmp()) break;

			// pre-check 
			UInt n = l.getNPieces();
			if(n<2) {
				cout << line << endl;
				Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
						"The parameter keyword in xcfunc should have its value defined");
				handleExcep(excep);
			}

			// now let's come up with a nesty if else code
			string keyword = l.findValue(0);
			w.capitalize(keyword);

			// if this is composite functional, then we should 
			// know its pieces - for exchange
			if (keyword == "EXCHANGE") {
				for(UInt i=1; i<l.getNPieces(); i++) {
					string funcname = l.findValue(i);
					w.capitalize(funcname);
					exchComponents.push_back(funcname);
					lenEx++;
				}
			}

			// factors for each exchange elementary piece
			if (keyword == "EXCHANGE_COEFFICIENTS") {
				defineExCoe = true;
				for(UInt i=1; i<l.getNPieces(); i++) {
					string c = l.findValue(i);
					Double coe = ZERO;
					if (! w.toDouble(c, coe)) {
						cout << line << endl;
						Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
								"Can not transform the coefficient into Double type.");
						handleExcep(excep);
					}else{
						exchFactors.push_back(coe);
					}
				}
			}

			// elementary correlation functional
			if (keyword == "CORRELATION") {
				for(UInt i=1; i<l.getNPieces(); i++) {
					string funcname = l.findValue(i);
					w.capitalize(funcname);
					corrComponents.push_back(funcname);
					lenEc++;
				}
			}

			if (keyword == "CORRELATION_COEFFICIENTS") {
				defineEcCoe = true;
				for(UInt i=1; i<l.getNPieces(); i++) {
					string c = l.findValue(i);
					Double coe = ZERO;
					if (! w.toDouble(c, coe)) {
						cout << line << endl;
						Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
								"Can not transform the coefficient into Double type.");
						handleExcep(excep);
					}else{
						corrFactors.push_back(coe);
					}
				}
			}

			// now let's go to the common part shared by exchange/correlation
			// variables for this functional
			if (keyword == "VAR_TYPE") {

				// normal type of functional
				if (n == 2) {
					string type = l.findValue(1);
					w.capitalize(type);
					if (type == "LDA") {
						hasRho = true;
					}else if (type == "GGA") {
						hasRho = true;
						hasGRho = true;
					}else if (type == "META-GGA") {
						hasRho = true;
						hasGRho = true;
						hasTau  = true;
					}else if (type == "META-GGA-WITH-LAP") {
						hasRho = true;
						hasGRho = true;
						hasTau  = true;
						hasLap  = true;
					}else if (type == "ORBITAL") {
						orbitalFunc  = true;
					}else{
						cout << line << endl;
						Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
								"Wrong type of var in given functional provided");
						handleExcep(excep);
					}
				}else{
					// in this case, user should list the variables explicitly
					for(UInt i=1; i<l.getNPieces(); i++) {
						string var = l.findValue(i);
						w.capitalize(var);
						if (var == "RHO") {
							hasRho = true;
						}else if (var == "GAMMA") {
							hasGRho = true;
						}else if (var == "TAU") {
							hasTau  = true;
						}else if (var == "LAP") {
							hasLap  = true;
						}else if (var == "ORBITAL") {
							orbitalFunc  = true;
						}else if (var == "EXRHO") {
							hasExchRho = true;
						}else{
							cout << line << endl;
							Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
									"In xcfunc.conf Wrong type of var in given functional provided");
							handleExcep(excep);
						}
					}
				}
			} 

			// linear or non-linear functional
			if (keyword == "LINEAR") {
				string state = l.findValue(1);
				w.capitalize(state);
				if (state == "FALSE" || state == "F") {
					islinear = false;
				}else{
					islinear = true;
				}
			}

			// linear or non-linear functional
			if (keyword == "FUNCTIONAL_DERIV") {
				string value = l.findValue(1);
				UInt tmp = 0;
				if (!w.toUInt(value,tmp)) {
					string infor = "the functional derivatives order setting is invalid, not an integer";
					Excep excep("XCFunc","setupXCFunc",EXCEPTION_CONVERT_TO_INT_FAILED,infor);
					handleExcep(excep);
				}
				maxFunDerivOrder = tmp;
			}

			// let's see whether the you have B13 correlation functional switch method defined?
			if (keyword == "B13_COOR_METHOD") {
				string c = l.findValue(1);
				UInt  method = -1;
				if (! w.toUInt(c, method)) {
					cout << line << endl;
					Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
							"Can not transform the value for B13_COOR_METHOD into UInt type.");
					handleExcep(excep);
				}
				B13CoorMethod = method;
			}

			// let's see whether the you have set choice for B05 nondynamic correlation for parallel spin.
			if (keyword == "B05_NDPAR_METHOD") {
				string c = l.findValue(1);
				UInt  method;
				if (! w.toUInt(c, method)) {
					cout << line << endl;
					Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
							"Can not transform the value for B05_NDPAR_METHOD into UInt type.");
					handleExcep(excep);
				}
				B05NDParMethod = method;
			}
		}
	}

	// close the file
	inf.close();

	// post processing
	// if the exchange functional is defined as a whole one (not composite),
	// then the name is just the exchange functional name
	if (lenEx == 0) {
		exchComponents.push_back(exName);
		exchFactors.push_back(ONE);
		lenEx++;
	}

	// so as the correlation
	if (ecName != "NONE" && lenEc == 0) {
		corrComponents.push_back(ecName);
		corrFactors.push_back(ONE);
		lenEc++;
	}

	// orbital functional may not be explictly given
	// but we can know it from the names
	if (isHybrid()) orbitalFunc = true;

	// we still need to check the coefficient arrays
	if (islinear) {
		if (lenEx > 0) {
			if (defineExCoe) {
				UInt n = exchFactors.size();
				if (n != lenEx) {
					Excep excep("XCFunc","setupXCFunc",INVALID_FUNCTIONAL_COEFFICIENTS,
							"exchange coef. array is not matching number of exchange functionals"); 
					handleExcep(excep);
				}
			}else{
				exchFactors.assign(lenEx, ONE);
			}
		}
		if (lenEc > 0) {
			if (defineEcCoe) {
				UInt n = corrFactors.size();
				if (n != lenEc) {
					Excep excep("XCFunc","setupXCFunc",INVALID_FUNCTIONAL_COEFFICIENTS,
							"exchange coef. array is not matching number of exchange functionals"); 
					handleExcep(excep);
				}
			}else{
				corrFactors.assign(lenEc, ONE);
			}
		}
	}

	// also we need to check the method switch for both exchange/correlation
	// the method must be 1 or 2
	// the default value of 0 is invalid
	if (useB13Coor()) {
		if (B13CoorMethod != 1 && B13CoorMethod!= 2) {
			cout << "For B13 correlation functional, please set B13_COOR_METHOD value in xcfunc.conf" << endl;
			cout << "B13_COOR_METHOD = 1 means you use lambda integrated form formula" << endl;
			cout << "B13_COOR_METHOD = 2 means you use lambda= 1 form formula" << endl;
			cout << "only the two values are allowed, for more detail please refer to the b13coor.f" << endl;
			Excep excep("XCFunc","setupXCFunc",EXCEPTION_INPUT_PARAMETER_INVALID,
					"the B13 correlation functional's method is not properly set");
			handleExcep(excep);
		}
	}
}

XCFunc::XCFunc(const string& name):B13CoorMethod(B13COORMETHOD_DEFAULT), B05NDParMethod(B05NDPARMETHOD_DEFAULT),
	becke05_p(BECKE05_P_DEFAULT),becke05_q(BECKE05_Q_DEFAULT),
	kp14_alpha(KP14_ALPHA_DEFAULT),kp14_b(KP14_B_DEFAULT),kp14_c_ndpar(KP14_C_NDPAR_DEFAULT), 
	kp14_c_ndpar_cap(KP14_C_NDPAR_CAP_DEFAULT),
	exName(name),ecName("NONE"),lenEx(0),lenEc(0),maxFunDerivOrder(1),
	islinear(true), hasRho(false), hasGRho(false),hasTau(false),hasLap(false),
	hasExchRho(false),orbitalFunc(false)
{
	setupXCFunc();
}

XCFunc::XCFunc(const string& ex, const string& ec):B13CoorMethod(B13COORMETHOD_DEFAULT), B05NDParMethod(B05NDPARMETHOD_DEFAULT),
	becke05_p(BECKE05_P_DEFAULT),becke05_q(BECKE05_Q_DEFAULT),
	kp14_alpha(KP14_ALPHA_DEFAULT),kp14_b(KP14_B_DEFAULT),kp14_c_ndpar(KP14_C_NDPAR_DEFAULT),
	kp14_c_ndpar_cap(KP14_C_NDPAR_CAP_DEFAULT),
	exName(ex),ecName(ec),lenEx(0),lenEc(0),maxFunDerivOrder(1),
	islinear(true), hasRho(false), hasGRho(false),hasTau(false),hasLap(false),
	hasExchRho(false),orbitalFunc(false)
{
	setupXCFunc();
}

XCFunc::XCFunc(const string& input, const UInt& section):B13CoorMethod(B13COORMETHOD_DEFAULT), B05NDParMethod(B05NDPARMETHOD_DEFAULT),
	becke05_p(BECKE05_P_DEFAULT),becke05_q(BECKE05_Q_DEFAULT),
	kp14_alpha(KP14_ALPHA_DEFAULT),kp14_b(KP14_B_DEFAULT),kp14_c_ndpar(KP14_C_NDPAR_DEFAULT),
	kp14_c_ndpar_cap(KP14_C_NDPAR_CAP_DEFAULT),
	exName("NONE"),ecName("NONE"),
	lenEx(0),lenEc(0),maxFunDerivOrder(1),islinear(true),hasRho(false),hasGRho(false),hasTau(false),
	hasLap(false),hasExchRho(false),orbitalFunc(false)
{
	// now parse the functional name from input file
	ParameterParsing pp(input,"xcfunc",section);
	if (! pp.hasAnyParameters()) {
		string infor = "xcfunc is misssing in the section: ";
		infor = infor + boost::lexical_cast<string>(section);
		Excep excep("XCFunc","constructor",EXCEPTION_MODULE_SEARCH_FAILED,infor);
		handleExcep(excep);
	}
	// now let's read in the keywords
	WordConvert w;
        //Catch the functional first so that it can be modified by other keywords.
        if (pp.hasKeyDefined("NAME"))
        {
			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  name
			/// group    xcfunc
			/// values   functional name
			/// 
			/// this keyword defines the exchange and/or correlation functional
			/// name. The values can be one string (for example, B3LYP represents
			/// a full name of a functional), or two strings (like Becke88 LYP);
			/// the first one represents the exchange functional and the second one
			/// is the correlation functional. No default value is set.
			///
			/// we will perform the functional name check according to the xcfunc.conf
			/// file in the setting folder to see whether the user input functional
			/// is defined in the code.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///

			string name = pp.getValue("NAME");
			w.capitalize(name);

			// parsing the name
			LineParse l(name);
			UInt nFunc = l.getNPieces();
			if (nFunc == 1) {

				// in this case, exchange and correlation are named together
				// or only one part exists
				exName = name;
				w.capitalize(exName);
				ecName ="NONE";

			}else if (nFunc == 2) {

				// in this case, exchange and correlation are named saperately
				exName = l.findValue(0);
				w.capitalize(exName);
				ecName = l.findValue(1);
				w.capitalize(ecName);

			}else{
				cout << "the correct way to use a functional in xcfunc is like this" << endl;
				cout << "\%xcfunc" << endl;
				cout << "name B3LYP (another way is \"name Slater VWN5\")" << endl;
				cout << "\%end" << endl;
				cout << "the functional name you input is " << name << endl;
				Excep excep("XCFunc","CONSTRUCTOR",INVALID_FUNCTIONAL_NAME,
						"you did not define the functional name properly");
				handleExcep(excep);
			}

			// now go to set up the functional information
			setupXCFunc();
	}
	else
	{
		Excep excep("XCFunc","CONSTRUCTOR",INVALID_FUNCTIONAL_NAME,
			"you did not define a functional.");
		handleExcep(excep);
	}
	//Now treat the rest of the keywords.
	UInt nKey = pp.getKeyNumber();
cout << "here nkey" << nKey << endl;
	for(UInt iKey=0; iKey<nKey; iKey++) {

		// get the key name
		string keyword = pp.getKey(iKey);

		// now let's read in
		if (w.compare("NAME",keyword)) { //Already processed.

		}else if (w.compare("KP14_ALPHA",keyword)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  KP14_alpha
			/// group    xcfunc
			/// values   the KP14 functional parameter of alpha
			///
			/// for the specific definition please see the functional
			/// file kp14ec.f. The default value is 1.355.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value kp14_alpha is invalid, must be a floating point number");
				handleExcep(excep);
			}
			kp14_alpha = c;

		}else if (w.compare("KP14_B",keyword)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  KP14_b
			/// group    xcfunc
			/// values   the KP14 functional parameter of b
			///
			/// for the specific definition please see the functional
			/// file kp14ec.f. The default value is 0.038. 
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value kp14_b is invalid, must be a floating point number");
				handleExcep(excep);
			}
			kp14_b = c;

		}else if (w.compare("KP14_C_NDPAR",keyword)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  KP14_C_NDPAR
			/// group    xcfunc
			/// values   the KP14 functional parameter of c_ndpar
			///
			/// for the specific definition please see the functional
			/// file kp14ec.f. The default value is 1.128.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value kp14_c_ndpar is invalid, must be a floating point number");
				handleExcep(excep);
			}
			kp14_c_ndpar = c;     

		}else if (w.compare("KP14_C_NDPAR_CAP",keyword)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  KP14_C_NDPAR_CAP
			/// group    xcfunc
			/// values   the KP14 functional parameter for capping the parallel NDC.  It is
			///          applicable to KP14_mNDPAR1, and KP14_mNDPAR2.  The default value
                        ///          is not good.  Roberto is still optimizing.
			/// 
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value kp14_c_ndpar_cap is invalid, must be a floating point number");
				handleExcep(excep);
			}
			kp14_c_ndpar_cap = c;     

		}else if (w.compare("BECKE05_P",keyword)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  BECKE05_P
			/// group    xcfunc
			/// values   the functional parameter of P value in becke05 functional
			///
			/// for the specific definition please see the functional
			/// file becke05ex.f. The default value is 115.0.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value BECKE05_P is invalid, must be a floating point number");
				handleExcep(excep);
			}
			becke05_p = c;     

		}else if (w.compare("BECKE05_Q",keyword)) {

			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  BECKE05_Q
			/// group    xcfunc
			/// values   the functional parameter of Q value in becke05 functional
			///
			/// for the specific definition please see the functional
			/// file becke05ex.f. The default value is 120.0.
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value BECKE05_Q is invalid, must be a floating point number");
				handleExcep(excep);
			}
			becke05_q = c;     

		}else if (w.compare("B13COOR_OPP",keyword)) {
			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  B13COOR_OPP
			/// group    xcfunc
			/// values   the coef for the B13COOR_OPP correlation component.
			///
			/// for the specific definition please see xcfunc.conf. The default is there. 
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value B13COOR_OPP is invalid, must be a floating point number");
				handleExcep(excep);
			}

			ptrdiff_t iCorr = find(corrComponents.begin(), corrComponents.end(), "B13COOR_OPP") 
                                         - corrComponents.begin();     
			corrFactors[iCorr] = c;

		}else if (w.compare("B13COOR_PAR",keyword)) {
			///
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			/// keyword  B13COOR_OPP
			/// group    xcfunc
			/// values   the coef for the B13COOR_OPP correlation component.
			///
			/// for the specific definition please see xcfunc.conf. The default is there. 
			/// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			///
			string v = pp.getValue(keyword);
			Double c = ZERO;
			if (! w.toDouble(v, c)) {
				cout << v << endl;
				Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
						"the keyword value B13COOR_PAR is invalid, must be a floating point number");
				handleExcep(excep);
			}

			ptrdiff_t iCorr = find(corrComponents.begin(), corrComponents.end(), "B13COOR_PAR") 
                                         - corrComponents.begin();     
			corrFactors[iCorr] = c;

		}else {
			cout << keyword << endl;
			Excep excep("XCFunc","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,
					"the keyword value in the input section can not be parsed");
			handleExcep(excep);
		}
	}
}

bool XCFunc::isHybrid() const {
	for(UInt i=0; i<getLenEx(); i++) {
		if (exchComponents[i] == "HF") {
			return true;
		}
	}
	return false;
}

Double XCFunc::HFCoeff() const {
	if (! isHybrid()) return ZERO;
	UInt pos = 0;
	for(UInt i=0; i<getLenEx(); i++) {
		if (exchComponents[i] == "HF") {
			pos = i;
			break;
		}
	}
	return exchFactors[pos];
}

bool XCFunc::hasExchangeHole() const {
	for(UInt i=0; i<getLenEx(); i++) {
		if (exchComponents[i] == "VDWBR89") {
			return true;
		}
	}
	return false;
}

bool XCFunc::useB05FormStaticCoor() const {
	for(UInt i=0; i<getLenEx(); i++) {
		if (exchComponents[i] == "B05_NDOP" || exchComponents[i] == "B05_NDPAR") {
			return true;
		}
	}
	return false;
}

bool XCFunc::useB13Coor() const {
	for(UInt i=0; i<getLenEc(); i++) {
		if (corrComponents[i] == "B13COOR_OPP" || corrComponents[i] == "B13COOR_PAR") {
			return true;
		}
	}
	return false;
}

bool XCFunc::withoutDFTFunc() const
{
	// whether this is pure HF calculation?
	// for exchange, only HF involves
	// on the other hand, no correlation exists
	if (getLenEx() == 1 && exchComponents[0] == "HF" && getLenEc() == 0) return true;

	// the rest of cases, we have dft functional
	return false;
}

void XCFunc::print() const
{

	cout << "*******************************************" << endl;
	cout << "*           XCFunc Class                  *" << endl;
	cout << "*******************************************" << endl;
	cout << "Exchange functional name: " <<  exName << endl;
	if (hasCorrelation()) {
		cout << "Correlation functional name: " <<  ecName << endl;
	}

	// var types
	cout << "Varaibles for this functional: ";
	if (hasRho) cout << "Rho ";
	if (hasGRho) cout << "GRho ";
	if (hasTau) cout << "Tau ";
	if (hasLap) cout << "Lap ";
	if (hasExchRho) cout << "Exchange-Density ";
	if (orbitalFunc) cout << "Orbital-Functional ";
	cout << endl;

	// detailed functional
	cout << "Elementary exchange functional details" <<  endl;
	for(UInt i=0; i<lenEx; i++) {
		cout << "Functional name: " << exchComponents[i] << " Coef. " << 
			exchFactors[i] << endl;
	}
	cout << endl;
	if (hasCorrelation()) {
		cout << "Elementary correlation functional details" <<  endl;
		for(UInt i=0; i<lenEc; i++) {
			cout << "Functional name: " << corrComponents[i] << " Coef. " << 
				corrFactors[i] << endl;
		}
		cout << endl;
	}
}
