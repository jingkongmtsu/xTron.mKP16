/**
 * \file    xcvarinfor.h
 * \brief   control center for the variable information
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef XCVARINFOR_H
#define XCVARINFOR_H
#include "libgen.h"
#include <string>
using namespace std;

namespace xcvarinfor {

	//
	//  In the following section, we define the variable informaiton.
	//  These information are constant for the DFT modules.
	//  User could modify it when they need more type of variables.
	//  Generally, the information is divided into two groups, one
	//  is related to detailed variable; the other is for variable
	//  type. Finally, we provide array and functions to link them
	//  together.
	//

	/**
	 * the current number of variable types supported in the program
	 */
	const UInt MAX_DFTVAR_TYPES = 5;

	// firstly, let's define the basic variables 
	// 0 here is used as wrong return
	const UInt RA  = 1;     ///< alpha rho
	const UInt RB  = 2;     ///< beta rho
	const UInt DAX = 3;     ///< alpha gradient rho (partial rho partial x)
	const UInt DAY = 4;     ///< alpha gradient rho (partial rho partial y)
	const UInt DAZ = 5;     ///< alpha gradient rho (partial rho partial z)
	const UInt DBX = 6;     ///< beta  gradient rho (partial rho partial x)
	const UInt DBY = 7;     ///< beta  gradient rho (partial rho partial y)
	const UInt DBZ = 8;     ///< beta  gradient rho (partial rho partial z)
	const UInt GAA = 9;     ///< gamma alpha alpha: DAX dot DAX etc.
	const UInt GAB = 10;    ///< gamma alpha beta : DAX dot DBX etc.
	const UInt GBB = 11;    ///< gamma beta  beta : DBX dot DBX etc.
	const UInt TA  = 12;    ///< alpha tau
	const UInt TB  = 13;    ///< beta  tau
	const UInt LA  = 14;    ///< alpha Laplacian
	const UInt LB  = 15;    ///< beta Laplacian
	const UInt EXA = 16;    ///< alpha exchange energy density
	const UInt EXB = 17;    ///< beta exchange energy density

	// now let's define the variable type, basically it's same
	// with the first variable in this type
	const UInt RHO   = 1;
	const UInt GRHO  = 3;
	const UInt GAMMA = 9;
	const UInt TAU   = 12;
	const UInt LAP   = 14;
	const UInt EXRHO = 16;

	// now define the additional variable related information
	const UInt ALPHA_VAR = 0;
	const UInt BETA_VAR  = 1;
	const UInt GRHO_X    = 0;
	const UInt GRHO_Y    = 1;
	const UInt GRHO_Z    = 2;

	// define const array to group the variables together
	const UInt RHO_ARRAY[ ]   = { RA, RB };
	const UInt GRHO_ARRAY[ ]  = { DAX, DAY, DAZ, DBX, DBY, DBZ };
	const UInt GAMMA_ARRAY[ ] = { GAA, GAB, GBB };
	const UInt TAU_ARRAY[ ]   = { TA, TB };
	const UInt LAP_ARRAY[ ]   = { LA, LB };
	const UInt EXRHO_ARRAY[ ] = { EXA, EXB };

	/**
	 * return the variable array according to the var type
	 */
	const UInt* getDFTVarArray (const UInt& var);

	/**
	 * return the variable array length for the given variable type
	 */
	UInt getDFTVarArrayLength (const UInt& var);

	/**
	 * judge whether the given variable is beta or not
	 */
	bool isBetaDFTVar(const UInt& var);

	/**
	 * if the given var is beta, then we will give its alpha counterpart
	 * for alpha we will given it's beta counterpart
	 */
	UInt getCounterpartVar(const UInt& var);

	/**
	 * this function is used to print out the variable name in string format
	 * this is used primarily for debug purpose
	 */
	string getVarName(const UInt& v);

	/**
	 * judge whether the given variable is Gamma variable
	 */
	bool isGammaVar(const UInt& var);

	/**
	 * judge that whether the given var is Drho var
	 */
	bool isDRhoVar(const UInt& var);

	/**
	 * judge whether the given gamma variable is matching the 
	 * given gradient variable. That means, for the gradient
	 * rho of dvar, whether it satisfy the relation below?
	 * F_dvar' = F_gvar'*some factor*some gradient rho 
	 * F is some functional, F_dvar' means the first order 
	 * functional derivatives with respect to the var
	 */
	bool doesDRhoVarMatchGammaVar(const UInt& gvar, const UInt& dvar);

	/**
	 * get the factor for the process that transform gamma into drho
	 * GAA = DAX * DAX + DAY * DAY + DAZ * DAZ ->
	 * F_DAX' = F_GAA'*2*DAX factor is 2
	 * so for GBB  
	 * GBB = DBX * DBX + DBY * DBY + DBZ * DBZ ->
	 * F_DBX' = F_GBB'*2*DBX factor is 2
	 * on the other hand, for F_DAX' etc. in terms of GAB:
	 * GAB = DAX * DBX + DAY * DBY + DAZ * DBZ ->
	 * F_DBX' = F_GAB'*DAX factor is 1
	 * F_DAX' = F_GAB'*DBX factor is 1
	 * F is some functional, F_DAX' means the first order 
	 * functional derivatives with respect to the DAX
	 */
	Double getFacGammaIntoDRho(const UInt& gammaVar); 

	/**
	 * get the corresponding var for the process that transform gamma into drho
	 * GAA = DAX * DAX + DAY * DAY + DAZ * DAZ ->
	 * F_DAX' = F_GAA'*2*DAX the corresponding var is DAX
	 * GAB = DAX * DBX + DAY * DBY + DAZ * DBZ ->
	 * F_DAX' = F_GAB'*DBX the corresponding var is DBX 
	 * F_DBX' = F_GAB'*DAX the corresponding var is DAX 
	 */
	UInt getVarGammaIntoDRho(const UInt& gammaVar, const UInt& dVar);

}

#endif

