/**
 * cpp file corresponding to shellprop.h
 */
#include "boost/lexical_cast.hpp"
#include "excep.h"
#include "functions.h"
#include "angmomlist.h"
#include "shellprop.h"
using namespace excep;
using namespace functions;
using namespace shellprop;

void shellprop::getlmn(const UInt& lmin, const UInt& lmax, const UInt& index, 
			UInt& l, UInt& m, UInt& n) 
{
	if (lmax<=MAX_ANG_GLOBAL_BASIS_SET_INDICES) {
		UInt k = lmin*(lmin+1)*(lmin+2)/6 + index; // always in cartesian form
		l = GLOBAL_BASIS_SET_INDICES[3*k+0];
		m = GLOBAL_BASIS_SET_INDICES[3*k+1];
		n = GLOBAL_BASIS_SET_INDICES[3*k+2];
	}else{

		// additional check
		if (lmin != lmax) {
			string infor = "lmin is: " + boost::lexical_cast<string>(lmin); 
			infor = infor + " lmax is: " + boost::lexical_cast<string>(lmax); 
			Excep excep("shellprop","getlmn",EXCEPTION_UNKNOWN_COMP_SHELL_TYPE,infor);
			handleExcep(excep);
		}

		// is the index outside range?
		if (index>=getCartBas(lmin,lmax)) {
			string infor = "index is: " + boost::lexical_cast<string>(index); 
			Excep excep("shellprop","getlmn",EXCEPTION_BASIS_SET_INDEX_OUTOF_RANGE,infor);
			handleExcep(excep);
		}

		// now we have to do it manually
		// this is how we calculate libint order
		UInt count = 0;
		UInt ang = lmin;
		for(UInt i=0; i<=ang; i++) {
			l = ang - i;
			for(UInt j=0; j<=i; j++) {
				m = i-j; 
				n = j;   
				if (count == index) {
					return;
				}
				count++;
			}
		}
	}
}

void shellprop::getL(const string& s, UInt& lmin, UInt& lmax) 
{
	//
	// now let's come to a list to get the angular momentum
	// we try to list all of known possible angular momentums
	// 
	if (s == "S" || s == "s") {
		lmin = 0;
		lmax = 0;
		return;
	}else if (s == "SP" || s == "sp") {
		lmin = 0;
		lmax = 1;
		return;
	}else if (s == "P" || s == "p") {
		lmin = 1;
		lmax = 1;
		return;
	}else if (s == "D" || s == "d") {
		lmin = 2;
		lmax = 2;
		return;
	}else if (s == "F" || s == "f") {
		lmin = 3;
		lmax = 3;
		return;
	}else if (s == "G" || s == "g") {
		lmin = 4;
		lmax = 4;
		return;
	}else if (s == "H" || s == "h") {
		lmin = 5;
		lmax = 5;
		return;
	}else if (s == "I" || s == "i") {
		lmin = 6;
		lmax = 6;
		return;
	}else if (s == "K" || s == "k") {
		lmin = 7;
		lmax = 7;
		return;
	}else if (s == "L" || s == "l") {
		lmin = 8;
		lmax = 8;
		return;
	}else if (s == "M" || s == "m") {
		lmin = 9;
		lmax = 9;
		return;
	}else if (s == "N" || s == "n") {
		lmin = 10;
		lmax = 10;
		return;
	}else if (s == "O" || s == "o") {
		lmin = 11;
		lmax = 11;
		return;
	}else if (s == "Q" || s == "q") {
		lmin = 12;
		lmax = 12;
		return;
	}else if (s == "R" || s == "r") {
		lmin = 13;
		lmax = 13;
		return;
	}else if (s == "T" || s == "t") {
		lmin = 14;
		lmax = 14;
		return;
	}else if (s == "U" || s == "u") {
		lmin = 15;
		lmax = 15;
		return;
	}else if (s == "V" || s == "v") {
		lmin = 16;
		lmax = 16;
		return;
	}else if (s == "W" || s == "w") {
		lmin = 17;
		lmax = 17;
		return;
	}else if (s == "X" || s == "x") {
		lmin = 18;
		lmax = 18;
		return;
	}else if (s == "Y" || s == "y") {
		lmin = 19;
		lmax = 19;
		return;
	}else if (s == "Z" || s == "z") {
		lmin = 20;
		lmax = 20;
		return;
	}else{
		string infor = "shell symbol is: " + s; 
		Excep excep("shellprop","getL",EXCEPTION_UNKNOWN_SHELL_SYMBOL,infor);
		handleExcep(excep);
	}
}

string shellprop::getShellName(const UInt& lmin, const UInt& lmax) 
{
	// is it SP shell?
	if (lmin == 0 && lmax == 1) return "SP";

	// additional check
	if (lmin != lmax) {
		string infor = "lmin is: " + boost::lexical_cast<string>(lmin); 
		infor = infor + " lmax is: " + boost::lexical_cast<string>(lmax); 
		Excep excep("shellprop","getShellName",EXCEPTION_UNKNOWN_COMP_SHELL_TYPE,infor);
		handleExcep(excep);
	}

	// now let's do switch
	UInt l = lmin;
	switch(l) {
		case 0:
			return "S";
			break;
		case 1:
			return "P";
			break;
		case 2:
			return "D";
			break;
		case 3:
			return "F";
			break;
		case 4:
			return "G";
			break;
		case 5:
			return "H";
			break;
		case 6:
			return "I";
			break;
		case 7:
			return "K";
			break;
		case 8:
			return "L";
			break;
		case 9:
			return "M";
			break;
		case 10:
			return "N";
			break;
		case 11:
			return "O";
			break;
		case 12:
			return "Q";
			break;
		case 13:
			return "R";
			break;
		case 14:
			return "T";
			break;
		case 15:
			return "U";
			break;
		case 16:
			return "V";
			break;
		case 17:
			return "W";
			break;
		case 18:
			return "X";
			break;
		case 19:
			return "Y";
			break;
		case 20:
			return "Z";
			break;
		default:
			string infor = "shell L is: " + boost::lexical_cast<string>(l); 
			Excep excep("shellprop","getShellName",EXCEPTION_UNKNOWN_SHELL_SYMBOL,infor);
			handleExcep(excep);
			break;
	}
	return "NONE";
}

UInt shellprop::codeShell(const UInt& section, const UInt& part) 
{
	// check the information
	if (section>=MAX_MOLECULE_SECTION_NUMBER) {
		string infor = "given section is: " + boost::lexical_cast<string>(section); 
		Excep excep("shellprop","codeShell",EXCEPTION_ILLEGAL_GEOM_SECTION_INDEX,infor);
		handleExcep(excep);
	}
	if (part>=MAX_MOLECULE_SECTION_NUMBER) {
		string infor = "aux shell index is : " + boost::lexical_cast<string>(part); 
		Excep excep("shellprop","codeShell",EXCEPTION_ILLEGAL_AUX_SHELL_INDEX,infor);
		handleExcep(excep);
	}
	return section*MAX_MOLECULE_SECTION_NUMBER + part;
}

void shellprop::decodeShell(const UInt& scode, UInt& section, UInt& part) 
{
	section = scode%MAX_MOLECULE_SECTION_NUMBER;
	part = scode - section*MAX_MOLECULE_SECTION_NUMBER;
}

UInt shellprop::getBasIndexForTheBasisSetOrder(const UInt& L, const UInt& index, const UInt& order)
{
	// is the index outside range?
	UInt nBas = getCartBas(L,L);
	if (index>=nBas) {
		string infor = "index is: " + boost::lexical_cast<string>(index); 
		Excep excep("shellprop","getBasIndexForTheBasisSetOrder",EXCEPTION_BASIS_SET_INDEX_OUTOF_RANGE,infor);
		handleExcep(excep);
	}

	// firstly get the l m n value for the input L and index
	UInt l = 0;
	UInt m = 0;
	UInt n = 0;
	getlmn(L,L,index,l,m,n); 

	// now let's see which index it is for the given order
	UInt newIndex = -1;
	if (order == AIMPAC_BASIS_SET_ORDER) {

		// double check whether the L is higher than the limit?
		if (L>MAX_ANG_AIMPAC_BASIS_SET_INDICES) {
			string infor = "for AIMPAC we only have basis set order information up to angular momentum: " + 
				boost::lexical_cast<string>(MAX_ANG_AIMPAC_BASIS_SET_INDICES); 
			Excep excep("shellprop","getBasIndexForTheBasisSetOrder",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
			handleExcep(excep);
		}

		// now let's check what is the offset for the given l m and n?
		UInt offset = 0;
		if (L>0) offset = getCartBas(0,L-1);
		for(UInt i=0; i<nBas; i++) {
			UInt k  = offset + i;
			UInt l1 = AIMPAC_BASIS_SET_INDICES[3*k+0];
			UInt m1 = AIMPAC_BASIS_SET_INDICES[3*k+1];
			UInt n1 = AIMPAC_BASIS_SET_INDICES[3*k+2];
			if (l1 == l && m1 == m && n1 == n) {
				newIndex = k;
				break;
			}
		}
	}else{

		// we do not have the given basis set order
		string infor = "the input basis set order is not recognized"; 
		Excep excep("shellprop","getBasIndexForTheBasisSetOrder",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	// check whether we have got the index
	if (newIndex == static_cast<UInt>(-1)) {
		string infor = "we failed to get a meaningful result basis set index for the input basis set order"; 
		Excep excep("shellprop","getBasIndexForTheBasisSetOrder",EXCEPTION_RESULT_IS_MEANINGLESS,infor);
		handleExcep(excep);
	}

	// finally return the new index
	return newIndex;
}

Double shellprop::normCartBasisSets(const UInt& l, const UInt& m, const UInt& n,
		const UInt& nPrim, const Double* c, const Double* e)
{
	Double sum = ZERO;
	for(UInt i=0; i<nPrim; i++) {
		for(UInt j=i; j<nPrim; j++) {
			Double exp2 = e[i]+e[j];
			Double coe  = c[i]*c[j];
			Double t    = ONE;
			if (exp2<THRESHOLD_MATH) {
				// in this case, exp2 is actually 0
				// so we do not do normalization factor calculation
				return ONE;
			}
			t = coe*gammaFunInt(2*l,exp2)*gammaFunInt(2*m,exp2)*gammaFunInt(2*n,exp2);
			if (j==i) {
				sum += t;
			} else {
				sum += TWO*t;
			}
		}
	}
	Double N = ONE/sqrt(sum);
	return N;
}

Double shellprop::normPureBasisSets(const UInt& iShell, const UInt& nPrim, const Double* c, 
	 const Double* e) {
	Double sum = ZERO;
	for(UInt i=0; i<nPrim; i++) {
		for(UInt j=i; j<nPrim; j++) {

			// coe and exp sum
			Double exp2 = e[i]+e[j];
			Double coe  = c[i]*c[j];
			if (exp2<THRESHOLD_MATH) {
				// in this case, exp2 is actually 0
				// so we do not do normalization factor calculation
				return ONE;
			}

			// t = (pi^(0.5)/4)*d1*dj/exp2^(1.5)
			Double t = (sqrt(PI)/FOUR)*coe/pow(exp2,THREE/TWO);

			// s = (2L+1)!!/2^L*exp2^L
			Double s = ONE;
			for(UInt ang=iShell; ang>0; ang--) {
				s = s*(TWO*ang+ONE)/(TWO*exp2);
			}

			// now doing sum
			if (j==i) {
				sum += s*t;
			} else {
				sum += TWO*s*t;
			}
		}
	}
	Double N = ONE/sqrt(sum);
	return N;
}

/*
void shellprop::scaleNormFactors(const Int& Lmin,const Int& Lmax,Double* N) 
{
	Int i = 0;
	for(Int iShell=Lmin;iShell<=Lmax;iShell++) {
		Int nbas= getCartBas(iShell,iShell);
		for(Int iBas=0;iBas<nbas;iBas++) {
			Int l = 0;
			Int m = 0;
			Int n = 0;
			getlmn(iShell,iShell,iBas,l,m,n);
			//Double n1 = sqrt(doubleFac(TWO*l-ONE)*doubleFac(TWO*m-ONE)*doubleFac(TWO*n-ONE));
			//Double n2 = sqrt(doubleFac(TWO*iShell-ONE));

			//
			// because we have to avoid overflow of double for large
			// L in double factor calculation, now we change into
			// new way
			//
			Double x = ONE;
			Int L = iShell;
			for(Int count=0; count<L; count++) {
				Int l0 = l-count;
				if (l0<1) l0 = 1;
				Int m0 = m-count;
				if (m0<1) m0 = 1;
				Int n0 = n-count;
				if (n0<1) n0 = 1;
				Int L0 = L-count;
				if (L0<1) L0 = 1;
				x = x*(TWO*l0-ONE)*(TWO*m0-ONE)*(TWO*n0-ONE)/(TWO*L0-ONE);
			}
			x = ONE/x;
			N[i] = sqrt(x);
			i++;
		}
	}
}

*/


