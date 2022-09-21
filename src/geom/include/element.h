/**
 * \file  element.h
 * \brief This file contains global constant data for the elements in periodic table
 * \author Fenglai Liu and Jing Kong
 */

#ifndef ELEMENT_H
#define ELEMENT_H
#include "libgen.h"
#include <string>

namespace element {

	using namespace std;

	///
	/// maximum number of atom type 
	///
	const UInt MAX_ATOMS_TYPE = 113; 

	/**
	 * \brief define the atomic symbols used in this program. The X is the ghost atom
	 */
	const string SYMBOLS [] =
	{  "X" ,
		"H" ,                                                       "He",
		"Li","Be",                         "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
		"Na","Mg",                         "Al","Si","P" ,"S" ,"Cl","Ar",
		"K" ,"Ca","Sc","Ti","V" ,"Cr","Mn",
		"Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
		"Rb","Sr","Y" ,"Zr","Nb","Mo","Tc",
		"Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
		"Cs","Ba","La",
		"Ce","Pr","Nd","Pm","Sm","Eu","Gd",
		"Tb","Dy","Ho","Er","Tm","Yb","Lu",
		"Hf","Ta","W" ,"Re",
		"Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
		"Fr","Ra","Ac",
		"Th","Pa","U" ,"Np","Pu","Am","Cm",
		"Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh",
		"Hs","Mt","Ds","Rg","Cn"
	};

	/**
	 * This is Bragg-Slater Radii data 
	 * unit is in angstrom
	 * We note that for the last few atoms, we do not have its data actually
	 * just use the 1.35
	 * the ghost atom of X is set to be same with hydrogen
	 */
	//                              X      H     He 
	const Double ATOM_RADII[ ] = {0.53E0,0.53E0,0.31E0,
		// Li   Be      B      C      N      O      F     Ne   
		1.45E0,1.05E0,0.85E0,0.70E0,0.65E0,0.60E0,0.50E0,0.38E0,
		// Na   Mg      Al     Si     P      S      Cl    Ar    
		1.80E0,1.50E0,1.25E0,1.10E0,1.00E0,1.00E0,1.00E0,0.71E0,
		// K    Ca      Sc     Ti     V      Cr     Mn     Fe     Co     Ni      Cu     Zn
		2.20E0,1.80E0,1.60E0,1.40E0,1.35E0,1.40E0,1.40E0,1.40E0,1.35E0,1.35E0,1.35E0,1.35E0,
		// Ga   Ge      As     Se     Br     Kr   
		1.30E0,1.25E0,1.15E0,1.15E0,1.15E0,0.88E0,
		// Rb   Sr      Y      Zr     Nb     Mo     Tc      Ru     Rh    Pd     Ag     Cd   
		2.35E0,2.00E0,1.80E0,1.55E0,1.45E0,1.45E0,1.35E0,1.30E0,1.35E0,1.40E0,1.60E0,1.55E0,
		// In   Sn      Sb     Te     I      Xe  
		1.55E0,1.45E0,1.45E0,1.40E0,1.40E0,1.08E0,
		// Cs   Ba      La     Ce     Pr     Nd     Pm      Sm     Eu    Gd      Tb    Dy   
		2.60E0,2.15E0,1.95E0,1.85E0,1.85E0,1.85E0,1.85E0,1.85E0,1.85E0,1.80E0,1.75E0,1.75E0,
		// Ho   Er      Tm     Yb     Lu     Hf     Ta      W      Re    Os      Ir    Pt  
		1.75E0,1.75E0,1.75E0,1.75E0,1.75E0,1.55E0,1.45E0,1.35E0,1.35E0,1.30E0,1.35E0,1.35E0,
		// Au   Hg      Tl     Pb     Bi     Po     At      Rn 
		1.35E0,1.50E0,1.90E0,1.80E0,1.60E0,1.90E0,1.27E0,1.20E0,
		// Fr   Ra      Ac     Th     Pa     U      Np      Pu     Am     Cm     Bk    Cf
		2.60E0,2.15E0,1.95E0,1.80E0,1.75E0,1.75E0,1.75E0,1.75E0,1.75E0,1.75E0,1.75E0,1.75E0,
		// Es   Fm	    Md     No     Lr     Rf     Db      Sg     Bh     Hs     Mt    Ds
		1.75E0,1.75E0,1.75E0,1.75E0,1.75E0,1.55E0,1.45E0,1.35E0,1.30E0,1.30E0,1.35E0,1.35E0,
	   //	Rg   Cn  
		1.35E0,1.35E0
	};

	/**
	 * This is spin multiplicity data for ground state atom defined in Periodic Table
	 * Obtained from the given term symbol
	 * The data below is obtained from http://www.webelements.com/
	 * We note that from the atom Lr, all of data are tentative
	 * the ghost atom of X is set to be 1 
	 */
	//                                             X   H     He 
	const UInt GROUND_STATE_ATOM_MULTIPLICITY[ ] = {1,  2,    1,
		// Li   Be      B      C      N      O      F     Ne   
		   2,   1,      2,     3,     4,     3,     2,    1,
		// Na   Mg      Al     Si     P      S      Cl    Ar    
		   2,   1,      2,     3,     4,     3,     2,    1,
		// K    Ca      Sc     Ti     V      Cr     Mn    Fe     Co     Ni   Cu     Zn
		   2,   1,      2,     3,     4,     7,     6,    5,     4,     3,   2,     1,
		// Ga   Ge      As     Se     Br     Kr   
		   2,   3,      4,     3,     2,     1,
		// Rb   Sr      Y      Zr     Nb     Mo     Tc    Ru     Rh    Pd    Ag    Cd   
		   2,   1,      2,     3,     6,     7,     6,    5,     4,     1,   2,     1,
		// In   Sn      Sb     Te     I      Xe  
		   2,   3,      4,     3,     2,     1,
		// Cs   Ba      La     Ce     Pr     Nd     Pm    Sm     Eu    Gd    Tb    Dy   
		   2,   1,      2,     1,     4,     5,     6,    7,     8,    9,    6,    5,  
		// Ho   Er      Tm     Yb     Lu     Hf     Ta    W      Re    Os    Ir    Pt  
		   4,   3,      2,     1,     2,     3,     4,    5,     6,    5,    4,    3,
		// Au   Hg      Tl     Pb     Bi     Po     At    Rn 
		   2,   1,      2,     3,     4,     3,     2,    1,
		// Fr   Ra      Ac     Th     Pa     U      Np    Pu     Am    Cm    Bk    Cf
		   2,   1,      2,     3,     4,     5,     6,    7,     8,    9,    6,    5,  
		// Es   Fm	    Md     No     Lr     Rf     Db    Sg     Bh    Hs    Mt    Ds
		   5,   3,      2,     1,     2,     3,     4,    5,     6,    5,    4,    3,
	   //	Rg   Cn 
		   2,   1
	};

	/**
	 * This is Free Atomic Polarizabilities 
	 * get it from CRC Handbook of Chemistry and Physics(94 Ed.), see the website below:
	 * http://www.hbcpnetbase.com/
	 * We do not have data after the element of Fr
	 * and the ghost atom is defined to have zero atomic polarizabilities
	 * untis are angstrom^3(1.0E-30 m^3), needed to be transformed to au (0.529...)^3.
	 */
	//                                                        X         H        He 
	const Double GROUND_STATE_ATOM_POLARIZABILITIES[ ] = { 0.0E0,   0.6668E0, 0.2051E0,  
		// Li        Be          B        C          N         O         F        Ne   
		24.3300E0, 5.6000E0, 3.0300E0, 1.6700E0, 1.1000E0, 0.8020E0, 0.5570E0, 0.3943E0,  
		// Na        Mg          Al       Si         P         S         Cl       Ar    
		24.1100E0, 10.6000E0,6.8000E0, 5.5300E0, 3.6300E0, 2.9000E0, 2.1800E0, 1.6411E0,
		// K         Ca          Sc       Ti         V         Cr        Mn       Fe     
      43.0600E0, 22.8000E0,17.8000E0,14.6000E0,12.4000E0,11.6000E0,9.4000E0, 8.4000E0,   
		// Co        Ni          Cu       Zn         Ga        Ge        As       Se     
		7.5000E0,  6.8000E0, 6.2000E0, 5.7500E0, 8.1200E0, 5.8400E0, 4.3100E0, 3.7700E0,   
		// Br        Kr   
		3.0500E0,  2.4844E0,
		// Rb        Sr          Y        Zr         Nb        Mo        Tc       Ru     
      47.2400E0, 27.6000E0,22.7000E0,17.9000E0,15.7000E0,12.8000E0,11.4000E0,9.6000E0,   
		// Rh        Pd          Ag       Cd         In        Sn        Sb       Te     
		8.6000E0,  4.8000E0, 7.2000E0, 7.3600E0, 10.2000E0,7.8400E0, 6.6000E0, 5.5000E0,  
		// I         Xe  
	  	5.3500E0,  4.0440E0,
		// Cs        Ba          La       Ce         Pr        Nd        Pm       Sm     
      59.4200E0,39.7000E0,31.1000E0, 29.6000E0,28.2000E0,31.4000E0,30.1000E0,28.8000E0,  
		// Eu        Gd          Tb       Dy         Ho        Er        Tm       Yb     
		27.7000E0,23.5000E0,25.5000E0, 24.5000E0,23.6000E0,22.7000E0,21.8000E0,21.0000E0,  
		// Lu        Hf          Ta       W          Re        Os        Ir       Pt  
		21.9000E0,16.2000E0,13.1000E0, 11.1000E0,9.7000E0, 8.5000E0, 7.6000E0, 6.5000E0,
		// Au        Hg          Tl       Pb         Bi        Po        At       Rn 
      5.8000E0, 5.0200E0, 7.6000E0,  7.0100E0, 7.4000E0, 6.8000E0, 6.0000E0, 5.3000E0,  
		
		// Fr        Ra          Ac       Th         Pa        U         Np       Pu     
		48.6000E0,36.4800E0,32.1000E0, 32.1000E0,25.4000E0,24.9000E0,24.8000E0,24.5000E0,  
		// Am        Cm          Bk       Cf         Es        Fm	     Md       No     
		23.3000E0,23.0000E0,22.7000E0, 20.5000E0,19.7000E0,23.8000E0,18.2000E0,16.4000E0,
		// Lr        Rf          Db       Sg         Bh        Hs        Mt       Ds
		-1.0E0,    -1.0E0,  -1.0E0,    -1.0E0,   -1.0E0,   -1.0E0,   -1.0E0,   -1.0E0,
	   //	Rg        Cn 
		-1.0E0,    -1.0E0
	};

	/**
	 *GetAtomicSymbolthis function is used to return the symbol of a given atomic number
	 */
	string getAtomicSymbol(UInt n); 

	/**
	 *getAtomicNumber is used to return the atomic number for a given atomic symbol
	 */
	UInt getAtomicNumber(string s);

	/**
	 * This function is used to determine whether the given string
	 * is a atomic symbol. At the same time, if it's the right symbol
	 * the s will be standardized.
	 */
	bool isAtomicSymbol(string& s);

	/**
	 * get the atomic radii in unit of au by given atomic number
	 */
	Double getAtomRadii(const UInt& Z);

	/**
	 * get the ground state atom multiplicity by given atomic number
	 */
	UInt getGroundStateAtomMult(const UInt& Z);

	/**
	 * get the atomic polarizability in unit of au^3 by given atomic number
	 */
	Double getAtomPolarizability(const UInt& Z);

	/**
	 * whether the given atom is ghost atom?
	 */
	bool isGhostAtom(const UInt& Z);

}

#endif
