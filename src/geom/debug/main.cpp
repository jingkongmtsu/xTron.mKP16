//
// testing log:
//
// Feb. 2014: 
// test the element.cpp and molecule.cpp, atom.h
//
// April 2014:
// after adding in stable sort and atom index to atom class;
// fully test molecule class and atom class
//
//
#include<cstdio>
#include<string>
#include<iostream>
#include<fstream>
#include "molecule.h"
#include "zmatrix.h"
#include "element.h"
using namespace molecule;
using namespace zmatrix;
using namespace element;

Int main(int argc, char* argv[])
{
	//
	// set input and output
	//
	string input = "input.in";
	string outf  = "test.out";
	freopen (outf.c_str(),"w",stdout);

	//
	// set test parameter
	//
	bool testMol  = false;
	bool testZ    = false;
	bool testEle  = false;
	for(Int i=1; i<argc; i++) {
		string com = argv[i];
		if (com == "molecule"  ) testMol  = true;
		if (com == "zmatrix"  )  testZ    = true;
		if (com == "element"  )  testEle  = true;
	}

	// test molecule
	if (testMol) {

		// test normal molecule reading
		Molecule mol1(input);
		mol1.print();
		printf("\n");

		Molecule mol2(input,2);
		mol2.print();
		printf("\n");

		// test one atom molecule
		Molecule mol3(1);
		mol3.print();
		printf("\n");

		// let's test that whether we can read in
		// zmatrix and transform the data into xyz
		Molecule mol4(input,3);
		mol4.print();
		printf("\n");
	}

	// test the zmatrix
	if (testZ) {

		// set the coordinates
		DoubleVec coord(3*100);
		ZMatrix z1(input,3);
		z1.print();
		z1.loadInToXYZ(coord);
		for(UInt i=0; i<z1.getNAtoms(); i++) {
			printf("for atom %d  %-12.6f %-12.6f  %-12.6f\n",(Int)i, coord[3*i+0],
					coord[3*i+1],coord[3*i+2]);
		}
		printf("\n");
		ZMatrix z2(input,4);
		z2.print();
		z2.loadInToXYZ(coord);
		for(UInt i=0; i<z2.getNAtoms(); i++) {
			printf("for atom %d %-12.6f  %-12.6f  %-12.6f\n",(Int)i, coord[3*i+0],
					coord[3*i+1],coord[3*i+2]);
		}
		printf("\n");
		ZMatrix z3(input,5);
		z3.print();
		z3.loadInToXYZ(coord);
		for(UInt i=0; i<z3.getNAtoms(); i++) {
			printf("for atom %d %-12.6f  %-12.6f  %-12.6f\n",(Int)i, coord[3*i+0],
					coord[3*i+1],coord[3*i+2]);
		}
	}

	// test element 
	if (testEle) {
		// output all of element radii
		for(UInt i=0; i<MAX_ATOMS_TYPE; i++) {
			UInt Z = i + 1;
			Double Radii = getAtomRadii(Z);
			printf("for atom %d radii is %-16.10f\n", (Int)Z, Radii);
		}
	}

	return 0;
}
