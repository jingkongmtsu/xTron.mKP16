//
// April 2014:
// here we only test the parallel timing
//
#include<iostream>
#include "molecule.h"
#include "shell.h"
#include "matrix.h"
#include "globalinfor.h"
#include "cptrans.h"
using namespace molecule;
using namespace shell;
using namespace matrix;
using namespace globalinfor;
using namespace cptrans;

void cptrans_para_test()
{
	string inf = "cptrans_test/dna.in";
	Molecule molecule(inf,1);
	MolShell ms(inf,molecule);
	GlobalInfor infor(inf);
	infor.enableMultiThreadsMode();
	UInt nbas = ms.getNBas();
	UInt ncarbas = ms.getNCarBas();
	std::cout << "totally number of basis sets: " << nbas << std::endl;
	std::cout << "totally number of Cartesian basis sets: " << ncarbas << std::endl;

	// 
	// here we need to do several tests
	// for CP transformation, we will test that:
	// 1  the job is done on shell pair
	// 2  on row shell then on column shell
	//
	bool doC2POnShellPair   = true;
	bool doC2POnRowColShell = true;

	// work on shell pair
	if (doC2POnShellPair) {

		// now create the input
		Mtrx M1(ncarbas,ncarbas);
		M1.set(ONE);

		// do the C2P trans
		CPTransBasisNorm trans(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW_COL);
		trans.transform(ms,ms,M1);
	}

	// work on row shell then on column shell
	if (doC2POnRowColShell) {

		// now create the input
		Mtrx M1(ncarbas,ncarbas);
		M1.set(ONE);

		// do the C2P trans on row
		CPTransBasisNorm trans1(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW);
		trans1.transform(ms,ms,M1);

		// do the C2P trans on col
		CPTransBasisNorm trans2(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_COL);
		trans2.transform(ms,ms,M1);
	}

}
