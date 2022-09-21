//
// April 2014:
// finish the cp transformation test as well as basis set scale test
// on both serial and parallel codes here
//
#include<iostream>
#include "molecule.h"
#include "shell.h"
#include "matrix.h"
#include "filerw.h"
#include "cptrans.h"
#include "globalinfor.h"
//#include "blockmatrixlist.h"
using namespace molecule;
using namespace shell;
using namespace matrix;
using namespace filerw;
using namespace globalinfor;
//using namespace blockmatrixlist;
using namespace cptrans;

void cptrans_test(bool withPara)
{
	string inf = "cptrans_test/input.in";
	Molecule molecule(inf,1);
	MolShell ms(inf,molecule);
	GlobalInfor infor(inf);
	if (withPara) {
		infor.enableMultiThreadsMode();
	}else{
		infor.disableMultiThreadsMode();
	}
	//ms.print(3);
	UInt nbas = ms.getNBas();
	UInt ncarbas = ms.getNCarBas();

	// 
	// here we need to do several tests
	// firstly, for CP transformation,
	// we will test that:
	// 1  the job is done on shell pair
	// 2  on row shell then on column shell
	// for basis set normalization, we test:
	// 1  the job is done on both row and column
	// 2  the job is done on row then col
	//
	// finally, we add a case that we originally
	// failed to do it correctly, this is the 
	// molecule 4
	//
	bool doC2POnShellPair   = true;
	bool doC2POnRowColShell = true;
	bool doNormOnShells     = true;
	bool doNormOnRowCol     = true;
	bool doMatrixTranspose  = true;
	bool specialCaseTest    = true;

	/////////////////////////////////////
	//  @@@@     C2P  work             //
	/////////////////////////////////////

	// work on shell pair
	if (doC2POnShellPair) {

		// set jobs
		bool doMatrix = true;
		bool doList   = true;

		// this is for matrix type of data
		if (doMatrix) {

			// now create the input
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// do the C2P trans
			CPTransBasisNorm trans(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW_COL);
			trans.transform(ms,ms,M1);

			// now read in output result
			Mtrx M2(nbas,nbas);
			string out = "cptrans_test/pure.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),nbas,nbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " C2P on shell pair with matrix data is failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " C2P on shell pair with matrix data is passed with accuracy of " << thresh << std::endl;
			}
		}

		// now test the list type of data
		/*
		if (doList) {

			// now create the input
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// next, let's creat the matrix list
			BlockMtrxList list;
			for(UInt jAtom=0; jAtom<ms.getNAtomShells(); jAtom++) {
				const AtomShell& cas = ms.getAtomShell(jAtom);
				UInt pCol = cas.getBasisStartIndex(TYPE_CART);
				UInt nCol = cas.getNCarBas();
				for(UInt iAtom=jAtom; iAtom<ms.getNAtomShells(); iAtom++) {
					const AtomShell& ras = ms.getAtomShell(iAtom);
					UInt pRow = ras.getBasisStartIndex(TYPE_CART);
					UInt nRow = ras.getNCarBas();

					// now let's build the block matrix
					// all of the lower-triangular part
					bool inTrans = false;
					BlockMtrx B(nRow,nCol,pRow,pCol,inTrans);

					// block matrix also need the proper index
					B.initPropIndex(iAtom,jAtom);

					// form B data and push it to list
					B.getData(M1);
					list.update(B,false);
				}
			}

			// now do transformation
			CPTransBasisNorm trans(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW_COL);
			trans.transform(ms,ms,list);

			// update list into the M1 matrix
			bool reform = true;
			M1.reset(nbas,nbas,reform);
			list.updateMatrix(M1);

			// now read in output result
			Mtrx M2(nbas,nbas);
			string out = "cptrans_test/pure.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),nbas,nbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " C2P on shell pair with matrix list is failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " C2P on shell pair with matrix list is passed with accuracy of " << thresh << std::endl;
			}
		}
		*/
	}

	// work on row shell then on column shell
	if (doC2POnRowColShell) {

		// set job
		bool onMatrix = true;
		bool onList = true;

		// work on matrix
		if (onMatrix) {

			// now create the input
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// do the C2P trans on row
			CPTransBasisNorm trans1(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW);
			trans1.transform(ms,ms,M1);

			// do the C2P trans on col
			CPTransBasisNorm trans2(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_COL);
			trans2.transform(ms,ms,M1);

			// now read in output result
			Mtrx M2(nbas,nbas);
			string out = "cptrans_test/pure.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),nbas,nbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " C2P on row/column shell with matrix data is failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " C2P on row/column shell with matrix data is passed with accuracy of " << thresh << std::endl;
			}
		}

		// on list
		/*
		if (onList) {

			// now create the input
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// next, let's creat the matrix list
			BlockMtrxList list;
			for(UInt jAtom=0; jAtom<ms.getNAtomShells(); jAtom++) {
				const AtomShell& cas = ms.getAtomShell(jAtom);
				UInt pCol = cas.getBasisStartIndex(TYPE_CART);
				UInt nCol = cas.getNCarBas();
				for(UInt iAtom=jAtom; iAtom<ms.getNAtomShells(); iAtom++) {
					const AtomShell& ras = ms.getAtomShell(iAtom);
					UInt pRow = ras.getBasisStartIndex(TYPE_CART);
					UInt nRow = ras.getNCarBas();

					// now let's build the block matrix
					// all of the lower-triangular part
					bool inTrans = false;
					BlockMtrx B(nRow,nCol,pRow,pCol,inTrans);

					// block matrix also need the proper index
					B.initPropIndex(iAtom,jAtom);

					// form B data and push it to list
					B.getData(M1);
					list.update(B,false);
				}
			}

			// do the C2P trans on row
			CPTransBasisNorm trans1(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW);
			trans1.transform(ms,ms,list);

			// do the C2P trans on col
			CPTransBasisNorm trans2(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_COL);
			trans2.transform(ms,ms,list);

			// update list into the M1 matrix
			bool reform = true;
			M1.reset(nbas,nbas,reform);
			list.updateMatrix(M1);

			// now read in output result
			Mtrx M2(nbas,nbas);
			string out = "cptrans_test/pure.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),nbas,nbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " C2P on row/column shell with matrix list is failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " C2P on row/column shell with matrix list is passed with accuracy of " << thresh << std::endl;
			}
		}
		*/
	}

	/////////////////////////////////////
	//  @@@@     Norm  work            //
	/////////////////////////////////////
	Molecule mol(inf,2);
	MolShell ms2(inf,mol);

	// work on shells 
	if (doNormOnShells) {

		// assign jobs
		bool onMatrix = true;
		bool onList = true;

		// work on matrix
		if (onMatrix) {

			// now create the input
			UInt ncarbas = ms2.getNCarBas();
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// do the normalization
			CPTransBasisNorm trans(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW_COL);
			trans.transform(ms2,ms2,M1);

			// now read in output result
			Mtrx M2(ncarbas,ncarbas);
			string out = "cptrans_test/cart2.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),ncarbas,ncarbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " normalization on shells with matrix is failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " normalization on shells with matrix is passed with accuracy of " << thresh << std::endl;
			}
		}

		// now on list
		/*
		if (onList) {

			// now create the input
			UInt ncarbas = ms2.getNCarBas();
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// next, let's creat the matrix list
			BlockMtrxList list;
			for(UInt jAtom=0; jAtom<ms2.getNAtomShells(); jAtom++) {
				const AtomShell& cas = ms2.getAtomShell(jAtom);
				UInt pCol = cas.getBasisStartIndex(TYPE_CART);
				UInt nCol = cas.getNCarBas();
				for(UInt iAtom=jAtom; iAtom<ms2.getNAtomShells(); iAtom++) {
					const AtomShell& ras = ms2.getAtomShell(iAtom);
					UInt pRow = ras.getBasisStartIndex(TYPE_CART);
					UInt nRow = ras.getNCarBas();

					// now let's build the block matrix
					// all of the lower-triangular part
					bool inTrans = false;
					BlockMtrx B(nRow,nCol,pRow,pCol,inTrans);

					// block matrix also need the proper index
					B.initPropIndex(iAtom,jAtom);

					// form B data and push it to list
					B.getData(M1);
					list.update(B,false);
				}
			}

			// do the normalization
			CPTransBasisNorm trans(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW_COL);
			trans.transform(ms2,ms2,list);

			// update list into the M1 matrix
			bool reform = true;
			M1.reset(ncarbas,ncarbas,reform);
			list.updateMatrix(M1);

			// now read in output result
			Mtrx M2(ncarbas,ncarbas);
			string out = "cptrans_test/cart2.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),ncarbas,ncarbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " normalization on shells with list is failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " normalization on shells with list is passed with accuracy of " << thresh << std::endl;
			}
		}
		*/
	}

	// work on row shell then on column shell
	if (doNormOnRowCol) {

		// assign jobs
		bool onMatrix = true;
		bool onList = true;

		// work on matrix
		if (onMatrix) {

			// now create the input
			UInt ncarbas = ms2.getNCarBas();
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// on row
			CPTransBasisNorm trans1(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW);
			trans1.transform(ms2,ms2,M1);

			// on col
			CPTransBasisNorm trans2(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_COL);
			trans2.transform(ms2,ms2,M1);

			// now read in output result
			Mtrx M2(ncarbas,ncarbas);
			string out = "cptrans_test/cart2.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),ncarbas,ncarbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " norm on row/column shell with matrix failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " norm on row/column shell with matrix passed with accuracy of " << thresh << std::endl;
			}
		}

		// now on list
		/*
		if (onList) {

			// now create the input
			UInt ncarbas = ms2.getNCarBas();
			Mtrx M1(ncarbas,ncarbas);
			string in = "cptrans_test/cart.txt";
			FileReadWrite filerw1(in);
			filerw1.readMatrixFromTextFile(M1.getPtr(),ncarbas,ncarbas);

			// next, let's creat the matrix list
			BlockMtrxList list;
			for(UInt jAtom=0; jAtom<ms2.getNAtomShells(); jAtom++) {
				const AtomShell& cas = ms2.getAtomShell(jAtom);
				UInt pCol = cas.getBasisStartIndex(TYPE_CART);
				UInt nCol = cas.getNCarBas();
				for(UInt iAtom=jAtom; iAtom<ms2.getNAtomShells(); iAtom++) {
					const AtomShell& ras = ms2.getAtomShell(iAtom);
					UInt pRow = ras.getBasisStartIndex(TYPE_CART);
					UInt nRow = ras.getNCarBas();

					// now let's build the block matrix
					// all of the lower-triangular part
					bool inTrans = false;
					BlockMtrx B(nRow,nCol,pRow,pCol,inTrans);

					// block matrix also need the proper index
					B.initPropIndex(iAtom,jAtom);

					// form B data and push it to list
					B.getData(M1);
					list.update(B,false);
				}
			}

			// on row
			CPTransBasisNorm trans1(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_ROW);
			trans1.transform(ms2,ms2,list);

			// on col
			CPTransBasisNorm trans2(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF,CP_WITH_COL);
			trans2.transform(ms2,ms2,list);

			// update list into the M1 matrix
			bool reform = true;
			M1.reset(ncarbas,ncarbas,reform);
			list.updateMatrix(M1);

			// now read in output result
			Mtrx M2(ncarbas,ncarbas);
			string out = "cptrans_test/cart2.txt";
			FileReadWrite filerw2(out);
			filerw2.readMatrixFromTextFile(M2.getPtr(),ncarbas,ncarbas);

			// now do compare work
			Double thresh = 1.0E-5;
			bool fail = M2.diffComp(M1,thresh,true);
			if (fail) {
				std::cout << " normalization on row/col with list is failed with accuracy of " << thresh << std::endl;
			}else{
				std::cout << " normalization on row/col with list is passed with accuracy of " << thresh << std::endl;
			}
		}
		*/
	}

	////////////////////////////////////////////
	//  @@@@  CP trans with matrix transpose  //
	////////////////////////////////////////////
	Molecule mol3(inf,3);
	MolShell ms3(inf,mol3);
	if (doMatrixTranspose) {

		// now create the input
		UInt nbas = ms3.getNBas();
		Mtrx M(nbas,nbas);
		M.set(ONE);

		// now do P2C, but with matrix transpose
		CPTransBasisNorm trans1(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW_COL);
		trans1.transform(ms3,ms3,M);

		// now let's do C2P, with matrix transpose too
		CPTransBasisNorm trans2(infor,P2C_WITH_L00,UNDO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW_COL);
		trans2.transform(ms3,ms3,M);

		// now do compare work
		Double thresh = 1.0E-7;
		Mtrx M2(nbas,nbas);
		M2.set(ONE);
		bool fail = M2.diffComp(M,thresh,true);
		if (fail) {
			std::cout << " testing with transpose data with matrix failed with accuracy of " << thresh << std::endl;
		}else{
			std::cout << " testing with transpose data with matrix passed with accuracy of " << thresh << std::endl;
		}
	}

	/////////////////////////////////////
	//  @@@@  C2P work on special case //
	/////////////////////////////////////

	if (specialCaseTest) {

		// set up shell and mol data
		Molecule mol4(inf,4);
		MolShell ms4(inf,mol4);
		UInt nbas     = ms4.getNBas();
		UInt ncarbas  = ms4.getNCarBas();

		std::cout << "////////////////////////////////////////" << std::endl;
		std::cout << "// this is a special case for testing //" << std::endl;
		std::cout << "// originally we failed this case     //" << std::endl;
		std::cout << "////////////////////////////////////////" << std::endl;

		// now create the input
		Mtrx M1(ncarbas,ncarbas);
		M1.init(nbas,nbas);
		string in = "cptrans_test/den_pure.txt";
		FileReadWrite filerw1(in);
		filerw1.readMatrixFromTextFile(M1.getPtr(),nbas,nbas);

		// do C2P transpose
		CPTransBasisNorm trans(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW_COL);
		trans.transform(ms4,ms4,M1);
		M1.copyLowToUpper();

		// now read in output result
		Mtrx M2(ncarbas,ncarbas);
		string out = "cptrans_test/den_cart.txt";
		FileReadWrite filerw2(out);
		filerw2.readMatrixFromTextFile(M2.getPtr(),ncarbas,ncarbas);

		// now do compare work
		Double thresh = 1.0E-5;
		bool fail = M2.diffComp(M1,thresh,true);
		if (fail) {
			std::cout << " C2P trans on shell pair with matrix data is failed with accuracy of " << thresh << std::endl;
		}else{
			std::cout << " C2P trans on shell pair with matrix data is passed with accuracy of " << thresh << std::endl;
		}

		// now reset
		M1.init(nbas,nbas);
		filerw1.readMatrixFromTextFile(M1.getPtr(),nbas,nbas);

		// do C2P transpose on row
		CPTransBasisNorm trans1(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW);
		trans1.transform(ms4,ms4,M1);

		// do C2P transpose on row
		CPTransBasisNorm trans2(infor,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_COL);
		trans2.transform(ms4,ms4,M1);

		// now do compare work
		fail = M2.diffComp(M1,thresh,true);
		if (fail) {
			std::cout << " C2P trans on row/col with matrix data is failed with accuracy of " << thresh << std::endl;
		}else{
			std::cout << " C2P trans on row/col with matrix data is passed with accuracy of " << thresh << std::endl;
		}
	}
}
