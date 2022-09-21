/**
 * \file    dftmatrix.cpp
 * \brief   produce matrix result based on signifacant basis set order
 * \author  Fenglai Liu and Jing Kong
 */
#include<iostream>
#include "blas.h"
#include "blas1.h"
#include "shell.h"
#include "excep.h"
#include "xcintsinfor.h"
#include "sigatombasis.h"
#include "dftmatrix.h"
using namespace blas;
using namespace shell;
using namespace excep;
using namespace xcintsinfor;
using namespace sigatombasis;
using namespace dftmatrix;
using namespace std;

DFTMatrix::DFTMatrix(const MolShell& rowShell, const SigAtomBasis& colSigList, 
		const UInt& nSpin):SpinMatrix(nSpin,rowShell.getNBas(),colSigList.getNSigBasis()),
	status(COL_IN_SIG_ORDER)
{ }

DFTMatrix::DFTMatrix(const SigAtomBasis& rowSigList, const MolShell& colShell,
		const UInt& nSpin):SpinMatrix(nSpin,rowSigList.getNSigBasis(),colShell.getNBas()),
	status(ROW_IN_SIG_ORDER)
{ }

DFTMatrix::DFTMatrix(const SigAtomBasis& rowSigList, const SigAtomBasis& colSigList,
		const UInt& nSpin):SpinMatrix(nSpin,rowSigList.getNSigBasis(),colSigList.getNSigBasis()),
	status(ROW_COL_IN_SIG_ORDER)
{ }

DFTMatrix::DFTMatrix( const SigAtomBasis& rowSigList, const SigAtomBasis& colSigList,
		const UInt& rowAtom, const UInt& colAtom, 
		const UInt& nSpin):SpinMatrix(nSpin,rowSigList.getNSigBasis(rowAtom),colSigList.getNSigBasis(colAtom)),
	status(ROW_COL_IN_SIG_ORDER)
{ }

void DFTMatrix::intoSigOrder(const UInt& iSpin, const Mtrx& source, 
		const SigAtomBasis& rowBasis, const SigAtomBasis& colBasis) 
{

	// get the sig basis index list
	// for the row matrix
	UInt nRowSigBasis = rowBasis.getNSigBasis();
	const UIntVec& sigRowBasisIndex = rowBasis.getGlobalSigBasisIndex();

	// now for the column part
	UInt nColSigBasis = colBasis.getNSigBasis();
	const UIntVec& sigColBasisIndex = colBasis.getGlobalSigBasisIndex();

	// additional checking for dimension
	if (iSpin>=mtrxList.size()) {
		string infor = "the given spin state is out of NSpin in DFT matrices";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	UInt row = mtrxList[iSpin].getRow();
	UInt col = mtrxList[iSpin].getCol();
	if (nRowSigBasis != row || nColSigBasis != col) {
		string infor = "row or col dimension does not match";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}

	// now doing real transformation of data
	Mtrx& goal = mtrxList[iSpin];
	for(UInt j=0; j<nColSigBasis; j++) {
		UInt colIndex = sigColBasisIndex[j];
		for(UInt i=0; i<nRowSigBasis; i++) {
			UInt rowIndex = sigRowBasisIndex[i];
			goal(i,j) += source.val(rowIndex,colIndex);
		}
	}
}

void DFTMatrix::intoNormOrder(const SigAtomBasis& rowBasis, 
		const SigAtomBasis& colBasis, const UInt& iSpin, Mtrx& goal) const
{

	// get the sig basis index list
	// for the row matrix
	UInt nRowSigBasis = rowBasis.getNSigBasis();
	const UIntVec& sigRowBasisIndex = rowBasis.getGlobalSigBasisIndex();

	// now for the column part
	UInt nColSigBasis = colBasis.getNSigBasis();
	const UIntVec& sigColBasisIndex = colBasis.getGlobalSigBasisIndex();

	if (iSpin>=mtrxList.size()) {
		string infor = "the given spin state is out of NSpin in DFT matrices";
		Excep excep("DFTMatrix","intoNormOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	UInt row = mtrxList[iSpin].getRow();
	UInt col = mtrxList[iSpin].getCol();
	if (nRowSigBasis != row || nColSigBasis != col) {
		string infor = "row or col dimension does not match";
		Excep excep("DFTMatrix","intoNormOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}

	const Mtrx& source = mtrxList[iSpin];
	for(UInt j=0; j<nColSigBasis; j++) {
		UInt colIndex = sigColBasisIndex[j];
		for(UInt i=0; i<nRowSigBasis; i++) {
			UInt rowIndex = sigRowBasisIndex[i];
			goal(rowIndex,colIndex) += source.val(i,j);
		}
	}
}

void DFTMatrix::intoSigOrder(const UInt& iSpin, const Mtrx& source, const SigAtomBasis& b) 
{

	// information
	UInt nSigBasis = b.getNSigBasis();
	const UIntVec& sigBasisIndex = b.getGlobalSigBasisIndex();

	// check the status
	if (iSpin>=mtrxList.size()) {
		string infor = "the given spin state is out of NSpin in DFT matrices";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	UInt row = mtrxList[iSpin].getRow();
	UInt col = mtrxList[iSpin].getCol();
	if (status == ROW_COL_IN_SIG_ORDER) {
		string infor = "wrong matrix status inside, we only do row or col conversion";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	if (status == ROW_IN_SIG_ORDER && nSigBasis != row && source.getCol() != col) {
		string infor = "row/col dimension does not match, row data converted";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	if (status == COL_IN_SIG_ORDER && nSigBasis != col && source.getRow() != row) {
		string infor = "row/col dimension does not match, col data converted";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}

	// matrix
	Mtrx& goal = mtrxList[iSpin];
	if(status == ROW_IN_SIG_ORDER) {
		for(UInt j=0; j<source.getCol(); j++) {
			for(UInt i=0; i<nSigBasis; i++) {
				UInt rowIndex = sigBasisIndex[i]; 
				goal(i,j) += source.val(rowIndex,j);
			}
		}
	}else {
		UInt len = source.getRow();
		for(UInt i=0; i<nSigBasis; i++) {
			UInt colIndex = sigBasisIndex[i]; 
			vaxpy(source.getPtr(0,colIndex),goal.getPtr(0,i),ONE,len);
		}
	}
}

void DFTMatrix::intoNormOrder(const SigAtomBasis& b, const UInt& iSpin, Mtrx& goal) const
{
	// information
	UInt nSigBasis = b.getNSigBasis();
	const UIntVec& sigBasisIndex = b.getGlobalSigBasisIndex();

	// check the status
	if (iSpin>=mtrxList.size()) {
		string infor = "the given spin state is out of NSpin in DFT matrices";
		Excep excep("DFTMatrix","intoNormOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	UInt row = mtrxList[iSpin].getRow();
	UInt col = mtrxList[iSpin].getCol();
	if (status == ROW_COL_IN_SIG_ORDER) {
		string infor = "wrong matrix status inside, we only do row or col conversion";
		Excep excep("DFTMatrix","intoNormOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	if (status == ROW_IN_SIG_ORDER && nSigBasis != row && goal.getCol() != col) {
		string infor = "row/col dimension does not match, row data converted";
		Excep excep("DFTMatrix","intoNormOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	if (status == COL_IN_SIG_ORDER && nSigBasis != col && goal.getRow() != row) {
		string infor = "row/col dimension does not match, col data converted";
		Excep excep("DFTMatrix","intoNormOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}

	// real work here
	const Mtrx& source = mtrxList[iSpin];
	if(status == ROW_IN_SIG_ORDER) {
		for(UInt j=0; j<source.getCol(); j++) {
			for(UInt i=0; i<nSigBasis; i++) {
				UInt rowIndex = sigBasisIndex[i]; 
				goal(rowIndex,j) += source.val(i,j);
			}
		}
	}else {
		UInt len = source.getRow();
		for(UInt i=0; i<nSigBasis; i++) {
			UInt colIndex = sigBasisIndex[i]; 
			vaxpy(source.getPtr(0,i),goal.getPtr(0,colIndex),ONE,len);
		}
	}
}

void DFTMatrix::intoSigOrder(const UInt& iSpin, const Mtrx& source, const AtomShell& rowAtomShell, 
		const AtomShell& colAtomShell, const SigAtomBasis& rowBasis, 
		const SigAtomBasis& colBasis, const UInt& iRowSigAtom,
		const UInt& iColSigAtom) 
{

	// get the sig basis dimension
	UInt nRowSigBasis = rowBasis.getNSigBasis(iRowSigAtom);
	UInt nColSigBasis = colBasis.getNSigBasis(iColSigAtom);

	// additional check
	if (iSpin>=mtrxList.size()) {
		string infor = "the given spin state is out of NSpin in DFT matrices";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}
	UInt row = mtrxList[iSpin].getRow();
	UInt col = mtrxList[iSpin].getCol();
	if (nRowSigBasis != row || nColSigBasis != col) {
		string infor = "row or col dimension does not match";
		Excep excep("DFTMatrix","intoSigOrder",EXCEPTION_XCINTS_DFTMATRIX_ERROR,infor);
		handleExcep(excep);
	}

	// real work
	Mtrx& goal = mtrxList[iSpin];
	for(UInt iCol=0; iCol<nColSigBasis; iCol++) {
		UInt colIndex = colBasis.getLocalAtomBasisIndex(colAtomShell,iColSigAtom,iCol);
		for(UInt iRow=0; iRow<nRowSigBasis; iRow++) {
			UInt rowIndex = rowBasis.getLocalAtomBasisIndex(rowAtomShell,iRowSigAtom,iRow);
			goal(iRow,iCol) += source.val(rowIndex,colIndex);
		}
	}
}

void DFTMatrix::print() const
{
	UInt nSpin = mtrxList.size();
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		string s;
		if (iSpin == 0) {
			s = "Alpha DFT matrix result";
		}else{
			s = "Beta DFT matrix result";
		}
		mtrxList[iSpin].print(s);
	}
}

DFTVect::DFTVect(const SigAtomBasis& sigList, const UInt& nSpin0):nSpin(nSpin0),
	nData(sigList.getNSigBasis()),vectList(nData*nSpin,ZERO)
{ }

void DFTVect::intoSigOrder(const DoubleVec& v, const SigAtomBasis& b)
{
	UInt nSigBasis = b.getNSigBasis();
	const UIntVec& sigBasisIndex = b.getGlobalSigBasisIndex();
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		UInt offset = iSpin*nSigBasis;
		for(UInt i=0; i<nSigBasis; i++) {
			UInt index = sigBasisIndex[i];
			vectList[i+offset] += v[index+offset];
		}
	}
}

void DFTVect::intoNormOrder(DoubleVec& v, const SigAtomBasis& b) const
{
	UInt nSigBasis = b.getNSigBasis();
	const UIntVec& sigBasisIndex = b.getGlobalSigBasisIndex();
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		UInt offset = iSpin*nSigBasis;
		for(UInt i=0; i<nSigBasis; i++) {
			UInt index = sigBasisIndex[i];
			v[index] += vectList[i+offset];
		}
	}
}

void DFTVect::print() const
{
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		if (iSpin == 0) {
			cout << "Alpha DFT Vector results:" << endl;
		}else{
			cout << "Beta  DFT Vector results:" << endl;
		}
		UInt offset = iSpin*nData;
		for(UInt i=0; i<nData; i++) {
			printf("%-12d  %-12.6f\n", (Int)i, vectList[i+offset]);
		}
	}
}
