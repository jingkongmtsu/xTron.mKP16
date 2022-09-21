/**
 * cpp file related to scfdiis.h
 * \author fenglai liu 
 */
#include<iostream>
#include<string>
#include<boost/lexical_cast.hpp>
#include "blas.h"
#include "blas1.h"
#include "excep.h"
#include "shell.h"
#include "lapack.h"
#include "denmtrx.h"
#include "scfparam.h"
#include "oneemtrx.h"
#include "spinmatrix.h"
#include "scfdiis.h"
using namespace blas;
using namespace excep;
using namespace lapack;
using namespace denmtrx;
using namespace scfparam;
using namespace oneemtrx;
using namespace spinmatrix;
using namespace scfdiis;
using namespace std;

SCFDIIS::SCFDIIS(const SCFParam& par):infor(par.getGlobalInfor()),error(ZERO),
	errMtrx(par.getMaxSCFCycles(),par.getMaxSCFCycles()),
	errVecList(par.getGlobalInfor(),"scfdiis",par.getSec(),par.getNSpin(),par.useFile()) 
{ }

void SCFDIIS::errorVec(const Mtrx& S, const Mtrx& O, const DenMtrx& den, 
		const SpinMatrix& fock, SpinMatrix& errVec) const
{
	//
	//error vector is in matrix form of E = F*P*S - S*P*F
	//finally, the result error vector must be orthogonalized, that is to say:
	// E' = O^{T}EO
	//we note, that fock, den and ov should be symmetrical matrix
	//
	if (! fock.isSquare()) {
		string infor = "fock matrix must be square in deriving the error vector"; 
		Excep excep("SCFDIIS","errorVec",EXCEPTION_SCFDIIS_ERROR,infor);
		handleExcep(excep);
	}
	if (! den.isSquare()) {
		string infor = "density matrix must be square in deriving the error vector"; 
		Excep excep("SCFDIIS","errorVec",EXCEPTION_SCFDIIS_ERROR,infor);
		handleExcep(excep);
	}
	if (! S.isSquare()) {
		string infor = "two body overlap matrix must be square in deriving the error vector"; 
		Excep excep("SCFDIIS","errorVec",EXCEPTION_SCFDIIS_ERROR,infor);
		handleExcep(excep);
	}

	// set up some tmp matrix
	UInt n = S.getRow();
	Mtrx tmp(n,n);

	// now work on each spin component matrix
	for(UInt iSpin=0; iSpin<fock.getNSpin(); iSpin++) {

		// get data matrix
		const Mtrx& F = fock.getMtrx(iSpin);
		const Mtrx& P = den.getMtrx(iSpin);
		Mtrx&       E = errVec.getMtrx(iSpin);

		// still check the dimension: F*P, (FP)*S and S
		if (F.getCol() != P.getRow()) {
			string infor = "dimension conflicts: Fock's col should equal to P's row"; 
			Excep excep("SCFDIIS","errorVec",EXCEPTION_SCFDIIS_ERROR,infor);
			handleExcep(excep);
		}
		if (S.getRow() != P.getCol()) {
			string infor = "dimension conflicts: (FP)'s col should equal to S's row"; 
			Excep excep("SCFDIIS","errorVec",EXCEPTION_SCFDIIS_ERROR,infor);
			handleExcep(excep);
		}

		// tmp = FP
		// just in case that only lower part of F is available
		tmp.symMatrixMult(F,P,'L',ONE,ZERO);
		//tmp.print("FP");

		// E = (FP)*S
		E.symMatrixMult(S,tmp,'R',ONE,ZERO);
		//E.print("FPS");

		// tmp = -(FPS)^{T} = -SPF
		// the square matrix transpose is easy
		tmp.copyMatrix(E);
		tmp.transpose(true);
		tmp.scale(MINUS_ONE);

		// E = FPS - SPF 
		E.add(tmp);
		//E.print("FPS-SPF");

		// now transform the error matrix into orthogonal basis set
		bool withRedundantBasis = false;
		if (O.getCol()!=O.getRow()) withRedundantBasis = true;

		// tmp = EO
		if (withRedundantBasis) {
			bool formData = false;
			UInt nMO  = O.getCol();
			UInt nBas = F.getRow();
			tmp.reset(nBas,nMO,formData);
		}
		tmp.mult(E,O,false,false,ONE,ZERO);
		//tmp.print("EO");

		// E = O^{T}EO
		if (withRedundantBasis) {
			bool formData = false;
			UInt nMO = O.getCol();
			E.reset(nMO,nMO,formData);
		}
		E.mult(O,tmp,true,false,ONE,ZERO);
		//E.print("final error matrix");

		// reset the tmp matrix dimension
		if (withRedundantBasis) {
			bool formData = true;
			tmp.reset(n,n,formData);
		}
	}
}

void SCFDIIS::updateError(const SpinMatrix& errVec) 
{
	//
	// here we compute the the maximum element
	// based on the errVec
	//
	// for close shell, this is just the maximum element
	// of errVec
	//
	// for open shell, we take the maximum value 
	// of max(fabs(alphaErrVec[i])+fabs(betaErrVec[i]))/2
	//
	UInt nSpin = errVec.getNSpin();
	error = ZERO;
	if (nSpin == 1) {
		const Mtrx& M = errVec.getMtrx(0);
		const DoubleVec& err = M.getVec();
		error = maxSearch(&err[0],err.size());
	}else{

		// get the data
		const Mtrx& M0 = errVec.getMtrx(0);
		const Mtrx& M1 = errVec.getMtrx(1);
		const DoubleVec& err0 = M0.getVec();
		const DoubleVec& err1 = M1.getVec();

		// check dimension
		if (err0.size() != err1.size()) {
			string infor = "dimension conflicts: the first error vector should have same dimension with second"; 
			Excep excep("SCFDIIS","updateError",EXCEPTION_SCFDIIS_ERROR,infor);
			handleExcep(excep);
		}

		// now compute error
		for(UInt i=0; i<err0.size(); i++) {
			Double val = (fabs(err0[i])+fabs(err1[i]))/TWO;
			if (val>error) error=val;
		}
	}
}

void SCFDIIS::formErrMtrx(const SpinMatrix& errVec) 
{
	// for the open shell, we do error matrix as
	// E = (E_alpha + E_beta)/2
	// so to take the average value 
	// for close shell, E = E_alpha
	UInt nSpin = errVec.getNSpin();
	Double fac = ONE;
	if (nSpin == 2) fac = HALF;

	// for first round of data, we do update
	if (errVecList.nData() == 0) {
		for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

			// form the error matrix
			const Mtrx& M = errVec.getMtrx(iSpin);
			const DoubleVec& err = M.getVec();
			Double errVal = vdot(&err.front(),&err.front(),err.size());
			errMtrx(0,0) += fac*errVal;

			// finally, store this error vec into the history
			errVecList.storeData(err,iSpin);
		}
		return;
	}

	// let's store the new error vector first
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		const Mtrx& M = errVec.getMtrx(iSpin);
		const DoubleVec& err = M.getVec();
		errVecList.storeData(err,iSpin);
	}

	// now let's do normal updating work - items 
	// related to the error matrix
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// input new error vector
		const Mtrx& M = errVec.getMtrx(iSpin);
		const DoubleVec& err = M.getVec();
		DoubleVec oldErr(err.size());

		// update the data inside error matrix
		UInt lastOne = errVecList.nData()-1;
		for(UInt i=0; i<errVecList.nData(); i++) {

			// get the old error vector and create the matrix element
			oldErr.assign(oldErr.size(),ZERO);
			errVecList.retrieveData(oldErr,i,iSpin);
			Double val = vdot(&oldErr.front(), &err.front(), err.size());

			// update matrix
			if (i != lastOne) {
				errMtrx(i,lastOne) += fac*val;
				errMtrx(lastOne,i) += fac*val;
			}else{
				errMtrx(lastOne,lastOne) += fac*val;
			}
		}
	}
}

void SCFDIIS::updateCoefs(const UIntVec& scfIndexArray, DoubleVec& coe) const 
{
	// for the first round of DIIS, that is to say;
	// the len is 1; then the coefficient is one
	// we do not need to solve it
	UInt len = scfIndexArray.size();
	coe.assign(len,ZERO);
	if (len == 1) {
		coe[0] = ONE;
		return;
	}

	// initilize the A to be error matrix
	// the result error matrix should have
	// the last column and last row as -1
	Mtrx A(len+1,len+1);
	for(UInt i=0; i<len+1; i++) {
		A(i,len) = MINUS_ONE;
		A(len,i) = MINUS_ONE;
	}
	A(len,len) = ZERO;

	// we double check the scf index array
	// it's length should be within the current data size
	UInt size = errVecList.nData();
	if (len>size) {
		string infor = "dimension conflicts: the input scf index array has size > error vector list size"; 
		Excep excep("SCFDIIS","updateCoefs",EXCEPTION_SCFDIIS_ERROR,infor);
		handleExcep(excep);
	}

	// now let's form the error matrix A
	// copy the content from error matrix
	for(UInt j=0; j<len; j++) {
		UInt errMtrxColIndex = scfIndexArray[j];

		// double check the index
		// because the row index will take data also from scfIndexArray
		// we do not need to check it twice
		if (errMtrxColIndex>=size) {
			cout << "errMtrxColIndex: " << errMtrxColIndex << endl;
			string infor = "the errMtrxColIndex from scfIndexArray is invalid"; 
			Excep excep("SCFDIIS","updateCoefs",EXCEPTION_SCFDIIS_ERROR,infor);
			handleExcep(excep);
		}

		// form data
		for(UInt i=0; i<len; i++) {
			UInt errMtrxRowIndex = scfIndexArray[i];
			A(i,j) = errMtrx.val(errMtrxRowIndex,errMtrxColIndex);
		}
	}

	// before solving the matrix of A, we need to check that 
	// whether A is singular
	// we will use the svd method to do that
	//if (A.rank(linearDepThresh) < len+1) {
	//}
	
	// now solve the matrix equation of ErrMtrx*C = [0,0,0.. -1]
	DoubleVec coefs(len+1,ZERO);
	coefs[len] = MINUS_ONE;
	symAXeqB('L',A.getPtr(),len+1,&coefs.front(),1);

	// let's write the result coefs back
	for(UInt i=0; i<len; i++) {
		coe[i] = coefs[i];
	}
}

void SCFDIIS::updateDIISInfor(const OneEMtrx& oneEMtrx, const DenMtrx& den, const SpinMatrix& fock)
{
	// let's firstly construct the error vector
	// because we are going to use the error vector for matrix
	// multiplication, therefore the error vector will be in
	// the same dimension of Fock matrix 
	// however, in deriving the error vectors it's dimension
	// will be changed into nDimErrVec later
	UInt nRowOrtho  = 0;
	UInt nDimErrVec = 0;
	oneEMtrx.getDim(ORTHOGONAL_MATRIX,nRowOrtho,nDimErrVec);
	SpinMatrix errVec(fock.getNSpin(),nRowOrtho,nRowOrtho);

	// we need to rebuild the overlap matrix as well as the orthogonal matrix
	bool generateErrVec = true;
	if (generateErrVec) {

		// get the overlap matrix
		UInt nRow = -1;
		UInt nCol = -1;
		oneEMtrx.getDim(TWO_BODY_OVERLAP,nRow,nCol);
		Mtrx S(nRow,nCol);
		oneEMtrx.getM(TWO_BODY_OVERLAP,S);

		// get the orthogonal matrix
		Mtrx O(nRowOrtho,nDimErrVec);
		oneEMtrx.getM(ORTHOGONAL_MATRIX,O);

		// produce the error vector
		errorVec(S,O,den,fock,errVec);
		//errVec.print("inside doDIIS, FPS-SPF");
	}

	// now compute the error
	updateError(errVec); 

	// form the error matrix
	formErrMtrx(errVec);
}

void SCFDIIS::print() const
{
	// print out the error matrix
	UInt n = errVecList.nData();
	if (n>0) {
		Mtrx tmp(n,n);
		tmp.copyFromMatrix(0,0,0,0,n,n,errMtrx);
		tmp.print("the diis error matrix");
	}
}

