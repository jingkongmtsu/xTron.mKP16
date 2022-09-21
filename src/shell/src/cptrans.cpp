/**
 * cpp file corresponding to cptrans.h
 * author fenglai liu
 */
#include <string>
#include <boost/lexical_cast.hpp>
#include "shellprop.h"
#include "shell.h"
#include "blas.h"
#include "blas1.h"
#include "omp_setting.h"
#include "excep.h"
#include "cptrans.h"
using namespace std;
using namespace shellprop;
using namespace shell;
using namespace blas;
using namespace omp_setting;
using namespace excep;
using namespace cptrans;

//////////////////////////////////
//      @@@@  CPTransData       //
//////////////////////////////////
CPTransData::CPTransData(const UInt& maxL0, const UInt& transWork, const UInt& scaleWork, 
		const UInt& matrixStatus):maxL(maxL0)
{
	//
	// transformation work
	//
	if (transWork != NO_TRANS && transWork != C2P_WITH_L00 && 
			transWork != C2P_WITH_XYZ && transWork != P2C_WITH_L00 && transWork != P2C_WITH_XYZ) {
		string infor = "invalid input parameters for transWork"; 
		Excep excep("CPTransData","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	//
	// scale work
	//
	if (scaleWork != NO_SCALE && scaleWork != DO_SCALE && scaleWork != UNDO_SCALE) { 
		string infor = "invalid input parameters for scaleWork"; 
		Excep excep("CPTransData","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	//
	// matrix status
	//
	if (matrixStatus != WITH_MATRIX_ITSELF && matrixStatus != WITH_MATRIX_TRANSPOSE) {
		string infor = "invalid input parameters for matrixStatus"; 
		Excep excep("CPTransData","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	// check maxL
	if (maxL>CPTRANS_MAX_SHELL) {
		string infor = "cptrans class only support maxL up to L = 6"; 
		Excep excep("CPTransData","constructor",EXCEPTION_UNSUPPORTED_ANGULAR_MOMENTUM,infor);
		handleExcep(excep);
	}

	// form the converting matrix
	if (transWork != NO_TRANS) {
		for(UInt L=2; L<=maxL; L++) {
			if (L == 2) {
				formConvertMatrix(L,transWork,matrixStatus,convertMtrx_D);
			}else if (L == 3) {
				formConvertMatrix(L,transWork,matrixStatus,convertMtrx_F);
			}else if (L == 4) {
				formConvertMatrix(L,transWork,matrixStatus,convertMtrx_G);
			}else if (L == 5) {
				formConvertMatrix(L,transWork,matrixStatus,convertMtrx_H);
			}else if (L == 6) {
				formConvertMatrix(L,transWork,matrixStatus,convertMtrx_I);
			}
		}
	}

	// form normalization vector if needed
	if (scaleWork != NO_SCALE) {
		for(UInt L=2; L<=maxL; L++) {
			if (L == 2) {
				formConvertVector(L,scaleWork,scaleVec_D);
			}else if (L == 3) {
				formConvertVector(L,scaleWork,scaleVec_F);
			}else if (L == 4) {
				formConvertVector(L,scaleWork,scaleVec_G);
			}else if (L == 5) {
				formConvertVector(L,scaleWork,scaleVec_H);
			}else if (L == 6) {
				formConvertVector(L,scaleWork,scaleVec_I);
			}
		}
	}
}

void CPTransData::initScaleData()
{
	UInt scaleWork = DO_SCALE;
	for(UInt L=2; L<=maxL; L++) {
		if (L == 2) {
			formConvertVector(L,scaleWork,scaleVec_D);
		}else if (L == 3) {
			formConvertVector(L,scaleWork,scaleVec_F);
		}else if (L == 4) {
			formConvertVector(L,scaleWork,scaleVec_G);
		}else if (L == 5) {
			formConvertVector(L,scaleWork,scaleVec_H);
		}else if (L == 6) {
			formConvertVector(L,scaleWork,scaleVec_I);
		}
	}
}

void CPTransData::formConvertMatrix(const UInt& L, const UInt& transWork, 
		const UInt& matrixStatus, Mtrx& convert)
{
	// set up the matrix
	UInt ncb = getCartBas(L,L);
	UInt npb = getPureBas(L,L);
	if (transWork == C2P_WITH_L00 || transWork == C2P_WITH_XYZ) {
		convert.init(ncb,npb);
	}else{
		convert.init(npb,ncb);
	}

	// now fill in data
	switch(transWork) {
		case C2P_WITH_L00:
			getC2PFromL00(L,convert.getPtr());
			break;
		case C2P_WITH_XYZ:
			getC2PFromLxLyLz(L,convert.getPtr()); 
			break;
		case P2C_WITH_L00:
			getP2CFromL00(L,convert.getPtr()); 
			break;
		case P2C_WITH_XYZ:
			getP2CFromLxLyLz(L,convert.getPtr()); 
			break;
	}

	// do matrix transpose?
	if (matrixStatus == WITH_MATRIX_TRANSPOSE) {
		bool inplace = true;
		convert.transpose(inplace);
	}
}

void CPTransData::formConvertVector(const UInt& L, const UInt& scaleWork, DoubleVec& convert)
{
	// set up the vector
	UInt ncb = getCartBas(L,L);
	convert.assign(ncb,ZERO);
	getNormScaleVector(L,&convert.front()); 

	//
	// if this is undo the work
	// we need to further do some work
	//
	if (scaleWork == UNDO_SCALE) {
		for(UInt i=0; i<ncb; i++) {
			convert[i] = ONE/convert[i];
		}
	}
}

const Mtrx& CPTransData::getConvertMatrix(const UInt& L) const 
{
	// check L
	if (L>maxL) {
		string infor = "C2P/P2C matrix is not available"; 
		Excep excep("CPTransData","constructor",EXCEPTION_CPTRANS_DATA_NOT_AVAILABLE,infor);
		handleExcep(excep);
	}

	// now get the matrix
	switch(L) {
		case 2:
			return convertMtrx_D;
			break;
		case 3:
			return convertMtrx_F;
			break;
		case 4:
			return convertMtrx_G;
			break;
		case 5:
			return convertMtrx_H;
			break;
		case 6:
			return convertMtrx_I;
			break;
		default:
			string infor = "input L:" + boost::lexical_cast<string>(L); 
			Excep excep("CPTransData","constructor",EXCEPTION_CPTRANS_DATA_NOT_AVAILABLE,infor);
			handleExcep(excep);
			break;
	}

	// default, we return something to avoid compiler warning
	return convertMtrx_D;
}

const DoubleVec& CPTransData::getConvertVec(const UInt& L) const 
{
	// check L
	if (L>maxL) {
		string infor = "basis set normalization vector is not available"; 
		Excep excep("CPTransData","constructor",EXCEPTION_CPTRANS_DATA_NOT_AVAILABLE,infor);
		handleExcep(excep);
	}

	// now get the matrix
	switch(L) {
		case 2:
			return scaleVec_D;
			break;
		case 3:
			return scaleVec_F;
			break;
		case 4:
			return scaleVec_G;
			break;
		case 5:
			return scaleVec_H;
			break;
		case 6:
			return scaleVec_I;
			break;
		default:
			string infor = "input L:" + boost::lexical_cast<string>(L); 
			Excep excep("CPTransData","constructor",EXCEPTION_CPTRANS_DATA_NOT_AVAILABLE,infor);
			handleExcep(excep);
			break;
	}

	// default, we return something
	return scaleVec_D;
}

///////////////////////////////////////////////
//      @@@@  CPTransAtomShellPair           //
///////////////////////////////////////////////
CPTransAtomShell::CPTransAtomShell(const UInt& maxL, const UInt& transWork, 
		const UInt& matrixStatus):sourceType(TYPE_NORM),
	targetType(TYPE_NORM),c2p(maxL,transWork,NO_SCALE,matrixStatus)
{
	//
	// determine the source and target type
	// we can derive it from the input work type and matrix status
	//
	if (transWork == C2P_WITH_L00 || transWork == C2P_WITH_XYZ) {
		sourceType = TYPE_CART;
	}else if (transWork == P2C_WITH_L00 || transWork == P2C_WITH_XYZ) {
		targetType = TYPE_CART;
	}

	//
	// if the matrix is actually used as transpose form, then
	// the source and target types must reverse
	//
	if (matrixStatus == WITH_MATRIX_TRANSPOSE) {
		UInt tmp = sourceType;
		sourceType = targetType;
		targetType = tmp;
	}
}

void CPTransAtomShell::cpTransformOnAtomShellPair(const AtomShell& rs, 
		const AtomShell& cs, const Mtrx& S, Mtrx& tmp, Mtrx& T) const
{
	//
	// firstly, let's copy the S/P/SP data block from S to T
	//
	for(UInt jShell=0; jShell<cs.getNShell(); jShell++) {
		const Shell& js = cs.getShell(jShell);
		UInt L2 = js.getLmax();
		if (L2 >= 2) continue;
		UInt sourceColPos  = js.getLocalBasisIndex(0,sourceType);
		UInt targetColPos  = js.getLocalBasisIndex(0,targetType);
		UInt nCol = js.getNBas();
		for(UInt iShell=0; iShell<rs.getNShell(); iShell++) {
			const Shell& is = rs.getShell(iShell);
			UInt L1 = is.getLmax();
			if (L1 >= 2) continue;
			UInt sourceRowPos  = is.getLocalBasisIndex(0,sourceType);
			UInt targetRowPos  = is.getLocalBasisIndex(0,targetType);
			UInt nRow = is.getNBas();
			S.copyToMatrix(sourceRowPos,sourceColPos,targetRowPos,targetColPos,nRow,nCol,T);
		}
	}

	//
	// next, let's deal with the highL-SP block 
	// highL is from rs, and SP block is from cs
	//
	for(UInt jShell=0; jShell<cs.getNShell(); jShell++) {
		const Shell& js = cs.getShell(jShell);
		UInt L2 = js.getLmax();
		if (L2 >= 2) continue;
		UInt sourceColPos  = js.getLocalBasisIndex(0,sourceType);
		UInt targetColPos  = js.getLocalBasisIndex(0,targetType);
		UInt nCol = js.getNBas();
		for(UInt iShell=0; iShell<rs.getNShell(); iShell++) {

			// get shell data
			const Shell& is = rs.getShell(iShell);
			UInt L1 = is.getLmax();
			if (L1 < 2) continue;

			// now let's perform the job for the two shells
			UInt sourceRowPos  = is.getLocalBasisIndex(0,sourceType);
			UInt targetRowPos  = is.getLocalBasisIndex(0,targetType);
			UInt nSourceRows   = is.getNBas(sourceType);

			//
			// here if i shell is pure, then the c2p/p2c operation is dorable
			// else we just copy the data block
			//
			if (is.isPure()) {

				// choose the convert matrix 
				// T = C^T*S
				// we call mmul directly
				const Mtrx& C = c2p.getConvertMatrix(L1);
				mmul(C.getPtr(),S.getPtr(sourceRowPos,sourceColPos),
						T.getPtr(targetRowPos,targetColPos),C.getRow(),C.getCol(),
						nSourceRows,nCol,C.getLd(),S.getLd(),T.getLd(),
						'T','N',ONE,ZERO);
			}else{
				S.copyToMatrix(sourceRowPos,sourceColPos,targetRowPos,targetColPos,nSourceRows,nCol,T);
			}
		}
	}

	// now do SP-highL block 
	for(UInt iShell=0; iShell<rs.getNShell(); iShell++) {

		// get shell data
		const Shell& is = rs.getShell(iShell);
		UInt L1 = is.getLmax();
		if (L1 >= 2) continue;
		UInt sourceRowPos  = is.getLocalBasisIndex(0,sourceType);
		UInt targetRowPos  = is.getLocalBasisIndex(0,targetType);
		UInt nRows   = is.getNBas();
		for(UInt jShell=0; jShell<cs.getNShell(); jShell++) {

			// get the shell data
			const Shell& js = cs.getShell(jShell);
			UInt L2 = js.getLmax();
			if (L2 < 2) continue;

			// now let's perform the job for the two shells
			UInt sourceColPos  = js.getLocalBasisIndex(0,sourceType);
			UInt targetColPos  = js.getLocalBasisIndex(0,targetType);
			UInt nSourceCols   = js.getNBas(sourceType);

			//
			// here if j shell is pure, then the c2p/p2c operation is dorable
			// else we just copy the data block
			//
			if (js.isPure()) {
				// T = S*C
				const Mtrx& C = c2p.getConvertMatrix(L2);
				mmul(S.getPtr(sourceRowPos,sourceColPos),C.getPtr(),
						T.getPtr(targetRowPos,targetColPos),nRows,nSourceCols,C.getRow(),
						C.getCol(),S.getLd(),C.getLd(),T.getLd(),
						'N','N',ONE,ZERO);
			}else{
				S.copyToMatrix(sourceRowPos,sourceColPos,targetRowPos,targetColPos,nRows,nSourceCols,T);
			}
		}
	}

	//
	// finally, let's process the highL-highL block
	// 
	for(UInt jShell=0; jShell<cs.getNShell(); jShell++) {
		const Shell& js = cs.getShell(jShell);
		if (! js.isHighL()) continue;

		// get the dimension information for js
		UInt sourceColPos  = js.getLocalBasisIndex(0,sourceType);
		UInt targetColPos  = js.getLocalBasisIndex(0,targetType);
		UInt nSourceCols   = js.getNBas(sourceType);

		for(UInt iShell=0; iShell<rs.getNShell(); iShell++) {

			// get the dimension information for row shell
			const Shell& is = rs.getShell(iShell);
			if (! is.isHighL()) continue;
			UInt sourceRowPos  = is.getLocalBasisIndex(0,sourceType);
			UInt targetRowPos  = is.getLocalBasisIndex(0,targetType);
			UInt nSourceRows   = is.getNBas(sourceType);

			// 
			// first case, if both of two shells are Cartesian
			// then we simply do copy work
			//
			if (! is.isPure() && ! js.isPure()) {
				S.copyToMatrix(sourceRowPos,sourceColPos,targetRowPos,targetColPos,
						nSourceRows,nSourceCols,T);
				continue;
			}

			//
			// second case, "is" is cart and "js" is pure
			//
			if (! is.isPure() && js.isPure()) {

				// T = S*C
				// we note, in this case the shell pair data block
				// can not be in the diagonal block
				UInt L = js.getLmax();
				const Mtrx& C = c2p.getConvertMatrix(L);
				mmul(S.getPtr(sourceRowPos,sourceColPos),C.getPtr(),
						T.getPtr(targetRowPos,targetColPos),nSourceRows,nSourceCols,C.getRow(),
						C.getCol(),S.getLd(),C.getLd(),T.getLd(),
						'N','N',ONE,ZERO);
				continue;
			}

			//
			// third case, "is" is pure and "js" is cart
			//
			if (is.isPure() && ! js.isPure()) {

				// T = C^T*S
				UInt L = is.getLmax();
				const Mtrx& C = c2p.getConvertMatrix(L);
				mmul(C.getPtr(),S.getPtr(sourceRowPos,sourceColPos),
						T.getPtr(targetRowPos,targetColPos),C.getRow(),C.getCol(),  
						nSourceRows,nSourceCols,C.getLd(),S.getLd(),T.getLd(),
						'T','N',ONE,ZERO);
				continue;
			}

			//
			// now finally, we need to do
			// T = C1^T*S*C2
			//

			//
			// TMP = S*C2
			//
			UInt L2 = js.getLmax();
			const Mtrx& C2 = c2p.getConvertMatrix(L2);
			tmp.init(nSourceRows,C2.getCol());
			mmul(S.getPtr(sourceRowPos,sourceColPos),C2.getPtr(),
					tmp.getPtr(),nSourceRows,nSourceCols,
					C2.getRow(),C2.getCol(),S.getLd(),C2.getLd(),tmp.getLd(),
					'N','N',ONE,ZERO);

			//
			// T = C1^T*TMP 
			//
			UInt L1 = is.getLmax();
			const Mtrx& C1 = c2p.getConvertMatrix(L1);
			mmul(C1.getPtr(),tmp.getPtr(),
					T.getPtr(targetRowPos,targetColPos),C1.getRow(),C1.getCol(),  
					tmp.getRow(),tmp.getCol(),C1.getLd(),tmp.getLd(),T.getLd(),
					'T','N',ONE,ZERO);
		}
	}
}

void CPTransAtomShell::cpTransformOnAtomShell(const UInt& work, 
		const AtomShell& s, const Mtrx& S, Mtrx& T) const
{
	//
	// firstly, let's copy the S/P/SP data block from S to T
	//
	if (work == CP_WITH_ROW) {
		for(UInt iShell=0; iShell<s.getNShell(); iShell++) {
			const Shell& is = s.getShell(iShell);
			UInt L1 = is.getLmax();
			if (L1 >= 2) continue;
			UInt sourceRowPos  = is.getLocalBasisIndex(0,sourceType);
			UInt targetRowPos  = is.getLocalBasisIndex(0,targetType);
			UInt nRow = is.getNBas();
			S.copyToMatrix(sourceRowPos,0,targetRowPos,0,nRow,S.getCol(),T);
		}
	}else{
		for(UInt jShell=0; jShell<s.getNShell(); jShell++) {
			const Shell& js = s.getShell(jShell);
			UInt L2 = js.getLmax();
			if (L2 >= 2) continue;
			UInt sourceColPos  = js.getLocalBasisIndex(0,sourceType);
			UInt targetColPos  = js.getLocalBasisIndex(0,targetType);
			UInt nCol = js.getNBas();
			S.copyToMatrix(0,sourceColPos,0,targetColPos,S.getRow(),nCol,T);
		}
	}

	//
	// next, let's deal with the highL block 
	//
	for(UInt iShell=0; iShell<s.getNShell(); iShell++) {

		// get shell data
		const Shell& is = s.getShell(iShell);
		UInt L = is.getLmax();
		if (L < 2) continue;

		// now let's perform the job for the two shells
		UInt sourceRowPos  = 0;
		UInt targetRowPos  = 0;
		UInt sourceColPos  = 0;
		UInt targetColPos  = 0;
		UInt nSourceRows   = 0;
		UInt nSourceCols   = 0;
		if (work == CP_WITH_ROW) {
			sourceRowPos  = is.getLocalBasisIndex(0,sourceType);
			targetRowPos  = is.getLocalBasisIndex(0,targetType);
			nSourceRows   = is.getNBas(sourceType);
			nSourceCols   = S.getCol();
		}else{
			sourceColPos  = is.getLocalBasisIndex(0,sourceType);
			targetColPos  = is.getLocalBasisIndex(0,targetType);
			nSourceRows   = S.getRow();
			nSourceCols   = is.getNBas(sourceType);
		}
		//cout << "sourceRowPos " << sourceRowPos << endl;
		//cout << "sourceColPos " << sourceColPos << endl;
		//cout << "targetRowPos " << targetRowPos << endl;
		//cout << "targetColPos " << targetColPos << endl;

		//
		// here if i shell is pure, then the c2p/p2c operation is dorable
		// else we just copy the data block
		//
		if (is.isPure()) {

			if (work == CP_WITH_ROW) {
				// T = C^T*S
				const Mtrx& C = c2p.getConvertMatrix(L);
				mmul(C.getPtr(),S.getPtr(sourceRowPos,sourceColPos),
						T.getPtr(targetRowPos,targetColPos),C.getRow(),C.getCol(),
						nSourceRows,nSourceCols,C.getLd(),S.getLd(),T.getLd(),
						'T','N',ONE,ZERO);
			}else{
				// T = S*C
				const Mtrx& C = c2p.getConvertMatrix(L);
				mmul(S.getPtr(sourceRowPos,sourceColPos),C.getPtr(),
						T.getPtr(targetRowPos,targetColPos),nSourceRows,nSourceCols,C.getRow(),
						C.getCol(),S.getLd(),C.getLd(),T.getLd(),
						'N','N',ONE,ZERO);
			}
		}else{
			S.copyToMatrix(sourceRowPos,sourceColPos,targetRowPos,targetColPos,
					nSourceRows,nSourceCols,T);
		}
	}
}

///////////////////////////////////////////////
//      @@@@     NormAtomShellData           //
///////////////////////////////////////////////
void NormAtomShellData::normRowData(const AtomShell& s, Mtrx& M) 
{
	// do we really need to do the work?
	if (s.getMaxL()<=1 || s.allPure()) return;

	// determine the matrix type for M
	// it's row maybe in Cartesian state, or normal state
	UInt matrixType = -1;
	if (M.getRow() == s.getNCarBas()) {
		matrixType = TYPE_CART;
	}else if (M.getRow() == s.getNBas()) {
		matrixType = TYPE_NORM;
	}else{
		string infor = "for the input matrix, we can not figure out it's shell dimension type"; 
		Excep excep("NormAtomShellData","normRowData", EXCEPTION_SHELL_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	// firstly, we need to form the row vector
	UInt nBas = s.getNBas(matrixType);
	rowVec.assign(nBas,ONE);
	UInt offset = 0;
	for(UInt iShell=0; iShell<s.getNShell(); iShell++) {

		// get shell data
		// we will omit these pure basis and S/P/SP shells
		const Shell& is = s.getShell(iShell);
		UInt L = is.getLmax();
		if (L < 2) {
			offset += is.getNBas(); 
			continue;
		}else if (is.isPure()) {
			offset += is.getNBas(matrixType); 
			continue;
		}

		// get the vector data
		const DoubleVec& vec = scale.getConvertVec(L);
		for(UInt i=0; i<vec.size(); i++) {
			Double N = vec[i];
			rowVec[offset] = N;
			offset++;
		}
	}

	// now do scale work
	for(UInt col=0; col<M.getCol(); col++) {
		vmul(M.getPtr(0,col),&rowVec[0],M.getPtr(0,col),M.getRow());
	}
}

void NormAtomShellData::normColData(const AtomShell& s, Mtrx& M) const
{
	// do we really need to do the work?
	if (s.getMaxL()<=1 || s.allPure()) return;

	// determine the matrix type for M
	// it's row maybe in Cartesian state, or normal state
	UInt matrixType = -1;
	if (M.getCol() == s.getNCarBas()) {
		matrixType = TYPE_CART;
	}else if (M.getCol() == s.getNBas()) {
		matrixType = TYPE_NORM;
	}else{
		string infor = "for the input matrix, we can not figure out it's shell dimension type"; 
		Excep excep("NormAtomShellData","normColData", EXCEPTION_SHELL_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	// for column part of data, we do not need rowVec
	for(UInt iShell=0; iShell<s.getNShell(); iShell++) {

		// get shell data
		// we will omit these pure basis and S/P/SP shells
		const Shell& is = s.getShell(iShell);
		UInt L = is.getLmax();
		if (L < 2 || is.isPure()) {
			continue;
		}

		// get the vector data
		UInt basisIndex = is.getLocalBasisIndex(0,matrixType);
		const DoubleVec& vec = scale.getConvertVec(L);
		for(UInt i=0; i<vec.size(); i++) {
			Double N = vec[i];
			vscal(M.getPtr(0,basisIndex+i),N,M.getRow());
		}
	}
}

///////////////////////////////////////////////
//      @@@@     TBB_CPTransMatrix           //
///////////////////////////////////////////////
TBB_CPTransMatrix::TBB_CPTransMatrix(const MolShell& rowShell, 
		const MolShell& colShell, const UInt& transWork, const UInt& matrixStatus,
		const Mtrx& source, Mtrx& target):rowColStatus(CP_WITH_ROW_COL),
	cptrans(max(rowShell.getMaxL(),colShell.getMaxL()),transWork,matrixStatus),
	rs(rowShell),cs(colShell),M0(source),M1(target)
{

}	

TBB_CPTransMatrix::TBB_CPTransMatrix(const MolShell& s, const UInt& transWork, 
		const UInt& matrixStatus, const UInt& rowColStatus0,
		const Mtrx& source, Mtrx& target):rowColStatus(rowColStatus0),
	cptrans(s.getMaxL(),transWork,matrixStatus),rs(s),cs(s),M0(source),M1(target)
{

}

void TBB_CPTransMatrix::doCPTransOnRowAndCol(const blocked_range<UInt>& r)  const
{
	// set up scratch matrix
	UInt maxRow = rs.getMaxBasis(TYPE_CART);
	UInt maxCol = cs.getMaxBasis(TYPE_CART);
	BlockMtrx S(maxRow,maxCol);
	BlockMtrx T(maxRow,maxCol);

	// set up tmp matrix
	UInt maxL = max(rs.getMaxL(),cs.getMaxL());
	UInt nMaxCarBas = getCartBas(maxL,maxL);
	Mtrx tmp(nMaxCarBas,nMaxCarBas);

	// let's see whether two shells are same?
	bool sameShell = true;
	if (! (rs == cs)) sameShell = false;

	// now loop over each atom shell pair
	for( UInt n=r.begin(); n!=r.end(); ++n ) {

		//
		//get the atom shell data
		//
		UInt rowShellIndex = -1;
		UInt colShellIndex = -1;
		UInt nAtoms = rs.getNAtomShells();
		if (sameShell) {
			getRowColIndexFromPackIndex(nAtoms,n,rowShellIndex,colShellIndex); 
		}else{
			colShellIndex = n/nAtoms;
			rowShellIndex = n-colShellIndex*nAtoms;
		}
		const AtomShell& rowShell = rs.getAtomShell(rowShellIndex);
		const AtomShell& colShell = cs.getAtomShell(colShellIndex);

		//
		//get the data block from source
		//here we do not need to set up the property index again
		//
		bool inTrans = false;
		UInt sourceType = cptrans.getSourceType();
		UInt nRowSource = rowShell.getNBas(sourceType);
		UInt nColSource = colShell.getNBas(sourceType);
		UInt oldRowPos  = rowShell.getBasisStartIndex(sourceType);
		UInt oldColPos  = colShell.getBasisStartIndex(sourceType);
		S.init(nRowSource,nColSource,oldRowPos,oldColPos);
		S.getData(M0,inTrans);

		//
		// set up local matrix to perserve the result
		//
		UInt targetType = cptrans.getTargetType();
		UInt nRowTarget = rowShell.getNBas(targetType);
		UInt nColTarget = colShell.getNBas(targetType);
		UInt newRowPos  = rowShell.getBasisStartIndex(targetType);
		UInt newColPos  = colShell.getBasisStartIndex(targetType);
		T.init(nRowTarget,nColTarget,newRowPos,newColPos);

		//
		// all cart situation we just do copy
		// else do CP transformation
		//
		if (rowShell.allCart() && colShell.allCart()){
			T.copyMatrix(S);
		}else{
			cptrans.cpTransformOnAtomShellPair(rowShell,colShell,S,tmp,T);
		};
		T.updateData(M1,inTrans);
	}
}

void TBB_CPTransMatrix::doCPTransOnRowOrCol(const blocked_range<UInt>& r) const
{
	// set up data
	UInt maxRowSource = M0.getRow();
	UInt maxColSource = M0.getCol();
	UInt maxRowTarget = M1.getRow();
	UInt maxColTarget = M1.getCol();
	if (rowColStatus == CP_WITH_ROW) {
		maxRowSource = rs.getMaxBasis(TYPE_CART);
		maxRowTarget = rs.getMaxBasis(TYPE_CART);
	}else{
		maxColSource = cs.getMaxBasis(TYPE_CART);
		maxColTarget = cs.getMaxBasis(TYPE_CART);
	}
	BlockMtrx S(maxRowSource,maxColSource);
	BlockMtrx T(maxRowTarget,maxColTarget);

	// loop over atom shell on row/col
	for( UInt n=r.begin(); n!=r.end(); ++n) {

		//
		//get the atom shell data
		//
		const AtomShell& s = rs.getAtomShell(n);

		//
		//get the data block from source
		//
		UInt nRowSource = M0.getRow();
		UInt nColSource = M0.getCol();
		UInt oldRowPos  = 0;
		UInt oldColPos  = 0;
		UInt sourceType = cptrans.getSourceType();
		if (rowColStatus == CP_WITH_ROW) {
			nRowSource = s.getNBas(sourceType);
			oldRowPos  = s.getBasisStartIndex(sourceType);
		}else{
			nColSource = s.getNBas(sourceType);
			oldColPos  = s.getBasisStartIndex(sourceType);
		}
		bool inTrans = false;
		S.init(nRowSource,nColSource,oldRowPos,oldColPos);
		S.getData(M0,inTrans);

		//
		// set up local matrix to perserve the result
		//
		UInt nRowTarget = M1.getRow();
		UInt nColTarget = M1.getCol();
		UInt newRowPos  = 0;
		UInt newColPos  = 0;
		UInt targetType = cptrans.getTargetType();
		if (rowColStatus == CP_WITH_ROW) {
			nRowTarget = s.getNBas(targetType);
			newRowPos  = s.getBasisStartIndex(targetType);
		}else{
			nColTarget = s.getNBas(targetType);
			newColPos  = s.getBasisStartIndex(targetType);
		}
		T.init(nRowTarget,nColTarget,newRowPos,newColPos);

		//
		// do CP transformation
		// if all cart situation we just do copy
		//
		if (s.allCart()) {
			T.copyMatrix(S);
		}else{
			cptrans.cpTransformOnAtomShell(rowColStatus,s,S,T);
		}
		T.updateData(M1,inTrans);
	}
}

///////////////////////////////////////////////
//      @@@@     TBB_NormBasis               //
///////////////////////////////////////////////
void TBB_NormBasisMatrix::operator()(const blocked_range<UInt>& r) const
{
	//
	// store a copy the reference
	//
	const Double* rowVec = &rowNormVec[0];
	const Double* colVec = &colNormVec[0];
	Mtrx& T = target;
	bool doRow = doRowWork;
	bool doCol = doColWork;

	//
	// we note that the loop over r is always done in column index
	//
	for( UInt n=r.begin(); n!=r.end(); ++n ) {

		//
		// first case, do column work
		// only scale this column
		//
		if (doCol) {
         Double N = colVec[n];
			vscal(T.getPtr(0,n),N,T.getRow());
		}

		//
		// second case, do row work
		//
		if (doRow) {
			vmul(T.getPtr(0,n),rowVec,T.getPtr(0,n),T.getRow());
		}
	}
}

///////////////////////////////////////////////
//      @@@@     top interface class         //
///////////////////////////////////////////////
CPTransBasisNorm::CPTransBasisNorm(const GlobalInfor& infor0, const UInt& transWork0, 
		const UInt& scaleWork0, const UInt& mtrxStatus0, 
		const UInt& rowCol):infor(infor0),transWork(transWork0),scaleWork(scaleWork0),
	matrixStatus(mtrxStatus0),rowColStatus(rowCol)
{
	//
	// transformation work
	//
	if (transWork != NO_TRANS && transWork != C2P_WITH_L00 && 
			transWork != C2P_WITH_XYZ && transWork != P2C_WITH_L00 && transWork != P2C_WITH_XYZ) {
		string infor = "invalid input parameters for transWork"; 
		Excep excep("CPTransBasisNorm","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	//
	// scale work
	//
	if (scaleWork != NO_SCALE && scaleWork != DO_SCALE && scaleWork != UNDO_SCALE) { 
		string infor = "invalid input parameters for scaleWork"; 
		Excep excep("CPTransBasisNorm","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	//
	// matrix status
	//
	if (matrixStatus != WITH_MATRIX_ITSELF && matrixStatus != WITH_MATRIX_TRANSPOSE) {
		string infor = "invalid input parameters for matrixStatus"; 
		Excep excep("CPTransBasisNorm","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

	// row col check
	if (rowColStatus != CP_WITH_ROW_COL && rowColStatus != CP_WITH_ROW && rowColStatus != CP_WITH_COL) {
		string infor = "invalid input parameters for rowColStatus"; 
		Excep excep("CPTransBasisNorm","constructor",EXCEPTION_INPUT_PARAMETER_INVALID,infor);
		handleExcep(excep);
	}

}

void CPTransBasisNorm::normBasis(const MolShell& rs, const MolShell& cs, Mtrx& M) const
{
	//
	// do we do the scale work?
	// two situations we just returned:
	// 1 all pure we do nothing inside
	// 2 for only S/P/SP shell, we do nothing 
	//
	if (scaleWork == NO_SCALE) return;
	if (rowColStatus == CP_WITH_ROW_COL) {
		if( rs.allPure() && cs.allPure()) return;
		if( rs.getMaxL()<2 && cs.getMaxL()<2) return;
	}else if (rowColStatus == CP_WITH_ROW) {
		if (rs.allPure() || rs.getMaxL()<2) return;
	}else if (rowColStatus == CP_WITH_COL) {
		if (cs.allPure() || cs.getMaxL()<2) return;
	}

	//
	// set work status
	//
	bool doRow = false;
	bool doCol = false;
	if (rowColStatus == CP_WITH_ROW_COL || rowColStatus == CP_WITH_ROW) doRow = true;
	if (rowColStatus == CP_WITH_ROW_COL || rowColStatus == CP_WITH_COL) doCol = true;

	//
	// build scale vector
	//
	UInt maxL = max(rs.getMaxL(),cs.getMaxL());
	CPTransData scale(maxL,NO_TRANS,scaleWork,WITH_MATRIX_ITSELF);

	//
	// build the row vec and col vec
	//
	DoubleVec rowVec(M.getRow(),ONE);
	DoubleVec colVec(M.getCol(),ONE);
	if (doRow) {

		// determine the matrix type
		UInt matrixType = -1;
		if (M.getRow() == rs.getNCarBas()) {
			matrixType = TYPE_CART;
		}else if (M.getRow() == rs.getNBas()) {
			matrixType = TYPE_NORM;
		}else{
			string infor = "for the given matrix, row does not match shell dimension"; 
			Excep excep("CPTransBasisNorm","normBasis",EXCEPTION_SHELL_DIMENSION_CONFLICT,infor);
			handleExcep(excep);
		}

		// now we build the row vector
		for(UInt iAtom=0; iAtom<rs.getNAtomShells(); iAtom++) {
			const AtomShell& rowShell = rs.getAtomShell(iAtom);
			for(UInt iShell=0; iShell<rowShell.getNShell(); iShell++) {

				// get shell data
				// we will omit these pure basis and S/P/SP shells
				const Shell& is = rowShell.getShell(iShell);
				UInt L = is.getLmax();
				if (L < 2 || is.isPure()) continue;

				// get the vector data
				const DoubleVec& vec = scale.getConvertVec(L);

				// do scaling
				UInt offset = is.getBasisIndex(0,matrixType);
				for(UInt i=0; i<vec.size(); i++) {
					Double N = vec[i];
					rowVec[offset+i] = N;
				}
			}
		}
	}
	if (doCol) {

		// for the situation that row is same with col shell,
		// we just copy
		if (rs == cs && doRow) {
			colVec = rowVec;
		}else{

			// determine the matrix type
			UInt matrixType = -1;
			if (M.getCol() == cs.getNCarBas()) {
				matrixType = TYPE_CART;
			}else if (M.getCol() == cs.getNBas()) {
				matrixType = TYPE_NORM;
			}else{
				string infor = "for the given matrix, col does not match shell dimension"; 
				Excep excep("CPTransBasisNorm","normBasis",EXCEPTION_SHELL_DIMENSION_CONFLICT,infor);
				handleExcep(excep);
			}

			// now we build the row vector
			for(UInt jAtom=0; jAtom<cs.getNAtomShells(); jAtom++) {
				const AtomShell& colShell = cs.getAtomShell(jAtom);
				for(UInt jShell=0; jShell<colShell.getNShell(); jShell++) {

					// get shell data
					// we will omit these pure basis and S/P/SP shells
					const Shell& js = colShell.getShell(jShell);
					UInt L = js.getLmax();
					if (L < 2 || js.isPure()) continue;

					// get the vector data
					const DoubleVec& vec = scale.getConvertVec(L);

					// do scaling
					UInt offset = js.getBasisIndex(0,matrixType);
					for(UInt i=0; i<vec.size(); i++) {
						Double N = vec[i];
						colVec[offset+i] = N;
					}
				}
			}
		}
	}
	
	// now it's time to set up the TBB object 
	TBB_NormBasisMatrix tbb_scale(doRow,doCol,rowVec,colVec,M);

	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	if (infor.useMultiThreads()) {

		init.initialize(infor.getNCPUThreads());
		
		// if with TBB, then we turn off the BLAS lib 
		// which may use the openmp mode
		omp_turnoff();
	}else{
		init.initialize(1);
	}

	// possible timing code
	tick_count t0 = tick_count::now();

	// this is the parallel body
	UInt len = M.getCol();
	parallel_for(blocked_range<UInt>(0,len), tbb_scale);

	// possible timing code
	//tick_count t1 = tick_count::now();
	//Double t = (t1-t0).seconds();
	//printf("%s  %-12.6f\n", "basis set scale with TBB threads, time in seconds ", t);
}

void CPTransBasisNorm::CPTransform(const MolShell& rs, const MolShell& cs, Mtrx& M) const
{
	// do we need to do the work?
	if (transWork == NO_TRANS) return;
	if (rowColStatus == CP_WITH_ROW_COL) {
		if(rs.allCart() && cs.allCart()) return;
		if(rs.getMaxL()<2 && cs.getMaxL()<2) return;
	}else if (rowColStatus == CP_WITH_ROW) {
		if(rs.allCart()) return;
		if(rs.getMaxL()<2) return;
	}else if (rowColStatus == CP_WITH_COL) {
		if(cs.allCart()) return;
		if(cs.getMaxL()<2) return;
	}

	// set up target type
	// the code is copied from  CPTransAtomShell
	// make sure if you modify one place, you need to modify 
	// the code here
	UInt targetType = TYPE_NORM;
	if (transWork == P2C_WITH_L00 || transWork == P2C_WITH_XYZ) {
		targetType = TYPE_CART;
	}

	//
	// if the matrix is actually used as transpose form, then
	// the source and target types must reverse
	//
	if (matrixStatus == WITH_MATRIX_TRANSPOSE) {
		UInt newType = -1;
		if (targetType == TYPE_CART) {
			newType = TYPE_NORM;
		}else{
			newType = TYPE_CART;
		}
		targetType = newType;
	}

	// set up additional matrix to do conversion
	UInt nRows = rs.getNBas(targetType);
	UInt nCols = cs.getNBas(targetType);
	if (rowColStatus == CP_WITH_ROW) {
		nCols = M.getCol();
	}else if (rowColStatus == CP_WITH_COL) {
		nRows = M.getRow();
	}
	Mtrx T(nRows,nCols);
	//cout << "rowcol for T" << nRows << " " << nCols << endl;
	//cout << "rowcol for M" << M.getRow() << " " << M.getCol() << endl;

	// set up the scheduler
	task_scheduler_init init(task_scheduler_init::deferred);
	if (infor.useMultiThreads()) {

		init.initialize(infor.getNCPUThreads());
		
		// if with TBB, then we turn off the BLAS lib 
		// which may use the openmp mode
		omp_turnoff();
	}else{
		init.initialize(1);
	}
	
	// possible timing code
	tick_count t0 = tick_count::now();

	// now the work is processed according to row/col
	if (rowColStatus == CP_WITH_ROW_COL) {
		TBB_CPTransMatrix tbb_trans(rs,cs,transWork,matrixStatus,M,T);
		UInt iAtoms = rs.getNAtomShells();
		UInt jAtoms = cs.getNAtomShells();
		UInt len = iAtoms*jAtoms;
		if (rs == cs) {
			len = iAtoms*(iAtoms+1)/2;
		}
		parallel_for(blocked_range<UInt>(0,len), tbb_trans);

		// if the transform work is only done for the lower triangular part
		// of matrix, here we will do additional copy to make it full
		if (rs == cs) {
			T.copyLowToUpper();
		}
	}else if (rowColStatus == CP_WITH_COL) {
		TBB_CPTransMatrix tbb_trans(cs,transWork,matrixStatus,CP_WITH_COL,M,T);
		UInt len = cs.getNAtomShells();
		parallel_for(blocked_range<UInt>(0,len), tbb_trans);
	}else{
		TBB_CPTransMatrix tbb_trans(rs,transWork,matrixStatus,CP_WITH_ROW,M,T);
		UInt len = rs.getNAtomShells();
		parallel_for(blocked_range<UInt>(0,len), tbb_trans);
	}

	// write the result back
	M = T;

	// possible timing code
	//tick_count t1 = tick_count::now();
	//Double t = (t1-t0).seconds();
	//printf("%s  %-12.6f\n", "C2P/P2C with TBB threads, time in seconds ", t);
}

