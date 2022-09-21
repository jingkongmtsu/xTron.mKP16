/**
 * cpp file related to oneemtrx.h
 */
#include <boost/lexical_cast.hpp>
#include "excep.h"
#include "shell.h"
#include "molecule.h"
#include "globalinfor.h"
#include "gintsinfor.h"
#include "gints2d.h"
#include "oneemtrx.h"
using namespace excep;
using namespace shell;
using namespace molecule;
using namespace globalinfor;
using namespace gintsinfor;
using namespace gints2d;
using namespace oneemtrx;

OneEMtrx::OneEMtrx(const GlobalInfor& infor, const UInt& section, 
		const UIntVec& matrixList, bool useFile):
	HistDataMan(infor,"oneemtrx",section,1,useFile),jobList(matrixList.size(),0)
{
	// we need to re-modify the joblist
	// so that the two body overlap and orthogonal matrix
	// is at the back
	UInt pos = 0;
	for(UInt iJob=0; iJob<matrixList.size(); iJob++) {
		UInt job = matrixList[iJob];
		if (job == TWO_BODY_OVERLAP || job == ORTHOGONAL_MATRIX) continue;
		jobList[pos] = job;
		pos++;
	}

	// now let's see whether we have ov and ortho
	bool hasOV = false;
	bool hasOrtho = false;
	for(UInt iJob=0; iJob<matrixList.size(); iJob++) {
		UInt job = matrixList[iJob];
		if (job == TWO_BODY_OVERLAP) hasOV = true;
		if (job == ORTHOGONAL_MATRIX) hasOrtho = true;
	}

	// now append them
	UInt n = matrixList.size();
	if (hasOV) {
		if (hasOrtho) {
			jobList[n-2] = TWO_BODY_OVERLAP;
			jobList[n-1] = ORTHOGONAL_MATRIX;
		}else{
			jobList[n-1] = TWO_BODY_OVERLAP;
		}
	}else{
		if (hasOrtho) {
			jobList[n-1] = ORTHOGONAL_MATRIX;
		}
	}
}	

OneEMtrx::OneEMtrx(const GlobalInfor& infor, const UInt& section, 
		bool useFile):HistDataMan(infor,"oneemtrx",section,1,useFile),jobList(4,0)
{
	// for regular SCF part, four matrices needed:
	// kineic
	// NAI
	// two body overlap
	// orthogonal matrix
	jobList[0] = KINETIC;
	jobList[1] = NUCLEAR_ATTRACTION;
	jobList[2] = TWO_BODY_OVERLAP;
	jobList[3] = ORTHOGONAL_MATRIX;
}

void OneEMtrx::formMtrx(const GIntsInfor& intInfor, const Molecule& mol, 
		const MolShell& rs, const MolShell& cs, bool printTiming, bool printMatrix)
{
	// get the dimension data first
	UInt nRow = rs.getNBas();
	UInt nCol = cs.getNBas();

	// now let's do the matrix forming
	bool hasOV = false;
	bool hasOrtho = false;
	UInt order = 0;
	Mtrx M(nRow,nCol);
	bool isSymm = false;
	if (rs == cs) isSymm = true;
	for(UInt iJob=0; iJob<jobList.size(); iJob++) {

		// get the job 
		// we will do overlap and orthogonal matrix later
		UInt job = jobList[iJob];
		if (job == TWO_BODY_OVERLAP) {
			hasOV = true;
			continue;
		}
		if (job == ORTHOGONAL_MATRIX) {
			hasOrtho = true;
			continue;
		}

		// reset the matrix
		if (iJob>0) {
			M.set(ZERO);
		}

		// form matrix and store it
		GInts2D gints(rs,cs,intInfor,job,order);
		gints.doMtrx(rs,cs,mol,M,printTiming);
		storeData(M,isSymm);
		if (printMatrix) {
			string title = getJobName(job);
			M.copyLowToUpper();
			M.print(title);
		}
	}

	// overlap and orthogonal
	if (hasOV || hasOrtho) {

		// overlap matrix
		M.set(ZERO);
		GInts2D gints(rs,cs,intInfor,TWO_BODY_OVERLAP,order);
		gints.doMtrx(rs,cs,mol,M,printTiming);
		if (hasOV) {
			storeData(M,isSymm);
			if (printMatrix) {
				M.copyLowToUpper();
				M.print("TWO BODY OVERLAP");
			}
		}

		// finally, form the orthogonal matrix
		// the orthogonal matrix will be stored in full
		// as a matrix, this is required by forming the 
		// mo (see mo.h and mo.cpp)
		// Therefore, we pretend that this matrix is 
		// assymetrical
		if (hasOrtho) {
			bool symmMatrix = false;
		   Mtrx ortho(nRow,nCol);
			ortho.formOrthoMatrix(M,intInfor.getLinearDepThresh());
			if (printMatrix) {
				M.print("ORTHOGONAL MATRIX(S^(-1/2))");
			}
			storeData(ortho,symmMatrix);
		}
	}
}	

void OneEMtrx::getM(const UInt& name, Mtrx& M, bool intoFullMatrix) const 
{
	// find the matrix position
	UInt pos = -1;
	bool findIt = false;
	for(UInt i=0; i<jobList.size(); i++) {
		if (jobList[i]	== name) {
			pos = i;
			findIt = true;
			break;
		}
	}

	// check whether we find it
	if (! findIt) {
		string infor = "the integral matrix name is: " + boost::lexical_cast<string>(name);
		Excep excep("OneEMtrx","getM",EXCEPTION_GINTS_FAIL_TO_FIND_MATRIX,infor);
		handleExcep(excep);
	}

	// check the dimension
	if (M.getRow() != rowDim(pos) || M.getCol() != colDim(pos)) {
		string infor = "the integral matrix name is: " + boost::lexical_cast<string>(name);
		Excep excep("OneEMtrx","getM",EXCEPTION_MATRIX_DIMENSION_CONFLICT,infor);
		handleExcep(excep);
	}

	// now we are safe to get the matrix
	retrieveData(M,pos);

	// do we also copy data to upper part?
	// if the matrix is symmetrical and also
	// you ask to make it into full matrix
	// we will do it
	if(intoFullMatrix) {
		bool isSymm = false;
		if (isSymmData(pos)) isSymm = true;
		if(isSymm) {
			M.copyLowToUpper();
		}
	}
}

void OneEMtrx::getDim(const UInt& name, UInt& nRow, UInt& nCol) const
{
	// find the matrix position
	UInt pos = -1;
	bool findIt = false;
	for(UInt i=0; i<jobList.size(); i++) {
		if (jobList[i]	== name) {
			pos = i;
			findIt = true;
			break;
		}
	}

	// check whether we find it
	if (! findIt) {
		string infor = "the integral matrix name is: " + boost::lexical_cast<string>(name);
		Excep excep("OneEMtrx","getDim",EXCEPTION_GINTS_FAIL_TO_FIND_MATRIX,infor);
		handleExcep(excep);
	}

	// row and col dimensions
	nRow = rowDim(pos);
	nCol = colDim(pos);
}
