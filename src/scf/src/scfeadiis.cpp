/**
 * cpp file related to scfeadiis.h
 * \author fenglai liu 
 * 
 */
#include<cmath>
#include<cstdio>
#include<iostream>
#include "excep.h"
#include "blas.h"
#include "blas1.h"
#include "denmtrx.h"
#include "spinmatrix.h"
#include "diiscontroller.h"
#include "scfconv.h"
#include "scfmacro.h"
#include "scfeadiis.h"
using namespace excep;
using namespace blas;
using namespace denmtrx;
using namespace spinmatrix;
using namespace scfconv;
using namespace scfmacro;
using namespace diiscontroller;
using namespace scfeadiis;
using namespace std;

SCFEADIIS::SCFEADIIS(const DIISController& cont):jobType(cont.getJob()),maxPrec(cont.getMaxEADIISPrec()),
	nTerms(cont.getNTerms()),params(getLenParams(),ZERO),coefs(nTerms,ZERO)
{ }

void SCFEADIIS::updateTerms(const DIISController& cont, const SCFConv& conv, 
		const DenMtrx& den, const SpinMatrix& fock)
{
	///
	/// now let's see how we can scale the error functional in terms 
	/// of different spin state.
	///
	/// For ADIIS, the original error functional is (equation 7)
	/// 
	/// \f$ f = E(D_{n}) + 2\sum_{i}^{n}c_{i}(D_{i}-D_{n}|F(D_{n})) + 
	/// \sum_{i=1}^{n}\sum_{j=1}^{n} c_{i}c_{j}(D_{i}-D_{n}|F(D_{j})-F(D_{n}))\f$
	///
	/// we note that we actually do not have the zero term E(D_{n}). Because it does not
	/// involves any varying coefficients, so for error functional calculation we can
	/// set it to 0 and does not affect the result.
	///
	/// both of the first order term (involving only c_{i}) and second order term(involving
	/// both c_{i} and c_{j}) are working on the spin dependent density matrix and fock
	/// matrix, therefore for any spin state we just use the original form. so the scaling
	/// factor below is all one for any cases.
	///
	/// For EDIIS, the error function is (equation 8):
	///
	/// \f$ f = \sum_{i}c_{i}E(D_{i}) - \frac{1}{2}
	/// \sum_{i=1}^{n}\sum_{j=1}^{n} c_{i}c_{j}(D_{i}-D_{j}|F(D_{i})-F(D_{j})) \f$
	///
	/// for close shell, the second order term only got the alpha part(since nSpin ==1 so
	/// it does not go with beta part), this brings unbalance between the energy (1st order term)
	/// and the second order term (involing D and F). so for close shell we have scaling factor
	/// as 2.
	///
	///
	UInt ADIIS_scale_t1 = ONE;
	UInt ADIIS_scale_t2 = ONE;
	UInt EDIIS_scale_t1 = ONE;
	UInt EDIIS_scale_t2 = TWO;
	if (conv.getNSpin() == 2) {
		EDIIS_scale_t2 = ONE;
	}

	//
	// the input data matrix (den and Fock) should be symmetrical
	// matrix, and all of history data is symmetrical,too
	//
	bool isSymm = true;
	if (! fock.isSquare()) {
		string infor = "fock matrix must be square in deriving the error functional"; 
		Excep excep("SCFEADIIS","updateTerms",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}
	if (! den.isSquare()) {
		string infor = "density matrix must be square in deriving the error functional"; 
		Excep excep("SCFEADIIS","updateTerms",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}

	// check the dimension
	const Mtrx& F = fock.getMtrx(0);
	const Mtrx& P = den.getMtrx(0);
	if (F.getRow() != P.getRow()) {
		string infor = "dimension conflicts: Fock's dimension should equal to P's dimension"; 
		Excep excep("SCFEADIIS","updateTerms",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}
	UInt n = F.getRow();

	// the first element in the params is the number of electrons
	//double NE = static_cast<double>(conv.getNEles());
	//params[0] = NE;
	
	// get the scf index array from diis controller
	const UIntVec& scfIndexArray = cont.getIndexArray();

	// now it's real work
	if (useADIIS(jobType)) {

		// term 1 and 2
		Mtrx Di(n,n);
		Mtrx Fj(n,n);
		const HistDataMan& oldDen = conv.getHistDenData();
		const HistDataMan& oldFock = conv.getHistFockData();
		for(UInt iSpin=0; iSpin<conv.getNSpin(); iSpin++) {

			// get the data
			const Mtrx& F = fock.getMtrx(iSpin);
			const Mtrx& D = den.getMtrx(iSpin);

			// term 1
			Double* t1 = getT1();
			for(UInt i=0; i<nTerms; i++) {

				// initilize the old density D_{i}
				Di.set(ZERO);

				// if it's the last term, the t1's contribution is zero
				if (i == nTerms-1) continue;

				// get the index for Di
				UInt index = scfIndexArray[i];

				// now retrieve data
				// and do Di = Di - Dn
				oldDen.retrieveData(Di,index,iSpin);
				Di.add(D,MINUS_ONE);

				// now the result term
				// we have a coefficient of 2 for t1 terms
				// see the eq. 7 on ADIIS paper
				t1[i] = TWO*ADIIS_scale_t1*Di.dotProduct(F,isSymm);
			}

			// now it's the second term
			Double* t2 = getT2();
			for(UInt j=0; j<nTerms; j++) {

				// initilize the old Fock F_{i}
				Fj.set(ZERO);

				// get the index for Fj
				UInt indexFj = scfIndexArray[j];

				// now retrieve data
				// and do Fj = Fj - Fn
				oldFock.retrieveData(Fj,indexFj,iSpin);
				Fj.add(F,MINUS_ONE);

				// now do the density part
				for(UInt i=0; i<nTerms; i++) {

					// initilize the old density D_{i}
					Di.set(ZERO);

					// get the index for Di
					UInt indexDi = scfIndexArray[i];

					// now retrieve data
					// and do Di = Di - Dn
					oldDen.retrieveData(Di,indexDi,iSpin);
					Di.add(D,MINUS_ONE);

					// now the result term
					t2[i+j*nTerms] = ADIIS_scale_t2*Di.dotProduct(Fj,isSymm);
				}
			}
		}

	}else if (useEDIIS(jobType)) {

		// term 1
		const DoubleVec& oldEnergy = conv.getHistEnergyData();
		Double* t1 = getT1();
		for(UInt i=0; i<nTerms; i++) {
			UInt index = scfIndexArray[i];
			t1[i] = EDIIS_scale_t1*oldEnergy[index];
		}

		// now it's the second term
		// this is a symmetrical matrix
		// originally for the T2 EDIIS paper has coefficients of 1/2
		// see equation of 8 in EDIIS paper
		// we will add it here
		Mtrx Di(n,n);
		Mtrx Dj(n,n);
		Mtrx Fi(n,n);
		Mtrx Fj(n,n);
		const HistDataMan& oldFock = conv.getHistFockData();
		const HistDataMan& oldDen  = conv.getHistDenData();
		Double* t2 = getT2();
		for(UInt iSpin=0; iSpin<conv.getNSpin(); iSpin++) {
			for(UInt j=0; j<nTerms; j++) {

				// form Dj and Fj
				Dj.set(ZERO);
				Fj.set(ZERO);
				UInt index1 = scfIndexArray[j];
				oldDen.retrieveData(Dj,index1,iSpin);
				oldFock.retrieveData(Fj,index1,iSpin);

				// now loop over i
				for(UInt i=j; i<nTerms; i++) {

					// for diagonal term, it's contribution is zero
					if (i == j) continue;

					// initilize the data
					Di.set(ZERO);
					Fi.set(ZERO);
					UInt index2 = scfIndexArray[i];
					oldDen.retrieveData(Di,index2,iSpin);
					oldFock.retrieveData(Fi,index2,iSpin);
					Fi.add(Fj,MINUS_ONE);
					Di.add(Dj,MINUS_ONE);

					// now the result term
					Double val = MINUS_HALF*Di.dotProduct(Fi,isSymm);
					t2[i+j*nTerms] = EDIIS_scale_t2*val;
					t2[j+i*nTerms] = EDIIS_scale_t2*val;
				}
			}
		}
	}
}

/*
void SCFEADIIS::updateIdemMtrx(const SCFConv& conv, const Mtrx& S)
{
	//
	// Note for the idempotent requirement for density matrix. Here for
	// both EDIIS/ADIIS the new density matrix (as a sum of the density
	// matrices over the coefficients) is required to be idempotent.
	// For achieving that purpose, we append the following term to the 
	// error functional:
	//
	// \f$ f = \sum_{i}\sum_{j}C_{i}C_{j}{-(D_{i}SD_{j}|S) + N} \f$
	//
	// When the new density matrix is:
	//
	// \f$ D = \sum_{i}C_{i} D_{i}\f$
	//
	// It's easy to see this is just DSD = D combined with tr(DS) = N
	// If the error function is minimized, we expected that the f 
	// above should be approaching to zero.
	//
	UInt scale = ONE;
	if (conv.getNSpin() == 1) scale = TWO;

	//
	// the input data matrix (den and Fock) should be symmetrical
	// matrix, and we assume that all of history data is symmetrical,too
	// we do not check it at this point
	//
	bool isSymm = true;
	UInt n = S.getRow();

	// finally, let's add in idempotent requirement
	Mtrx Di(n,n);
	Mtrx Dj(n,n);
	Mtrx tmp(n,n);
	const HistDataMan& oldDen  = conv.getHistDenData();
	double* idem = getIdemMatrix();
	for(UInt iSpin=0; iSpin<conv.getNSpin(); iSpin++) {
		for(UInt j=beginIndex; j<=endIndex; j++) {

			// form Dj 
			Dj.set(ZERO);
			oldDen.retrieveData(Dj,j,iSpin);

			// now loop over i
			for(UInt i=j; i<=endIndex; i++) {

				// initilize the data
				Di.set(ZERO);
				oldDen.retrieveData(Di,i,iSpin);

				// tmp = Di*S
				tmp.symMatrixMult(Di,S,'L',ONE,ZERO);
				//Di.print("Di");
				//S.print("S");
				//tmp.print("Di*S");

				// Di = (Di*S)*Dj
				Di.symMatrixMult(Dj,tmp,'R',ONE,ZERO);
				//Dj.print("Dj");
				//Di.print("Di*S*Dj");

				// now the result term, which is (DiSDj|S)
				Double val = Di.dotProduct(S,isSymm);
				val = scale*val;
				idem[i+j*nTerms] += ZERO;

				// making (j,i) term
				if (i != j) {

					// doing transpose
					// DjSDi = (DiSDj)^T
					bool inPlace = true;
					Di.transpose(inPlace);

					// now the result term, which is (DjSDi|S)
					Double val = Di.dotProduct(S,isSymm);
					val = scale*val;
					idem[j+i*nTerms] += ZERO;
				}
			}
		}
	}
}
*/

void SCFEADIIS::minimization(const DIISController& cont, const SCFConv& conv, 
		const DenMtrx& den, const SpinMatrix& Fock)
{
	// if nTerms is 1, then we do not need to do anything
	if (nTerms == 1) {
		coefs[0] = ONE;
		return;
	}

	// firstly, update the normal terms
	updateTerms(cont,conv,den,Fock);

	// now update idempotent matrix
	/*
	bool update = true;
	if (update) {

		// get the overlap matrix
		// make sure that overlap matrix should be a full matrix
		UInt nRow = -1;
		UInt nCol = -1;
		oneEMtrx.getDim(TWO_BODY_OVERLAP,nRow,nCol);
		Mtrx S(nRow,nCol);
		bool withFull = true;
		oneEMtrx.getM(TWO_BODY_OVERLAP,S,withFull);

		// do the computation
		updateIdemMtrx(conv,S);
	}
	*/

	// derive the real coefficients
	// let's see the space to determine which way we will go
	// for the exact solution, for each precision digit the
	// loop number will take (2*9+1)^(nTerms-1). So for 
	// nTerms = 7 the cost is 47045881 ~ 5.0E+8
	//
	// for the approximate solution, the cost is 2*10^(nTerms-1)
	// so limit is 10 is reasonable.
	//
   if (nTerms<=SCF_EADIIS_SPACE_LIMIT) {
		findCoefs(cont.withMorePreciseCoef());
	}else {
		string infor = "the number of SCF vector is larger than the limit";
		Excep excep("SCFEADIIS","minimization",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}	

	// finally check the result
	checkCoefs(); 
}

void SCFEADIIS::findCoefs(bool withExpensiveSolver)
{
	//
	// this function try to get the result set of coefficients
	// for a given error functional:
	//
	// \f$ f = f_{0} + \sum_{i} c_{i}T1[i] + \sum_{i}\sum_{j} c_{i}c_{j}T2[i,j] \f$
	//
	// T1 and T2 are the components for the error functional f, by varing
	// the coefficients of c in constraints:
	//
	// \f$ c_{i} >=0 and \sum_{i} c_{i} = 1\f$
	//
	// we are looking for the coefficients set which minimizes the f.
	//
	// this searching algorithm is divided into two round:
	//
	// the first round it searchs the special cases, that for each i
	// we set c_{i} to be 1 and the rest of them are 0. This will
	// generate the first trial result for coefficients.
	//
	// The second round is divided into two steps. The first step is to
	// do an exact optimal coefficients search. If the coefficient space
	// is small, then we will generate the answer as accurate as 0.001.
	// else we will generate the answer as accurate as 0.01.
	//
	// the second step will refine the answer for the first step until
	// the given precision is met. It will update the c_{i} for both left
	// and right intervals until the global minimum is found.
	//
	// the search algorithm itself is exhaustive search. Because sum of coefficients
	// is 1, so the varying space number is nTerms-1. We use loop_identifier to 
	// store all of possible combinations of the coefficients, you can commit out
	// the debugging code to see how to search progress.
	//

	// set the number of coefficients 
	UInt nCoefSpace = nTerms;
	//cout << "number of coefs space: " << nCoefSpace << endl;

	// set the number of precisions
	UInt nPreNum = maxPrec;
	//cout << "number of coefs precision: " << nPreNum << endl;
	
	// set up an constant for comparing two numbers
	// because we only support limitted precision case,
	// for example; 1.0E-5 probably is the most accurate one
	// therefore in this precision two coefficients will have
	// difference > 1.0E-5
	// so we can set a constant of 1.0E-10 so that to compare
	// that whether two coefficients are same
	const Double SMALLEST_NUM = 1.0E-10;
	if (nPreNum>=5) {
		string infor = "the result coefficients rerquirs too high precision, please make it less than 1.0E-5";
		Excep excep("SCFEADIIS","findCoefs",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}

	// set up a very big number used to initialize the error functional
	const Double BIGGEST_NUM = 1.0E+10;

	// set the terms of error functional
	const Double* T1 = getT1();
	const Double* T2 = getT2();
	Double T0 = ZERO;

	// set up the optimum coefficents list
	// that's the result for the 1st round of 
	// search
	DoubleVec optCoefList1(nCoefSpace);

	// now let's get the initialization of 
	// the coef list
	// because \sum c_{i} = 1, and c_{i} >=0
	// therefore the first trial search,
	// is that for each i the corresponding c_{i} is 1,
	// and others are 0
	Double f1 = BIGGEST_NUM;
	UInt choice = 0;
	for(UInt i=0; i<nCoefSpace; i++) {

		// now compute the error functional
		Double f = T0;

		// first order term
		f += T1[i];

		// second order term
		f += T2[i+i*nCoefSpace];

		// now let's get the smallest 
		if (f<f1) {
			choice = i;
			f1 = f;
		}
	}
	optCoefList1.assign(nCoefSpace,ZERO);
	optCoefList1[choice] = ONE;
	//printf("choice is %d and f1 is %f\n",(Int)choice,f1);

	//
	// now let's search the optimum solution
	// the loop_identifier
	// because \sum_{i}c_{i} = 1, the last coefficients
	// is depending on the rest of the others
	// so the length of loop_identifier is nCoefSpace-1
	//
	// here different from the one used in findCoefs function,
	// the value in loop_identifier starts from -9 to 9; therefore
	// it's IntVec, and we do not have sign array
	//
	UInt nTerm = nCoefSpace-1;
	IntVec loop_identifier(nTerm, 0);

	// 
	// this one stores the second round of coefficients exploring
	// optCoefList2 is the result coefficients, and tempCoefList
	// stores the temp best coefficients
	//
	DoubleVec optCoefList2(nCoefSpace);
	DoubleVec tempCoefList(nCoefSpace);
	optCoefList2.assign(nCoefSpace,ZERO);

	//
	// set up the trial coefficients list
	//
	DoubleVec trailCoefList(nCoefSpace);

	//
	// now let's set the initial value for the search
	// if the total space number is less than the given limit,
	// then we can carry on the search exact to 0.001
	//
	// the limitThresh is 99 for np = 2, 999 for np = 3 (0.01*99 = 0.99,
	// this is the limit for np=2; so as np=3)
	//
	// testingSum is related to the requirement that sum_{i}C_{i} = 1.
	// if sum of varying coefficients is larger than the value, then
	// with base multiplied (np=2 is 0.01, 3 is 0.001) you will see
	// that the sum of coefficients will be larger than 1
	//
	UInt initialNP  = 2;
	Int testingSum  = 100;
	Int limitThresh = 99;
	if (withExpensiveSolver) {
		initialNP   = 3;
		testingSum  = 1000;
		limitThresh = 999;
	}

	// now we start the opt position search
	// f2 stores the minimum error functional value
	Double f2 = BIGGEST_NUM;
	for(UInt np=initialNP; np<=nPreNum; np++) {

		// set the base of inc value
		// from 0.01 to 0.001 etc.
		Double base = 0.01E0;
		if (np == 3) base = 0.001E0;
		if (np > 3) {
			base = pow(0.1E0,np);
		}

		//
		// for the first round of search, the value in
		// loop_identifier starts from 1 to 99 for np = 2,
		// or from 1 to 999 for np = 3. 
		//
		// for the rest of search, the values starts from -9 to 9
		// to sample the values of left interval and right interval
		//
		if (np==initialNP) {
			loop_identifier.assign(nTerm,0); 
		}else{
			loop_identifier.assign(nTerm,-9); 
		}

		// now do the search
		UInt curPos = 0;     // current working position in loop_identifier
		while(true) {

			// test that whether we omit this calculation
			// we test the first cycle, since this is 
			// the most expensive one
			bool goodToGo = true;
			if (np == initialNP) {
				Int sum = 0;
				for(UInt i=0; i<nTerm; i++) sum += loop_identifier[i];
				if (sum>testingSum) goodToGo = false;
			}

			// may I proceed?
			if (goodToGo) {

				//
				// first step, deriving the trail coefficients
				//
				Double sum = ZERO;
				bool positiveCoef = true;
				for(UInt i=0; i<nTerm; i++) {
					Double c = optCoefList2[i]+base*loop_identifier[i];
					if (c<ZERO) {
						positiveCoef = false;
						break;
					}
					trailCoefList[i] = c;
					sum += c;
				}

				//
				// whether the coefficients are valid?
				// for negative coefficient, just go to next round
				// we only perform the error functional calculation for 
				// the case that 1-sum >= 0
				//
				if (positiveCoef) {
					if (sum<ONE || fabs(ONE-sum)<SMALLEST_NUM) {
						trailCoefList[nTerm] = ONE-sum;

						// now compute the error functional
						Double f = T0;

						// first order term
						Double fx1 = ZERO;
						for(UInt i=0; i<nCoefSpace; i++) {
							fx1 += trailCoefList[i]*T1[i];
						}
						f += fx1;

						// second order term
						Double fx2 = ZERO;
						for(UInt i=0; i<nCoefSpace; i++) {
							Double ci = trailCoefList[i];
							for(UInt j=0; j<nCoefSpace; j++) {
								Double cj = trailCoefList[j];
								fx2 += ci*cj*T2[j+i*nCoefSpace];
							}
						}
						f += fx2;

						// now let's get the smallest 
						if (f<f2) {
							f2 = f;
							tempCoefList = trailCoefList;
						}

						//
						// debug
						// print out the loop_identifier
						//
						/*
							for(int i=0; i<nTerm; i++) {
							printf("%-2d ", loop_identifier[i]);
							}
							printf("\n");
							for(int i=0; i<nCoefSpace; i++) {
							printf("%-4f ", trailCoefList[i]);
							}
							printf("\n");
							printf("fx1 is %-16.8f, fx2 is %-16.8f, f is %-16.8f, f2 is %-16.8f\n", fx1, fx2, f, f2);
							*/
					}
				}
			}

			// increase position
			// each position should be between -9 to 9 (or 0 to 99/999 for first round)
			Int limit = 9;
			if (np == initialNP) limit = limitThresh;
			if (loop_identifier[curPos] < limit) {
				loop_identifier[curPos] += 1;
			}else{

				// let's go to see whether we have loop over all possible arragements
				if (loop_identifier[nTerm-1] == limit) {
					bool finish = true;

					// here we need to be careful
					// the nTerms is possibilily less than 2
					// for example, the second SCF cycle etc.
					if (nTerm>=2) {
						for(Int i=(Int)(nTerm-2); i>=0; i--) {
							if (loop_identifier[i] != limit) {
								finish = false;
								break;
							}
						}
					}
					if (finish) break;
				}

				// carry the position in the loop to the next possible one
				UInt iPos = curPos+1;
				while(iPos<nTerm) {
					if (loop_identifier[iPos] < limit) {
						loop_identifier[iPos] += 1;
						break;
					}
					iPos++;
				}

				// clear all of the position in front of iPos
				// so that to start a new round of searching
				Int val = -9;
				if (np==initialNP) val = 0;
				for(UInt i=0; i<iPos; i++) {
					loop_identifier[i] = val;
				}

				// set current position to the original
				curPos = 0;
			}
		}

		// update the opt coefficients
		// if temp is different from the opt one
		bool same = true;
		for(UInt i=0; i<nCoefSpace; i++) {
			if(fabs(optCoefList2[i]-tempCoefList[i])>SMALLEST_NUM) {
				same = false;
				break;
			}
		}
		if (!same) {
			optCoefList2 = tempCoefList;
		}

		// debug
		// print out final result
		/*
		printf("for np %d f2 is %f\n", (Int)np, f2);
		for(UInt i=0; i<nCoefSpace; i++) {
			printf("%-10.6f  ", optCoefList2[i]);
		}
		printf("\n");
		*/
	}

	//
	// finally, let's compare the two rounds of search
	// optCoefList2 and optCoefList1
	//
	if (f2<f1) {
		coefs = optCoefList2;
	}else{
		coefs = optCoefList1;
	}
}

void SCFEADIIS::checkCoefs() const 
{
	// firstly, see whether each term >= 0
	bool pass = true;
	for(UInt i=0; i<nTerms; i++) {
		if (coefs[i]<-THRESHOLD_MATH) {
			pass = false;
			break;
		}
	}
	if (!pass) {
		for(UInt i=0; i<nTerms; i++) {
			printf("%-15.8f  ", coefs[i]);
		}
		printf("\n");
		string infor = "the result coefficients have elements less than zero, which is invalid";
		Excep excep("SCFEADIIS","checkCoefs",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}

	// now let's check the sum
	Double sum = ZERO;
	for(UInt i=0; i<nTerms; i++) {
		sum += coefs[i];
	}
	if (fabs(sum-ONE)>0.00001E0) {
		string infor = "the sum of result coefficients is not 1, which is invalid";
		Excep excep("SCFEADIIS","checkCoefs",EXCEPTION_SCFEADIIS_ERROR,infor);
		handleExcep(excep);
	}
}

/*
bool SCFEADIIS::checkIdempotent() const 
{
	UInt offset = 1+nTerms+nTerms*nTerms;
	if (useADIIS()) offset += 1;
	const double* idemMtrx = &params[offset];
	double f = ZERO;
	for(UInt j=0; j<nTerms; j++) {
		for(UInt i=0; i<nTerms; i++) {
			f += coefs[i]*coefs[j]*idemMtrx[i+j*nTerms];
		}
	}
	Double NE = params[0];
	if (fabs(f-NE)>tol) {
		printf("idempotent condition value is %-12.8f, and NE is %-12.8f\n", f, NE);
		return false;
		//string infor = "the idempotent requirement is not satisfied";
		//Excep excep("SCFEADIIS","checkCoefs",EXCEPTION_SCFEADIIS_ERROR,infor);
		//handleExcep(excep);
	}
	return true;
}
*/

void SCFEADIIS::print() const 
{
	cout << "***************************************" << endl;
	cout << "debug information for SCF EDIIS/ADIIS  " << endl;
	cout << "***************************************" << endl;

	// general information
	if (useADIIS(jobType)) {
		cout << "the method we currently use is: ADIIS" << endl;
	}else if (useEDIIS(jobType)) {
		cout << "the method we currently use is: EDIIS" << endl;
	}
	cout << "totally " << nTerms << " terms involved " << endl;
	cout << "the result mo is with precision of " << pow(0.1E0,maxPrec) << endl;
	cout << endl;

	// parameter array
	cout << "error functional components:" << endl;
	//cout << "number of total electrons: " << params[0] << endl;
	cout << "first order term: " << endl;
	UInt offset = 0;
	const Double* t1 = &params[offset];
	for(UInt i=0; i<nTerms; i++) {
		cout << i << " " << t1[i] << endl;
	}
	offset += nTerms;
	Mtrx t2(nTerms,nTerms);
	vcopy(&params[offset],t2.getPtr(),nTerms*nTerms);
	t2.print("second order term is a square matrix:");
	cout << endl;

	// idempotent matrix 
	/*
	t2.set(ZERO);
	offset += nTerms*nTerms;
	vcopy(&params[offset],t2.getPtr(),nTerms*nTerms);
	t2.print("idempotent matrix");
	cout << endl;
	*/

	// result
	cout << "result coefficients is: " << endl;
	for(UInt i=0; i<nTerms; i++) {
		cout << i << " " << coefs[i] << endl;
	}
}
