/**
 * cpp file for generating the matrix for pure basis set (angular part with
 * solid hamonics) and Cartesian basis set transformation. Basically, the 
 * C2P and P2C transformation.
 *
 * On the other hand, it also generates the scale vectors for Cartesian
 * Basis sets(normalization from Lx=L,Ly=0,Lz=0 to Lx,Ly,Lz). 
 *
 * the implementation is made according to the reference below:
 * "Transformation between Cartesian and pure spherical harmonic Gaussians"
 * H. Bernhard Schlegel and Michael J. Frisch
 * international Journal of Quantum Chemistry
 * Volume 54, Issue 2, pages 83–87
 *
 * This program is used to generate the matrix element for us to use
 * in the real programming. Since the matrix (C2P and P2C) are constant
 * for all kinds of basis sets, therefore once they are calculated;
 * then they are constants therefore we do not need to repeatly calculated
 * them again and again.
 *
 * We may have deep concerns for the fatorial calculation. We did some
 * test for it in the test_large_l folder. Please see the README for 
 * more details.
 *
 * For the basis sets, we need to be careful about the normalizations.
 * In the paper above, they are actually talking about the "normalized"
 * Cartesian basis set and "normalized" spherical basis sets. However,
 * in real programming we usually normalized the basis sets with 
 * (lx=L,ly=0,lz=0) rather than (lx,ly,lz). This is because all of basis sets
 * information are compressed together in the shell. Therefore, below
 * we have a function to get the C2P transformation matrix (cart to pure)
 * from the basis set normalized with (L,0,0) to the normalized pure basis sets.
 * and for C2P case, this is really the one we need to use in practical
 * programming.
 *
 */
#include "general.h"
#include "angmomlist.h"
#include "functions.h"
#include "boost/lexical_cast.hpp"
using boost::lexical_cast;
using namespace functions;

#define  WITH_L00       1
#define  WITH_LX_LY_LZ  2
#define  WITH_UNNORM_L  3


// this function is c2p from normalized cart basis to normalized 
// pure basis
Double getElementTransFormationMatrixToPure(const Int& L, const Int& M, 
		const Int& lx, const Int& ly, const Int& lz);

// this function is c2p from normalized cart basis with (lx=L, ly=0,lz=0)
// to normalized  pure basis
Double getElementTransFormationMatrixToPureFromL00(const Int& L, const Int& M, 
		const Int& lx, const Int& ly, const Int& lz);

// this function is c2p from un-normalized cart basis to un-normalized  pure basis
// which means, basically this is the orignial equation for an arbitrary basis set 
// function
Double getElementTransFormationMatrixToPureFromUnNormlized(const Int& L, const Int& M, 
		const Int& lx, const Int& ly, const Int& lz);

//
// this function is used to generate the S matrix used in P2C matrix generation
// The S matrix here is for normalized Cartesian basis
//
void getSMatrix(const Int& Lmin,const Int& Lmax, Double* S);

//
// this function is used to generate the S matrix used in P2C matrix generation
// The S matrix here is for normalized Cartesian basis, but the normalization
// factor is for Lx=L,Ly=0,Lz=0
//
void getSMatrixFromL00(const Int& Lmin,const Int& Lmax,Double* S);

//
// this is the top routine to form C2P matrix
// state see above
//
void formMatrixCartToPure(Double* toPure, const Int& Lmin, const Int& Lmax, const Int& state);

//
// this is the top routine to form P2C matrix
//
void formMatrixPureToCart(Double* toCart, const Int& Lmin, const Int& Lmax, const Int& state); 

//
// this is the top routine to produce the normalized factor ratio
// ratio = N(Lx,Ly,Lz)/N(Lx=L,Ly=0,Lz=0)
// so that we can scale the normalization factor from N(Lx=L,Ly=0,Lz=0)
// to N(Lx,Ly,Lz)
//
void scaleNormFactors(Double* N, const Int& Lmin, const Int& Lmax);

//
// print the data in matrix form
//
void print(string title, Int rows, Int column, Int cols, Double* M); 

//
// code printing
//
void printCode(Int rows, Int column, Double* M); 

//
// find the largest data for the given data array
//
Double largestData(Int rows, Int column, Double* M);

///////////////////////////////////////////////////////////////////////////
//  @@@@ functions to generate C2P matrix element etc.
///////////////////////////////////////////////////////////////////////////

Double getElementTransFormationMatrixToPure(const Int& L, const Int& M, 
		const Int& lx, const Int& ly, const Int& lz)
{

	// get the j
	Int j;
	Int absm = abs(M);
	if ((L-absm-lz)%2 == 0 && L-absm-lz >= 0) {
		j = (L-absm-lz)/2;
	}else{
		return ZERO;
	}

	// calculate the prefactor
	// this is for converting the normalized Cartesian basis set to normalized 
	// spherical hamonic functions
	// for large L, we will make sure that it's correct
	Double f1 = getFactorial(2*lx)*getFactorial(2*ly)*getFactorial(2*lz)*getFactorial(L);
	Double f2 = getFactorial(lx)*getFactorial(ly)*getFactorial(lz)*getFactorial(2*L);
	Double f3 = getFactorial(L-absm)/getFactorial(L+absm);
	Double C  = sqrt((f1/f2)*f3);	
	C        *= ONE/(getFactorial(L)*pow(TWO,L));
	//cout << "C " << C << endl;

	// get the result over summation
	Int ilimit = (L-absm)/2;
	Double s = ZERO;
	for(Int i=0; i<=ilimit; i++) {
		if (j>i) continue;
		Double c1 = getFactorial(2*(L-i))*pow(MINUS_ONE,i)/getFactorial(L-absm-2*i);
		c1  *=  bionomialCoe(L,i)*bionomialCoe(i,j);
		for (Int k=0; k<=j; k++) {
			if (lx-2*k >= 0  && lx-2*k <= absm) {
				Double c2 = bionomialCoe(absm,lx-2*k)*bionomialCoe(j,k);
				Double e  = ZERO;
				if (M == 0 && lx%2 == 0) {
					Int p = k-lx/2;
					e = pow(MINUS_ONE,p);
				}
				if (M>0 && abs(absm-lx)%2 == 0) {
					Int p = (2*k+absm-lx)/2;
					e = sqrt(TWO)*pow(MINUS_ONE,p);
				}
				if (M<0 && abs(absm-lx)%2 == 1) {
					Int p = (2*k+absm-lx)/2;
					e = sqrt(TWO)*pow(MINUS_ONE,p);
				}
				s += c1*c2*e;
				//cout << "c1 " << c1 << endl;
				//cout << "c2 " << c2 << endl;
				//cout << "e  " << e << endl;
			}
		}
	}
	Double result = C*s;
	//cout << "lx , ly, lz L M : " << lx << " " << ly << " " << lz << " " << L << " " << M << endl;
	//cout << "value " << result << endl;
	return result;
}

Double getElementTransFormationMatrixToPureFromL00(const Int& L, const Int& M, 
		const Int& lx, const Int& ly, const Int& lz)
{

	// get the j
	Int j;
	Int absm = abs(M);
	if ((L-absm-lz)%2 == 0 && L-absm-lz >= 0) {
		j = (L-absm-lz)/2;
	}else{
		return ZERO;
	}

	// calculate the prefactor
	Double f = getFactorial(L-absm)/getFactorial(L+absm);
	Double C = sqrt(f);	
	C        *= ONE/(getFactorial(L)*pow(TWO,L));
	//cout << "C " << C << endl;

	// get the result over summation
	Int ilimit = (L-absm)/2;
	Double s = ZERO;
	for(Int i=0; i<=ilimit; i++) {
		if (j>i) continue;
		Double c1 = getFactorial(2*(L-i))*pow(MINUS_ONE,i)/getFactorial(L-absm-2*i);
		c1  *=  bionomialCoe(L,i)*bionomialCoe(i,j);
		for (Int k=0; k<=j; k++) {
			if (lx-2*k >= 0  && lx-2*k <= absm) {
				Double c2 = bionomialCoe(absm,lx-2*k)*bionomialCoe(j,k);
				Double e  = ZERO;
				if (M == 0 && lx%2 == 0) {
					Int p = k-lx/2;
					e = pow(MINUS_ONE,p);
				}
				if (M>0 && abs(absm-lx)%2 == 0) {
					Int p = (2*k+absm-lx)/2;
					e = sqrt(TWO)*pow(MINUS_ONE,p);
				}
				if (M<0 && abs(absm-lx)%2 == 1) {
					Int p = (2*k+absm-lx)/2;
					e = sqrt(TWO)*pow(MINUS_ONE,p);
				}
				s += c1*c2*e;
				//cout << "c1 " << c1 << endl;
				//cout << "c2 " << c2 << endl;
				//cout << "e  " << e << endl;
			}
		}
	}
	Double result = C*s;
	//cout << "lx , ly, lz L M : " << lx << " " << ly << " " << lz << " " << L << " " << M << endl;
	//cout << "value " << result << endl;
	return result;
}

Double getElementTransFormationMatrixToPureFromUnNormlized(const Int& L, const Int& M, 
		const Int& lx, const Int& ly, const Int& lz)
{

	// get the j
	Int j;
	Int absm = abs(M);
	if ((L-absm-lz)%2 == 0 && L-absm-lz >= 0) {
		j = (L-absm-lz)/2;
	}else{
		return ZERO;
	}

	// calculate the prefactor
	Double f = getFactorial(L-absm)/getFactorial(L+absm)*(TWO*L+ONE)/(FOUR*PI);
	Double C = sqrt(f);	
	C        *= ONE/(getFactorial(L)*pow(TWO,L));
	//cout << "C " << C << endl;

	// get the result over summation
	Int ilimit = (L-absm)/2;
	Double s = ZERO;
	for(Int i=0; i<=ilimit; i++) {
		if (j>i) continue;
		Double c1 = getFactorial(2*(L-i))*pow(MINUS_ONE,i)/getFactorial(L-absm-2*i);
		c1  *=  bionomialCoe(L,i)*bionomialCoe(i,j);
		for (Int k=0; k<=j; k++) {
			if (lx-2*k >= 0  && lx-2*k <= absm) {
				Double c2 = bionomialCoe(absm,lx-2*k)*bionomialCoe(j,k);
				Double e  = ZERO;
				if (M == 0 && lx%2 == 0) {
					Int p = k-lx/2;
					e = pow(MINUS_ONE,p);
				}
				if (M>0 && abs(absm-lx)%2 == 0) {
					Int p = (2*k+absm-lx)/2;
					e = sqrt(TWO)*pow(MINUS_ONE,p);
				}
				if (M<0 && abs(absm-lx)%2 == 1) {
					Int p = (2*k+absm-lx)/2;
					e = sqrt(TWO)*pow(MINUS_ONE,p);
				}
				s += c1*c2*e;
				//cout << "c1 " << c1 << endl;
				//cout << "c2 " << c2 << endl;
				//cout << "e  " << e << endl;
			}
		}
	}
	Double result = C*s;
	//cout << "lx , ly, lz L M : " << lx << " " << ly << " " << lz << " " << L << " " << M << endl;
	//cout << "value " << result << endl;
	return result;
}

void formMatrixCartToPure(Double* toPure, const Int& Lmin, const Int& Lmax, const Int& state) 
{

	// get the dimension
	//Int nGenBas = getPureBas(Lmin,Lmax);
	Int nCarBas = getCartBas(Lmin,Lmax);

	/*
	// if it's S, P or SP shell, then we just return a unitary matrix
	// if this function is called
	if (Lmax <= 1) {
		for(Int i=0; i<nGenBas; i++) {
			for(Int j=0; j<nCarBas; j++) {
				if (j==i) {
					toPure[j+i*nCarBas] = ONE;
				}else{
					toPure[j+i*nCarBas] = ZERO;
				}
			}
		}
		return;
	}
	*/

	// build the matrix
	// here we note that the whole matrix will be like this(row in Cartesian basis index, 
	// column in pure basis set index):
	// A1
	//   A2
	//     A3
	//       ...
	// each L has its own matrix A, and the block between different L
	// is zero, which is obvious.       
	Int posr    = 0;  //the row starting index for the submatrix with respect to each L
	Int posc    = 0;  //the col starting index for the submatrix with respect to each L
	for(Int L=Lmin; L<=Lmax; L++) {
		Int nbascar = getCartBas(L,L);  // number of column for each submatrix
		Int col = posc;               
		for (Int M=-L; M <= L; M++) {
			for(Int ibas=0; ibas < nbascar; ibas++) {

				// calculate the matrix element
				// stuff the transformation matrix, first is "+" then the "-"
				Int l,m,n;
				getlmn(L,L,ibas,l,m,n); 
				Double x;
				if (state == WITH_L00) {
					x = getElementTransFormationMatrixToPureFromL00(L,M,l,m,n);
				}else if (state == WITH_LX_LY_LZ) {
					x = getElementTransFormationMatrixToPure(L,M,l,m,n);
				}else if (state == WITH_UNNORM_L) {
					x = getElementTransFormationMatrixToPureFromUnNormlized(L,M,l,m,n);
				}else{
					cout << "wrong place we entered " << endl;
					exit(1);
				}
				Int row = ibas+posr;
				toPure[row+col*nCarBas] = x;
			}
			col += 1;
		}
		posr += nbascar;
		posc += 2*L+1;  
	}
}

///////////////////////////////////////////////////////////////////////////
//  @@@@ functions to generate P2C matrix 
///////////////////////////////////////////////////////////////////////////
void getSMatrix(const Int& Lmin, const Int& Lmax, Double* S) 
{
	Int posr    = 0;  //the row starting index for the submatrix with respect to each L
	Int posc    = 0;  //the col starting index for the submatrix with respect to each L
	Int nCarBas = getCartBas(Lmin,Lmax);
	for(Int L=Lmin; L<=Lmax; L++) {
		Int nbascar = getCartBas(L,L);  // number of column for each submatrix
		for(Int jbas=0; jbas < nbascar; jbas++) {

			// angular momentum: corresponding to the column
			Int lj = 0;
			Int mj = 0;
			Int nj = 0;
			getlmn(L,L,jbas,lj,mj,nj); 
			Int col = jbas+posc;               

			// we only build the lower-tri-angular part since S is symmetrical
			for(Int ibas=jbas; ibas < nbascar; ibas++) {

				// angular momentum: corresponding to the row
				Int li = 0;
				Int mi = 0;
				Int ni = 0;
				getlmn(L,L,ibas,li,mi,ni); 
				Int row = ibas+posr;

				// now let's calculate the element
				Double x = ZERO;
				Int l = li+lj;
				Int m = mi+mj;
				Int n = ni+nj;
				if (l%2 == 0 && m%2 == 0 && n%2 == 0) {

					//numerator
					//x  = getFactorial(l)*getFactorial(m)*getFactorial(n);
					//x  = x/(getFactorial(l/2)*getFactorial(m/2)*getFactorial(n/2)); 
					x = getFactorial(l)/getFactorial(l/2);
					x = x*getFactorial(m)/getFactorial(m/2);
					x = x*getFactorial(n)/getFactorial(n/2);


					//denominator
					//Double d1 = getFactorial(2*li)*getFactorial(2*mi)*getFactorial(2*ni);
					//Double d2 = getFactorial(2*lj)*getFactorial(2*mj)*getFactorial(2*nj);
					//Double d3 = getFactorial(li)*getFactorial(mi)*getFactorial(ni);
					//Double d4 = getFactorial(lj)*getFactorial(mj)*getFactorial(nj);
					//Double d  = sqrt(d3*d4/(d1*d2));
					Double di1 = sqrt(getFactorial(li)/getFactorial(2*li));
					Double di2 = sqrt(getFactorial(mi)/getFactorial(2*mi));
					Double di3 = sqrt(getFactorial(ni)/getFactorial(2*ni));
					Double dj1 = sqrt(getFactorial(lj)/getFactorial(2*lj));
					Double dj2 = sqrt(getFactorial(mj)/getFactorial(2*mj));
					Double dj3 = sqrt(getFactorial(nj)/getFactorial(2*nj));
					Double d   = di1*di2*di3*dj1*dj2*dj3;
					x *= d; 
				}

				// stuff the S
				S[row+col*nCarBas] = x;
				S[col+row*nCarBas] = x;
			}
		}
		posr += nbascar;
		posc += nbascar;
	}
}

void getSMatrixFromL00(const Int& Lmin, const Int& Lmax, Double* S)
{
	Int posr    = 0;  //the row starting index for the submatrix with respect to each L
	Int posc    = 0;  //the col starting index for the submatrix with respect to each L
	Int nCarBas = getCartBas(Lmin,Lmax);
	for(Int L=Lmin; L<=Lmax; L++) {
		Int nbascar = getCartBas(L,L);  // number of column for each submatrix
		for(Int jbas=0; jbas < nbascar; jbas++) {

			// angular momentum: corresponding to the column
			Int lj = 0;
			Int mj = 0;
			Int nj = 0;
			getlmn(L,L,jbas,lj,mj,nj); 
			Int col = jbas+posc;               

			// we only build the lower-tri-angular part since S is symmetrical
			for(Int ibas=jbas; ibas < nbascar; ibas++) {

				// angular momentum: corresponding to the row
				Int li = 0;
				Int mi = 0;
				Int ni = 0;
				getlmn(L,L,ibas,li,mi,ni); 
				Int row = ibas+posr;

				// now let's calculate the element
				Double x = ZERO;
				Int l = li+lj;
				Int m = mi+mj;
				Int n = ni+nj;
				if (l%2 == 0 && m%2 == 0 && n%2 == 0) {
					x = doubleFac(l-1)/doubleFac(2*L-1);
					x = x*doubleFac(m-1);
					x = x*doubleFac(n-1);
				}

				// stuff the S
				S[row+col*nCarBas] = x;
				S[col+row*nCarBas] = x;
			}
		}
		posr += nbascar;
		posc += nbascar;
	}
}

void formMatrixPureToCart(Double* toCart, const Int& Lmin, const Int& Lmax, const Int& state) 
{
	// get the dimension
	Int nGenBas = getPureBas(Lmin,Lmax);
	Int nCarBas = getCartBas(Lmin,Lmax);

	/*
	// if it's S, P or SP shell, then we just return a unitary matrix
	// if this function is called
	if (Lmax <= 1) {
		for(Int i=0; i<nGenBas; i++) {
			for(Int j=0; j<nCarBas; j++) {
				if (j==i) {
					toCart[j+i*nCarBas] = ONE;
				}else{
					toCart[j+i*nCarBas] = ZERO;
				}
			}
		}
		return;
	}
	*/

	// firstly, let build the Cart->Pure matrix first
	vector<Double> toPure(nCarBas*nGenBas,ZERO);
	formMatrixCartToPure(&toPure.front(),Lmin,Lmax,state);

	// as we know, that we have toCart = toPure^T*S
	// S is the overlap matrix 
	// there could be two choice, that we form the matrix
	// either from normalized cartesian(Lx=L,Ly=0,Lz=0) to normalized pure
	// or from normalized cartesian (Lx,Ly,Lz) to normalized pure
	// the unormalized situation is invalid
	// since we can not produce the S matrix for unormalized situation
	vector<Double>S(nCarBas*nCarBas,ZERO);
	if (state == WITH_L00) {
		getSMatrixFromL00(Lmin,Lmax,&S.front());
	}else if (state == WITH_LX_LY_LZ) {
		getSMatrix(Lmin,Lmax,&S.front());
	}else{
		cout << "the state in P2C is invalid, can not construct the proper S matrix " << endl;
		exit(1);
	}

	// finally, combine them together to get toCart
	//mmul(&toPure.front(),&S.front(),toCart,nCarBas,nGenBas,nCarBas,nCarBas,'T','N'); 
	for(Int iCol2P=0; iCol2P<nGenBas; iCol2P++) {

		// get the column of toPure matrix
		vector<Double> col2Pure(nCarBas);
		for(Int i=0; i<nCarBas; i++) {
			col2Pure[i] = toPure[iCol2P*nCarBas+i];
		}

		for(Int iColS=0; iColS<nCarBas; iColS++) {

			// get the column of S matrix
			vector<Double> colS(nCarBas);
			for(Int i=0; i<nCarBas; i++) {
				colS[i] = S[iColS*nCarBas+i];
			}

			// now we ultiply them together
			Double val = ZERO;
			for(Int i=0; i<nCarBas; i++) {
				val += col2Pure[i]*colS[i];
			}
			toCart[iCol2P+iColS*nGenBas] = val;
		}
	}
	//print("tocart", nGenBas,nCarBas,6,&toCart[0]);

	// we need to see that whether they could give the identity matrix
	// that is to say, toCart*toPure = I
	vector<Double> I(nGenBas*nGenBas,ZERO);
	for(Int iRow=0; iRow<nGenBas; iRow++) {

		// get the row of to cart matrix
		vector<Double> rowdata(nCarBas);
		for(Int i=0; i<nCarBas; i++) {
			rowdata[i] = toCart[iRow+nGenBas*i];
		}

		// get the column of toPure matrix
		for(Int iCol=0; iCol<nGenBas; iCol++) {
			vector<Double> coldata(nCarBas);
			for(Int i=0; i<nCarBas; i++) {
				coldata[i] = toPure[iCol*nCarBas+i];
			}

			// now we ultiply them together
			Double val = ZERO;
			for(Int i=0; i<nCarBas; i++) {
				val += rowdata[i]*coldata[i];
			}
			I[iRow+iCol*nGenBas] = val;
		}
	}
	/*
	cout << "should be I matrix" << endl;
	for(Int i=0; i<nGenBas; i++) {
		for(Int j=0; j<nGenBas; j++) {
			printf("%-12.7Lf", I[j+i*nGenBas]);
		}
		printf("\n");
	}
	*/

	// test that whether I is idenity matrix?
	bool passit = true;
	for(Int i=0; i<nGenBas; i++) {
		for(Int j=0; j<nGenBas; j++) {
			Double val = I[j+i*nGenBas];
			if (i==j) {
				if (fabs(val-1.0E0)>0.00000001E0) {
					cout << "toCart matrix forming has problem" << endl;
					cout << "position " << j << " " << i << endl;
					passit = false;
				}
			}else{
				if (fabs(val-0.0E0)>0.00000001E0) {
					cout << val << endl;
					cout << "toCart matrix forming has problem" << endl;
					cout << "position " << j << " " << i << endl;
					passit = false;
				}
			}
		}
	}
	if (! passit) {
		exit(1);
	}
}

///////////////////////////////////////////////////////////////////////////
//  @@@@ functions to scale the normalization factors
///////////////////////////////////////////////////////////////////////////
void scaleNormFactors(Double* N, const Int& Lmin, const Int& Lmax) 
{
	Int i = 0;
	for(Int iShell=Lmin;iShell<=Lmax;iShell++) {
		Int nbas= getCartBas(iShell,iShell);
		for(Int iBas=0;iBas<nbas;iBas++) {
			Int l = 0;
			Int m = 0;
			Int n = 0;
			getlmn(iShell,iShell,iBas,l,m,n);
			Double n1 = sqrt(doubleFac(TWO*l-ONE)*doubleFac(TWO*m-ONE)*doubleFac(TWO*n-ONE));
			Double n2 = sqrt(doubleFac(TWO*iShell-ONE));
			N[i] = n2/n1;
			i++;
		}
	}
}

///////////////////////////////////////////////////////////////////////////
//  @@@@ functions to print etc.
///////////////////////////////////////////////////////////////////////////
void print(string title, Int rows, Int column, Int cols, Double* M) 
{
	// check the number of rows
	Int iwidth = 5; // this is related to the order of rows
	string space = "     ";
	if (rows > 99999 && rows <= 9999999) {
		iwidth = 7;
		cols  = 4;
		space = "       ";
	}

	// choose the printing precision in accordance with the column lines
	Int width = -1; 
	Int ndigits = -1; 
	switch (cols) {
		case 3:
			width   = 20; 
			ndigits = 10;
			break;
		case 4:
			width   = 16; 
			ndigits = 10;
			break;
		case 5:
			width   = 14; 
			ndigits = 8;
			break;
			break;
		case 6:
			width   = 12; 
			ndigits = 7;
			break;
		case 7:
			width   = 10; 
			ndigits = 5;
			break;
		default:
			break;
	}
	cout.precision(ndigits);

	// real printing work
	Int start = 0;
	Int colIndex = 1;
	Int end = 0;
	cout << title << endl;
	while(start<column) {

		// how much data we will print each time
		if (start+cols <= column) {
			end = start+cols;
		}else{
			end = column;
		}

		// print out the line number for the data
		cout << space << " ";
		for(Int i=start; i<end; i++) {
			cout << setw(width) << colIndex << " ";
			colIndex++;
		}
		cout << endl;

		// print out the data
		for(Int i=0; i<rows; i++) {
			cout << setw(iwidth) << left << i+1 << " ";
			for(Int j=start; j<end; j++) {
				Double v = M[i+j*rows];
				cout<<setw(width)<<showpoint<<fixed<<right<< v << " ";
			}
			cout << endl;
		}
		start = end;
	}
}

void printCode(Int rows, Int column, Double* M) 
{
	// matrix will be starting from each column to another 
	for(Int j=0; j<column; j++) {
		for(Int i=0; i<rows; i++) {
			Int index = i + j*rows;
			string lhs = "M[" + lexical_cast<string>(index) + "] = ";
			printf("%s %-25.18Lf%s\n", lhs.c_str(),M[index], ";");
		}
	}
}

Double largestData(Int rows, Int column, Double* M) 
{
	// matrix will be starting from each column to another 
	Double data = ZERO;
	for(Int j=0; j<column; j++) {
		for(Int i=0; i<rows; i++) {
			Int index = i + j*rows;
			if (fabs(M[index]) > data) data = fabs(M[index]);
		}
	}
	return data;
}

///////////////////////////////////////////////////////////////////////////
//  @@@@  main function
///////////////////////////////////////////////////////////////////////////
Int main()
{
	Int maxL = 10;

	// this is to transform the basis set function 
	// from L,0,0 to normalized pure
	for(Int L=0; L<=maxL; L++) {

		cout << "// C2P matrix from cartesian basis set normalized with (L,0,0) " << endl;
		cout << "// to normailized pure basis set functions." << endl;
		cout << "// For the angular momentum L: " << L << endl;
		cout << "if (L ==  " << L  << ") {" << endl;

		// to pure matrix
		Int ncb = getCartBas(L,L);
		Int npb = getPureBas(L,L);
		vector<Double> M(ncb*npb,ZERO);
		Int state = WITH_L00;
		formMatrixCartToPure(&M[0], L, L, state);
		//scaleNormFactorsForC2PMtrx(L,L,ncb,npb,&M[0]); 
		//print("topure", ncb,npb,6,&M[0]);
		printCode(ncb,npb,&M[0]);
		cout << "}" << endl;
		cout << endl;
		//Double maxdata = largestData(ncb,npb,&M[0]);
		//cout << "max data for C2P with L00 is " << maxdata << endl;
	}

	// this is to transform the basis set function 
	// from normalized cartesian to normalized pure
	for(Int L=0; L<=maxL; L++) {

		cout << "// C2P matrix from cartesian basis set normalized with (lx,ly,lz) " << endl;
		cout << "// to normailized pure basis set functions." << endl;
		cout << "// For the angular momentum L: " << L << endl;
		cout << "if (L ==  " << L  << ") {" << endl;

		// to pure matrix
		Int ncb = getCartBas(L,L);
		Int npb = getPureBas(L,L);
		vector<Double> M(ncb*npb,ZERO);
		Int state = WITH_LX_LY_LZ;
		formMatrixCartToPure(&M[0], L, L, state);
		//scaleNormFactorsForC2PMtrx(L,L,ncb,npb,&M[0]); 
		//print("topure", ncb,npb,6,&M[0]);
		printCode(ncb,npb,&M[0]);
		cout << "}" << endl;
		cout << endl;
		//Double maxdata = largestData(ncb,npb,&M[0]);
		//cout << "max data for C2P is " << maxdata << endl;
	}

	// this is to transform the basis set function 
	// from arbitrary cartesian to arbitrary pure
	for(Int L=0; L<=maxL; L++) {

		cout << "// C2P matrix from arbitrary cartesian basis set " << endl;
		cout << "// to arbitrary pure basis set functions." << endl;
		cout << "// For the angular momentum L: " << L << endl;
		cout << "if (L ==  " << L  << ") {" << endl;

		// to pure matrix
		Int ncb = getCartBas(L,L);
		Int npb = getPureBas(L,L);
		vector<Double> M(ncb*npb,ZERO);
		Int state = WITH_UNNORM_L;
		formMatrixCartToPure(&M[0], L, L, state);
		//scaleNormFactorsForC2PMtrx(L,L,ncb,npb,&M[0]); 
		//print("topure", ncb,npb,6,&M[0]);
		printCode(ncb,npb,&M[0]);
		cout << "}" << endl;
		cout << endl;
		//Double maxdata = largestData(ncb,npb,&M[0]);
		//cout << "max data for C2P is " << maxdata << endl;
	}

	// this is to transform the basis set function 
	// from normalized pure to normalized cart
	for(Int L=0; L<=maxL; L++) {

		cout << "// P2C matrix from normailized pure basis set functions" << endl;
		cout << "// to cartesian basis set normalized with (lx=L,ly=0,lz=0) " << endl;
		cout << "// For the angular momentum L: " << L << endl;
		cout << "if (L ==  " << L  << ") {" << endl;
		Int ncb = getCartBas(L,L);
		Int npb = getPureBas(L,L);
		vector<Double> M(ncb*npb);
		M.assign(ncb*npb,ZERO);
		Int state = WITH_L00;
		formMatrixPureToCart(&M[0], L, L,state);
		//print("tocart", npb,ncb,6,&M[0]);
		printCode(npb,ncb,&M[0]);
		cout << "}" << endl;
		cout << endl;
		//Double maxdata = largestData(npb,ncb,&M[0]);
		//cout << "max data for P2C is " << maxdata << endl;
	}

	// this is to transform the basis set function 
	// from normalized pure to normalized cart
	for(Int L=0; L<=maxL; L++) {

		cout << "// P2C matrix from normailized pure basis set functions" << endl;
		cout << "// to cartesian basis set normalized with (lx,ly,lz) " << endl;
		cout << "// For the angular momentum L: " << L << endl;
		cout << "if (L ==  " << L  << ") {" << endl;
		Int ncb = getCartBas(L,L);
		Int npb = getPureBas(L,L);
		vector<Double> M(ncb*npb);
		M.assign(ncb*npb,ZERO);
		Int state = WITH_LX_LY_LZ;
		formMatrixPureToCart(&M[0], L, L, state);
		//print("tocart", npb,ncb,6,&M[0]);
		printCode(npb,ncb,&M[0]);
		cout << "}" << endl;
		cout << endl;
		//Double maxdata = largestData(npb,ncb,&M[0]);
		//cout << "max data for P2C is " << maxdata << endl;
	}

	// 
	// this is to scale the basis set function 
	// from normalized cartesian basis set (lx=L,ly=0,lz=0)
	// into the normalized basis set (lx,ly,lz)
	// we do not have it for S and P case
	for(Int L=2; L<=maxL; L++) {

		cout << "// scale vectors for normailized Cartesian basis set functions" << endl;
		cout << "// from (lx=L,ly=0,lz=0) to (lx,ly,lz)" << endl;
		cout << "// For the angular momentum L: " << L << endl;
		cout << "if (L ==  " << L  << ") {" << endl;
		Int ncb = getCartBas(L,L);
		vector<Double> M(ncb);
		M.assign(ncb,ZERO);
		scaleNormFactors(&M[0],L,L); 
		printCode(ncb,1,&M[0]);
		cout << "}" << endl;
		cout << endl;
		//Double maxdata = largestData(ncb,1,&M[0]);
		//cout << "max data for scale data is " << maxdata << endl;
	}

	return 0;
}

