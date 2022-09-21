/**
 * CPP files corresponding to the batchbasis.h
 * \author  Fenglai Liu and Jing Kong
 */
#include <iostream>
#include "excep.h"
#include "xcvar.h"
#include "shellprop.h"
#include "batchgrid.h"
#include "shell.h"
#include "cptrans.h"
#include "blas.h"
#include "blas1.h"
#include "matrix.h"
#include "dftderivinfor.h"
#include "sigatombasis.h"
#include "xcintsinfor.h"
#include "denmtrx.h"
#include "dftmatrix.h"
#include "batchbasis.h"
using namespace excep;
using namespace xcvar;
using namespace batchgrid;
using namespace shell;
using namespace cptrans;
using namespace shellprop;
using namespace blas;
using namespace matrix;
using namespace sigatombasis;
using namespace xcintsinfor;
using namespace denmtrx;
using namespace dftmatrix;
using namespace batchbasis;
using namespace std;

void BatchBasis::radDFTBasis(const UInt& np, const Double* coe, const Double* e, 
		const Double* pts, const Double* xyz, DoubleVec& rad) const
{
	for(UInt i=0; i<nGrids; i++) {

		// get the r2: distance between the center and grid points
		Double PCX = pts[i*3  ]  - xyz[0];  // X
		Double PCY = pts[i*3+1]  - xyz[1];  // Y
		Double PCZ = pts[i*3+2]  - xyz[2];  // Z
		Double r2  = PCX*PCX+PCY*PCY+PCZ*PCZ;

		// initilize the derivatives etc.
		Double d0r = ZERO;
		Double d1r = ZERO;
		Double d2r = ZERO;
		Double d3r = ZERO;
		Double d4r = ZERO;

		// radial part
		for(UInt p=0; p<np; p++) {
			Double radial  = coe[p]*exp(MINUS_ONE*e[p]*r2);
			Double ep2     = TWO*e[p]; // minus sign will be added in the derivatives code
			d0r += radial;
			d1r += ep2*radial;
			d2r += ep2*ep2*radial;
			d3r += ep2*ep2*ep2*radial;
			d4r += ep2*ep2*ep2*ep2*radial;
		}
		rad[i] = d0r;
		if (maxDerivOrder >= 1) rad[i+nGrids]   = d1r;
		if (maxDerivOrder >= 2) rad[i+2*nGrids] = d2r;
		if (maxDerivOrder >= 3) rad[i+3*nGrids] = d3r;
		if (maxDerivOrder >= 4) rad[i+4*nGrids] = d4r;
	}
}

void BatchBasis::angDFTBasis(const UInt& lmax, const Double* pts, const Double* c, DoubleVec& ang) const
{
   for(UInt i=0; i<nGrids; i++) {

		// set the x y z
      Double GCX = pts[i*3  ]  - c[0];  // X
      Double GCY = pts[i*3+1]  - c[1];  // Y
      Double GCZ = pts[i*3+2]  - c[2];  // Z

      // this is to evaluate total number of basis sets, L from 0 to lmax
      UInt nTolBas = (lmax+1)*(lmax+2)*(lmax+3)/6; 

		// this is made in libint order
		// you may need to change it if you change the basis set order
      for(UInt L=0; L<= lmax; L++) {
         if(L == 0) {
            ang[0+i*nTolBas] = ONE;
         } else if(L == 1) {
            ang[1+i*nTolBas] = GCX;
            ang[2+i*nTolBas] = GCY;
            ang[3+i*nTolBas] = GCZ;
         } else if(L == 2) {
            ang[4+i*nTolBas] = GCX*GCX;
            ang[5+i*nTolBas] = GCX*GCY;
            ang[6+i*nTolBas] = GCX*GCZ;
            ang[7+i*nTolBas] = GCY*GCY;
            ang[8+i*nTolBas] = GCY*GCZ;
            ang[9+i*nTolBas] = GCZ*GCZ;
         } else if(L == 3) {
            ang[10+i*nTolBas] = GCX*GCX*GCX;
            ang[11+i*nTolBas] = GCX*GCX*GCY;
            ang[12+i*nTolBas] = GCX*GCX*GCZ;
            ang[13+i*nTolBas] = GCX*GCY*GCY;
            ang[14+i*nTolBas] = GCX*GCY*GCZ;
            ang[15+i*nTolBas] = GCX*GCZ*GCZ;
            ang[16+i*nTolBas] = GCY*GCY*GCY;
            ang[17+i*nTolBas] = GCY*GCY*GCZ;
            ang[18+i*nTolBas] = GCY*GCZ*GCZ;
            ang[19+i*nTolBas] = GCZ*GCZ*GCZ;
         } else if(L == 4) {
            ang[20+i*nTolBas] = GCX*GCX*GCX*GCX;
            ang[21+i*nTolBas] = GCX*GCX*GCX*GCY;
            ang[22+i*nTolBas] = GCX*GCX*GCX*GCZ;
            ang[23+i*nTolBas] = GCX*GCX*GCY*GCY;
            ang[24+i*nTolBas] = GCX*GCX*GCY*GCZ;
            ang[25+i*nTolBas] = GCX*GCX*GCZ*GCZ;
            ang[26+i*nTolBas] = GCX*GCY*GCY*GCY;
            ang[27+i*nTolBas] = GCX*GCY*GCY*GCZ;
            ang[28+i*nTolBas] = GCX*GCY*GCZ*GCZ;
            ang[29+i*nTolBas] = GCX*GCZ*GCZ*GCZ;
            ang[30+i*nTolBas] = GCY*GCY*GCY*GCY;
            ang[31+i*nTolBas] = GCY*GCY*GCY*GCZ;
            ang[32+i*nTolBas] = GCY*GCY*GCZ*GCZ;
            ang[33+i*nTolBas] = GCY*GCZ*GCZ*GCZ;
            ang[34+i*nTolBas] = GCZ*GCZ*GCZ*GCZ;
         } else if(L == 5) {
            ang[35+i*nTolBas] = GCX*GCX*GCX*GCX*GCX;
            ang[36+i*nTolBas] = GCX*GCX*GCX*GCX*GCY;
            ang[37+i*nTolBas] = GCX*GCX*GCX*GCX*GCZ;
            ang[38+i*nTolBas] = GCX*GCX*GCX*GCY*GCY;
            ang[39+i*nTolBas] = GCX*GCX*GCX*GCY*GCZ;
            ang[40+i*nTolBas] = GCX*GCX*GCX*GCZ*GCZ;
            ang[41+i*nTolBas] = GCX*GCX*GCY*GCY*GCY;
            ang[42+i*nTolBas] = GCX*GCX*GCY*GCY*GCZ;
            ang[43+i*nTolBas] = GCX*GCX*GCY*GCZ*GCZ;
            ang[44+i*nTolBas] = GCX*GCX*GCZ*GCZ*GCZ;
            ang[45+i*nTolBas] = GCX*GCY*GCY*GCY*GCY;
            ang[46+i*nTolBas] = GCX*GCY*GCY*GCY*GCZ;
            ang[47+i*nTolBas] = GCX*GCY*GCY*GCZ*GCZ;
            ang[48+i*nTolBas] = GCX*GCY*GCZ*GCZ*GCZ;
            ang[49+i*nTolBas] = GCX*GCZ*GCZ*GCZ*GCZ;
            ang[50+i*nTolBas] = GCY*GCY*GCY*GCY*GCY;
            ang[51+i*nTolBas] = GCY*GCY*GCY*GCY*GCZ;
            ang[52+i*nTolBas] = GCY*GCY*GCY*GCZ*GCZ;
            ang[53+i*nTolBas] = GCY*GCY*GCZ*GCZ*GCZ;
            ang[54+i*nTolBas] = GCY*GCZ*GCZ*GCZ*GCZ;
            ang[55+i*nTolBas] = GCZ*GCZ*GCZ*GCZ*GCZ;
         } else if(L == 6) {
            ang[56+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX;
            ang[57+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY;
            ang[58+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCZ;
            ang[59+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY;
            ang[60+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCZ;
            ang[61+i*nTolBas] = GCX*GCX*GCX*GCX*GCZ*GCZ;
            ang[62+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY;
            ang[63+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCZ;
            ang[64+i*nTolBas] = GCX*GCX*GCX*GCY*GCZ*GCZ;
            ang[65+i*nTolBas] = GCX*GCX*GCX*GCZ*GCZ*GCZ;
            ang[66+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY;
            ang[67+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCZ;
            ang[68+i*nTolBas] = GCX*GCX*GCY*GCY*GCZ*GCZ;
            ang[69+i*nTolBas] = GCX*GCX*GCY*GCZ*GCZ*GCZ;
            ang[70+i*nTolBas] = GCX*GCX*GCZ*GCZ*GCZ*GCZ;
            ang[71+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY;
            ang[72+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCZ;
            ang[73+i*nTolBas] = GCX*GCY*GCY*GCY*GCZ*GCZ;
            ang[74+i*nTolBas] = GCX*GCY*GCY*GCZ*GCZ*GCZ;
            ang[75+i*nTolBas] = GCX*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[76+i*nTolBas] = GCX*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[77+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY;
            ang[78+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[79+i*nTolBas] = GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[80+i*nTolBas] = GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[81+i*nTolBas] = GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[82+i*nTolBas] = GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[83+i*nTolBas] = GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
         } else if(L == 7) {
            ang[84+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX;
            ang[85+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY;
            ang[86+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCZ;
            ang[87+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY;
            ang[88+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCZ;
            ang[89+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCZ*GCZ;
            ang[90+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY;
            ang[91+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCZ;
            ang[92+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCZ*GCZ;
            ang[93+i*nTolBas] = GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ;
            ang[94+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY;
            ang[95+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCZ;
            ang[96+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCZ*GCZ;
            ang[97+i*nTolBas] = GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ;
            ang[98+i*nTolBas] = GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ;
            ang[99+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY;
            ang[100+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCZ;
            ang[101+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCZ*GCZ;
            ang[102+i*nTolBas] = GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ;
            ang[103+i*nTolBas] = GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[104+i*nTolBas] = GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[105+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[106+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[107+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[108+i*nTolBas] = GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[109+i*nTolBas] = GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[110+i*nTolBas] = GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[111+i*nTolBas] = GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[112+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[113+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[114+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[115+i*nTolBas] = GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[116+i*nTolBas] = GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[117+i*nTolBas] = GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[118+i*nTolBas] = GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[119+i*nTolBas] = GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
         } else if(L == 8) {
            ang[120+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX;
            ang[121+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY;
            ang[122+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCZ;
            ang[123+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY;
            ang[124+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCZ;
            ang[125+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCZ*GCZ;
            ang[126+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY;
            ang[127+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCZ;
            ang[128+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCZ*GCZ;
            ang[129+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ;
            ang[130+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY;
            ang[131+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCZ;
            ang[132+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCZ*GCZ;
            ang[133+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ;
            ang[134+i*nTolBas] = GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ;
            ang[135+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY;
            ang[136+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCZ;
            ang[137+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCZ*GCZ;
            ang[138+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ;
            ang[139+i*nTolBas] = GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[140+i*nTolBas] = GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[141+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[142+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[143+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[144+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[145+i*nTolBas] = GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[146+i*nTolBas] = GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[147+i*nTolBas] = GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[148+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[149+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[150+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[151+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[152+i*nTolBas] = GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[153+i*nTolBas] = GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[154+i*nTolBas] = GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[155+i*nTolBas] = GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[156+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[157+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[158+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[159+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[160+i*nTolBas] = GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[161+i*nTolBas] = GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[162+i*nTolBas] = GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[163+i*nTolBas] = GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[164+i*nTolBas] = GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
         } else if(L == 9) {
            ang[165+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX;
            ang[166+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY;
            ang[167+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCZ;
            ang[168+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY;
            ang[169+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCZ;
            ang[170+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCZ*GCZ;
            ang[171+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY;
            ang[172+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCZ;
            ang[173+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCZ*GCZ;
            ang[174+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ;
            ang[175+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY;
            ang[176+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCZ;
            ang[177+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCZ*GCZ;
            ang[178+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ;
            ang[179+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ;
            ang[180+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY;
            ang[181+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCZ;
            ang[182+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCZ*GCZ;
            ang[183+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ;
            ang[184+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[185+i*nTolBas] = GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[186+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[187+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[188+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[189+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[190+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[191+i*nTolBas] = GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[192+i*nTolBas] = GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[193+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[194+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[195+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[196+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[197+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[198+i*nTolBas] = GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[199+i*nTolBas] = GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[200+i*nTolBas] = GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[201+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[202+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[203+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[204+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[205+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[206+i*nTolBas] = GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[207+i*nTolBas] = GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[208+i*nTolBas] = GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[209+i*nTolBas] = GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[210+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[211+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[212+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[213+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[214+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[215+i*nTolBas] = GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[216+i*nTolBas] = GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[217+i*nTolBas] = GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[218+i*nTolBas] = GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[219+i*nTolBas] = GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
         } else if(L == 10) {
            ang[220+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX;
            ang[221+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY;
            ang[222+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCZ;
            ang[223+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY;
            ang[224+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCZ;
            ang[225+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCZ*GCZ;
            ang[226+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY;
            ang[227+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCZ;
            ang[228+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCZ*GCZ;
            ang[229+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ;
            ang[230+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY;
            ang[231+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCZ;
            ang[232+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCZ*GCZ;
            ang[233+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ;
            ang[234+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ;
            ang[235+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY;
            ang[236+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCZ;
            ang[237+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCZ*GCZ;
            ang[238+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ;
            ang[239+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[240+i*nTolBas] = GCX*GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[241+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[242+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[243+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[244+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[245+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[246+i*nTolBas] = GCX*GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[247+i*nTolBas] = GCX*GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[248+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[249+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[250+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[251+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[252+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[253+i*nTolBas] = GCX*GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[254+i*nTolBas] = GCX*GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[255+i*nTolBas] = GCX*GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[256+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[257+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[258+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[259+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[260+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[261+i*nTolBas] = GCX*GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[262+i*nTolBas] = GCX*GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[263+i*nTolBas] = GCX*GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[264+i*nTolBas] = GCX*GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[265+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[266+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[267+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[268+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[269+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[270+i*nTolBas] = GCX*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[271+i*nTolBas] = GCX*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[272+i*nTolBas] = GCX*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[273+i*nTolBas] = GCX*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[274+i*nTolBas] = GCX*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[275+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY;
            ang[276+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ;
            ang[277+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ;
            ang[278+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ;
            ang[279+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ;
            ang[280+i*nTolBas] = GCY*GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[281+i*nTolBas] = GCY*GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[282+i*nTolBas] = GCY*GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[283+i*nTolBas] = GCY*GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[284+i*nTolBas] = GCY*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
            ang[285+i*nTolBas] = GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ*GCZ;
         }
      }
   }
}

void BatchBasis::dftBasis0(const UInt& L, const UInt& nTotalBas, const DoubleVec& ang, 
		const DoubleVec& rad, Double* basis) const 
{
	//
	// here the basis set value is arranged in the order of (nBas,nGrids)
	// such order will be reversed later 
	//
	if(L == 0) {
		vcopy(&rad.front(),basis,nGrids);
	}else{
		UInt nCarBas = getCartBas(L,L);
		UInt nCarBasOffset = getCartBas(0,L-1); // this is referring offset of ang array
		for(UInt iG=0; iG<nGrids; iG++) {
			Double radial = rad[iG];
			for(UInt iBas=0; iBas<nCarBas; iBas++) {
				UInt basOffset = iBas + nCarBasOffset;
				basis[iBas+iG*nCarBas] = ang[basOffset+iG*nTotalBas]*radial;
			}
		}
	}
}

void BatchBasis::setupPhi(const BatchGrid& grid, const MolShell& ms, const SigAtomBasis& sigList)
{
	//
	// initilize the radial array - they have fixed length for any
	// atom and basis set, only related to the maxOrder
	//
	DoubleVec rad(nGrids*(maxDerivOrder+1));

	//
	// initilize the ang array
	//
	UInt minL = 0;
	UInt maxL = ms.getMaxL() + maxDerivOrder;
	UInt nTotalBas = getCartBas(minL,maxL);
	DoubleVec ang(nGrids*nTotalBas);

	//
	// set up the tmp phi vector
	// which is in cartesian order
	// the array's order is:
	// (nCarBas,nGrids,nDeriv)
	//
	UInt nDeriv = 0;
	for(UInt i=1; i<=maxDerivOrder; i++) {
		if (i == 1) nDeriv += N_DERIV_1;
		if (i == 2) nDeriv += N_DERIV_2;
		if (i == 3) nDeriv += N_DERIV_3;
		if (i == 4) nDeriv += N_DERIV_4;
	}
	UInt nMaxCarBas = getCartBas(ms.getMaxL(),ms.getMaxL());
	UInt totalDataLen = nGrids*(nDeriv+1)*nMaxCarBas;
	DoubleVec cartPhi(totalDataLen);

	//
	// we need to do shell transformation and scale work
	// either one is needed
	//
	maxL = ms.getMaxL();
	CPTransData c2p(maxL,C2P_WITH_L00,DO_SCALE,WITH_MATRIX_ITSELF);

	// now let's generate the basis set value and its derivatives
	const UIntVec& sigAtomList = sigList.getSigAtoms();
	for(UInt iAtom=0; iAtom<sigList.getNSigAtoms(); iAtom++) {

		// get the atom shell 
		UInt atomIndex = sigAtomList[iAtom];
		const AtomShell& atomShl = ms.getAtomShell(atomIndex);

		// let's generate the angular part information;
		// they are shared by all of basis sets inside 
		// this atom
		UInt maxAtomL = atomShl.getMaxL() + maxDerivOrder;
		UInt nTolBas  = getCartBas(minL,maxAtomL);
		const Double* pts = grid.getGridCoord();
		const Double* xyz = atomShl.getXYZ();
		angDFTBasis(maxAtomL,pts,xyz,ang);

		// now it's the loop for the shell data
		for(UInt iShell=0; iShell<atomShl.getNShell(); iShell++) {
			const Shell& s = atomShl.getShell(iShell);
			UInt basOffset = s.getBasisIndex(0,TYPE_NORM);

			// do it only for significant shell
			if (! sigList.isSigShell(basOffset)) continue;

			// basic information for this shell
			UInt lmin = s.getLmin();
			UInt lmax = s.getLmax();
			UInt np   = s.getNPrim();
			const Double* e = s.getExp();

			// loop over this shell
			for(UInt L=lmin; L<=lmax; L++) {

				// get the dimension of basis sets in this sub-shell
				UInt ncb  = getCartBas(L,L);
				UInt nb   = getPureBas(L,L);

				// get radial information 
				const Double* c = s.getCoe(L);
				radDFTBasis(np,c,e,pts,xyz,rad);

				// zero order basis set values
				Double* phi0 = &cartPhi[0];
				dftBasis0(L,nTolBas,ang,rad,phi0);

				// first derivatives
				if (maxDerivOrder >= 1) {
					Double* phi1 = phi0 + nGrids*ncb;
					dftbasisderiv1(nGrids,L,nTolBas,&ang.front(),&rad.front(),phi1);
				}

				// second derivatives
				if (maxDerivOrder >= 2) {
					Double* phi2 = phi0 + (N_DERIV_1+1)*nGrids*ncb;
					dftbasisderiv2(nGrids,L,nTolBas,&ang.front(),&rad.front(),phi2);
				}

				// third derivatives
				if (maxDerivOrder >= 3) {
					Double* phi3 = phi0 + (N_DERIV_2+N_DERIV_1+1)*nGrids*ncb;
					dftbasisderiv3(nGrids,L,nTolBas,&ang.front(),&rad.front(),phi3);
				}

				// fourth derivatives
				if (maxDerivOrder >= 4) {
					Double* phi4 = phi0 + (N_DERIV_3+N_DERIV_2+N_DERIV_1+1)*nGrids*ncb;
					dftbasisderiv4(nGrids,L,nTolBas,&ang.front(),&rad.front(),phi4);
				}

				// for composite shell quartet case, we need to count in the 
				// basis set offset
				// for example, P shell in the SP shell; it has an offset of 1
				// since S shell is ahead
				UInt offsetInShell = 0;
				if (lmin != lmax && L>lmin) {
					if (s.isPure()) {
						offsetInShell = getPureBas(lmin,L-1);
					}else{
						offsetInShell = getCartBas(lmin,L-1);
					}
				}
				basOffset += offsetInShell;

				// get the significant basis set index for this ishell
				UInt iSigBasisIndex = sigList.getSigBasisIndex(basOffset);

				// scale and transformation work
				if (s.isPure()) {
					const Mtrx& toPure = c2p.getConvertMatrix(L);
					for(UInt iDeriv=0; iDeriv<nDeriv+1; iDeriv++) {
						const Double* carPhi = &cartPhi[iDeriv*nGrids*ncb];
						Mtrx& result = phi[iDeriv];
						mmul(carPhi,toPure.getPtr(),result.getPtr(0,iSigBasisIndex),ncb,nGrids,
								ncb,nb,ncb,toPure.getLd(),result.getLd(),'T','N',ONE,ZERO);
					}
				}else{

					// firstly, let's see whether we need to scale the basis set
					if (L >=2) {
						const DoubleVec& convert = c2p.getConvertVec(L);
						for(UInt iDeriv=0; iDeriv<nDeriv+1; iDeriv++) {
							Double* carPhi = &cartPhi[iDeriv*nGrids*ncb];
							for(UInt iG=0; iG<nGrids; iG++) {
								vmul(&convert.front(),&carPhi[iG*ncb],&carPhi[iG*ncb],ncb);
							}
						}
					}

					// here we need to convert the data into (nGrids,nBas) type
					// that is to say, reverse the cartPhi's data order and write
					// it back to result
					// we note, that here the basis set dimension is always ncb
					// this is true for S/P/SP shell, as well as cartesian D, F etc.
					for(UInt iDeriv=0; iDeriv<nDeriv+1; iDeriv++) {
						const Double* carPhi = &cartPhi[iDeriv*nGrids*ncb];
						Mtrx& result = phi[iDeriv];
						for(UInt iBas=0; iBas<ncb; iBas++) {
							UInt offset = iSigBasisIndex + iBas;
							for(UInt iG=0; iG<nGrids; iG++) {
								result(iG,offset) = carPhi[iBas+iG*ncb];
							}
						}
					}
				}
			}
		}
	}
}

BatchBasis::BatchBasis(const XCVar& xcvar, const BatchGrid& grid, const MolShell& molShl, 
		const SigAtomBasis& sigList):maxDerivOrder(xcvar.getMaxBasisDerivOrder()),
	nGrids(grid.getNGrids()),nSigBasis(sigList.getNSigBasis())
{
	UInt nDeriv = 0;
	for(UInt i=1; i<=maxDerivOrder; i++) {
		if (i == 1) nDeriv += N_DERIV_1;
		if (i == 2) nDeriv += N_DERIV_2;
		if (i == 3) nDeriv += N_DERIV_3;
		if (i == 4) nDeriv += N_DERIV_4;
	}
	phi.reserve(nDeriv+1);
	Mtrx tmp(nGrids,nSigBasis);
	for(UInt i=0; i<nDeriv+1; i++) {
		phi.push_back(tmp);
	}
}

BatchBasis::BatchBasis(const UInt& nSpin, const UInt& maxOrder, const UInt& ng, 
		const UInt& nb):maxDerivOrder(maxOrder),nGrids(ng),nSigBasis(nb)
{
	// just do initilization here
	UInt nDeriv = 0;
	for(UInt i=1; i<=maxDerivOrder; i++) {
		if (i == 1) nDeriv += N_DERIV_1;
		if (i == 2) nDeriv += N_DERIV_2;
		if (i == 3) nDeriv += N_DERIV_3;
		if (i == 4) nDeriv += N_DERIV_4;
	}
	phi.reserve(nDeriv+1);
	Mtrx tmp(nGrids,nSigBasis);
	for(UInt i=0; i<nDeriv+1; i++) {
		phi.push_back(tmp);
	}
}

void BatchBasis::batchBasisSplit(const BatchBasis& bbasis, 
		const SigAtomBasis& sigList, const UInt& iSigAtom)
{
	// count how many orders it has
	// in default we have value to copy
	UInt nOrder = 1;
	for(UInt i=1; i<=maxDerivOrder; i++) {
		if (i == 1) nOrder += N_DERIV_1;
		if (i == 2) nOrder += N_DERIV_2;
		if (i == 3) nOrder += N_DERIV_3;
		if (i == 4) nOrder += N_DERIV_4;
	}

	// copy the given batch basis set value for each order
	UInt len  = nGrids*nSigBasis;
	UInt startIndex  = sigList.getSigBasisBeginIndex(iSigAtom);
	for(UInt order=0; order<nOrder; order++) {
		const Mtrx& p0 = bbasis.getPhi(order);
		Mtrx& p1       = phi[order];
		vcopy(p0.getPtr(0,startIndex),p1.getPtr(),len);
	}
}

void BatchBasis::print() const {

	cout << "*******************************************" << endl;
	cout << "*         BatchBasis Results              *" << endl;
	cout << "*******************************************" << endl;

	// print out the phi values etc.
	const Mtrx& phiValue = getPhi(0);
	phiValue.print("Phi value for this batch");

	if (maxDerivOrder >=1) {
		for(UInt i=0; i<N_DERIV_1; i++) {
			UInt order = XC_DERIV_ORDER_1[i];
			const Mtrx& phiMatrix = getPhi(order);
			cout << "For the section " << order;
			phiMatrix.print(" Phi Derivatives value for this batch");
		}
	}
	cout << endl;

	if (maxDerivOrder >=2) {
		for(UInt i=0; i<N_DERIV_2; i++) {
			UInt order = XC_DERIV_ORDER_2[i];
			const Mtrx& phiMatrix = getPhi(order);
			cout << "For the section " << order;
			phiMatrix.print(" Phi Derivatives value for this batch");
		}
	}
	cout << endl;

	if (maxDerivOrder >=3) {
		for(UInt i=0; i<N_DERIV_3; i++) {
			UInt order = XC_DERIV_ORDER_3[i];
			const Mtrx& phiMatrix = getPhi(order);
			cout << "For the section " << order;
			phiMatrix.print(" Phi Derivatives value for this batch");
		}
	}
	cout << endl;

	if (maxDerivOrder >=4) {
		for(UInt i=0; i<N_DERIV_4; i++) {
			UInt order = XC_DERIV_ORDER_4[i];
			const Mtrx& phiMatrix = getPhi(order);
			cout << "For the section " << order;
			phiMatrix.print(" Phi Derivatives value for this batch");
		}
	}
	cout << endl;
}

