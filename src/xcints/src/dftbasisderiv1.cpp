/**
 * This function is used to generate 1 derivatives for basis set 
 * The basis set derivatives are evaluated for the given shell which 
 * is characterized by the L(no composite shell!). Generally, by given the 
 * derivative order (for exmaple, X, Y Z or XX, YY or XYY etc.)
 * for an arbitrary shell we could combine the radial part and 
 * the angular part together so to form the result.
 * The result is arranged as: (nBas, ng, nDerivOrder)
 * nBas is the number of Cartesian type basis set for shell with L
 * \param ng         number of grid points 
 * \param L          angular momentum of the shell 
 * \param nTolCarBas number of Cartesian basis set in the ang array 
 * \param ang        angular part of the basis set values(nTolCarBas,ng) 
 * \param rad        radial part of the basis set values 
 * \return basis     derivatives of basis set values for the given order
 * \author Fenglai Liu and Jing Kong 
 */
#include"libgen.h"
#include"batchbasis.h"
using namespace batchbasis;

void BatchBasis::dftbasisderiv1(const UInt& ng, const UInt& L, const UInt& nTolCarBas, const Double* ang, const Double* rad, Double* basis) const 
{

   // now we set up the nBas for the computation
   UInt nBas = (L+1)*(L+2)/2;

   // now we do derivatives for the given basis set to X
   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[1];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*angArray[0]-rad[ip+ng]*angArray[4];
         bas[1] = -rad[ip+ng]*angArray[5];
         bas[2] = -rad[ip+ng]*angArray[6];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*2*angArray[1]-rad[ip+ng]*angArray[10];
         bas[1] = rad[ip]*angArray[2]-rad[ip+ng]*angArray[11];
         bas[2] = rad[ip]*angArray[3]-rad[ip+ng]*angArray[12];
         bas[3] = -rad[ip+ng]*angArray[13];
         bas[4] = -rad[ip+ng]*angArray[14];
         bas[5] = -rad[ip+ng]*angArray[15];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*3*angArray[4]-rad[ip+ng]*angArray[20];
         bas[1] = rad[ip]*2*angArray[5]-rad[ip+ng]*angArray[21];
         bas[2] = rad[ip]*2*angArray[6]-rad[ip+ng]*angArray[22];
         bas[3] = rad[ip]*angArray[7]-rad[ip+ng]*angArray[23];
         bas[4] = rad[ip]*angArray[8]-rad[ip+ng]*angArray[24];
         bas[5] = rad[ip]*angArray[9]-rad[ip+ng]*angArray[25];
         bas[6] = -rad[ip+ng]*angArray[26];
         bas[7] = -rad[ip+ng]*angArray[27];
         bas[8] = -rad[ip+ng]*angArray[28];
         bas[9] = -rad[ip+ng]*angArray[29];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*4*angArray[10]-rad[ip+ng]*angArray[35];
         bas[1] = rad[ip]*3*angArray[11]-rad[ip+ng]*angArray[36];
         bas[2] = rad[ip]*3*angArray[12]-rad[ip+ng]*angArray[37];
         bas[3] = rad[ip]*2*angArray[13]-rad[ip+ng]*angArray[38];
         bas[4] = rad[ip]*2*angArray[14]-rad[ip+ng]*angArray[39];
         bas[5] = rad[ip]*2*angArray[15]-rad[ip+ng]*angArray[40];
         bas[6] = rad[ip]*angArray[16]-rad[ip+ng]*angArray[41];
         bas[7] = rad[ip]*angArray[17]-rad[ip+ng]*angArray[42];
         bas[8] = rad[ip]*angArray[18]-rad[ip+ng]*angArray[43];
         bas[9] = rad[ip]*angArray[19]-rad[ip+ng]*angArray[44];
         bas[10] = -rad[ip+ng]*angArray[45];
         bas[11] = -rad[ip+ng]*angArray[46];
         bas[12] = -rad[ip+ng]*angArray[47];
         bas[13] = -rad[ip+ng]*angArray[48];
         bas[14] = -rad[ip+ng]*angArray[49];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*5*angArray[20]-rad[ip+ng]*angArray[56];
         bas[1] = rad[ip]*4*angArray[21]-rad[ip+ng]*angArray[57];
         bas[2] = rad[ip]*4*angArray[22]-rad[ip+ng]*angArray[58];
         bas[3] = rad[ip]*3*angArray[23]-rad[ip+ng]*angArray[59];
         bas[4] = rad[ip]*3*angArray[24]-rad[ip+ng]*angArray[60];
         bas[5] = rad[ip]*3*angArray[25]-rad[ip+ng]*angArray[61];
         bas[6] = rad[ip]*2*angArray[26]-rad[ip+ng]*angArray[62];
         bas[7] = rad[ip]*2*angArray[27]-rad[ip+ng]*angArray[63];
         bas[8] = rad[ip]*2*angArray[28]-rad[ip+ng]*angArray[64];
         bas[9] = rad[ip]*2*angArray[29]-rad[ip+ng]*angArray[65];
         bas[10] = rad[ip]*angArray[30]-rad[ip+ng]*angArray[66];
         bas[11] = rad[ip]*angArray[31]-rad[ip+ng]*angArray[67];
         bas[12] = rad[ip]*angArray[32]-rad[ip+ng]*angArray[68];
         bas[13] = rad[ip]*angArray[33]-rad[ip+ng]*angArray[69];
         bas[14] = rad[ip]*angArray[34]-rad[ip+ng]*angArray[70];
         bas[15] = -rad[ip+ng]*angArray[71];
         bas[16] = -rad[ip+ng]*angArray[72];
         bas[17] = -rad[ip+ng]*angArray[73];
         bas[18] = -rad[ip+ng]*angArray[74];
         bas[19] = -rad[ip+ng]*angArray[75];
         bas[20] = -rad[ip+ng]*angArray[76];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*6*angArray[35]-rad[ip+ng]*angArray[84];
         bas[1] = rad[ip]*5*angArray[36]-rad[ip+ng]*angArray[85];
         bas[2] = rad[ip]*5*angArray[37]-rad[ip+ng]*angArray[86];
         bas[3] = rad[ip]*4*angArray[38]-rad[ip+ng]*angArray[87];
         bas[4] = rad[ip]*4*angArray[39]-rad[ip+ng]*angArray[88];
         bas[5] = rad[ip]*4*angArray[40]-rad[ip+ng]*angArray[89];
         bas[6] = rad[ip]*3*angArray[41]-rad[ip+ng]*angArray[90];
         bas[7] = rad[ip]*3*angArray[42]-rad[ip+ng]*angArray[91];
         bas[8] = rad[ip]*3*angArray[43]-rad[ip+ng]*angArray[92];
         bas[9] = rad[ip]*3*angArray[44]-rad[ip+ng]*angArray[93];
         bas[10] = rad[ip]*2*angArray[45]-rad[ip+ng]*angArray[94];
         bas[11] = rad[ip]*2*angArray[46]-rad[ip+ng]*angArray[95];
         bas[12] = rad[ip]*2*angArray[47]-rad[ip+ng]*angArray[96];
         bas[13] = rad[ip]*2*angArray[48]-rad[ip+ng]*angArray[97];
         bas[14] = rad[ip]*2*angArray[49]-rad[ip+ng]*angArray[98];
         bas[15] = rad[ip]*angArray[50]-rad[ip+ng]*angArray[99];
         bas[16] = rad[ip]*angArray[51]-rad[ip+ng]*angArray[100];
         bas[17] = rad[ip]*angArray[52]-rad[ip+ng]*angArray[101];
         bas[18] = rad[ip]*angArray[53]-rad[ip+ng]*angArray[102];
         bas[19] = rad[ip]*angArray[54]-rad[ip+ng]*angArray[103];
         bas[20] = rad[ip]*angArray[55]-rad[ip+ng]*angArray[104];
         bas[21] = -rad[ip+ng]*angArray[105];
         bas[22] = -rad[ip+ng]*angArray[106];
         bas[23] = -rad[ip+ng]*angArray[107];
         bas[24] = -rad[ip+ng]*angArray[108];
         bas[25] = -rad[ip+ng]*angArray[109];
         bas[26] = -rad[ip+ng]*angArray[110];
         bas[27] = -rad[ip+ng]*angArray[111];
      }

   }


   // now we do derivatives for the given basis set to Y
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[2];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[5];
         bas[1] = rad[ip]*angArray[0]-rad[ip+ng]*angArray[7];
         bas[2] = -rad[ip+ng]*angArray[8];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[11];
         bas[1] = rad[ip]*angArray[1]-rad[ip+ng]*angArray[13];
         bas[2] = -rad[ip+ng]*angArray[14];
         bas[3] = rad[ip]*2*angArray[2]-rad[ip+ng]*angArray[16];
         bas[4] = rad[ip]*angArray[3]-rad[ip+ng]*angArray[17];
         bas[5] = -rad[ip+ng]*angArray[18];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[21];
         bas[1] = rad[ip]*angArray[4]-rad[ip+ng]*angArray[23];
         bas[2] = -rad[ip+ng]*angArray[24];
         bas[3] = rad[ip]*2*angArray[5]-rad[ip+ng]*angArray[26];
         bas[4] = rad[ip]*angArray[6]-rad[ip+ng]*angArray[27];
         bas[5] = -rad[ip+ng]*angArray[28];
         bas[6] = rad[ip]*3*angArray[7]-rad[ip+ng]*angArray[30];
         bas[7] = rad[ip]*2*angArray[8]-rad[ip+ng]*angArray[31];
         bas[8] = rad[ip]*angArray[9]-rad[ip+ng]*angArray[32];
         bas[9] = -rad[ip+ng]*angArray[33];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[36];
         bas[1] = rad[ip]*angArray[10]-rad[ip+ng]*angArray[38];
         bas[2] = -rad[ip+ng]*angArray[39];
         bas[3] = rad[ip]*2*angArray[11]-rad[ip+ng]*angArray[41];
         bas[4] = rad[ip]*angArray[12]-rad[ip+ng]*angArray[42];
         bas[5] = -rad[ip+ng]*angArray[43];
         bas[6] = rad[ip]*3*angArray[13]-rad[ip+ng]*angArray[45];
         bas[7] = rad[ip]*2*angArray[14]-rad[ip+ng]*angArray[46];
         bas[8] = rad[ip]*angArray[15]-rad[ip+ng]*angArray[47];
         bas[9] = -rad[ip+ng]*angArray[48];
         bas[10] = rad[ip]*4*angArray[16]-rad[ip+ng]*angArray[50];
         bas[11] = rad[ip]*3*angArray[17]-rad[ip+ng]*angArray[51];
         bas[12] = rad[ip]*2*angArray[18]-rad[ip+ng]*angArray[52];
         bas[13] = rad[ip]*angArray[19]-rad[ip+ng]*angArray[53];
         bas[14] = -rad[ip+ng]*angArray[54];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[57];
         bas[1] = rad[ip]*angArray[20]-rad[ip+ng]*angArray[59];
         bas[2] = -rad[ip+ng]*angArray[60];
         bas[3] = rad[ip]*2*angArray[21]-rad[ip+ng]*angArray[62];
         bas[4] = rad[ip]*angArray[22]-rad[ip+ng]*angArray[63];
         bas[5] = -rad[ip+ng]*angArray[64];
         bas[6] = rad[ip]*3*angArray[23]-rad[ip+ng]*angArray[66];
         bas[7] = rad[ip]*2*angArray[24]-rad[ip+ng]*angArray[67];
         bas[8] = rad[ip]*angArray[25]-rad[ip+ng]*angArray[68];
         bas[9] = -rad[ip+ng]*angArray[69];
         bas[10] = rad[ip]*4*angArray[26]-rad[ip+ng]*angArray[71];
         bas[11] = rad[ip]*3*angArray[27]-rad[ip+ng]*angArray[72];
         bas[12] = rad[ip]*2*angArray[28]-rad[ip+ng]*angArray[73];
         bas[13] = rad[ip]*angArray[29]-rad[ip+ng]*angArray[74];
         bas[14] = -rad[ip+ng]*angArray[75];
         bas[15] = rad[ip]*5*angArray[30]-rad[ip+ng]*angArray[77];
         bas[16] = rad[ip]*4*angArray[31]-rad[ip+ng]*angArray[78];
         bas[17] = rad[ip]*3*angArray[32]-rad[ip+ng]*angArray[79];
         bas[18] = rad[ip]*2*angArray[33]-rad[ip+ng]*angArray[80];
         bas[19] = rad[ip]*angArray[34]-rad[ip+ng]*angArray[81];
         bas[20] = -rad[ip+ng]*angArray[82];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[85];
         bas[1] = rad[ip]*angArray[35]-rad[ip+ng]*angArray[87];
         bas[2] = -rad[ip+ng]*angArray[88];
         bas[3] = rad[ip]*2*angArray[36]-rad[ip+ng]*angArray[90];
         bas[4] = rad[ip]*angArray[37]-rad[ip+ng]*angArray[91];
         bas[5] = -rad[ip+ng]*angArray[92];
         bas[6] = rad[ip]*3*angArray[38]-rad[ip+ng]*angArray[94];
         bas[7] = rad[ip]*2*angArray[39]-rad[ip+ng]*angArray[95];
         bas[8] = rad[ip]*angArray[40]-rad[ip+ng]*angArray[96];
         bas[9] = -rad[ip+ng]*angArray[97];
         bas[10] = rad[ip]*4*angArray[41]-rad[ip+ng]*angArray[99];
         bas[11] = rad[ip]*3*angArray[42]-rad[ip+ng]*angArray[100];
         bas[12] = rad[ip]*2*angArray[43]-rad[ip+ng]*angArray[101];
         bas[13] = rad[ip]*angArray[44]-rad[ip+ng]*angArray[102];
         bas[14] = -rad[ip+ng]*angArray[103];
         bas[15] = rad[ip]*5*angArray[45]-rad[ip+ng]*angArray[105];
         bas[16] = rad[ip]*4*angArray[46]-rad[ip+ng]*angArray[106];
         bas[17] = rad[ip]*3*angArray[47]-rad[ip+ng]*angArray[107];
         bas[18] = rad[ip]*2*angArray[48]-rad[ip+ng]*angArray[108];
         bas[19] = rad[ip]*angArray[49]-rad[ip+ng]*angArray[109];
         bas[20] = -rad[ip+ng]*angArray[110];
         bas[21] = rad[ip]*6*angArray[50]-rad[ip+ng]*angArray[112];
         bas[22] = rad[ip]*5*angArray[51]-rad[ip+ng]*angArray[113];
         bas[23] = rad[ip]*4*angArray[52]-rad[ip+ng]*angArray[114];
         bas[24] = rad[ip]*3*angArray[53]-rad[ip+ng]*angArray[115];
         bas[25] = rad[ip]*2*angArray[54]-rad[ip+ng]*angArray[116];
         bas[26] = rad[ip]*angArray[55]-rad[ip+ng]*angArray[117];
         bas[27] = -rad[ip+ng]*angArray[118];
      }

   }


   // now we do derivatives for the given basis set to Z
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[3];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[6];
         bas[1] = -rad[ip+ng]*angArray[8];
         bas[2] = rad[ip]*angArray[0]-rad[ip+ng]*angArray[9];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[12];
         bas[1] = -rad[ip+ng]*angArray[14];
         bas[2] = rad[ip]*angArray[1]-rad[ip+ng]*angArray[15];
         bas[3] = -rad[ip+ng]*angArray[17];
         bas[4] = rad[ip]*angArray[2]-rad[ip+ng]*angArray[18];
         bas[5] = rad[ip]*2*angArray[3]-rad[ip+ng]*angArray[19];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[22];
         bas[1] = -rad[ip+ng]*angArray[24];
         bas[2] = rad[ip]*angArray[4]-rad[ip+ng]*angArray[25];
         bas[3] = -rad[ip+ng]*angArray[27];
         bas[4] = rad[ip]*angArray[5]-rad[ip+ng]*angArray[28];
         bas[5] = rad[ip]*2*angArray[6]-rad[ip+ng]*angArray[29];
         bas[6] = -rad[ip+ng]*angArray[31];
         bas[7] = rad[ip]*angArray[7]-rad[ip+ng]*angArray[32];
         bas[8] = rad[ip]*2*angArray[8]-rad[ip+ng]*angArray[33];
         bas[9] = rad[ip]*3*angArray[9]-rad[ip+ng]*angArray[34];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[37];
         bas[1] = -rad[ip+ng]*angArray[39];
         bas[2] = rad[ip]*angArray[10]-rad[ip+ng]*angArray[40];
         bas[3] = -rad[ip+ng]*angArray[42];
         bas[4] = rad[ip]*angArray[11]-rad[ip+ng]*angArray[43];
         bas[5] = rad[ip]*2*angArray[12]-rad[ip+ng]*angArray[44];
         bas[6] = -rad[ip+ng]*angArray[46];
         bas[7] = rad[ip]*angArray[13]-rad[ip+ng]*angArray[47];
         bas[8] = rad[ip]*2*angArray[14]-rad[ip+ng]*angArray[48];
         bas[9] = rad[ip]*3*angArray[15]-rad[ip+ng]*angArray[49];
         bas[10] = -rad[ip+ng]*angArray[51];
         bas[11] = rad[ip]*angArray[16]-rad[ip+ng]*angArray[52];
         bas[12] = rad[ip]*2*angArray[17]-rad[ip+ng]*angArray[53];
         bas[13] = rad[ip]*3*angArray[18]-rad[ip+ng]*angArray[54];
         bas[14] = rad[ip]*4*angArray[19]-rad[ip+ng]*angArray[55];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[58];
         bas[1] = -rad[ip+ng]*angArray[60];
         bas[2] = rad[ip]*angArray[20]-rad[ip+ng]*angArray[61];
         bas[3] = -rad[ip+ng]*angArray[63];
         bas[4] = rad[ip]*angArray[21]-rad[ip+ng]*angArray[64];
         bas[5] = rad[ip]*2*angArray[22]-rad[ip+ng]*angArray[65];
         bas[6] = -rad[ip+ng]*angArray[67];
         bas[7] = rad[ip]*angArray[23]-rad[ip+ng]*angArray[68];
         bas[8] = rad[ip]*2*angArray[24]-rad[ip+ng]*angArray[69];
         bas[9] = rad[ip]*3*angArray[25]-rad[ip+ng]*angArray[70];
         bas[10] = -rad[ip+ng]*angArray[72];
         bas[11] = rad[ip]*angArray[26]-rad[ip+ng]*angArray[73];
         bas[12] = rad[ip]*2*angArray[27]-rad[ip+ng]*angArray[74];
         bas[13] = rad[ip]*3*angArray[28]-rad[ip+ng]*angArray[75];
         bas[14] = rad[ip]*4*angArray[29]-rad[ip+ng]*angArray[76];
         bas[15] = -rad[ip+ng]*angArray[78];
         bas[16] = rad[ip]*angArray[30]-rad[ip+ng]*angArray[79];
         bas[17] = rad[ip]*2*angArray[31]-rad[ip+ng]*angArray[80];
         bas[18] = rad[ip]*3*angArray[32]-rad[ip+ng]*angArray[81];
         bas[19] = rad[ip]*4*angArray[33]-rad[ip+ng]*angArray[82];
         bas[20] = rad[ip]*5*angArray[34]-rad[ip+ng]*angArray[83];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[86];
         bas[1] = -rad[ip+ng]*angArray[88];
         bas[2] = rad[ip]*angArray[35]-rad[ip+ng]*angArray[89];
         bas[3] = -rad[ip+ng]*angArray[91];
         bas[4] = rad[ip]*angArray[36]-rad[ip+ng]*angArray[92];
         bas[5] = rad[ip]*2*angArray[37]-rad[ip+ng]*angArray[93];
         bas[6] = -rad[ip+ng]*angArray[95];
         bas[7] = rad[ip]*angArray[38]-rad[ip+ng]*angArray[96];
         bas[8] = rad[ip]*2*angArray[39]-rad[ip+ng]*angArray[97];
         bas[9] = rad[ip]*3*angArray[40]-rad[ip+ng]*angArray[98];
         bas[10] = -rad[ip+ng]*angArray[100];
         bas[11] = rad[ip]*angArray[41]-rad[ip+ng]*angArray[101];
         bas[12] = rad[ip]*2*angArray[42]-rad[ip+ng]*angArray[102];
         bas[13] = rad[ip]*3*angArray[43]-rad[ip+ng]*angArray[103];
         bas[14] = rad[ip]*4*angArray[44]-rad[ip+ng]*angArray[104];
         bas[15] = -rad[ip+ng]*angArray[106];
         bas[16] = rad[ip]*angArray[45]-rad[ip+ng]*angArray[107];
         bas[17] = rad[ip]*2*angArray[46]-rad[ip+ng]*angArray[108];
         bas[18] = rad[ip]*3*angArray[47]-rad[ip+ng]*angArray[109];
         bas[19] = rad[ip]*4*angArray[48]-rad[ip+ng]*angArray[110];
         bas[20] = rad[ip]*5*angArray[49]-rad[ip+ng]*angArray[111];
         bas[21] = -rad[ip+ng]*angArray[113];
         bas[22] = rad[ip]*angArray[50]-rad[ip+ng]*angArray[114];
         bas[23] = rad[ip]*2*angArray[51]-rad[ip+ng]*angArray[115];
         bas[24] = rad[ip]*3*angArray[52]-rad[ip+ng]*angArray[116];
         bas[25] = rad[ip]*4*angArray[53]-rad[ip+ng]*angArray[117];
         bas[26] = rad[ip]*5*angArray[54]-rad[ip+ng]*angArray[118];
         bas[27] = rad[ip]*6*angArray[55]-rad[ip+ng]*angArray[119];
      }

   }


}


