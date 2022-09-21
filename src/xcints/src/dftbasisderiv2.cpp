/**
 * This function is used to generate 2 derivatives for basis set 
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

void BatchBasis::dftbasisderiv2(const UInt& ng, const UInt& L, const UInt& nTolCarBas, const Double* ang, const Double* rad, Double* basis) const 
{

   // now we set up the nBas for the computation
   UInt nBas = (L+1)*(L+2)/2;

   // now we do derivatives for the given basis set to XX
   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*angArray[4];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*angArray[10];
         bas[1] = -rad[ip+ng]*angArray[2]+rad[ip+2*ng]*angArray[11];
         bas[2] = -rad[ip+ng]*angArray[3]+rad[ip+2*ng]*angArray[12];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*2*angArray[0]-rad[ip+ng]*5*angArray[4]+rad[ip+2*ng]*angArray[20];
         bas[1] = -rad[ip+ng]*3*angArray[5]+rad[ip+2*ng]*angArray[21];
         bas[2] = -rad[ip+ng]*3*angArray[6]+rad[ip+2*ng]*angArray[22];
         bas[3] = -rad[ip+ng]*angArray[7]+rad[ip+2*ng]*angArray[23];
         bas[4] = -rad[ip+ng]*angArray[8]+rad[ip+2*ng]*angArray[24];
         bas[5] = -rad[ip+ng]*angArray[9]+rad[ip+2*ng]*angArray[25];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*6*angArray[1]-rad[ip+ng]*7*angArray[10]+rad[ip+2*ng]*angArray[35];
         bas[1] = rad[ip]*2*angArray[2]-rad[ip+ng]*5*angArray[11]+rad[ip+2*ng]*angArray[36];
         bas[2] = rad[ip]*2*angArray[3]-rad[ip+ng]*5*angArray[12]+rad[ip+2*ng]*angArray[37];
         bas[3] = -rad[ip+ng]*3*angArray[13]+rad[ip+2*ng]*angArray[38];
         bas[4] = -rad[ip+ng]*3*angArray[14]+rad[ip+2*ng]*angArray[39];
         bas[5] = -rad[ip+ng]*3*angArray[15]+rad[ip+2*ng]*angArray[40];
         bas[6] = -rad[ip+ng]*angArray[16]+rad[ip+2*ng]*angArray[41];
         bas[7] = -rad[ip+ng]*angArray[17]+rad[ip+2*ng]*angArray[42];
         bas[8] = -rad[ip+ng]*angArray[18]+rad[ip+2*ng]*angArray[43];
         bas[9] = -rad[ip+ng]*angArray[19]+rad[ip+2*ng]*angArray[44];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*12*angArray[4]-rad[ip+ng]*9*angArray[20]+rad[ip+2*ng]*angArray[56];
         bas[1] = rad[ip]*6*angArray[5]-rad[ip+ng]*7*angArray[21]+rad[ip+2*ng]*angArray[57];
         bas[2] = rad[ip]*6*angArray[6]-rad[ip+ng]*7*angArray[22]+rad[ip+2*ng]*angArray[58];
         bas[3] = rad[ip]*2*angArray[7]-rad[ip+ng]*5*angArray[23]+rad[ip+2*ng]*angArray[59];
         bas[4] = rad[ip]*2*angArray[8]-rad[ip+ng]*5*angArray[24]+rad[ip+2*ng]*angArray[60];
         bas[5] = rad[ip]*2*angArray[9]-rad[ip+ng]*5*angArray[25]+rad[ip+2*ng]*angArray[61];
         bas[6] = -rad[ip+ng]*3*angArray[26]+rad[ip+2*ng]*angArray[62];
         bas[7] = -rad[ip+ng]*3*angArray[27]+rad[ip+2*ng]*angArray[63];
         bas[8] = -rad[ip+ng]*3*angArray[28]+rad[ip+2*ng]*angArray[64];
         bas[9] = -rad[ip+ng]*3*angArray[29]+rad[ip+2*ng]*angArray[65];
         bas[10] = -rad[ip+ng]*angArray[30]+rad[ip+2*ng]*angArray[66];
         bas[11] = -rad[ip+ng]*angArray[31]+rad[ip+2*ng]*angArray[67];
         bas[12] = -rad[ip+ng]*angArray[32]+rad[ip+2*ng]*angArray[68];
         bas[13] = -rad[ip+ng]*angArray[33]+rad[ip+2*ng]*angArray[69];
         bas[14] = -rad[ip+ng]*angArray[34]+rad[ip+2*ng]*angArray[70];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*20*angArray[10]-rad[ip+ng]*11*angArray[35]+rad[ip+2*ng]*angArray[84];
         bas[1] = rad[ip]*12*angArray[11]-rad[ip+ng]*9*angArray[36]+rad[ip+2*ng]*angArray[85];
         bas[2] = rad[ip]*12*angArray[12]-rad[ip+ng]*9*angArray[37]+rad[ip+2*ng]*angArray[86];
         bas[3] = rad[ip]*6*angArray[13]-rad[ip+ng]*7*angArray[38]+rad[ip+2*ng]*angArray[87];
         bas[4] = rad[ip]*6*angArray[14]-rad[ip+ng]*7*angArray[39]+rad[ip+2*ng]*angArray[88];
         bas[5] = rad[ip]*6*angArray[15]-rad[ip+ng]*7*angArray[40]+rad[ip+2*ng]*angArray[89];
         bas[6] = rad[ip]*2*angArray[16]-rad[ip+ng]*5*angArray[41]+rad[ip+2*ng]*angArray[90];
         bas[7] = rad[ip]*2*angArray[17]-rad[ip+ng]*5*angArray[42]+rad[ip+2*ng]*angArray[91];
         bas[8] = rad[ip]*2*angArray[18]-rad[ip+ng]*5*angArray[43]+rad[ip+2*ng]*angArray[92];
         bas[9] = rad[ip]*2*angArray[19]-rad[ip+ng]*5*angArray[44]+rad[ip+2*ng]*angArray[93];
         bas[10] = -rad[ip+ng]*3*angArray[45]+rad[ip+2*ng]*angArray[94];
         bas[11] = -rad[ip+ng]*3*angArray[46]+rad[ip+2*ng]*angArray[95];
         bas[12] = -rad[ip+ng]*3*angArray[47]+rad[ip+2*ng]*angArray[96];
         bas[13] = -rad[ip+ng]*3*angArray[48]+rad[ip+2*ng]*angArray[97];
         bas[14] = -rad[ip+ng]*3*angArray[49]+rad[ip+2*ng]*angArray[98];
         bas[15] = -rad[ip+ng]*angArray[50]+rad[ip+2*ng]*angArray[99];
         bas[16] = -rad[ip+ng]*angArray[51]+rad[ip+2*ng]*angArray[100];
         bas[17] = -rad[ip+ng]*angArray[52]+rad[ip+2*ng]*angArray[101];
         bas[18] = -rad[ip+ng]*angArray[53]+rad[ip+2*ng]*angArray[102];
         bas[19] = -rad[ip+ng]*angArray[54]+rad[ip+2*ng]*angArray[103];
         bas[20] = -rad[ip+ng]*angArray[55]+rad[ip+2*ng]*angArray[104];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*30*angArray[20]-rad[ip+ng]*13*angArray[56]+rad[ip+2*ng]*angArray[120];
         bas[1] = rad[ip]*20*angArray[21]-rad[ip+ng]*11*angArray[57]+rad[ip+2*ng]*angArray[121];
         bas[2] = rad[ip]*20*angArray[22]-rad[ip+ng]*11*angArray[58]+rad[ip+2*ng]*angArray[122];
         bas[3] = rad[ip]*12*angArray[23]-rad[ip+ng]*9*angArray[59]+rad[ip+2*ng]*angArray[123];
         bas[4] = rad[ip]*12*angArray[24]-rad[ip+ng]*9*angArray[60]+rad[ip+2*ng]*angArray[124];
         bas[5] = rad[ip]*12*angArray[25]-rad[ip+ng]*9*angArray[61]+rad[ip+2*ng]*angArray[125];
         bas[6] = rad[ip]*6*angArray[26]-rad[ip+ng]*7*angArray[62]+rad[ip+2*ng]*angArray[126];
         bas[7] = rad[ip]*6*angArray[27]-rad[ip+ng]*7*angArray[63]+rad[ip+2*ng]*angArray[127];
         bas[8] = rad[ip]*6*angArray[28]-rad[ip+ng]*7*angArray[64]+rad[ip+2*ng]*angArray[128];
         bas[9] = rad[ip]*6*angArray[29]-rad[ip+ng]*7*angArray[65]+rad[ip+2*ng]*angArray[129];
         bas[10] = rad[ip]*2*angArray[30]-rad[ip+ng]*5*angArray[66]+rad[ip+2*ng]*angArray[130];
         bas[11] = rad[ip]*2*angArray[31]-rad[ip+ng]*5*angArray[67]+rad[ip+2*ng]*angArray[131];
         bas[12] = rad[ip]*2*angArray[32]-rad[ip+ng]*5*angArray[68]+rad[ip+2*ng]*angArray[132];
         bas[13] = rad[ip]*2*angArray[33]-rad[ip+ng]*5*angArray[69]+rad[ip+2*ng]*angArray[133];
         bas[14] = rad[ip]*2*angArray[34]-rad[ip+ng]*5*angArray[70]+rad[ip+2*ng]*angArray[134];
         bas[15] = -rad[ip+ng]*3*angArray[71]+rad[ip+2*ng]*angArray[135];
         bas[16] = -rad[ip+ng]*3*angArray[72]+rad[ip+2*ng]*angArray[136];
         bas[17] = -rad[ip+ng]*3*angArray[73]+rad[ip+2*ng]*angArray[137];
         bas[18] = -rad[ip+ng]*3*angArray[74]+rad[ip+2*ng]*angArray[138];
         bas[19] = -rad[ip+ng]*3*angArray[75]+rad[ip+2*ng]*angArray[139];
         bas[20] = -rad[ip+ng]*3*angArray[76]+rad[ip+2*ng]*angArray[140];
         bas[21] = -rad[ip+ng]*angArray[77]+rad[ip+2*ng]*angArray[141];
         bas[22] = -rad[ip+ng]*angArray[78]+rad[ip+2*ng]*angArray[142];
         bas[23] = -rad[ip+ng]*angArray[79]+rad[ip+2*ng]*angArray[143];
         bas[24] = -rad[ip+ng]*angArray[80]+rad[ip+2*ng]*angArray[144];
         bas[25] = -rad[ip+ng]*angArray[81]+rad[ip+2*ng]*angArray[145];
         bas[26] = -rad[ip+ng]*angArray[82]+rad[ip+2*ng]*angArray[146];
         bas[27] = -rad[ip+ng]*angArray[83]+rad[ip+2*ng]*angArray[147];
      }

   }


   // now we do derivatives for the given basis set to XY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[5];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[2]+rad[ip+2*ng]*angArray[11];
         bas[1] = -rad[ip+ng]*angArray[1]+rad[ip+2*ng]*angArray[13];
         bas[2] = rad[ip+2*ng]*angArray[14];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*angArray[21];
         bas[1] = rad[ip]*angArray[0]-rad[ip+ng]*(angArray[4]+angArray[7])+rad[ip+2*ng]*angArray[23];
         bas[2] = -rad[ip+ng]*angArray[8]+rad[ip+2*ng]*angArray[24];
         bas[3] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*angArray[26];
         bas[4] = -rad[ip+ng]*angArray[6]+rad[ip+2*ng]*angArray[27];
         bas[5] = rad[ip+2*ng]*angArray[28];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*3*angArray[11]+rad[ip+2*ng]*angArray[36];
         bas[1] = rad[ip]*2*angArray[1]-rad[ip+ng]*(2*angArray[13]+angArray[10])+rad[ip+2*ng]*angArray[38];
         bas[2] = -rad[ip+ng]*2*angArray[14]+rad[ip+2*ng]*angArray[39];
         bas[3] = rad[ip]*2*angArray[2]-rad[ip+ng]*(2*angArray[11]+angArray[16])+rad[ip+2*ng]*angArray[41];
         bas[4] = rad[ip]*angArray[3]-rad[ip+ng]*(angArray[12]+angArray[17])+rad[ip+2*ng]*angArray[42];
         bas[5] = -rad[ip+ng]*angArray[18]+rad[ip+2*ng]*angArray[43];
         bas[6] = -rad[ip+ng]*3*angArray[13]+rad[ip+2*ng]*angArray[45];
         bas[7] = -rad[ip+ng]*2*angArray[14]+rad[ip+2*ng]*angArray[46];
         bas[8] = -rad[ip+ng]*angArray[15]+rad[ip+2*ng]*angArray[47];
         bas[9] = rad[ip+2*ng]*angArray[48];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*4*angArray[21]+rad[ip+2*ng]*angArray[57];
         bas[1] = rad[ip]*3*angArray[4]-rad[ip+ng]*(angArray[20]+3*angArray[23])+rad[ip+2*ng]*angArray[59];
         bas[2] = -rad[ip+ng]*3*angArray[24]+rad[ip+2*ng]*angArray[60];
         bas[3] = rad[ip]*4*angArray[5]-rad[ip+ng]*(2*angArray[26]+2*angArray[21])+rad[ip+2*ng]*angArray[62];
         bas[4] = rad[ip]*2*angArray[6]-rad[ip+ng]*(2*angArray[27]+angArray[22])+rad[ip+2*ng]*angArray[63];
         bas[5] = -rad[ip+ng]*2*angArray[28]+rad[ip+2*ng]*angArray[64];
         bas[6] = rad[ip]*3*angArray[7]-rad[ip+ng]*(3*angArray[23]+angArray[30])+rad[ip+2*ng]*angArray[66];
         bas[7] = rad[ip]*2*angArray[8]-rad[ip+ng]*(2*angArray[24]+angArray[31])+rad[ip+2*ng]*angArray[67];
         bas[8] = rad[ip]*angArray[9]-rad[ip+ng]*(angArray[25]+angArray[32])+rad[ip+2*ng]*angArray[68];
         bas[9] = -rad[ip+ng]*angArray[33]+rad[ip+2*ng]*angArray[69];
         bas[10] = -rad[ip+ng]*4*angArray[26]+rad[ip+2*ng]*angArray[71];
         bas[11] = -rad[ip+ng]*3*angArray[27]+rad[ip+2*ng]*angArray[72];
         bas[12] = -rad[ip+ng]*2*angArray[28]+rad[ip+2*ng]*angArray[73];
         bas[13] = -rad[ip+ng]*angArray[29]+rad[ip+2*ng]*angArray[74];
         bas[14] = rad[ip+2*ng]*angArray[75];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*5*angArray[36]+rad[ip+2*ng]*angArray[85];
         bas[1] = rad[ip]*4*angArray[10]-rad[ip+ng]*(4*angArray[38]+angArray[35])+rad[ip+2*ng]*angArray[87];
         bas[2] = -rad[ip+ng]*4*angArray[39]+rad[ip+2*ng]*angArray[88];
         bas[3] = rad[ip]*6*angArray[11]-rad[ip+ng]*(2*angArray[36]+3*angArray[41])+rad[ip+2*ng]*angArray[90];
         bas[4] = rad[ip]*3*angArray[12]-rad[ip+ng]*(angArray[37]+3*angArray[42])+rad[ip+2*ng]*angArray[91];
         bas[5] = -rad[ip+ng]*3*angArray[43]+rad[ip+2*ng]*angArray[92];
         bas[6] = rad[ip]*6*angArray[13]-rad[ip+ng]*(2*angArray[45]+3*angArray[38])+rad[ip+2*ng]*angArray[94];
         bas[7] = rad[ip]*4*angArray[14]-rad[ip+ng]*(2*angArray[46]+2*angArray[39])+rad[ip+2*ng]*angArray[95];
         bas[8] = rad[ip]*2*angArray[15]-rad[ip+ng]*(2*angArray[47]+angArray[40])+rad[ip+2*ng]*angArray[96];
         bas[9] = -rad[ip+ng]*2*angArray[48]+rad[ip+2*ng]*angArray[97];
         bas[10] = rad[ip]*4*angArray[16]-rad[ip+ng]*(4*angArray[41]+angArray[50])+rad[ip+2*ng]*angArray[99];
         bas[11] = rad[ip]*3*angArray[17]-rad[ip+ng]*(3*angArray[42]+angArray[51])+rad[ip+2*ng]*angArray[100];
         bas[12] = rad[ip]*2*angArray[18]-rad[ip+ng]*(2*angArray[43]+angArray[52])+rad[ip+2*ng]*angArray[101];
         bas[13] = rad[ip]*angArray[19]-rad[ip+ng]*(angArray[44]+angArray[53])+rad[ip+2*ng]*angArray[102];
         bas[14] = -rad[ip+ng]*angArray[54]+rad[ip+2*ng]*angArray[103];
         bas[15] = -rad[ip+ng]*5*angArray[45]+rad[ip+2*ng]*angArray[105];
         bas[16] = -rad[ip+ng]*4*angArray[46]+rad[ip+2*ng]*angArray[106];
         bas[17] = -rad[ip+ng]*3*angArray[47]+rad[ip+2*ng]*angArray[107];
         bas[18] = -rad[ip+ng]*2*angArray[48]+rad[ip+2*ng]*angArray[108];
         bas[19] = -rad[ip+ng]*angArray[49]+rad[ip+2*ng]*angArray[109];
         bas[20] = rad[ip+2*ng]*angArray[110];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[57]+rad[ip+2*ng]*angArray[121];
         bas[1] = rad[ip]*5*angArray[20]-rad[ip+ng]*(5*angArray[59]+angArray[56])+rad[ip+2*ng]*angArray[123];
         bas[2] = -rad[ip+ng]*5*angArray[60]+rad[ip+2*ng]*angArray[124];
         bas[3] = rad[ip]*8*angArray[21]-rad[ip+ng]*(2*angArray[57]+4*angArray[62])+rad[ip+2*ng]*angArray[126];
         bas[4] = rad[ip]*4*angArray[22]-rad[ip+ng]*(4*angArray[63]+angArray[58])+rad[ip+2*ng]*angArray[127];
         bas[5] = -rad[ip+ng]*4*angArray[64]+rad[ip+2*ng]*angArray[128];
         bas[6] = rad[ip]*9*angArray[23]-rad[ip+ng]*(3*angArray[59]+3*angArray[66])+rad[ip+2*ng]*angArray[130];
         bas[7] = rad[ip]*6*angArray[24]-rad[ip+ng]*(2*angArray[60]+3*angArray[67])+rad[ip+2*ng]*angArray[131];
         bas[8] = rad[ip]*3*angArray[25]-rad[ip+ng]*(angArray[61]+3*angArray[68])+rad[ip+2*ng]*angArray[132];
         bas[9] = -rad[ip+ng]*3*angArray[69]+rad[ip+2*ng]*angArray[133];
         bas[10] = rad[ip]*8*angArray[26]-rad[ip+ng]*(4*angArray[62]+2*angArray[71])+rad[ip+2*ng]*angArray[135];
         bas[11] = rad[ip]*6*angArray[27]-rad[ip+ng]*(2*angArray[72]+3*angArray[63])+rad[ip+2*ng]*angArray[136];
         bas[12] = rad[ip]*4*angArray[28]-rad[ip+ng]*(2*angArray[73]+2*angArray[64])+rad[ip+2*ng]*angArray[137];
         bas[13] = rad[ip]*2*angArray[29]-rad[ip+ng]*(2*angArray[74]+angArray[65])+rad[ip+2*ng]*angArray[138];
         bas[14] = -rad[ip+ng]*2*angArray[75]+rad[ip+2*ng]*angArray[139];
         bas[15] = rad[ip]*5*angArray[30]-rad[ip+ng]*(5*angArray[66]+angArray[77])+rad[ip+2*ng]*angArray[141];
         bas[16] = rad[ip]*4*angArray[31]-rad[ip+ng]*(4*angArray[67]+angArray[78])+rad[ip+2*ng]*angArray[142];
         bas[17] = rad[ip]*3*angArray[32]-rad[ip+ng]*(3*angArray[68]+angArray[79])+rad[ip+2*ng]*angArray[143];
         bas[18] = rad[ip]*2*angArray[33]-rad[ip+ng]*(2*angArray[69]+angArray[80])+rad[ip+2*ng]*angArray[144];
         bas[19] = rad[ip]*angArray[34]-rad[ip+ng]*(angArray[70]+angArray[81])+rad[ip+2*ng]*angArray[145];
         bas[20] = -rad[ip+ng]*angArray[82]+rad[ip+2*ng]*angArray[146];
         bas[21] = -rad[ip+ng]*6*angArray[71]+rad[ip+2*ng]*angArray[148];
         bas[22] = -rad[ip+ng]*5*angArray[72]+rad[ip+2*ng]*angArray[149];
         bas[23] = -rad[ip+ng]*4*angArray[73]+rad[ip+2*ng]*angArray[150];
         bas[24] = -rad[ip+ng]*3*angArray[74]+rad[ip+2*ng]*angArray[151];
         bas[25] = -rad[ip+ng]*2*angArray[75]+rad[ip+2*ng]*angArray[152];
         bas[26] = -rad[ip+ng]*angArray[76]+rad[ip+2*ng]*angArray[153];
         bas[27] = rad[ip+2*ng]*angArray[154];
      }

   }


   // now we do derivatives for the given basis set to YY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*angArray[7];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[1]+rad[ip+2*ng]*angArray[13];
         bas[1] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*angArray[16];
         bas[2] = -rad[ip+ng]*angArray[3]+rad[ip+2*ng]*angArray[17];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[4]+rad[ip+2*ng]*angArray[23];
         bas[1] = -rad[ip+ng]*3*angArray[5]+rad[ip+2*ng]*angArray[26];
         bas[2] = -rad[ip+ng]*angArray[6]+rad[ip+2*ng]*angArray[27];
         bas[3] = rad[ip]*2*angArray[0]-rad[ip+ng]*5*angArray[7]+rad[ip+2*ng]*angArray[30];
         bas[4] = -rad[ip+ng]*3*angArray[8]+rad[ip+2*ng]*angArray[31];
         bas[5] = -rad[ip+ng]*angArray[9]+rad[ip+2*ng]*angArray[32];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[10]+rad[ip+2*ng]*angArray[38];
         bas[1] = -rad[ip+ng]*3*angArray[11]+rad[ip+2*ng]*angArray[41];
         bas[2] = -rad[ip+ng]*angArray[12]+rad[ip+2*ng]*angArray[42];
         bas[3] = rad[ip]*2*angArray[1]-rad[ip+ng]*5*angArray[13]+rad[ip+2*ng]*angArray[45];
         bas[4] = -rad[ip+ng]*3*angArray[14]+rad[ip+2*ng]*angArray[46];
         bas[5] = -rad[ip+ng]*angArray[15]+rad[ip+2*ng]*angArray[47];
         bas[6] = rad[ip]*6*angArray[2]-rad[ip+ng]*7*angArray[16]+rad[ip+2*ng]*angArray[50];
         bas[7] = rad[ip]*2*angArray[3]-rad[ip+ng]*5*angArray[17]+rad[ip+2*ng]*angArray[51];
         bas[8] = -rad[ip+ng]*3*angArray[18]+rad[ip+2*ng]*angArray[52];
         bas[9] = -rad[ip+ng]*angArray[19]+rad[ip+2*ng]*angArray[53];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[20]+rad[ip+2*ng]*angArray[59];
         bas[1] = -rad[ip+ng]*3*angArray[21]+rad[ip+2*ng]*angArray[62];
         bas[2] = -rad[ip+ng]*angArray[22]+rad[ip+2*ng]*angArray[63];
         bas[3] = rad[ip]*2*angArray[4]-rad[ip+ng]*5*angArray[23]+rad[ip+2*ng]*angArray[66];
         bas[4] = -rad[ip+ng]*3*angArray[24]+rad[ip+2*ng]*angArray[67];
         bas[5] = -rad[ip+ng]*angArray[25]+rad[ip+2*ng]*angArray[68];
         bas[6] = rad[ip]*6*angArray[5]-rad[ip+ng]*7*angArray[26]+rad[ip+2*ng]*angArray[71];
         bas[7] = rad[ip]*2*angArray[6]-rad[ip+ng]*5*angArray[27]+rad[ip+2*ng]*angArray[72];
         bas[8] = -rad[ip+ng]*3*angArray[28]+rad[ip+2*ng]*angArray[73];
         bas[9] = -rad[ip+ng]*angArray[29]+rad[ip+2*ng]*angArray[74];
         bas[10] = rad[ip]*12*angArray[7]-rad[ip+ng]*9*angArray[30]+rad[ip+2*ng]*angArray[77];
         bas[11] = rad[ip]*6*angArray[8]-rad[ip+ng]*7*angArray[31]+rad[ip+2*ng]*angArray[78];
         bas[12] = rad[ip]*2*angArray[9]-rad[ip+ng]*5*angArray[32]+rad[ip+2*ng]*angArray[79];
         bas[13] = -rad[ip+ng]*3*angArray[33]+rad[ip+2*ng]*angArray[80];
         bas[14] = -rad[ip+ng]*angArray[34]+rad[ip+2*ng]*angArray[81];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[35]+rad[ip+2*ng]*angArray[87];
         bas[1] = -rad[ip+ng]*3*angArray[36]+rad[ip+2*ng]*angArray[90];
         bas[2] = -rad[ip+ng]*angArray[37]+rad[ip+2*ng]*angArray[91];
         bas[3] = rad[ip]*2*angArray[10]-rad[ip+ng]*5*angArray[38]+rad[ip+2*ng]*angArray[94];
         bas[4] = -rad[ip+ng]*3*angArray[39]+rad[ip+2*ng]*angArray[95];
         bas[5] = -rad[ip+ng]*angArray[40]+rad[ip+2*ng]*angArray[96];
         bas[6] = rad[ip]*6*angArray[11]-rad[ip+ng]*7*angArray[41]+rad[ip+2*ng]*angArray[99];
         bas[7] = rad[ip]*2*angArray[12]-rad[ip+ng]*5*angArray[42]+rad[ip+2*ng]*angArray[100];
         bas[8] = -rad[ip+ng]*3*angArray[43]+rad[ip+2*ng]*angArray[101];
         bas[9] = -rad[ip+ng]*angArray[44]+rad[ip+2*ng]*angArray[102];
         bas[10] = rad[ip]*12*angArray[13]-rad[ip+ng]*9*angArray[45]+rad[ip+2*ng]*angArray[105];
         bas[11] = rad[ip]*6*angArray[14]-rad[ip+ng]*7*angArray[46]+rad[ip+2*ng]*angArray[106];
         bas[12] = rad[ip]*2*angArray[15]-rad[ip+ng]*5*angArray[47]+rad[ip+2*ng]*angArray[107];
         bas[13] = -rad[ip+ng]*3*angArray[48]+rad[ip+2*ng]*angArray[108];
         bas[14] = -rad[ip+ng]*angArray[49]+rad[ip+2*ng]*angArray[109];
         bas[15] = rad[ip]*20*angArray[16]-rad[ip+ng]*11*angArray[50]+rad[ip+2*ng]*angArray[112];
         bas[16] = rad[ip]*12*angArray[17]-rad[ip+ng]*9*angArray[51]+rad[ip+2*ng]*angArray[113];
         bas[17] = rad[ip]*6*angArray[18]-rad[ip+ng]*7*angArray[52]+rad[ip+2*ng]*angArray[114];
         bas[18] = rad[ip]*2*angArray[19]-rad[ip+ng]*5*angArray[53]+rad[ip+2*ng]*angArray[115];
         bas[19] = -rad[ip+ng]*3*angArray[54]+rad[ip+2*ng]*angArray[116];
         bas[20] = -rad[ip+ng]*angArray[55]+rad[ip+2*ng]*angArray[117];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[56]+rad[ip+2*ng]*angArray[123];
         bas[1] = -rad[ip+ng]*3*angArray[57]+rad[ip+2*ng]*angArray[126];
         bas[2] = -rad[ip+ng]*angArray[58]+rad[ip+2*ng]*angArray[127];
         bas[3] = rad[ip]*2*angArray[20]-rad[ip+ng]*5*angArray[59]+rad[ip+2*ng]*angArray[130];
         bas[4] = -rad[ip+ng]*3*angArray[60]+rad[ip+2*ng]*angArray[131];
         bas[5] = -rad[ip+ng]*angArray[61]+rad[ip+2*ng]*angArray[132];
         bas[6] = rad[ip]*6*angArray[21]-rad[ip+ng]*7*angArray[62]+rad[ip+2*ng]*angArray[135];
         bas[7] = rad[ip]*2*angArray[22]-rad[ip+ng]*5*angArray[63]+rad[ip+2*ng]*angArray[136];
         bas[8] = -rad[ip+ng]*3*angArray[64]+rad[ip+2*ng]*angArray[137];
         bas[9] = -rad[ip+ng]*angArray[65]+rad[ip+2*ng]*angArray[138];
         bas[10] = rad[ip]*12*angArray[23]-rad[ip+ng]*9*angArray[66]+rad[ip+2*ng]*angArray[141];
         bas[11] = rad[ip]*6*angArray[24]-rad[ip+ng]*7*angArray[67]+rad[ip+2*ng]*angArray[142];
         bas[12] = rad[ip]*2*angArray[25]-rad[ip+ng]*5*angArray[68]+rad[ip+2*ng]*angArray[143];
         bas[13] = -rad[ip+ng]*3*angArray[69]+rad[ip+2*ng]*angArray[144];
         bas[14] = -rad[ip+ng]*angArray[70]+rad[ip+2*ng]*angArray[145];
         bas[15] = rad[ip]*20*angArray[26]-rad[ip+ng]*11*angArray[71]+rad[ip+2*ng]*angArray[148];
         bas[16] = rad[ip]*12*angArray[27]-rad[ip+ng]*9*angArray[72]+rad[ip+2*ng]*angArray[149];
         bas[17] = rad[ip]*6*angArray[28]-rad[ip+ng]*7*angArray[73]+rad[ip+2*ng]*angArray[150];
         bas[18] = rad[ip]*2*angArray[29]-rad[ip+ng]*5*angArray[74]+rad[ip+2*ng]*angArray[151];
         bas[19] = -rad[ip+ng]*3*angArray[75]+rad[ip+2*ng]*angArray[152];
         bas[20] = -rad[ip+ng]*angArray[76]+rad[ip+2*ng]*angArray[153];
         bas[21] = rad[ip]*30*angArray[30]-rad[ip+ng]*13*angArray[77]+rad[ip+2*ng]*angArray[156];
         bas[22] = rad[ip]*20*angArray[31]-rad[ip+ng]*11*angArray[78]+rad[ip+2*ng]*angArray[157];
         bas[23] = rad[ip]*12*angArray[32]-rad[ip+ng]*9*angArray[79]+rad[ip+2*ng]*angArray[158];
         bas[24] = rad[ip]*6*angArray[33]-rad[ip+ng]*7*angArray[80]+rad[ip+2*ng]*angArray[159];
         bas[25] = rad[ip]*2*angArray[34]-rad[ip+ng]*5*angArray[81]+rad[ip+2*ng]*angArray[160];
         bas[26] = -rad[ip+ng]*3*angArray[82]+rad[ip+2*ng]*angArray[161];
         bas[27] = -rad[ip+ng]*angArray[83]+rad[ip+2*ng]*angArray[162];
      }

   }


   // now we do derivatives for the given basis set to XZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[6];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[3]+rad[ip+2*ng]*angArray[12];
         bas[1] = rad[ip+2*ng]*angArray[14];
         bas[2] = -rad[ip+ng]*angArray[1]+rad[ip+2*ng]*angArray[15];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*angArray[22];
         bas[1] = -rad[ip+ng]*angArray[8]+rad[ip+2*ng]*angArray[24];
         bas[2] = rad[ip]*angArray[0]-rad[ip+ng]*(angArray[4]+angArray[9])+rad[ip+2*ng]*angArray[25];
         bas[3] = rad[ip+2*ng]*angArray[27];
         bas[4] = -rad[ip+ng]*angArray[5]+rad[ip+2*ng]*angArray[28];
         bas[5] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*angArray[29];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*3*angArray[12]+rad[ip+2*ng]*angArray[37];
         bas[1] = -rad[ip+ng]*2*angArray[14]+rad[ip+2*ng]*angArray[39];
         bas[2] = rad[ip]*2*angArray[1]-rad[ip+ng]*(2*angArray[15]+angArray[10])+rad[ip+2*ng]*angArray[40];
         bas[3] = -rad[ip+ng]*angArray[17]+rad[ip+2*ng]*angArray[42];
         bas[4] = rad[ip]*angArray[2]-rad[ip+ng]*(angArray[11]+angArray[18])+rad[ip+2*ng]*angArray[43];
         bas[5] = rad[ip]*2*angArray[3]-rad[ip+ng]*(2*angArray[12]+angArray[19])+rad[ip+2*ng]*angArray[44];
         bas[6] = rad[ip+2*ng]*angArray[46];
         bas[7] = -rad[ip+ng]*angArray[13]+rad[ip+2*ng]*angArray[47];
         bas[8] = -rad[ip+ng]*2*angArray[14]+rad[ip+2*ng]*angArray[48];
         bas[9] = -rad[ip+ng]*3*angArray[15]+rad[ip+2*ng]*angArray[49];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*4*angArray[22]+rad[ip+2*ng]*angArray[58];
         bas[1] = -rad[ip+ng]*3*angArray[24]+rad[ip+2*ng]*angArray[60];
         bas[2] = rad[ip]*3*angArray[4]-rad[ip+ng]*(angArray[20]+3*angArray[25])+rad[ip+2*ng]*angArray[61];
         bas[3] = -rad[ip+ng]*2*angArray[27]+rad[ip+2*ng]*angArray[63];
         bas[4] = rad[ip]*2*angArray[5]-rad[ip+ng]*(2*angArray[28]+angArray[21])+rad[ip+2*ng]*angArray[64];
         bas[5] = rad[ip]*4*angArray[6]-rad[ip+ng]*(2*angArray[29]+2*angArray[22])+rad[ip+2*ng]*angArray[65];
         bas[6] = -rad[ip+ng]*angArray[31]+rad[ip+2*ng]*angArray[67];
         bas[7] = rad[ip]*angArray[7]-rad[ip+ng]*(angArray[23]+angArray[32])+rad[ip+2*ng]*angArray[68];
         bas[8] = rad[ip]*2*angArray[8]-rad[ip+ng]*(2*angArray[24]+angArray[33])+rad[ip+2*ng]*angArray[69];
         bas[9] = rad[ip]*3*angArray[9]-rad[ip+ng]*(3*angArray[25]+angArray[34])+rad[ip+2*ng]*angArray[70];
         bas[10] = rad[ip+2*ng]*angArray[72];
         bas[11] = -rad[ip+ng]*angArray[26]+rad[ip+2*ng]*angArray[73];
         bas[12] = -rad[ip+ng]*2*angArray[27]+rad[ip+2*ng]*angArray[74];
         bas[13] = -rad[ip+ng]*3*angArray[28]+rad[ip+2*ng]*angArray[75];
         bas[14] = -rad[ip+ng]*4*angArray[29]+rad[ip+2*ng]*angArray[76];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*5*angArray[37]+rad[ip+2*ng]*angArray[86];
         bas[1] = -rad[ip+ng]*4*angArray[39]+rad[ip+2*ng]*angArray[88];
         bas[2] = rad[ip]*4*angArray[10]-rad[ip+ng]*(4*angArray[40]+angArray[35])+rad[ip+2*ng]*angArray[89];
         bas[3] = -rad[ip+ng]*3*angArray[42]+rad[ip+2*ng]*angArray[91];
         bas[4] = rad[ip]*3*angArray[11]-rad[ip+ng]*(angArray[36]+3*angArray[43])+rad[ip+2*ng]*angArray[92];
         bas[5] = rad[ip]*6*angArray[12]-rad[ip+ng]*(2*angArray[37]+3*angArray[44])+rad[ip+2*ng]*angArray[93];
         bas[6] = -rad[ip+ng]*2*angArray[46]+rad[ip+2*ng]*angArray[95];
         bas[7] = rad[ip]*2*angArray[13]-rad[ip+ng]*(2*angArray[47]+angArray[38])+rad[ip+2*ng]*angArray[96];
         bas[8] = rad[ip]*4*angArray[14]-rad[ip+ng]*(2*angArray[48]+2*angArray[39])+rad[ip+2*ng]*angArray[97];
         bas[9] = rad[ip]*6*angArray[15]-rad[ip+ng]*(2*angArray[49]+3*angArray[40])+rad[ip+2*ng]*angArray[98];
         bas[10] = -rad[ip+ng]*angArray[51]+rad[ip+2*ng]*angArray[100];
         bas[11] = rad[ip]*angArray[16]-rad[ip+ng]*(angArray[41]+angArray[52])+rad[ip+2*ng]*angArray[101];
         bas[12] = rad[ip]*2*angArray[17]-rad[ip+ng]*(2*angArray[42]+angArray[53])+rad[ip+2*ng]*angArray[102];
         bas[13] = rad[ip]*3*angArray[18]-rad[ip+ng]*(3*angArray[43]+angArray[54])+rad[ip+2*ng]*angArray[103];
         bas[14] = rad[ip]*4*angArray[19]-rad[ip+ng]*(4*angArray[44]+angArray[55])+rad[ip+2*ng]*angArray[104];
         bas[15] = rad[ip+2*ng]*angArray[106];
         bas[16] = -rad[ip+ng]*angArray[45]+rad[ip+2*ng]*angArray[107];
         bas[17] = -rad[ip+ng]*2*angArray[46]+rad[ip+2*ng]*angArray[108];
         bas[18] = -rad[ip+ng]*3*angArray[47]+rad[ip+2*ng]*angArray[109];
         bas[19] = -rad[ip+ng]*4*angArray[48]+rad[ip+2*ng]*angArray[110];
         bas[20] = -rad[ip+ng]*5*angArray[49]+rad[ip+2*ng]*angArray[111];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[58]+rad[ip+2*ng]*angArray[122];
         bas[1] = -rad[ip+ng]*5*angArray[60]+rad[ip+2*ng]*angArray[124];
         bas[2] = rad[ip]*5*angArray[20]-rad[ip+ng]*(5*angArray[61]+angArray[56])+rad[ip+2*ng]*angArray[125];
         bas[3] = -rad[ip+ng]*4*angArray[63]+rad[ip+2*ng]*angArray[127];
         bas[4] = rad[ip]*4*angArray[21]-rad[ip+ng]*(4*angArray[64]+angArray[57])+rad[ip+2*ng]*angArray[128];
         bas[5] = rad[ip]*8*angArray[22]-rad[ip+ng]*(2*angArray[58]+4*angArray[65])+rad[ip+2*ng]*angArray[129];
         bas[6] = -rad[ip+ng]*3*angArray[67]+rad[ip+2*ng]*angArray[131];
         bas[7] = rad[ip]*3*angArray[23]-rad[ip+ng]*(angArray[59]+3*angArray[68])+rad[ip+2*ng]*angArray[132];
         bas[8] = rad[ip]*6*angArray[24]-rad[ip+ng]*(2*angArray[60]+3*angArray[69])+rad[ip+2*ng]*angArray[133];
         bas[9] = rad[ip]*9*angArray[25]-rad[ip+ng]*(3*angArray[61]+3*angArray[70])+rad[ip+2*ng]*angArray[134];
         bas[10] = -rad[ip+ng]*2*angArray[72]+rad[ip+2*ng]*angArray[136];
         bas[11] = rad[ip]*2*angArray[26]-rad[ip+ng]*(2*angArray[73]+angArray[62])+rad[ip+2*ng]*angArray[137];
         bas[12] = rad[ip]*4*angArray[27]-rad[ip+ng]*(2*angArray[74]+2*angArray[63])+rad[ip+2*ng]*angArray[138];
         bas[13] = rad[ip]*6*angArray[28]-rad[ip+ng]*(2*angArray[75]+3*angArray[64])+rad[ip+2*ng]*angArray[139];
         bas[14] = rad[ip]*8*angArray[29]-rad[ip+ng]*(4*angArray[65]+2*angArray[76])+rad[ip+2*ng]*angArray[140];
         bas[15] = -rad[ip+ng]*angArray[78]+rad[ip+2*ng]*angArray[142];
         bas[16] = rad[ip]*angArray[30]-rad[ip+ng]*(angArray[66]+angArray[79])+rad[ip+2*ng]*angArray[143];
         bas[17] = rad[ip]*2*angArray[31]-rad[ip+ng]*(2*angArray[67]+angArray[80])+rad[ip+2*ng]*angArray[144];
         bas[18] = rad[ip]*3*angArray[32]-rad[ip+ng]*(3*angArray[68]+angArray[81])+rad[ip+2*ng]*angArray[145];
         bas[19] = rad[ip]*4*angArray[33]-rad[ip+ng]*(4*angArray[69]+angArray[82])+rad[ip+2*ng]*angArray[146];
         bas[20] = rad[ip]*5*angArray[34]-rad[ip+ng]*(5*angArray[70]+angArray[83])+rad[ip+2*ng]*angArray[147];
         bas[21] = rad[ip+2*ng]*angArray[149];
         bas[22] = -rad[ip+ng]*angArray[71]+rad[ip+2*ng]*angArray[150];
         bas[23] = -rad[ip+ng]*2*angArray[72]+rad[ip+2*ng]*angArray[151];
         bas[24] = -rad[ip+ng]*3*angArray[73]+rad[ip+2*ng]*angArray[152];
         bas[25] = -rad[ip+ng]*4*angArray[74]+rad[ip+2*ng]*angArray[153];
         bas[26] = -rad[ip+ng]*5*angArray[75]+rad[ip+2*ng]*angArray[154];
         bas[27] = -rad[ip+ng]*6*angArray[76]+rad[ip+2*ng]*angArray[155];
      }

   }


   // now we do derivatives for the given basis set to YZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[8];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[14];
         bas[1] = -rad[ip+ng]*angArray[3]+rad[ip+2*ng]*angArray[17];
         bas[2] = -rad[ip+ng]*angArray[2]+rad[ip+2*ng]*angArray[18];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[24];
         bas[1] = -rad[ip+ng]*angArray[6]+rad[ip+2*ng]*angArray[27];
         bas[2] = -rad[ip+ng]*angArray[5]+rad[ip+2*ng]*angArray[28];
         bas[3] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*angArray[31];
         bas[4] = rad[ip]*angArray[0]-rad[ip+ng]*(angArray[7]+angArray[9])+rad[ip+2*ng]*angArray[32];
         bas[5] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*angArray[33];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[39];
         bas[1] = -rad[ip+ng]*angArray[12]+rad[ip+2*ng]*angArray[42];
         bas[2] = -rad[ip+ng]*angArray[11]+rad[ip+2*ng]*angArray[43];
         bas[3] = -rad[ip+ng]*2*angArray[14]+rad[ip+2*ng]*angArray[46];
         bas[4] = rad[ip]*angArray[1]-rad[ip+ng]*(angArray[13]+angArray[15])+rad[ip+2*ng]*angArray[47];
         bas[5] = -rad[ip+ng]*2*angArray[14]+rad[ip+2*ng]*angArray[48];
         bas[6] = -rad[ip+ng]*3*angArray[17]+rad[ip+2*ng]*angArray[51];
         bas[7] = rad[ip]*2*angArray[2]-rad[ip+ng]*(2*angArray[18]+angArray[16])+rad[ip+2*ng]*angArray[52];
         bas[8] = rad[ip]*2*angArray[3]-rad[ip+ng]*(2*angArray[17]+angArray[19])+rad[ip+2*ng]*angArray[53];
         bas[9] = -rad[ip+ng]*3*angArray[18]+rad[ip+2*ng]*angArray[54];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[60];
         bas[1] = -rad[ip+ng]*angArray[22]+rad[ip+2*ng]*angArray[63];
         bas[2] = -rad[ip+ng]*angArray[21]+rad[ip+2*ng]*angArray[64];
         bas[3] = -rad[ip+ng]*2*angArray[24]+rad[ip+2*ng]*angArray[67];
         bas[4] = rad[ip]*angArray[4]-rad[ip+ng]*(angArray[23]+angArray[25])+rad[ip+2*ng]*angArray[68];
         bas[5] = -rad[ip+ng]*2*angArray[24]+rad[ip+2*ng]*angArray[69];
         bas[6] = -rad[ip+ng]*3*angArray[27]+rad[ip+2*ng]*angArray[72];
         bas[7] = rad[ip]*2*angArray[5]-rad[ip+ng]*(2*angArray[28]+angArray[26])+rad[ip+2*ng]*angArray[73];
         bas[8] = rad[ip]*2*angArray[6]-rad[ip+ng]*(2*angArray[27]+angArray[29])+rad[ip+2*ng]*angArray[74];
         bas[9] = -rad[ip+ng]*3*angArray[28]+rad[ip+2*ng]*angArray[75];
         bas[10] = -rad[ip+ng]*4*angArray[31]+rad[ip+2*ng]*angArray[78];
         bas[11] = rad[ip]*3*angArray[7]-rad[ip+ng]*(angArray[30]+3*angArray[32])+rad[ip+2*ng]*angArray[79];
         bas[12] = rad[ip]*4*angArray[8]-rad[ip+ng]*(2*angArray[33]+2*angArray[31])+rad[ip+2*ng]*angArray[80];
         bas[13] = rad[ip]*3*angArray[9]-rad[ip+ng]*(3*angArray[32]+angArray[34])+rad[ip+2*ng]*angArray[81];
         bas[14] = -rad[ip+ng]*4*angArray[33]+rad[ip+2*ng]*angArray[82];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[88];
         bas[1] = -rad[ip+ng]*angArray[37]+rad[ip+2*ng]*angArray[91];
         bas[2] = -rad[ip+ng]*angArray[36]+rad[ip+2*ng]*angArray[92];
         bas[3] = -rad[ip+ng]*2*angArray[39]+rad[ip+2*ng]*angArray[95];
         bas[4] = rad[ip]*angArray[10]-rad[ip+ng]*(angArray[38]+angArray[40])+rad[ip+2*ng]*angArray[96];
         bas[5] = -rad[ip+ng]*2*angArray[39]+rad[ip+2*ng]*angArray[97];
         bas[6] = -rad[ip+ng]*3*angArray[42]+rad[ip+2*ng]*angArray[100];
         bas[7] = rad[ip]*2*angArray[11]-rad[ip+ng]*(2*angArray[43]+angArray[41])+rad[ip+2*ng]*angArray[101];
         bas[8] = rad[ip]*2*angArray[12]-rad[ip+ng]*(2*angArray[42]+angArray[44])+rad[ip+2*ng]*angArray[102];
         bas[9] = -rad[ip+ng]*3*angArray[43]+rad[ip+2*ng]*angArray[103];
         bas[10] = -rad[ip+ng]*4*angArray[46]+rad[ip+2*ng]*angArray[106];
         bas[11] = rad[ip]*3*angArray[13]-rad[ip+ng]*(angArray[45]+3*angArray[47])+rad[ip+2*ng]*angArray[107];
         bas[12] = rad[ip]*4*angArray[14]-rad[ip+ng]*(2*angArray[48]+2*angArray[46])+rad[ip+2*ng]*angArray[108];
         bas[13] = rad[ip]*3*angArray[15]-rad[ip+ng]*(3*angArray[47]+angArray[49])+rad[ip+2*ng]*angArray[109];
         bas[14] = -rad[ip+ng]*4*angArray[48]+rad[ip+2*ng]*angArray[110];
         bas[15] = -rad[ip+ng]*5*angArray[51]+rad[ip+2*ng]*angArray[113];
         bas[16] = rad[ip]*4*angArray[16]-rad[ip+ng]*(4*angArray[52]+angArray[50])+rad[ip+2*ng]*angArray[114];
         bas[17] = rad[ip]*6*angArray[17]-rad[ip+ng]*(2*angArray[51]+3*angArray[53])+rad[ip+2*ng]*angArray[115];
         bas[18] = rad[ip]*6*angArray[18]-rad[ip+ng]*(2*angArray[54]+3*angArray[52])+rad[ip+2*ng]*angArray[116];
         bas[19] = rad[ip]*4*angArray[19]-rad[ip+ng]*(4*angArray[53]+angArray[55])+rad[ip+2*ng]*angArray[117];
         bas[20] = -rad[ip+ng]*5*angArray[54]+rad[ip+2*ng]*angArray[118];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[124];
         bas[1] = -rad[ip+ng]*angArray[58]+rad[ip+2*ng]*angArray[127];
         bas[2] = -rad[ip+ng]*angArray[57]+rad[ip+2*ng]*angArray[128];
         bas[3] = -rad[ip+ng]*2*angArray[60]+rad[ip+2*ng]*angArray[131];
         bas[4] = rad[ip]*angArray[20]-rad[ip+ng]*(angArray[59]+angArray[61])+rad[ip+2*ng]*angArray[132];
         bas[5] = -rad[ip+ng]*2*angArray[60]+rad[ip+2*ng]*angArray[133];
         bas[6] = -rad[ip+ng]*3*angArray[63]+rad[ip+2*ng]*angArray[136];
         bas[7] = rad[ip]*2*angArray[21]-rad[ip+ng]*(2*angArray[64]+angArray[62])+rad[ip+2*ng]*angArray[137];
         bas[8] = rad[ip]*2*angArray[22]-rad[ip+ng]*(2*angArray[63]+angArray[65])+rad[ip+2*ng]*angArray[138];
         bas[9] = -rad[ip+ng]*3*angArray[64]+rad[ip+2*ng]*angArray[139];
         bas[10] = -rad[ip+ng]*4*angArray[67]+rad[ip+2*ng]*angArray[142];
         bas[11] = rad[ip]*3*angArray[23]-rad[ip+ng]*(angArray[66]+3*angArray[68])+rad[ip+2*ng]*angArray[143];
         bas[12] = rad[ip]*4*angArray[24]-rad[ip+ng]*(2*angArray[69]+2*angArray[67])+rad[ip+2*ng]*angArray[144];
         bas[13] = rad[ip]*3*angArray[25]-rad[ip+ng]*(3*angArray[68]+angArray[70])+rad[ip+2*ng]*angArray[145];
         bas[14] = -rad[ip+ng]*4*angArray[69]+rad[ip+2*ng]*angArray[146];
         bas[15] = -rad[ip+ng]*5*angArray[72]+rad[ip+2*ng]*angArray[149];
         bas[16] = rad[ip]*4*angArray[26]-rad[ip+ng]*(4*angArray[73]+angArray[71])+rad[ip+2*ng]*angArray[150];
         bas[17] = rad[ip]*6*angArray[27]-rad[ip+ng]*(2*angArray[72]+3*angArray[74])+rad[ip+2*ng]*angArray[151];
         bas[18] = rad[ip]*6*angArray[28]-rad[ip+ng]*(2*angArray[75]+3*angArray[73])+rad[ip+2*ng]*angArray[152];
         bas[19] = rad[ip]*4*angArray[29]-rad[ip+ng]*(4*angArray[74]+angArray[76])+rad[ip+2*ng]*angArray[153];
         bas[20] = -rad[ip+ng]*5*angArray[75]+rad[ip+2*ng]*angArray[154];
         bas[21] = -rad[ip+ng]*6*angArray[78]+rad[ip+2*ng]*angArray[157];
         bas[22] = rad[ip]*5*angArray[30]-rad[ip+ng]*(5*angArray[79]+angArray[77])+rad[ip+2*ng]*angArray[158];
         bas[23] = rad[ip]*8*angArray[31]-rad[ip+ng]*(2*angArray[78]+4*angArray[80])+rad[ip+2*ng]*angArray[159];
         bas[24] = rad[ip]*9*angArray[32]-rad[ip+ng]*(3*angArray[79]+3*angArray[81])+rad[ip+2*ng]*angArray[160];
         bas[25] = rad[ip]*8*angArray[33]-rad[ip+ng]*(4*angArray[80]+2*angArray[82])+rad[ip+2*ng]*angArray[161];
         bas[26] = rad[ip]*5*angArray[34]-rad[ip+ng]*(5*angArray[81]+angArray[83])+rad[ip+2*ng]*angArray[162];
         bas[27] = -rad[ip+ng]*6*angArray[82]+rad[ip+2*ng]*angArray[163];
      }

   }


   // now we do derivatives for the given basis set to ZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*angArray[9];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[1]+rad[ip+2*ng]*angArray[15];
         bas[1] = -rad[ip+ng]*angArray[2]+rad[ip+2*ng]*angArray[18];
         bas[2] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*angArray[19];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[4]+rad[ip+2*ng]*angArray[25];
         bas[1] = -rad[ip+ng]*angArray[5]+rad[ip+2*ng]*angArray[28];
         bas[2] = -rad[ip+ng]*3*angArray[6]+rad[ip+2*ng]*angArray[29];
         bas[3] = -rad[ip+ng]*angArray[7]+rad[ip+2*ng]*angArray[32];
         bas[4] = -rad[ip+ng]*3*angArray[8]+rad[ip+2*ng]*angArray[33];
         bas[5] = rad[ip]*2*angArray[0]-rad[ip+ng]*5*angArray[9]+rad[ip+2*ng]*angArray[34];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[10]+rad[ip+2*ng]*angArray[40];
         bas[1] = -rad[ip+ng]*angArray[11]+rad[ip+2*ng]*angArray[43];
         bas[2] = -rad[ip+ng]*3*angArray[12]+rad[ip+2*ng]*angArray[44];
         bas[3] = -rad[ip+ng]*angArray[13]+rad[ip+2*ng]*angArray[47];
         bas[4] = -rad[ip+ng]*3*angArray[14]+rad[ip+2*ng]*angArray[48];
         bas[5] = rad[ip]*2*angArray[1]-rad[ip+ng]*5*angArray[15]+rad[ip+2*ng]*angArray[49];
         bas[6] = -rad[ip+ng]*angArray[16]+rad[ip+2*ng]*angArray[52];
         bas[7] = -rad[ip+ng]*3*angArray[17]+rad[ip+2*ng]*angArray[53];
         bas[8] = rad[ip]*2*angArray[2]-rad[ip+ng]*5*angArray[18]+rad[ip+2*ng]*angArray[54];
         bas[9] = rad[ip]*6*angArray[3]-rad[ip+ng]*7*angArray[19]+rad[ip+2*ng]*angArray[55];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[20]+rad[ip+2*ng]*angArray[61];
         bas[1] = -rad[ip+ng]*angArray[21]+rad[ip+2*ng]*angArray[64];
         bas[2] = -rad[ip+ng]*3*angArray[22]+rad[ip+2*ng]*angArray[65];
         bas[3] = -rad[ip+ng]*angArray[23]+rad[ip+2*ng]*angArray[68];
         bas[4] = -rad[ip+ng]*3*angArray[24]+rad[ip+2*ng]*angArray[69];
         bas[5] = rad[ip]*2*angArray[4]-rad[ip+ng]*5*angArray[25]+rad[ip+2*ng]*angArray[70];
         bas[6] = -rad[ip+ng]*angArray[26]+rad[ip+2*ng]*angArray[73];
         bas[7] = -rad[ip+ng]*3*angArray[27]+rad[ip+2*ng]*angArray[74];
         bas[8] = rad[ip]*2*angArray[5]-rad[ip+ng]*5*angArray[28]+rad[ip+2*ng]*angArray[75];
         bas[9] = rad[ip]*6*angArray[6]-rad[ip+ng]*7*angArray[29]+rad[ip+2*ng]*angArray[76];
         bas[10] = -rad[ip+ng]*angArray[30]+rad[ip+2*ng]*angArray[79];
         bas[11] = -rad[ip+ng]*3*angArray[31]+rad[ip+2*ng]*angArray[80];
         bas[12] = rad[ip]*2*angArray[7]-rad[ip+ng]*5*angArray[32]+rad[ip+2*ng]*angArray[81];
         bas[13] = rad[ip]*6*angArray[8]-rad[ip+ng]*7*angArray[33]+rad[ip+2*ng]*angArray[82];
         bas[14] = rad[ip]*12*angArray[9]-rad[ip+ng]*9*angArray[34]+rad[ip+2*ng]*angArray[83];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[35]+rad[ip+2*ng]*angArray[89];
         bas[1] = -rad[ip+ng]*angArray[36]+rad[ip+2*ng]*angArray[92];
         bas[2] = -rad[ip+ng]*3*angArray[37]+rad[ip+2*ng]*angArray[93];
         bas[3] = -rad[ip+ng]*angArray[38]+rad[ip+2*ng]*angArray[96];
         bas[4] = -rad[ip+ng]*3*angArray[39]+rad[ip+2*ng]*angArray[97];
         bas[5] = rad[ip]*2*angArray[10]-rad[ip+ng]*5*angArray[40]+rad[ip+2*ng]*angArray[98];
         bas[6] = -rad[ip+ng]*angArray[41]+rad[ip+2*ng]*angArray[101];
         bas[7] = -rad[ip+ng]*3*angArray[42]+rad[ip+2*ng]*angArray[102];
         bas[8] = rad[ip]*2*angArray[11]-rad[ip+ng]*5*angArray[43]+rad[ip+2*ng]*angArray[103];
         bas[9] = rad[ip]*6*angArray[12]-rad[ip+ng]*7*angArray[44]+rad[ip+2*ng]*angArray[104];
         bas[10] = -rad[ip+ng]*angArray[45]+rad[ip+2*ng]*angArray[107];
         bas[11] = -rad[ip+ng]*3*angArray[46]+rad[ip+2*ng]*angArray[108];
         bas[12] = rad[ip]*2*angArray[13]-rad[ip+ng]*5*angArray[47]+rad[ip+2*ng]*angArray[109];
         bas[13] = rad[ip]*6*angArray[14]-rad[ip+ng]*7*angArray[48]+rad[ip+2*ng]*angArray[110];
         bas[14] = rad[ip]*12*angArray[15]-rad[ip+ng]*9*angArray[49]+rad[ip+2*ng]*angArray[111];
         bas[15] = -rad[ip+ng]*angArray[50]+rad[ip+2*ng]*angArray[114];
         bas[16] = -rad[ip+ng]*3*angArray[51]+rad[ip+2*ng]*angArray[115];
         bas[17] = rad[ip]*2*angArray[16]-rad[ip+ng]*5*angArray[52]+rad[ip+2*ng]*angArray[116];
         bas[18] = rad[ip]*6*angArray[17]-rad[ip+ng]*7*angArray[53]+rad[ip+2*ng]*angArray[117];
         bas[19] = rad[ip]*12*angArray[18]-rad[ip+ng]*9*angArray[54]+rad[ip+2*ng]*angArray[118];
         bas[20] = rad[ip]*20*angArray[19]-rad[ip+ng]*11*angArray[55]+rad[ip+2*ng]*angArray[119];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[56]+rad[ip+2*ng]*angArray[125];
         bas[1] = -rad[ip+ng]*angArray[57]+rad[ip+2*ng]*angArray[128];
         bas[2] = -rad[ip+ng]*3*angArray[58]+rad[ip+2*ng]*angArray[129];
         bas[3] = -rad[ip+ng]*angArray[59]+rad[ip+2*ng]*angArray[132];
         bas[4] = -rad[ip+ng]*3*angArray[60]+rad[ip+2*ng]*angArray[133];
         bas[5] = rad[ip]*2*angArray[20]-rad[ip+ng]*5*angArray[61]+rad[ip+2*ng]*angArray[134];
         bas[6] = -rad[ip+ng]*angArray[62]+rad[ip+2*ng]*angArray[137];
         bas[7] = -rad[ip+ng]*3*angArray[63]+rad[ip+2*ng]*angArray[138];
         bas[8] = rad[ip]*2*angArray[21]-rad[ip+ng]*5*angArray[64]+rad[ip+2*ng]*angArray[139];
         bas[9] = rad[ip]*6*angArray[22]-rad[ip+ng]*7*angArray[65]+rad[ip+2*ng]*angArray[140];
         bas[10] = -rad[ip+ng]*angArray[66]+rad[ip+2*ng]*angArray[143];
         bas[11] = -rad[ip+ng]*3*angArray[67]+rad[ip+2*ng]*angArray[144];
         bas[12] = rad[ip]*2*angArray[23]-rad[ip+ng]*5*angArray[68]+rad[ip+2*ng]*angArray[145];
         bas[13] = rad[ip]*6*angArray[24]-rad[ip+ng]*7*angArray[69]+rad[ip+2*ng]*angArray[146];
         bas[14] = rad[ip]*12*angArray[25]-rad[ip+ng]*9*angArray[70]+rad[ip+2*ng]*angArray[147];
         bas[15] = -rad[ip+ng]*angArray[71]+rad[ip+2*ng]*angArray[150];
         bas[16] = -rad[ip+ng]*3*angArray[72]+rad[ip+2*ng]*angArray[151];
         bas[17] = rad[ip]*2*angArray[26]-rad[ip+ng]*5*angArray[73]+rad[ip+2*ng]*angArray[152];
         bas[18] = rad[ip]*6*angArray[27]-rad[ip+ng]*7*angArray[74]+rad[ip+2*ng]*angArray[153];
         bas[19] = rad[ip]*12*angArray[28]-rad[ip+ng]*9*angArray[75]+rad[ip+2*ng]*angArray[154];
         bas[20] = rad[ip]*20*angArray[29]-rad[ip+ng]*11*angArray[76]+rad[ip+2*ng]*angArray[155];
         bas[21] = -rad[ip+ng]*angArray[77]+rad[ip+2*ng]*angArray[158];
         bas[22] = -rad[ip+ng]*3*angArray[78]+rad[ip+2*ng]*angArray[159];
         bas[23] = rad[ip]*2*angArray[30]-rad[ip+ng]*5*angArray[79]+rad[ip+2*ng]*angArray[160];
         bas[24] = rad[ip]*6*angArray[31]-rad[ip+ng]*7*angArray[80]+rad[ip+2*ng]*angArray[161];
         bas[25] = rad[ip]*12*angArray[32]-rad[ip+ng]*9*angArray[81]+rad[ip+2*ng]*angArray[162];
         bas[26] = rad[ip]*20*angArray[33]-rad[ip+ng]*11*angArray[82]+rad[ip+2*ng]*angArray[163];
         bas[27] = rad[ip]*30*angArray[34]-rad[ip+ng]*13*angArray[83]+rad[ip+2*ng]*angArray[164];
      }

   }


}


