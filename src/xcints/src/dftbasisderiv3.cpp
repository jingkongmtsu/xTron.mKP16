/**
 * This function is used to generate 3 derivatives for basis set 
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

void BatchBasis::dftbasisderiv3(const UInt& ng, const UInt& L, const UInt& nTolCarBas, const Double* ang, const Double* rad, Double* basis) const 
{

   // now we set up the nBas for the computation
   UInt nBas = (L+1)*(L+2)/2;

   // now we do derivatives for the given basis set to XXX
   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*angArray[10];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*6*angArray[4]-rad[ip+3*ng]*angArray[20];
         bas[1] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*angArray[21];
         bas[2] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*angArray[22];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*12*angArray[1]+rad[ip+2*ng]*9*angArray[10]-rad[ip+3*ng]*angArray[35];
         bas[1] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*6*angArray[11]-rad[ip+3*ng]*angArray[36];
         bas[2] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*6*angArray[12]-rad[ip+3*ng]*angArray[37];
         bas[3] = rad[ip+2*ng]*3*angArray[13]-rad[ip+3*ng]*angArray[38];
         bas[4] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[39];
         bas[5] = rad[ip+2*ng]*3*angArray[15]-rad[ip+3*ng]*angArray[40];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*6*angArray[0]-rad[ip+ng]*27*angArray[4]+rad[ip+2*ng]*12*angArray[20]-rad[ip+3*ng]*angArray[56];
         bas[1] = -rad[ip+ng]*12*angArray[5]+rad[ip+2*ng]*9*angArray[21]-rad[ip+3*ng]*angArray[57];
         bas[2] = -rad[ip+ng]*12*angArray[6]+rad[ip+2*ng]*9*angArray[22]-rad[ip+3*ng]*angArray[58];
         bas[3] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*6*angArray[23]-rad[ip+3*ng]*angArray[59];
         bas[4] = -rad[ip+ng]*3*angArray[8]+rad[ip+2*ng]*6*angArray[24]-rad[ip+3*ng]*angArray[60];
         bas[5] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*6*angArray[25]-rad[ip+3*ng]*angArray[61];
         bas[6] = rad[ip+2*ng]*3*angArray[26]-rad[ip+3*ng]*angArray[62];
         bas[7] = rad[ip+2*ng]*3*angArray[27]-rad[ip+3*ng]*angArray[63];
         bas[8] = rad[ip+2*ng]*3*angArray[28]-rad[ip+3*ng]*angArray[64];
         bas[9] = rad[ip+2*ng]*3*angArray[29]-rad[ip+3*ng]*angArray[65];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*24*angArray[1]-rad[ip+ng]*48*angArray[10]+rad[ip+2*ng]*15*angArray[35]-rad[ip+3*ng]*angArray[84];
         bas[1] = rad[ip]*6*angArray[2]-rad[ip+ng]*27*angArray[11]+rad[ip+2*ng]*12*angArray[36]-rad[ip+3*ng]*angArray[85];
         bas[2] = rad[ip]*6*angArray[3]-rad[ip+ng]*27*angArray[12]+rad[ip+2*ng]*12*angArray[37]-rad[ip+3*ng]*angArray[86];
         bas[3] = -rad[ip+ng]*12*angArray[13]+rad[ip+2*ng]*9*angArray[38]-rad[ip+3*ng]*angArray[87];
         bas[4] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*9*angArray[39]-rad[ip+3*ng]*angArray[88];
         bas[5] = -rad[ip+ng]*12*angArray[15]+rad[ip+2*ng]*9*angArray[40]-rad[ip+3*ng]*angArray[89];
         bas[6] = -rad[ip+ng]*3*angArray[16]+rad[ip+2*ng]*6*angArray[41]-rad[ip+3*ng]*angArray[90];
         bas[7] = -rad[ip+ng]*3*angArray[17]+rad[ip+2*ng]*6*angArray[42]-rad[ip+3*ng]*angArray[91];
         bas[8] = -rad[ip+ng]*3*angArray[18]+rad[ip+2*ng]*6*angArray[43]-rad[ip+3*ng]*angArray[92];
         bas[9] = -rad[ip+ng]*3*angArray[19]+rad[ip+2*ng]*6*angArray[44]-rad[ip+3*ng]*angArray[93];
         bas[10] = rad[ip+2*ng]*3*angArray[45]-rad[ip+3*ng]*angArray[94];
         bas[11] = rad[ip+2*ng]*3*angArray[46]-rad[ip+3*ng]*angArray[95];
         bas[12] = rad[ip+2*ng]*3*angArray[47]-rad[ip+3*ng]*angArray[96];
         bas[13] = rad[ip+2*ng]*3*angArray[48]-rad[ip+3*ng]*angArray[97];
         bas[14] = rad[ip+2*ng]*3*angArray[49]-rad[ip+3*ng]*angArray[98];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*60*angArray[4]-rad[ip+ng]*75*angArray[20]+rad[ip+2*ng]*18*angArray[56]-rad[ip+3*ng]*angArray[120];
         bas[1] = rad[ip]*24*angArray[5]-rad[ip+ng]*48*angArray[21]+rad[ip+2*ng]*15*angArray[57]-rad[ip+3*ng]*angArray[121];
         bas[2] = rad[ip]*24*angArray[6]-rad[ip+ng]*48*angArray[22]+rad[ip+2*ng]*15*angArray[58]-rad[ip+3*ng]*angArray[122];
         bas[3] = rad[ip]*6*angArray[7]-rad[ip+ng]*27*angArray[23]+rad[ip+2*ng]*12*angArray[59]-rad[ip+3*ng]*angArray[123];
         bas[4] = rad[ip]*6*angArray[8]-rad[ip+ng]*27*angArray[24]+rad[ip+2*ng]*12*angArray[60]-rad[ip+3*ng]*angArray[124];
         bas[5] = rad[ip]*6*angArray[9]-rad[ip+ng]*27*angArray[25]+rad[ip+2*ng]*12*angArray[61]-rad[ip+3*ng]*angArray[125];
         bas[6] = -rad[ip+ng]*12*angArray[26]+rad[ip+2*ng]*9*angArray[62]-rad[ip+3*ng]*angArray[126];
         bas[7] = -rad[ip+ng]*12*angArray[27]+rad[ip+2*ng]*9*angArray[63]-rad[ip+3*ng]*angArray[127];
         bas[8] = -rad[ip+ng]*12*angArray[28]+rad[ip+2*ng]*9*angArray[64]-rad[ip+3*ng]*angArray[128];
         bas[9] = -rad[ip+ng]*12*angArray[29]+rad[ip+2*ng]*9*angArray[65]-rad[ip+3*ng]*angArray[129];
         bas[10] = -rad[ip+ng]*3*angArray[30]+rad[ip+2*ng]*6*angArray[66]-rad[ip+3*ng]*angArray[130];
         bas[11] = -rad[ip+ng]*3*angArray[31]+rad[ip+2*ng]*6*angArray[67]-rad[ip+3*ng]*angArray[131];
         bas[12] = -rad[ip+ng]*3*angArray[32]+rad[ip+2*ng]*6*angArray[68]-rad[ip+3*ng]*angArray[132];
         bas[13] = -rad[ip+ng]*3*angArray[33]+rad[ip+2*ng]*6*angArray[69]-rad[ip+3*ng]*angArray[133];
         bas[14] = -rad[ip+ng]*3*angArray[34]+rad[ip+2*ng]*6*angArray[70]-rad[ip+3*ng]*angArray[134];
         bas[15] = rad[ip+2*ng]*3*angArray[71]-rad[ip+3*ng]*angArray[135];
         bas[16] = rad[ip+2*ng]*3*angArray[72]-rad[ip+3*ng]*angArray[136];
         bas[17] = rad[ip+2*ng]*3*angArray[73]-rad[ip+3*ng]*angArray[137];
         bas[18] = rad[ip+2*ng]*3*angArray[74]-rad[ip+3*ng]*angArray[138];
         bas[19] = rad[ip+2*ng]*3*angArray[75]-rad[ip+3*ng]*angArray[139];
         bas[20] = rad[ip+2*ng]*3*angArray[76]-rad[ip+3*ng]*angArray[140];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*120*angArray[10]-rad[ip+ng]*108*angArray[35]+rad[ip+2*ng]*21*angArray[84]-rad[ip+3*ng]*angArray[165];
         bas[1] = rad[ip]*60*angArray[11]-rad[ip+ng]*75*angArray[36]+rad[ip+2*ng]*18*angArray[85]-rad[ip+3*ng]*angArray[166];
         bas[2] = rad[ip]*60*angArray[12]-rad[ip+ng]*75*angArray[37]+rad[ip+2*ng]*18*angArray[86]-rad[ip+3*ng]*angArray[167];
         bas[3] = rad[ip]*24*angArray[13]-rad[ip+ng]*48*angArray[38]+rad[ip+2*ng]*15*angArray[87]-rad[ip+3*ng]*angArray[168];
         bas[4] = rad[ip]*24*angArray[14]-rad[ip+ng]*48*angArray[39]+rad[ip+2*ng]*15*angArray[88]-rad[ip+3*ng]*angArray[169];
         bas[5] = rad[ip]*24*angArray[15]-rad[ip+ng]*48*angArray[40]+rad[ip+2*ng]*15*angArray[89]-rad[ip+3*ng]*angArray[170];
         bas[6] = rad[ip]*6*angArray[16]-rad[ip+ng]*27*angArray[41]+rad[ip+2*ng]*12*angArray[90]-rad[ip+3*ng]*angArray[171];
         bas[7] = rad[ip]*6*angArray[17]-rad[ip+ng]*27*angArray[42]+rad[ip+2*ng]*12*angArray[91]-rad[ip+3*ng]*angArray[172];
         bas[8] = rad[ip]*6*angArray[18]-rad[ip+ng]*27*angArray[43]+rad[ip+2*ng]*12*angArray[92]-rad[ip+3*ng]*angArray[173];
         bas[9] = rad[ip]*6*angArray[19]-rad[ip+ng]*27*angArray[44]+rad[ip+2*ng]*12*angArray[93]-rad[ip+3*ng]*angArray[174];
         bas[10] = -rad[ip+ng]*12*angArray[45]+rad[ip+2*ng]*9*angArray[94]-rad[ip+3*ng]*angArray[175];
         bas[11] = -rad[ip+ng]*12*angArray[46]+rad[ip+2*ng]*9*angArray[95]-rad[ip+3*ng]*angArray[176];
         bas[12] = -rad[ip+ng]*12*angArray[47]+rad[ip+2*ng]*9*angArray[96]-rad[ip+3*ng]*angArray[177];
         bas[13] = -rad[ip+ng]*12*angArray[48]+rad[ip+2*ng]*9*angArray[97]-rad[ip+3*ng]*angArray[178];
         bas[14] = -rad[ip+ng]*12*angArray[49]+rad[ip+2*ng]*9*angArray[98]-rad[ip+3*ng]*angArray[179];
         bas[15] = -rad[ip+ng]*3*angArray[50]+rad[ip+2*ng]*6*angArray[99]-rad[ip+3*ng]*angArray[180];
         bas[16] = -rad[ip+ng]*3*angArray[51]+rad[ip+2*ng]*6*angArray[100]-rad[ip+3*ng]*angArray[181];
         bas[17] = -rad[ip+ng]*3*angArray[52]+rad[ip+2*ng]*6*angArray[101]-rad[ip+3*ng]*angArray[182];
         bas[18] = -rad[ip+ng]*3*angArray[53]+rad[ip+2*ng]*6*angArray[102]-rad[ip+3*ng]*angArray[183];
         bas[19] = -rad[ip+ng]*3*angArray[54]+rad[ip+2*ng]*6*angArray[103]-rad[ip+3*ng]*angArray[184];
         bas[20] = -rad[ip+ng]*3*angArray[55]+rad[ip+2*ng]*6*angArray[104]-rad[ip+3*ng]*angArray[185];
         bas[21] = rad[ip+2*ng]*3*angArray[105]-rad[ip+3*ng]*angArray[186];
         bas[22] = rad[ip+2*ng]*3*angArray[106]-rad[ip+3*ng]*angArray[187];
         bas[23] = rad[ip+2*ng]*3*angArray[107]-rad[ip+3*ng]*angArray[188];
         bas[24] = rad[ip+2*ng]*3*angArray[108]-rad[ip+3*ng]*angArray[189];
         bas[25] = rad[ip+2*ng]*3*angArray[109]-rad[ip+3*ng]*angArray[190];
         bas[26] = rad[ip+2*ng]*3*angArray[110]-rad[ip+3*ng]*angArray[191];
         bas[27] = rad[ip+2*ng]*3*angArray[111]-rad[ip+3*ng]*angArray[192];
      }

   }


   // now we do derivatives for the given basis set to XXY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[2]-rad[ip+3*ng]*angArray[11];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*angArray[21];
         bas[1] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[7]+angArray[4])-rad[ip+3*ng]*angArray[23];
         bas[2] = rad[ip+2*ng]*angArray[8]-rad[ip+3*ng]*angArray[24];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*5*angArray[11]-rad[ip+3*ng]*angArray[36];
         bas[1] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*(3*angArray[13]+angArray[10])-rad[ip+3*ng]*angArray[38];
         bas[2] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[39];
         bas[3] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(angArray[16]+2*angArray[11])-rad[ip+3*ng]*angArray[41];
         bas[4] = -rad[ip+ng]*angArray[3]+rad[ip+2*ng]*(angArray[17]+angArray[12])-rad[ip+3*ng]*angArray[42];
         bas[5] = rad[ip+2*ng]*angArray[18]-rad[ip+3*ng]*angArray[43];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*7*angArray[21]-rad[ip+3*ng]*angArray[57];
         bas[1] = rad[ip]*2*angArray[0]-rad[ip+ng]*(5*angArray[4]+2*angArray[7])+rad[ip+2*ng]*(5*angArray[23]+angArray[20])-rad[ip+3*ng]*angArray[59];
         bas[2] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*5*angArray[24]-rad[ip+3*ng]*angArray[60];
         bas[3] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(3*angArray[26]+2*angArray[21])-rad[ip+3*ng]*angArray[62];
         bas[4] = -rad[ip+ng]*3*angArray[6]+rad[ip+2*ng]*(3*angArray[27]+angArray[22])-rad[ip+3*ng]*angArray[63];
         bas[5] = rad[ip+2*ng]*3*angArray[28]-rad[ip+3*ng]*angArray[64];
         bas[6] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*(angArray[30]+3*angArray[23])-rad[ip+3*ng]*angArray[66];
         bas[7] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*(angArray[31]+2*angArray[24])-rad[ip+3*ng]*angArray[67];
         bas[8] = -rad[ip+ng]*angArray[9]+rad[ip+2*ng]*(angArray[32]+angArray[25])-rad[ip+3*ng]*angArray[68];
         bas[9] = rad[ip+2*ng]*angArray[33]-rad[ip+3*ng]*angArray[69];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*12*angArray[11]+rad[ip+2*ng]*9*angArray[36]-rad[ip+3*ng]*angArray[85];
         bas[1] = rad[ip]*6*angArray[1]-rad[ip+ng]*(7*angArray[10]+6*angArray[13])+rad[ip+2*ng]*(7*angArray[38]+angArray[35])-rad[ip+3*ng]*angArray[87];
         bas[2] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*7*angArray[39]-rad[ip+3*ng]*angArray[88];
         bas[3] = rad[ip]*4*angArray[2]-rad[ip+ng]*(10*angArray[11]+2*angArray[16])+rad[ip+2*ng]*(5*angArray[41]+2*angArray[36])-rad[ip+3*ng]*angArray[90];
         bas[4] = rad[ip]*2*angArray[3]-rad[ip+ng]*(5*angArray[12]+2*angArray[17])+rad[ip+2*ng]*(5*angArray[42]+angArray[37])-rad[ip+3*ng]*angArray[91];
         bas[5] = -rad[ip+ng]*2*angArray[18]+rad[ip+2*ng]*5*angArray[43]-rad[ip+3*ng]*angArray[92];
         bas[6] = -rad[ip+ng]*9*angArray[13]+rad[ip+2*ng]*(3*angArray[45]+3*angArray[38])-rad[ip+3*ng]*angArray[94];
         bas[7] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*(3*angArray[46]+2*angArray[39])-rad[ip+3*ng]*angArray[95];
         bas[8] = -rad[ip+ng]*3*angArray[15]+rad[ip+2*ng]*(3*angArray[47]+angArray[40])-rad[ip+3*ng]*angArray[96];
         bas[9] = rad[ip+2*ng]*3*angArray[48]-rad[ip+3*ng]*angArray[97];
         bas[10] = -rad[ip+ng]*4*angArray[16]+rad[ip+2*ng]*(4*angArray[41]+angArray[50])-rad[ip+3*ng]*angArray[99];
         bas[11] = -rad[ip+ng]*3*angArray[17]+rad[ip+2*ng]*(angArray[51]+3*angArray[42])-rad[ip+3*ng]*angArray[100];
         bas[12] = -rad[ip+ng]*2*angArray[18]+rad[ip+2*ng]*(angArray[52]+2*angArray[43])-rad[ip+3*ng]*angArray[101];
         bas[13] = -rad[ip+ng]*angArray[19]+rad[ip+2*ng]*(angArray[53]+angArray[44])-rad[ip+3*ng]*angArray[102];
         bas[14] = rad[ip+2*ng]*angArray[54]-rad[ip+3*ng]*angArray[103];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*20*angArray[21]+rad[ip+2*ng]*11*angArray[57]-rad[ip+3*ng]*angArray[121];
         bas[1] = rad[ip]*12*angArray[4]-rad[ip+ng]*(9*angArray[20]+12*angArray[23])+rad[ip+2*ng]*(9*angArray[59]+angArray[56])-rad[ip+3*ng]*angArray[123];
         bas[2] = -rad[ip+ng]*12*angArray[24]+rad[ip+2*ng]*9*angArray[60]-rad[ip+3*ng]*angArray[124];
         bas[3] = rad[ip]*12*angArray[5]-rad[ip+ng]*(14*angArray[21]+6*angArray[26])+rad[ip+2*ng]*(7*angArray[62]+2*angArray[57])-rad[ip+3*ng]*angArray[126];
         bas[4] = rad[ip]*6*angArray[6]-rad[ip+ng]*(7*angArray[22]+6*angArray[27])+rad[ip+2*ng]*(7*angArray[63]+angArray[58])-rad[ip+3*ng]*angArray[127];
         bas[5] = -rad[ip+ng]*6*angArray[28]+rad[ip+2*ng]*7*angArray[64]-rad[ip+3*ng]*angArray[128];
         bas[6] = rad[ip]*6*angArray[7]-rad[ip+ng]*(15*angArray[23]+2*angArray[30])+rad[ip+2*ng]*(5*angArray[66]+3*angArray[59])-rad[ip+3*ng]*angArray[130];
         bas[7] = rad[ip]*4*angArray[8]-rad[ip+ng]*(10*angArray[24]+2*angArray[31])+rad[ip+2*ng]*(5*angArray[67]+2*angArray[60])-rad[ip+3*ng]*angArray[131];
         bas[8] = rad[ip]*2*angArray[9]-rad[ip+ng]*(5*angArray[25]+2*angArray[32])+rad[ip+2*ng]*(5*angArray[68]+angArray[61])-rad[ip+3*ng]*angArray[132];
         bas[9] = -rad[ip+ng]*2*angArray[33]+rad[ip+2*ng]*5*angArray[69]-rad[ip+3*ng]*angArray[133];
         bas[10] = -rad[ip+ng]*12*angArray[26]+rad[ip+2*ng]*(3*angArray[71]+4*angArray[62])-rad[ip+3*ng]*angArray[135];
         bas[11] = -rad[ip+ng]*9*angArray[27]+rad[ip+2*ng]*(3*angArray[72]+3*angArray[63])-rad[ip+3*ng]*angArray[136];
         bas[12] = -rad[ip+ng]*6*angArray[28]+rad[ip+2*ng]*(3*angArray[73]+2*angArray[64])-rad[ip+3*ng]*angArray[137];
         bas[13] = -rad[ip+ng]*3*angArray[29]+rad[ip+2*ng]*(3*angArray[74]+angArray[65])-rad[ip+3*ng]*angArray[138];
         bas[14] = rad[ip+2*ng]*3*angArray[75]-rad[ip+3*ng]*angArray[139];
         bas[15] = -rad[ip+ng]*5*angArray[30]+rad[ip+2*ng]*(5*angArray[66]+angArray[77])-rad[ip+3*ng]*angArray[141];
         bas[16] = -rad[ip+ng]*4*angArray[31]+rad[ip+2*ng]*(4*angArray[67]+angArray[78])-rad[ip+3*ng]*angArray[142];
         bas[17] = -rad[ip+ng]*3*angArray[32]+rad[ip+2*ng]*(angArray[79]+3*angArray[68])-rad[ip+3*ng]*angArray[143];
         bas[18] = -rad[ip+ng]*2*angArray[33]+rad[ip+2*ng]*(angArray[80]+2*angArray[69])-rad[ip+3*ng]*angArray[144];
         bas[19] = -rad[ip+ng]*angArray[34]+rad[ip+2*ng]*(angArray[81]+angArray[70])-rad[ip+3*ng]*angArray[145];
         bas[20] = rad[ip+2*ng]*angArray[82]-rad[ip+3*ng]*angArray[146];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*30*angArray[36]+rad[ip+2*ng]*13*angArray[85]-rad[ip+3*ng]*angArray[166];
         bas[1] = rad[ip]*20*angArray[10]-rad[ip+ng]*(11*angArray[35]+20*angArray[38])+rad[ip+2*ng]*(11*angArray[87]+angArray[84])-rad[ip+3*ng]*angArray[168];
         bas[2] = -rad[ip+ng]*20*angArray[39]+rad[ip+2*ng]*11*angArray[88]-rad[ip+3*ng]*angArray[169];
         bas[3] = rad[ip]*24*angArray[11]-rad[ip+ng]*(12*angArray[41]+18*angArray[36])+rad[ip+2*ng]*(2*angArray[85]+9*angArray[90])-rad[ip+3*ng]*angArray[171];
         bas[4] = rad[ip]*12*angArray[12]-rad[ip+ng]*(9*angArray[37]+12*angArray[42])+rad[ip+2*ng]*(9*angArray[91]+angArray[86])-rad[ip+3*ng]*angArray[172];
         bas[5] = -rad[ip+ng]*12*angArray[43]+rad[ip+2*ng]*9*angArray[92]-rad[ip+3*ng]*angArray[173];
         bas[6] = rad[ip]*18*angArray[13]-rad[ip+ng]*(21*angArray[38]+6*angArray[45])+rad[ip+2*ng]*(7*angArray[94]+3*angArray[87])-rad[ip+3*ng]*angArray[175];
         bas[7] = rad[ip]*12*angArray[14]-rad[ip+ng]*(14*angArray[39]+6*angArray[46])+rad[ip+2*ng]*(7*angArray[95]+2*angArray[88])-rad[ip+3*ng]*angArray[176];
         bas[8] = rad[ip]*6*angArray[15]-rad[ip+ng]*(7*angArray[40]+6*angArray[47])+rad[ip+2*ng]*(7*angArray[96]+angArray[89])-rad[ip+3*ng]*angArray[177];
         bas[9] = -rad[ip+ng]*6*angArray[48]+rad[ip+2*ng]*7*angArray[97]-rad[ip+3*ng]*angArray[178];
         bas[10] = rad[ip]*8*angArray[16]-rad[ip+ng]*(20*angArray[41]+2*angArray[50])+rad[ip+2*ng]*(5*angArray[99]+4*angArray[90])-rad[ip+3*ng]*angArray[180];
         bas[11] = rad[ip]*6*angArray[17]-rad[ip+ng]*(15*angArray[42]+2*angArray[51])+rad[ip+2*ng]*(5*angArray[100]+3*angArray[91])-rad[ip+3*ng]*angArray[181];
         bas[12] = rad[ip]*4*angArray[18]-rad[ip+ng]*(10*angArray[43]+2*angArray[52])+rad[ip+2*ng]*(5*angArray[101]+2*angArray[92])-rad[ip+3*ng]*angArray[182];
         bas[13] = rad[ip]*2*angArray[19]-rad[ip+ng]*(5*angArray[44]+2*angArray[53])+rad[ip+2*ng]*(5*angArray[102]+angArray[93])-rad[ip+3*ng]*angArray[183];
         bas[14] = -rad[ip+ng]*2*angArray[54]+rad[ip+2*ng]*5*angArray[103]-rad[ip+3*ng]*angArray[184];
         bas[15] = -rad[ip+ng]*15*angArray[45]+rad[ip+2*ng]*(3*angArray[105]+5*angArray[94])-rad[ip+3*ng]*angArray[186];
         bas[16] = -rad[ip+ng]*12*angArray[46]+rad[ip+2*ng]*(3*angArray[106]+4*angArray[95])-rad[ip+3*ng]*angArray[187];
         bas[17] = -rad[ip+ng]*9*angArray[47]+rad[ip+2*ng]*(3*angArray[107]+3*angArray[96])-rad[ip+3*ng]*angArray[188];
         bas[18] = -rad[ip+ng]*6*angArray[48]+rad[ip+2*ng]*(3*angArray[108]+2*angArray[97])-rad[ip+3*ng]*angArray[189];
         bas[19] = -rad[ip+ng]*3*angArray[49]+rad[ip+2*ng]*(3*angArray[109]+angArray[98])-rad[ip+3*ng]*angArray[190];
         bas[20] = rad[ip+2*ng]*3*angArray[110]-rad[ip+3*ng]*angArray[191];
         bas[21] = -rad[ip+ng]*6*angArray[50]+rad[ip+2*ng]*(6*angArray[99]+angArray[112])-rad[ip+3*ng]*angArray[193];
         bas[22] = -rad[ip+ng]*5*angArray[51]+rad[ip+2*ng]*(5*angArray[100]+angArray[113])-rad[ip+3*ng]*angArray[194];
         bas[23] = -rad[ip+ng]*4*angArray[52]+rad[ip+2*ng]*(4*angArray[101]+angArray[114])-rad[ip+3*ng]*angArray[195];
         bas[24] = -rad[ip+ng]*3*angArray[53]+rad[ip+2*ng]*(angArray[115]+3*angArray[102])-rad[ip+3*ng]*angArray[196];
         bas[25] = -rad[ip+ng]*2*angArray[54]+rad[ip+2*ng]*(angArray[116]+2*angArray[103])-rad[ip+3*ng]*angArray[197];
         bas[26] = -rad[ip+ng]*angArray[55]+rad[ip+2*ng]*(angArray[117]+angArray[104])-rad[ip+3*ng]*angArray[198];
         bas[27] = rad[ip+2*ng]*angArray[118]-rad[ip+3*ng]*angArray[199];
      }

   }


   // now we do derivatives for the given basis set to XYY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[1]-rad[ip+3*ng]*angArray[13];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[7]+angArray[4])-rad[ip+3*ng]*angArray[23];
         bas[1] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*angArray[26];
         bas[2] = rad[ip+2*ng]*angArray[6]-rad[ip+3*ng]*angArray[27];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(2*angArray[13]+angArray[10])-rad[ip+3*ng]*angArray[38];
         bas[1] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*(angArray[16]+3*angArray[11])-rad[ip+3*ng]*angArray[41];
         bas[2] = -rad[ip+ng]*angArray[3]+rad[ip+2*ng]*(angArray[17]+angArray[12])-rad[ip+3*ng]*angArray[42];
         bas[3] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*5*angArray[13]-rad[ip+3*ng]*angArray[45];
         bas[4] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[46];
         bas[5] = rad[ip+2*ng]*angArray[15]-rad[ip+3*ng]*angArray[47];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*(3*angArray[23]+angArray[20])-rad[ip+3*ng]*angArray[59];
         bas[1] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(2*angArray[26]+3*angArray[21])-rad[ip+3*ng]*angArray[62];
         bas[2] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*(2*angArray[27]+angArray[22])-rad[ip+3*ng]*angArray[63];
         bas[3] = rad[ip]*2*angArray[0]-rad[ip+ng]*(5*angArray[7]+2*angArray[4])+rad[ip+2*ng]*(5*angArray[23]+angArray[30])-rad[ip+3*ng]*angArray[66];
         bas[4] = -rad[ip+ng]*3*angArray[8]+rad[ip+2*ng]*(angArray[31]+3*angArray[24])-rad[ip+3*ng]*angArray[67];
         bas[5] = -rad[ip+ng]*angArray[9]+rad[ip+2*ng]*(angArray[32]+angArray[25])-rad[ip+3*ng]*angArray[68];
         bas[6] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*7*angArray[26]-rad[ip+3*ng]*angArray[71];
         bas[7] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*5*angArray[27]-rad[ip+3*ng]*angArray[72];
         bas[8] = rad[ip+2*ng]*3*angArray[28]-rad[ip+3*ng]*angArray[73];
         bas[9] = rad[ip+2*ng]*angArray[29]-rad[ip+3*ng]*angArray[74];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*4*angArray[10]+rad[ip+2*ng]*(4*angArray[38]+angArray[35])-rad[ip+3*ng]*angArray[87];
         bas[1] = -rad[ip+ng]*9*angArray[11]+rad[ip+2*ng]*(3*angArray[41]+3*angArray[36])-rad[ip+3*ng]*angArray[90];
         bas[2] = -rad[ip+ng]*3*angArray[12]+rad[ip+2*ng]*(3*angArray[42]+angArray[37])-rad[ip+3*ng]*angArray[91];
         bas[3] = rad[ip]*4*angArray[1]-rad[ip+ng]*(10*angArray[13]+2*angArray[10])+rad[ip+2*ng]*(2*angArray[45]+5*angArray[38])-rad[ip+3*ng]*angArray[94];
         bas[4] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*(2*angArray[46]+3*angArray[39])-rad[ip+3*ng]*angArray[95];
         bas[5] = -rad[ip+ng]*2*angArray[15]+rad[ip+2*ng]*(2*angArray[47]+angArray[40])-rad[ip+3*ng]*angArray[96];
         bas[6] = rad[ip]*6*angArray[2]-rad[ip+ng]*(6*angArray[11]+7*angArray[16])+rad[ip+2*ng]*(angArray[50]+7*angArray[41])-rad[ip+3*ng]*angArray[99];
         bas[7] = rad[ip]*2*angArray[3]-rad[ip+ng]*(5*angArray[17]+2*angArray[12])+rad[ip+2*ng]*(5*angArray[42]+angArray[51])-rad[ip+3*ng]*angArray[100];
         bas[8] = -rad[ip+ng]*3*angArray[18]+rad[ip+2*ng]*(angArray[52]+3*angArray[43])-rad[ip+3*ng]*angArray[101];
         bas[9] = -rad[ip+ng]*angArray[19]+rad[ip+2*ng]*(angArray[53]+angArray[44])-rad[ip+3*ng]*angArray[102];
         bas[10] = -rad[ip+ng]*12*angArray[13]+rad[ip+2*ng]*9*angArray[45]-rad[ip+3*ng]*angArray[105];
         bas[11] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*7*angArray[46]-rad[ip+3*ng]*angArray[106];
         bas[12] = -rad[ip+ng]*2*angArray[15]+rad[ip+2*ng]*5*angArray[47]-rad[ip+3*ng]*angArray[107];
         bas[13] = rad[ip+2*ng]*3*angArray[48]-rad[ip+3*ng]*angArray[108];
         bas[14] = rad[ip+2*ng]*angArray[49]-rad[ip+3*ng]*angArray[109];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*5*angArray[20]+rad[ip+2*ng]*(angArray[56]+5*angArray[59])-rad[ip+3*ng]*angArray[123];
         bas[1] = -rad[ip+ng]*12*angArray[21]+rad[ip+2*ng]*(4*angArray[62]+3*angArray[57])-rad[ip+3*ng]*angArray[126];
         bas[2] = -rad[ip+ng]*4*angArray[22]+rad[ip+2*ng]*(4*angArray[63]+angArray[58])-rad[ip+3*ng]*angArray[127];
         bas[3] = rad[ip]*6*angArray[4]-rad[ip+ng]*(15*angArray[23]+2*angArray[20])+rad[ip+2*ng]*(5*angArray[59]+3*angArray[66])-rad[ip+3*ng]*angArray[130];
         bas[4] = -rad[ip+ng]*9*angArray[24]+rad[ip+2*ng]*(3*angArray[67]+3*angArray[60])-rad[ip+3*ng]*angArray[131];
         bas[5] = -rad[ip+ng]*3*angArray[25]+rad[ip+2*ng]*(3*angArray[68]+angArray[61])-rad[ip+3*ng]*angArray[132];
         bas[6] = rad[ip]*12*angArray[5]-rad[ip+ng]*(6*angArray[21]+14*angArray[26])+rad[ip+2*ng]*(2*angArray[71]+7*angArray[62])-rad[ip+3*ng]*angArray[135];
         bas[7] = rad[ip]*4*angArray[6]-rad[ip+ng]*(10*angArray[27]+2*angArray[22])+rad[ip+2*ng]*(2*angArray[72]+5*angArray[63])-rad[ip+3*ng]*angArray[136];
         bas[8] = -rad[ip+ng]*6*angArray[28]+rad[ip+2*ng]*(2*angArray[73]+3*angArray[64])-rad[ip+3*ng]*angArray[137];
         bas[9] = -rad[ip+ng]*2*angArray[29]+rad[ip+2*ng]*(2*angArray[74]+angArray[65])-rad[ip+3*ng]*angArray[138];
         bas[10] = rad[ip]*12*angArray[7]-rad[ip+ng]*(12*angArray[23]+9*angArray[30])+rad[ip+2*ng]*(9*angArray[66]+angArray[77])-rad[ip+3*ng]*angArray[141];
         bas[11] = rad[ip]*6*angArray[8]-rad[ip+ng]*(6*angArray[24]+7*angArray[31])+rad[ip+2*ng]*(angArray[78]+7*angArray[67])-rad[ip+3*ng]*angArray[142];
         bas[12] = rad[ip]*2*angArray[9]-rad[ip+ng]*(5*angArray[32]+2*angArray[25])+rad[ip+2*ng]*(5*angArray[68]+angArray[79])-rad[ip+3*ng]*angArray[143];
         bas[13] = -rad[ip+ng]*3*angArray[33]+rad[ip+2*ng]*(angArray[80]+3*angArray[69])-rad[ip+3*ng]*angArray[144];
         bas[14] = -rad[ip+ng]*angArray[34]+rad[ip+2*ng]*(angArray[81]+angArray[70])-rad[ip+3*ng]*angArray[145];
         bas[15] = -rad[ip+ng]*20*angArray[26]+rad[ip+2*ng]*11*angArray[71]-rad[ip+3*ng]*angArray[148];
         bas[16] = -rad[ip+ng]*12*angArray[27]+rad[ip+2*ng]*9*angArray[72]-rad[ip+3*ng]*angArray[149];
         bas[17] = -rad[ip+ng]*6*angArray[28]+rad[ip+2*ng]*7*angArray[73]-rad[ip+3*ng]*angArray[150];
         bas[18] = -rad[ip+ng]*2*angArray[29]+rad[ip+2*ng]*5*angArray[74]-rad[ip+3*ng]*angArray[151];
         bas[19] = rad[ip+2*ng]*3*angArray[75]-rad[ip+3*ng]*angArray[152];
         bas[20] = rad[ip+2*ng]*angArray[76]-rad[ip+3*ng]*angArray[153];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[35]+rad[ip+2*ng]*(6*angArray[87]+angArray[84])-rad[ip+3*ng]*angArray[168];
         bas[1] = -rad[ip+ng]*15*angArray[36]+rad[ip+2*ng]*(5*angArray[90]+3*angArray[85])-rad[ip+3*ng]*angArray[171];
         bas[2] = -rad[ip+ng]*5*angArray[37]+rad[ip+2*ng]*(angArray[86]+5*angArray[91])-rad[ip+3*ng]*angArray[172];
         bas[3] = rad[ip]*8*angArray[10]-rad[ip+ng]*(20*angArray[38]+2*angArray[35])+rad[ip+2*ng]*(4*angArray[94]+5*angArray[87])-rad[ip+3*ng]*angArray[175];
         bas[4] = -rad[ip+ng]*12*angArray[39]+rad[ip+2*ng]*(4*angArray[95]+3*angArray[88])-rad[ip+3*ng]*angArray[176];
         bas[5] = -rad[ip+ng]*4*angArray[40]+rad[ip+2*ng]*(4*angArray[96]+angArray[89])-rad[ip+3*ng]*angArray[177];
         bas[6] = rad[ip]*18*angArray[11]-rad[ip+ng]*(6*angArray[36]+21*angArray[41])+rad[ip+2*ng]*(3*angArray[99]+7*angArray[90])-rad[ip+3*ng]*angArray[180];
         bas[7] = rad[ip]*6*angArray[12]-rad[ip+ng]*(15*angArray[42]+2*angArray[37])+rad[ip+2*ng]*(5*angArray[91]+3*angArray[100])-rad[ip+3*ng]*angArray[181];
         bas[8] = -rad[ip+ng]*9*angArray[43]+rad[ip+2*ng]*(3*angArray[101]+3*angArray[92])-rad[ip+3*ng]*angArray[182];
         bas[9] = -rad[ip+ng]*3*angArray[44]+rad[ip+2*ng]*(3*angArray[102]+angArray[93])-rad[ip+3*ng]*angArray[183];
         bas[10] = rad[ip]*24*angArray[13]-rad[ip+ng]*(12*angArray[38]+18*angArray[45])+rad[ip+2*ng]*(9*angArray[94]+2*angArray[105])-rad[ip+3*ng]*angArray[186];
         bas[11] = rad[ip]*12*angArray[14]-rad[ip+ng]*(6*angArray[39]+14*angArray[46])+rad[ip+2*ng]*(2*angArray[106]+7*angArray[95])-rad[ip+3*ng]*angArray[187];
         bas[12] = rad[ip]*4*angArray[15]-rad[ip+ng]*(10*angArray[47]+2*angArray[40])+rad[ip+2*ng]*(2*angArray[107]+5*angArray[96])-rad[ip+3*ng]*angArray[188];
         bas[13] = -rad[ip+ng]*6*angArray[48]+rad[ip+2*ng]*(2*angArray[108]+3*angArray[97])-rad[ip+3*ng]*angArray[189];
         bas[14] = -rad[ip+ng]*2*angArray[49]+rad[ip+2*ng]*(2*angArray[109]+angArray[98])-rad[ip+3*ng]*angArray[190];
         bas[15] = rad[ip]*20*angArray[16]-rad[ip+ng]*(11*angArray[50]+20*angArray[41])+rad[ip+2*ng]*(angArray[112]+11*angArray[99])-rad[ip+3*ng]*angArray[193];
         bas[16] = rad[ip]*12*angArray[17]-rad[ip+ng]*(12*angArray[42]+9*angArray[51])+rad[ip+2*ng]*(9*angArray[100]+angArray[113])-rad[ip+3*ng]*angArray[194];
         bas[17] = rad[ip]*6*angArray[18]-rad[ip+ng]*(6*angArray[43]+7*angArray[52])+rad[ip+2*ng]*(angArray[114]+7*angArray[101])-rad[ip+3*ng]*angArray[195];
         bas[18] = rad[ip]*2*angArray[19]-rad[ip+ng]*(5*angArray[53]+2*angArray[44])+rad[ip+2*ng]*(5*angArray[102]+angArray[115])-rad[ip+3*ng]*angArray[196];
         bas[19] = -rad[ip+ng]*3*angArray[54]+rad[ip+2*ng]*(angArray[116]+3*angArray[103])-rad[ip+3*ng]*angArray[197];
         bas[20] = -rad[ip+ng]*angArray[55]+rad[ip+2*ng]*(angArray[117]+angArray[104])-rad[ip+3*ng]*angArray[198];
         bas[21] = -rad[ip+ng]*30*angArray[45]+rad[ip+2*ng]*13*angArray[105]-rad[ip+3*ng]*angArray[201];
         bas[22] = -rad[ip+ng]*20*angArray[46]+rad[ip+2*ng]*11*angArray[106]-rad[ip+3*ng]*angArray[202];
         bas[23] = -rad[ip+ng]*12*angArray[47]+rad[ip+2*ng]*9*angArray[107]-rad[ip+3*ng]*angArray[203];
         bas[24] = -rad[ip+ng]*6*angArray[48]+rad[ip+2*ng]*7*angArray[108]-rad[ip+3*ng]*angArray[204];
         bas[25] = -rad[ip+ng]*2*angArray[49]+rad[ip+2*ng]*5*angArray[109]-rad[ip+3*ng]*angArray[205];
         bas[26] = rad[ip+2*ng]*3*angArray[110]-rad[ip+3*ng]*angArray[206];
         bas[27] = rad[ip+2*ng]*angArray[111]-rad[ip+3*ng]*angArray[207];
      }

   }


   // now we do derivatives for the given basis set to YYY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*angArray[16];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*angArray[26];
         bas[1] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*6*angArray[7]-rad[ip+3*ng]*angArray[30];
         bas[2] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*angArray[31];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[11]-rad[ip+3*ng]*angArray[41];
         bas[1] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*6*angArray[13]-rad[ip+3*ng]*angArray[45];
         bas[2] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[46];
         bas[3] = -rad[ip+ng]*12*angArray[2]+rad[ip+2*ng]*9*angArray[16]-rad[ip+3*ng]*angArray[50];
         bas[4] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*6*angArray[17]-rad[ip+3*ng]*angArray[51];
         bas[5] = rad[ip+2*ng]*3*angArray[18]-rad[ip+3*ng]*angArray[52];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[21]-rad[ip+3*ng]*angArray[62];
         bas[1] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*6*angArray[23]-rad[ip+3*ng]*angArray[66];
         bas[2] = rad[ip+2*ng]*3*angArray[24]-rad[ip+3*ng]*angArray[67];
         bas[3] = -rad[ip+ng]*12*angArray[5]+rad[ip+2*ng]*9*angArray[26]-rad[ip+3*ng]*angArray[71];
         bas[4] = -rad[ip+ng]*3*angArray[6]+rad[ip+2*ng]*6*angArray[27]-rad[ip+3*ng]*angArray[72];
         bas[5] = rad[ip+2*ng]*3*angArray[28]-rad[ip+3*ng]*angArray[73];
         bas[6] = rad[ip]*6*angArray[0]-rad[ip+ng]*27*angArray[7]+rad[ip+2*ng]*12*angArray[30]-rad[ip+3*ng]*angArray[77];
         bas[7] = -rad[ip+ng]*12*angArray[8]+rad[ip+2*ng]*9*angArray[31]-rad[ip+3*ng]*angArray[78];
         bas[8] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*6*angArray[32]-rad[ip+3*ng]*angArray[79];
         bas[9] = rad[ip+2*ng]*3*angArray[33]-rad[ip+3*ng]*angArray[80];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[36]-rad[ip+3*ng]*angArray[90];
         bas[1] = -rad[ip+ng]*3*angArray[10]+rad[ip+2*ng]*6*angArray[38]-rad[ip+3*ng]*angArray[94];
         bas[2] = rad[ip+2*ng]*3*angArray[39]-rad[ip+3*ng]*angArray[95];
         bas[3] = -rad[ip+ng]*12*angArray[11]+rad[ip+2*ng]*9*angArray[41]-rad[ip+3*ng]*angArray[99];
         bas[4] = -rad[ip+ng]*3*angArray[12]+rad[ip+2*ng]*6*angArray[42]-rad[ip+3*ng]*angArray[100];
         bas[5] = rad[ip+2*ng]*3*angArray[43]-rad[ip+3*ng]*angArray[101];
         bas[6] = rad[ip]*6*angArray[1]-rad[ip+ng]*27*angArray[13]+rad[ip+2*ng]*12*angArray[45]-rad[ip+3*ng]*angArray[105];
         bas[7] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*9*angArray[46]-rad[ip+3*ng]*angArray[106];
         bas[8] = -rad[ip+ng]*3*angArray[15]+rad[ip+2*ng]*6*angArray[47]-rad[ip+3*ng]*angArray[107];
         bas[9] = rad[ip+2*ng]*3*angArray[48]-rad[ip+3*ng]*angArray[108];
         bas[10] = rad[ip]*24*angArray[2]-rad[ip+ng]*48*angArray[16]+rad[ip+2*ng]*15*angArray[50]-rad[ip+3*ng]*angArray[112];
         bas[11] = rad[ip]*6*angArray[3]-rad[ip+ng]*27*angArray[17]+rad[ip+2*ng]*12*angArray[51]-rad[ip+3*ng]*angArray[113];
         bas[12] = -rad[ip+ng]*12*angArray[18]+rad[ip+2*ng]*9*angArray[52]-rad[ip+3*ng]*angArray[114];
         bas[13] = -rad[ip+ng]*3*angArray[19]+rad[ip+2*ng]*6*angArray[53]-rad[ip+3*ng]*angArray[115];
         bas[14] = rad[ip+2*ng]*3*angArray[54]-rad[ip+3*ng]*angArray[116];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[57]-rad[ip+3*ng]*angArray[126];
         bas[1] = -rad[ip+ng]*3*angArray[20]+rad[ip+2*ng]*6*angArray[59]-rad[ip+3*ng]*angArray[130];
         bas[2] = rad[ip+2*ng]*3*angArray[60]-rad[ip+3*ng]*angArray[131];
         bas[3] = -rad[ip+ng]*12*angArray[21]+rad[ip+2*ng]*9*angArray[62]-rad[ip+3*ng]*angArray[135];
         bas[4] = -rad[ip+ng]*3*angArray[22]+rad[ip+2*ng]*6*angArray[63]-rad[ip+3*ng]*angArray[136];
         bas[5] = rad[ip+2*ng]*3*angArray[64]-rad[ip+3*ng]*angArray[137];
         bas[6] = rad[ip]*6*angArray[4]-rad[ip+ng]*27*angArray[23]+rad[ip+2*ng]*12*angArray[66]-rad[ip+3*ng]*angArray[141];
         bas[7] = -rad[ip+ng]*12*angArray[24]+rad[ip+2*ng]*9*angArray[67]-rad[ip+3*ng]*angArray[142];
         bas[8] = -rad[ip+ng]*3*angArray[25]+rad[ip+2*ng]*6*angArray[68]-rad[ip+3*ng]*angArray[143];
         bas[9] = rad[ip+2*ng]*3*angArray[69]-rad[ip+3*ng]*angArray[144];
         bas[10] = rad[ip]*24*angArray[5]-rad[ip+ng]*48*angArray[26]+rad[ip+2*ng]*15*angArray[71]-rad[ip+3*ng]*angArray[148];
         bas[11] = rad[ip]*6*angArray[6]-rad[ip+ng]*27*angArray[27]+rad[ip+2*ng]*12*angArray[72]-rad[ip+3*ng]*angArray[149];
         bas[12] = -rad[ip+ng]*12*angArray[28]+rad[ip+2*ng]*9*angArray[73]-rad[ip+3*ng]*angArray[150];
         bas[13] = -rad[ip+ng]*3*angArray[29]+rad[ip+2*ng]*6*angArray[74]-rad[ip+3*ng]*angArray[151];
         bas[14] = rad[ip+2*ng]*3*angArray[75]-rad[ip+3*ng]*angArray[152];
         bas[15] = rad[ip]*60*angArray[7]-rad[ip+ng]*75*angArray[30]+rad[ip+2*ng]*18*angArray[77]-rad[ip+3*ng]*angArray[156];
         bas[16] = rad[ip]*24*angArray[8]-rad[ip+ng]*48*angArray[31]+rad[ip+2*ng]*15*angArray[78]-rad[ip+3*ng]*angArray[157];
         bas[17] = rad[ip]*6*angArray[9]-rad[ip+ng]*27*angArray[32]+rad[ip+2*ng]*12*angArray[79]-rad[ip+3*ng]*angArray[158];
         bas[18] = -rad[ip+ng]*12*angArray[33]+rad[ip+2*ng]*9*angArray[80]-rad[ip+3*ng]*angArray[159];
         bas[19] = -rad[ip+ng]*3*angArray[34]+rad[ip+2*ng]*6*angArray[81]-rad[ip+3*ng]*angArray[160];
         bas[20] = rad[ip+2*ng]*3*angArray[82]-rad[ip+3*ng]*angArray[161];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[85]-rad[ip+3*ng]*angArray[171];
         bas[1] = -rad[ip+ng]*3*angArray[35]+rad[ip+2*ng]*6*angArray[87]-rad[ip+3*ng]*angArray[175];
         bas[2] = rad[ip+2*ng]*3*angArray[88]-rad[ip+3*ng]*angArray[176];
         bas[3] = -rad[ip+ng]*12*angArray[36]+rad[ip+2*ng]*9*angArray[90]-rad[ip+3*ng]*angArray[180];
         bas[4] = -rad[ip+ng]*3*angArray[37]+rad[ip+2*ng]*6*angArray[91]-rad[ip+3*ng]*angArray[181];
         bas[5] = rad[ip+2*ng]*3*angArray[92]-rad[ip+3*ng]*angArray[182];
         bas[6] = rad[ip]*6*angArray[10]-rad[ip+ng]*27*angArray[38]+rad[ip+2*ng]*12*angArray[94]-rad[ip+3*ng]*angArray[186];
         bas[7] = -rad[ip+ng]*12*angArray[39]+rad[ip+2*ng]*9*angArray[95]-rad[ip+3*ng]*angArray[187];
         bas[8] = -rad[ip+ng]*3*angArray[40]+rad[ip+2*ng]*6*angArray[96]-rad[ip+3*ng]*angArray[188];
         bas[9] = rad[ip+2*ng]*3*angArray[97]-rad[ip+3*ng]*angArray[189];
         bas[10] = rad[ip]*24*angArray[11]-rad[ip+ng]*48*angArray[41]+rad[ip+2*ng]*15*angArray[99]-rad[ip+3*ng]*angArray[193];
         bas[11] = rad[ip]*6*angArray[12]-rad[ip+ng]*27*angArray[42]+rad[ip+2*ng]*12*angArray[100]-rad[ip+3*ng]*angArray[194];
         bas[12] = -rad[ip+ng]*12*angArray[43]+rad[ip+2*ng]*9*angArray[101]-rad[ip+3*ng]*angArray[195];
         bas[13] = -rad[ip+ng]*3*angArray[44]+rad[ip+2*ng]*6*angArray[102]-rad[ip+3*ng]*angArray[196];
         bas[14] = rad[ip+2*ng]*3*angArray[103]-rad[ip+3*ng]*angArray[197];
         bas[15] = rad[ip]*60*angArray[13]-rad[ip+ng]*75*angArray[45]+rad[ip+2*ng]*18*angArray[105]-rad[ip+3*ng]*angArray[201];
         bas[16] = rad[ip]*24*angArray[14]-rad[ip+ng]*48*angArray[46]+rad[ip+2*ng]*15*angArray[106]-rad[ip+3*ng]*angArray[202];
         bas[17] = rad[ip]*6*angArray[15]-rad[ip+ng]*27*angArray[47]+rad[ip+2*ng]*12*angArray[107]-rad[ip+3*ng]*angArray[203];
         bas[18] = -rad[ip+ng]*12*angArray[48]+rad[ip+2*ng]*9*angArray[108]-rad[ip+3*ng]*angArray[204];
         bas[19] = -rad[ip+ng]*3*angArray[49]+rad[ip+2*ng]*6*angArray[109]-rad[ip+3*ng]*angArray[205];
         bas[20] = rad[ip+2*ng]*3*angArray[110]-rad[ip+3*ng]*angArray[206];
         bas[21] = rad[ip]*120*angArray[16]-rad[ip+ng]*108*angArray[50]+rad[ip+2*ng]*21*angArray[112]-rad[ip+3*ng]*angArray[210];
         bas[22] = rad[ip]*60*angArray[17]-rad[ip+ng]*75*angArray[51]+rad[ip+2*ng]*18*angArray[113]-rad[ip+3*ng]*angArray[211];
         bas[23] = rad[ip]*24*angArray[18]-rad[ip+ng]*48*angArray[52]+rad[ip+2*ng]*15*angArray[114]-rad[ip+3*ng]*angArray[212];
         bas[24] = rad[ip]*6*angArray[19]-rad[ip+ng]*27*angArray[53]+rad[ip+2*ng]*12*angArray[115]-rad[ip+3*ng]*angArray[213];
         bas[25] = -rad[ip+ng]*12*angArray[54]+rad[ip+2*ng]*9*angArray[116]-rad[ip+3*ng]*angArray[214];
         bas[26] = -rad[ip+ng]*3*angArray[55]+rad[ip+2*ng]*6*angArray[117]-rad[ip+3*ng]*angArray[215];
         bas[27] = rad[ip+2*ng]*3*angArray[118]-rad[ip+3*ng]*angArray[216];
      }

   }


   // now we do derivatives for the given basis set to XXZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[3]-rad[ip+3*ng]*angArray[12];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*angArray[22];
         bas[1] = rad[ip+2*ng]*angArray[8]-rad[ip+3*ng]*angArray[24];
         bas[2] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[9]+angArray[4])-rad[ip+3*ng]*angArray[25];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*5*angArray[12]-rad[ip+3*ng]*angArray[37];
         bas[1] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[39];
         bas[2] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*(3*angArray[15]+angArray[10])-rad[ip+3*ng]*angArray[40];
         bas[3] = rad[ip+2*ng]*angArray[17]-rad[ip+3*ng]*angArray[42];
         bas[4] = -rad[ip+ng]*angArray[2]+rad[ip+2*ng]*(angArray[18]+angArray[11])-rad[ip+3*ng]*angArray[43];
         bas[5] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(angArray[19]+2*angArray[12])-rad[ip+3*ng]*angArray[44];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*7*angArray[22]-rad[ip+3*ng]*angArray[58];
         bas[1] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*5*angArray[24]-rad[ip+3*ng]*angArray[60];
         bas[2] = rad[ip]*2*angArray[0]-rad[ip+ng]*(5*angArray[4]+2*angArray[9])+rad[ip+2*ng]*(5*angArray[25]+angArray[20])-rad[ip+3*ng]*angArray[61];
         bas[3] = rad[ip+2*ng]*3*angArray[27]-rad[ip+3*ng]*angArray[63];
         bas[4] = -rad[ip+ng]*3*angArray[5]+rad[ip+2*ng]*(3*angArray[28]+angArray[21])-rad[ip+3*ng]*angArray[64];
         bas[5] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(3*angArray[29]+2*angArray[22])-rad[ip+3*ng]*angArray[65];
         bas[6] = rad[ip+2*ng]*angArray[31]-rad[ip+3*ng]*angArray[67];
         bas[7] = -rad[ip+ng]*angArray[7]+rad[ip+2*ng]*(angArray[32]+angArray[23])-rad[ip+3*ng]*angArray[68];
         bas[8] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*(angArray[33]+2*angArray[24])-rad[ip+3*ng]*angArray[69];
         bas[9] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*(angArray[34]+3*angArray[25])-rad[ip+3*ng]*angArray[70];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*12*angArray[12]+rad[ip+2*ng]*9*angArray[37]-rad[ip+3*ng]*angArray[86];
         bas[1] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*7*angArray[39]-rad[ip+3*ng]*angArray[88];
         bas[2] = rad[ip]*6*angArray[1]-rad[ip+ng]*(7*angArray[10]+6*angArray[15])+rad[ip+2*ng]*(7*angArray[40]+angArray[35])-rad[ip+3*ng]*angArray[89];
         bas[3] = -rad[ip+ng]*2*angArray[17]+rad[ip+2*ng]*5*angArray[42]-rad[ip+3*ng]*angArray[91];
         bas[4] = rad[ip]*2*angArray[2]-rad[ip+ng]*(5*angArray[11]+2*angArray[18])+rad[ip+2*ng]*(5*angArray[43]+angArray[36])-rad[ip+3*ng]*angArray[92];
         bas[5] = rad[ip]*4*angArray[3]-rad[ip+ng]*(10*angArray[12]+2*angArray[19])+rad[ip+2*ng]*(5*angArray[44]+2*angArray[37])-rad[ip+3*ng]*angArray[93];
         bas[6] = rad[ip+2*ng]*3*angArray[46]-rad[ip+3*ng]*angArray[95];
         bas[7] = -rad[ip+ng]*3*angArray[13]+rad[ip+2*ng]*(3*angArray[47]+angArray[38])-rad[ip+3*ng]*angArray[96];
         bas[8] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*(3*angArray[48]+2*angArray[39])-rad[ip+3*ng]*angArray[97];
         bas[9] = -rad[ip+ng]*9*angArray[15]+rad[ip+2*ng]*(3*angArray[49]+3*angArray[40])-rad[ip+3*ng]*angArray[98];
         bas[10] = rad[ip+2*ng]*angArray[51]-rad[ip+3*ng]*angArray[100];
         bas[11] = -rad[ip+ng]*angArray[16]+rad[ip+2*ng]*(angArray[52]+angArray[41])-rad[ip+3*ng]*angArray[101];
         bas[12] = -rad[ip+ng]*2*angArray[17]+rad[ip+2*ng]*(angArray[53]+2*angArray[42])-rad[ip+3*ng]*angArray[102];
         bas[13] = -rad[ip+ng]*3*angArray[18]+rad[ip+2*ng]*(angArray[54]+3*angArray[43])-rad[ip+3*ng]*angArray[103];
         bas[14] = -rad[ip+ng]*4*angArray[19]+rad[ip+2*ng]*(4*angArray[44]+angArray[55])-rad[ip+3*ng]*angArray[104];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*20*angArray[22]+rad[ip+2*ng]*11*angArray[58]-rad[ip+3*ng]*angArray[122];
         bas[1] = -rad[ip+ng]*12*angArray[24]+rad[ip+2*ng]*9*angArray[60]-rad[ip+3*ng]*angArray[124];
         bas[2] = rad[ip]*12*angArray[4]-rad[ip+ng]*(9*angArray[20]+12*angArray[25])+rad[ip+2*ng]*(9*angArray[61]+angArray[56])-rad[ip+3*ng]*angArray[125];
         bas[3] = -rad[ip+ng]*6*angArray[27]+rad[ip+2*ng]*7*angArray[63]-rad[ip+3*ng]*angArray[127];
         bas[4] = rad[ip]*6*angArray[5]-rad[ip+ng]*(7*angArray[21]+6*angArray[28])+rad[ip+2*ng]*(7*angArray[64]+angArray[57])-rad[ip+3*ng]*angArray[128];
         bas[5] = rad[ip]*12*angArray[6]-rad[ip+ng]*(14*angArray[22]+6*angArray[29])+rad[ip+2*ng]*(7*angArray[65]+2*angArray[58])-rad[ip+3*ng]*angArray[129];
         bas[6] = -rad[ip+ng]*2*angArray[31]+rad[ip+2*ng]*5*angArray[67]-rad[ip+3*ng]*angArray[131];
         bas[7] = rad[ip]*2*angArray[7]-rad[ip+ng]*(5*angArray[23]+2*angArray[32])+rad[ip+2*ng]*(5*angArray[68]+angArray[59])-rad[ip+3*ng]*angArray[132];
         bas[8] = rad[ip]*4*angArray[8]-rad[ip+ng]*(10*angArray[24]+2*angArray[33])+rad[ip+2*ng]*(5*angArray[69]+2*angArray[60])-rad[ip+3*ng]*angArray[133];
         bas[9] = rad[ip]*6*angArray[9]-rad[ip+ng]*(15*angArray[25]+2*angArray[34])+rad[ip+2*ng]*(5*angArray[70]+3*angArray[61])-rad[ip+3*ng]*angArray[134];
         bas[10] = rad[ip+2*ng]*3*angArray[72]-rad[ip+3*ng]*angArray[136];
         bas[11] = -rad[ip+ng]*3*angArray[26]+rad[ip+2*ng]*(3*angArray[73]+angArray[62])-rad[ip+3*ng]*angArray[137];
         bas[12] = -rad[ip+ng]*6*angArray[27]+rad[ip+2*ng]*(3*angArray[74]+2*angArray[63])-rad[ip+3*ng]*angArray[138];
         bas[13] = -rad[ip+ng]*9*angArray[28]+rad[ip+2*ng]*(3*angArray[75]+3*angArray[64])-rad[ip+3*ng]*angArray[139];
         bas[14] = -rad[ip+ng]*12*angArray[29]+rad[ip+2*ng]*(3*angArray[76]+4*angArray[65])-rad[ip+3*ng]*angArray[140];
         bas[15] = rad[ip+2*ng]*angArray[78]-rad[ip+3*ng]*angArray[142];
         bas[16] = -rad[ip+ng]*angArray[30]+rad[ip+2*ng]*(angArray[79]+angArray[66])-rad[ip+3*ng]*angArray[143];
         bas[17] = -rad[ip+ng]*2*angArray[31]+rad[ip+2*ng]*(angArray[80]+2*angArray[67])-rad[ip+3*ng]*angArray[144];
         bas[18] = -rad[ip+ng]*3*angArray[32]+rad[ip+2*ng]*(angArray[81]+3*angArray[68])-rad[ip+3*ng]*angArray[145];
         bas[19] = -rad[ip+ng]*4*angArray[33]+rad[ip+2*ng]*(4*angArray[69]+angArray[82])-rad[ip+3*ng]*angArray[146];
         bas[20] = -rad[ip+ng]*5*angArray[34]+rad[ip+2*ng]*(5*angArray[70]+angArray[83])-rad[ip+3*ng]*angArray[147];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*30*angArray[37]+rad[ip+2*ng]*13*angArray[86]-rad[ip+3*ng]*angArray[167];
         bas[1] = -rad[ip+ng]*20*angArray[39]+rad[ip+2*ng]*11*angArray[88]-rad[ip+3*ng]*angArray[169];
         bas[2] = rad[ip]*20*angArray[10]-rad[ip+ng]*(11*angArray[35]+20*angArray[40])+rad[ip+2*ng]*(11*angArray[89]+angArray[84])-rad[ip+3*ng]*angArray[170];
         bas[3] = -rad[ip+ng]*12*angArray[42]+rad[ip+2*ng]*9*angArray[91]-rad[ip+3*ng]*angArray[172];
         bas[4] = rad[ip]*12*angArray[11]-rad[ip+ng]*(9*angArray[36]+12*angArray[43])+rad[ip+2*ng]*(9*angArray[92]+angArray[85])-rad[ip+3*ng]*angArray[173];
         bas[5] = rad[ip]*24*angArray[12]-rad[ip+ng]*(12*angArray[44]+18*angArray[37])+rad[ip+2*ng]*(2*angArray[86]+9*angArray[93])-rad[ip+3*ng]*angArray[174];
         bas[6] = -rad[ip+ng]*6*angArray[46]+rad[ip+2*ng]*7*angArray[95]-rad[ip+3*ng]*angArray[176];
         bas[7] = rad[ip]*6*angArray[13]-rad[ip+ng]*(7*angArray[38]+6*angArray[47])+rad[ip+2*ng]*(7*angArray[96]+angArray[87])-rad[ip+3*ng]*angArray[177];
         bas[8] = rad[ip]*12*angArray[14]-rad[ip+ng]*(14*angArray[39]+6*angArray[48])+rad[ip+2*ng]*(7*angArray[97]+2*angArray[88])-rad[ip+3*ng]*angArray[178];
         bas[9] = rad[ip]*18*angArray[15]-rad[ip+ng]*(21*angArray[40]+6*angArray[49])+rad[ip+2*ng]*(7*angArray[98]+3*angArray[89])-rad[ip+3*ng]*angArray[179];
         bas[10] = -rad[ip+ng]*2*angArray[51]+rad[ip+2*ng]*5*angArray[100]-rad[ip+3*ng]*angArray[181];
         bas[11] = rad[ip]*2*angArray[16]-rad[ip+ng]*(5*angArray[41]+2*angArray[52])+rad[ip+2*ng]*(5*angArray[101]+angArray[90])-rad[ip+3*ng]*angArray[182];
         bas[12] = rad[ip]*4*angArray[17]-rad[ip+ng]*(10*angArray[42]+2*angArray[53])+rad[ip+2*ng]*(5*angArray[102]+2*angArray[91])-rad[ip+3*ng]*angArray[183];
         bas[13] = rad[ip]*6*angArray[18]-rad[ip+ng]*(15*angArray[43]+2*angArray[54])+rad[ip+2*ng]*(5*angArray[103]+3*angArray[92])-rad[ip+3*ng]*angArray[184];
         bas[14] = rad[ip]*8*angArray[19]-rad[ip+ng]*(20*angArray[44]+2*angArray[55])+rad[ip+2*ng]*(5*angArray[104]+4*angArray[93])-rad[ip+3*ng]*angArray[185];
         bas[15] = rad[ip+2*ng]*3*angArray[106]-rad[ip+3*ng]*angArray[187];
         bas[16] = -rad[ip+ng]*3*angArray[45]+rad[ip+2*ng]*(3*angArray[107]+angArray[94])-rad[ip+3*ng]*angArray[188];
         bas[17] = -rad[ip+ng]*6*angArray[46]+rad[ip+2*ng]*(3*angArray[108]+2*angArray[95])-rad[ip+3*ng]*angArray[189];
         bas[18] = -rad[ip+ng]*9*angArray[47]+rad[ip+2*ng]*(3*angArray[109]+3*angArray[96])-rad[ip+3*ng]*angArray[190];
         bas[19] = -rad[ip+ng]*12*angArray[48]+rad[ip+2*ng]*(3*angArray[110]+4*angArray[97])-rad[ip+3*ng]*angArray[191];
         bas[20] = -rad[ip+ng]*15*angArray[49]+rad[ip+2*ng]*(3*angArray[111]+5*angArray[98])-rad[ip+3*ng]*angArray[192];
         bas[21] = rad[ip+2*ng]*angArray[113]-rad[ip+3*ng]*angArray[194];
         bas[22] = -rad[ip+ng]*angArray[50]+rad[ip+2*ng]*(angArray[114]+angArray[99])-rad[ip+3*ng]*angArray[195];
         bas[23] = -rad[ip+ng]*2*angArray[51]+rad[ip+2*ng]*(angArray[115]+2*angArray[100])-rad[ip+3*ng]*angArray[196];
         bas[24] = -rad[ip+ng]*3*angArray[52]+rad[ip+2*ng]*(angArray[116]+3*angArray[101])-rad[ip+3*ng]*angArray[197];
         bas[25] = -rad[ip+ng]*4*angArray[53]+rad[ip+2*ng]*(4*angArray[102]+angArray[117])-rad[ip+3*ng]*angArray[198];
         bas[26] = -rad[ip+ng]*5*angArray[54]+rad[ip+2*ng]*(5*angArray[103]+angArray[118])-rad[ip+3*ng]*angArray[199];
         bas[27] = -rad[ip+ng]*6*angArray[55]+rad[ip+2*ng]*(6*angArray[104]+angArray[119])-rad[ip+3*ng]*angArray[200];
      }

   }


   // now we do derivatives for the given basis set to XYZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*angArray[14];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[8]-rad[ip+3*ng]*angArray[24];
         bas[1] = rad[ip+2*ng]*angArray[6]-rad[ip+3*ng]*angArray[27];
         bas[2] = rad[ip+2*ng]*angArray[5]-rad[ip+3*ng]*angArray[28];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*2*angArray[14]-rad[ip+3*ng]*angArray[39];
         bas[1] = -rad[ip+ng]*angArray[3]+rad[ip+2*ng]*(angArray[12]+angArray[17])-rad[ip+3*ng]*angArray[42];
         bas[2] = -rad[ip+ng]*angArray[2]+rad[ip+2*ng]*(angArray[18]+angArray[11])-rad[ip+3*ng]*angArray[43];
         bas[3] = rad[ip+2*ng]*2*angArray[14]-rad[ip+3*ng]*angArray[46];
         bas[4] = -rad[ip+ng]*angArray[1]+rad[ip+2*ng]*(angArray[15]+angArray[13])-rad[ip+3*ng]*angArray[47];
         bas[5] = rad[ip+2*ng]*2*angArray[14]-rad[ip+3*ng]*angArray[48];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[24]-rad[ip+3*ng]*angArray[60];
         bas[1] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*(2*angArray[27]+angArray[22])-rad[ip+3*ng]*angArray[63];
         bas[2] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*(2*angArray[28]+angArray[21])-rad[ip+3*ng]*angArray[64];
         bas[3] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*(2*angArray[24]+angArray[31])-rad[ip+3*ng]*angArray[67];
         bas[4] = rad[ip]*angArray[0]-rad[ip+ng]*(angArray[4]+angArray[9]+angArray[7])+rad[ip+2*ng]*(angArray[25]+angArray[32]+angArray[23])-rad[ip+3*ng]*angArray[68];
         bas[5] = -rad[ip+ng]*2*angArray[8]+rad[ip+2*ng]*(angArray[33]+2*angArray[24])-rad[ip+3*ng]*angArray[69];
         bas[6] = rad[ip+2*ng]*3*angArray[27]-rad[ip+3*ng]*angArray[72];
         bas[7] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*(2*angArray[28]+angArray[26])-rad[ip+3*ng]*angArray[73];
         bas[8] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*(angArray[29]+2*angArray[27])-rad[ip+3*ng]*angArray[74];
         bas[9] = rad[ip+2*ng]*3*angArray[28]-rad[ip+3*ng]*angArray[75];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*4*angArray[39]-rad[ip+3*ng]*angArray[88];
         bas[1] = -rad[ip+ng]*3*angArray[12]+rad[ip+2*ng]*(angArray[37]+3*angArray[42])-rad[ip+3*ng]*angArray[91];
         bas[2] = -rad[ip+ng]*3*angArray[11]+rad[ip+2*ng]*(3*angArray[43]+angArray[36])-rad[ip+3*ng]*angArray[92];
         bas[3] = -rad[ip+ng]*4*angArray[14]+rad[ip+2*ng]*(2*angArray[46]+2*angArray[39])-rad[ip+3*ng]*angArray[95];
         bas[4] = rad[ip]*2*angArray[1]-rad[ip+ng]*(2*angArray[13]+angArray[10]+2*angArray[15])+rad[ip+2*ng]*(2*angArray[47]+angArray[40]+angArray[38])-rad[ip+3*ng]*angArray[96];
         bas[5] = -rad[ip+ng]*4*angArray[14]+rad[ip+2*ng]*(2*angArray[48]+2*angArray[39])-rad[ip+3*ng]*angArray[97];
         bas[6] = -rad[ip+ng]*3*angArray[17]+rad[ip+2*ng]*(3*angArray[42]+angArray[51])-rad[ip+3*ng]*angArray[100];
         bas[7] = rad[ip]*2*angArray[2]-rad[ip+ng]*(2*angArray[18]+2*angArray[11]+angArray[16])+rad[ip+2*ng]*(2*angArray[43]+angArray[52]+angArray[41])-rad[ip+3*ng]*angArray[101];
         bas[8] = rad[ip]*2*angArray[3]-rad[ip+ng]*(angArray[19]+2*angArray[12]+2*angArray[17])+rad[ip+2*ng]*(angArray[53]+angArray[44]+2*angArray[42])-rad[ip+3*ng]*angArray[102];
         bas[9] = -rad[ip+ng]*3*angArray[18]+rad[ip+2*ng]*(angArray[54]+3*angArray[43])-rad[ip+3*ng]*angArray[103];
         bas[10] = rad[ip+2*ng]*4*angArray[46]-rad[ip+3*ng]*angArray[106];
         bas[11] = -rad[ip+ng]*3*angArray[13]+rad[ip+2*ng]*(3*angArray[47]+angArray[45])-rad[ip+3*ng]*angArray[107];
         bas[12] = -rad[ip+ng]*4*angArray[14]+rad[ip+2*ng]*(2*angArray[48]+2*angArray[46])-rad[ip+3*ng]*angArray[108];
         bas[13] = -rad[ip+ng]*3*angArray[15]+rad[ip+2*ng]*(angArray[49]+3*angArray[47])-rad[ip+3*ng]*angArray[109];
         bas[14] = rad[ip+2*ng]*4*angArray[48]-rad[ip+3*ng]*angArray[110];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*5*angArray[60]-rad[ip+3*ng]*angArray[124];
         bas[1] = -rad[ip+ng]*4*angArray[22]+rad[ip+2*ng]*(4*angArray[63]+angArray[58])-rad[ip+3*ng]*angArray[127];
         bas[2] = -rad[ip+ng]*4*angArray[21]+rad[ip+2*ng]*(4*angArray[64]+angArray[57])-rad[ip+3*ng]*angArray[128];
         bas[3] = -rad[ip+ng]*6*angArray[24]+rad[ip+2*ng]*(2*angArray[60]+3*angArray[67])-rad[ip+3*ng]*angArray[131];
         bas[4] = rad[ip]*3*angArray[4]-rad[ip+ng]*(angArray[20]+3*angArray[23]+3*angArray[25])+rad[ip+2*ng]*(3*angArray[68]+angArray[61]+angArray[59])-rad[ip+3*ng]*angArray[132];
         bas[5] = -rad[ip+ng]*6*angArray[24]+rad[ip+2*ng]*(2*angArray[60]+3*angArray[69])-rad[ip+3*ng]*angArray[133];
         bas[6] = -rad[ip+ng]*6*angArray[27]+rad[ip+2*ng]*(2*angArray[72]+3*angArray[63])-rad[ip+3*ng]*angArray[136];
         bas[7] = rad[ip]*4*angArray[5]-rad[ip+ng]*(2*angArray[26]+4*angArray[28]+2*angArray[21])+rad[ip+2*ng]*(2*angArray[73]+2*angArray[64]+angArray[62])-rad[ip+3*ng]*angArray[137];
         bas[8] = rad[ip]*4*angArray[6]-rad[ip+ng]*(4*angArray[27]+2*angArray[22]+2*angArray[29])+rad[ip+2*ng]*(2*angArray[74]+2*angArray[63]+angArray[65])-rad[ip+3*ng]*angArray[138];
         bas[9] = -rad[ip+ng]*6*angArray[28]+rad[ip+2*ng]*(2*angArray[75]+3*angArray[64])-rad[ip+3*ng]*angArray[139];
         bas[10] = -rad[ip+ng]*4*angArray[31]+rad[ip+2*ng]*(4*angArray[67]+angArray[78])-rad[ip+3*ng]*angArray[142];
         bas[11] = rad[ip]*3*angArray[7]-rad[ip+ng]*(3*angArray[23]+3*angArray[32]+angArray[30])+rad[ip+2*ng]*(angArray[79]+3*angArray[68]+angArray[66])-rad[ip+3*ng]*angArray[143];
         bas[12] = rad[ip]*4*angArray[8]-rad[ip+ng]*(2*angArray[33]+4*angArray[24]+2*angArray[31])+rad[ip+2*ng]*(2*angArray[69]+angArray[80]+2*angArray[67])-rad[ip+3*ng]*angArray[144];
         bas[13] = rad[ip]*3*angArray[9]-rad[ip+ng]*(angArray[34]+3*angArray[25]+3*angArray[32])+rad[ip+2*ng]*(angArray[81]+angArray[70]+3*angArray[68])-rad[ip+3*ng]*angArray[145];
         bas[14] = -rad[ip+ng]*4*angArray[33]+rad[ip+2*ng]*(4*angArray[69]+angArray[82])-rad[ip+3*ng]*angArray[146];
         bas[15] = rad[ip+2*ng]*5*angArray[72]-rad[ip+3*ng]*angArray[149];
         bas[16] = -rad[ip+ng]*4*angArray[26]+rad[ip+2*ng]*(4*angArray[73]+angArray[71])-rad[ip+3*ng]*angArray[150];
         bas[17] = -rad[ip+ng]*6*angArray[27]+rad[ip+2*ng]*(3*angArray[74]+2*angArray[72])-rad[ip+3*ng]*angArray[151];
         bas[18] = -rad[ip+ng]*6*angArray[28]+rad[ip+2*ng]*(2*angArray[75]+3*angArray[73])-rad[ip+3*ng]*angArray[152];
         bas[19] = -rad[ip+ng]*4*angArray[29]+rad[ip+2*ng]*(4*angArray[74]+angArray[76])-rad[ip+3*ng]*angArray[153];
         bas[20] = rad[ip+2*ng]*5*angArray[75]-rad[ip+3*ng]*angArray[154];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*6*angArray[88]-rad[ip+3*ng]*angArray[169];
         bas[1] = -rad[ip+ng]*5*angArray[37]+rad[ip+2*ng]*(angArray[86]+5*angArray[91])-rad[ip+3*ng]*angArray[172];
         bas[2] = -rad[ip+ng]*5*angArray[36]+rad[ip+2*ng]*(angArray[85]+5*angArray[92])-rad[ip+3*ng]*angArray[173];
         bas[3] = -rad[ip+ng]*8*angArray[39]+rad[ip+2*ng]*(4*angArray[95]+2*angArray[88])-rad[ip+3*ng]*angArray[176];
         bas[4] = rad[ip]*4*angArray[10]-rad[ip+ng]*(angArray[35]+4*angArray[40]+4*angArray[38])+rad[ip+2*ng]*(4*angArray[96]+angArray[89]+angArray[87])-rad[ip+3*ng]*angArray[177];
         bas[5] = -rad[ip+ng]*8*angArray[39]+rad[ip+2*ng]*(4*angArray[97]+2*angArray[88])-rad[ip+3*ng]*angArray[178];
         bas[6] = -rad[ip+ng]*9*angArray[42]+rad[ip+2*ng]*(3*angArray[100]+3*angArray[91])-rad[ip+3*ng]*angArray[181];
         bas[7] = rad[ip]*6*angArray[11]-rad[ip+ng]*(6*angArray[43]+2*angArray[36]+3*angArray[41])+rad[ip+2*ng]*(2*angArray[92]+3*angArray[101]+angArray[90])-rad[ip+3*ng]*angArray[182];
         bas[8] = rad[ip]*6*angArray[12]-rad[ip+ng]*(6*angArray[42]+2*angArray[37]+3*angArray[44])+rad[ip+2*ng]*(3*angArray[102]+angArray[93]+2*angArray[91])-rad[ip+3*ng]*angArray[183];
         bas[9] = -rad[ip+ng]*9*angArray[43]+rad[ip+2*ng]*(3*angArray[92]+3*angArray[103])-rad[ip+3*ng]*angArray[184];
         bas[10] = -rad[ip+ng]*8*angArray[46]+rad[ip+2*ng]*(4*angArray[95]+2*angArray[106])-rad[ip+3*ng]*angArray[187];
         bas[11] = rad[ip]*6*angArray[13]-rad[ip+ng]*(2*angArray[45]+3*angArray[38]+6*angArray[47])+rad[ip+2*ng]*(2*angArray[107]+3*angArray[96]+angArray[94])-rad[ip+3*ng]*angArray[188];
         bas[12] = rad[ip]*8*angArray[14]-rad[ip+ng]*(4*angArray[48]+4*angArray[46]+4*angArray[39])+rad[ip+2*ng]*(2*angArray[108]+2*angArray[95]+2*angArray[97])-rad[ip+3*ng]*angArray[189];
         bas[13] = rad[ip]*6*angArray[15]-rad[ip+ng]*(6*angArray[47]+3*angArray[40]+2*angArray[49])+rad[ip+2*ng]*(2*angArray[109]+angArray[98]+3*angArray[96])-rad[ip+3*ng]*angArray[190];
         bas[14] = -rad[ip+ng]*8*angArray[48]+rad[ip+2*ng]*(2*angArray[110]+4*angArray[97])-rad[ip+3*ng]*angArray[191];
         bas[15] = -rad[ip+ng]*5*angArray[51]+rad[ip+2*ng]*(angArray[113]+5*angArray[100])-rad[ip+3*ng]*angArray[194];
         bas[16] = rad[ip]*4*angArray[16]-rad[ip+ng]*(4*angArray[41]+4*angArray[52]+angArray[50])+rad[ip+2*ng]*(4*angArray[101]+angArray[114]+angArray[99])-rad[ip+3*ng]*angArray[195];
         bas[17] = rad[ip]*6*angArray[17]-rad[ip+ng]*(6*angArray[42]+3*angArray[53]+2*angArray[51])+rad[ip+2*ng]*(angArray[115]+3*angArray[102]+2*angArray[100])-rad[ip+3*ng]*angArray[196];
         bas[18] = rad[ip]*6*angArray[18]-rad[ip+ng]*(2*angArray[54]+6*angArray[43]+3*angArray[52])+rad[ip+2*ng]*(2*angArray[103]+angArray[116]+3*angArray[101])-rad[ip+3*ng]*angArray[197];
         bas[19] = rad[ip]*4*angArray[19]-rad[ip+ng]*(4*angArray[44]+angArray[55]+4*angArray[53])+rad[ip+2*ng]*(angArray[117]+4*angArray[102]+angArray[104])-rad[ip+3*ng]*angArray[198];
         bas[20] = -rad[ip+ng]*5*angArray[54]+rad[ip+2*ng]*(5*angArray[103]+angArray[118])-rad[ip+3*ng]*angArray[199];
         bas[21] = rad[ip+2*ng]*6*angArray[106]-rad[ip+3*ng]*angArray[202];
         bas[22] = -rad[ip+ng]*5*angArray[45]+rad[ip+2*ng]*(angArray[105]+5*angArray[107])-rad[ip+3*ng]*angArray[203];
         bas[23] = -rad[ip+ng]*8*angArray[46]+rad[ip+2*ng]*(4*angArray[108]+2*angArray[106])-rad[ip+3*ng]*angArray[204];
         bas[24] = -rad[ip+ng]*9*angArray[47]+rad[ip+2*ng]*(3*angArray[107]+3*angArray[109])-rad[ip+3*ng]*angArray[205];
         bas[25] = -rad[ip+ng]*8*angArray[48]+rad[ip+2*ng]*(2*angArray[110]+4*angArray[108])-rad[ip+3*ng]*angArray[206];
         bas[26] = -rad[ip+ng]*5*angArray[49]+rad[ip+2*ng]*(5*angArray[109]+angArray[111])-rad[ip+3*ng]*angArray[207];
         bas[27] = rad[ip+2*ng]*6*angArray[110]-rad[ip+3*ng]*angArray[208];
      }

   }


   // now we do derivatives for the given basis set to YYZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[3]-rad[ip+3*ng]*angArray[17];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[6]-rad[ip+3*ng]*angArray[27];
         bas[1] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*angArray[31];
         bas[2] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[9]+angArray[7])-rad[ip+3*ng]*angArray[32];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[12]-rad[ip+3*ng]*angArray[42];
         bas[1] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[46];
         bas[2] = -rad[ip+ng]*angArray[1]+rad[ip+2*ng]*(angArray[15]+angArray[13])-rad[ip+3*ng]*angArray[47];
         bas[3] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*5*angArray[17]-rad[ip+3*ng]*angArray[51];
         bas[4] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*(3*angArray[18]+angArray[16])-rad[ip+3*ng]*angArray[52];
         bas[5] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(angArray[19]+2*angArray[17])-rad[ip+3*ng]*angArray[53];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[22]-rad[ip+3*ng]*angArray[63];
         bas[1] = rad[ip+2*ng]*3*angArray[24]-rad[ip+3*ng]*angArray[67];
         bas[2] = -rad[ip+ng]*angArray[4]+rad[ip+2*ng]*(angArray[25]+angArray[23])-rad[ip+3*ng]*angArray[68];
         bas[3] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*5*angArray[27]-rad[ip+3*ng]*angArray[72];
         bas[4] = -rad[ip+ng]*3*angArray[5]+rad[ip+2*ng]*(3*angArray[28]+angArray[26])-rad[ip+3*ng]*angArray[73];
         bas[5] = -rad[ip+ng]*2*angArray[6]+rad[ip+2*ng]*(angArray[29]+2*angArray[27])-rad[ip+3*ng]*angArray[74];
         bas[6] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*7*angArray[31]-rad[ip+3*ng]*angArray[78];
         bas[7] = rad[ip]*2*angArray[0]-rad[ip+ng]*(5*angArray[7]+2*angArray[9])+rad[ip+2*ng]*(5*angArray[32]+angArray[30])-rad[ip+3*ng]*angArray[79];
         bas[8] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(3*angArray[33]+2*angArray[31])-rad[ip+3*ng]*angArray[80];
         bas[9] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*(angArray[34]+3*angArray[32])-rad[ip+3*ng]*angArray[81];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[37]-rad[ip+3*ng]*angArray[91];
         bas[1] = rad[ip+2*ng]*3*angArray[39]-rad[ip+3*ng]*angArray[95];
         bas[2] = -rad[ip+ng]*angArray[10]+rad[ip+2*ng]*(angArray[40]+angArray[38])-rad[ip+3*ng]*angArray[96];
         bas[3] = -rad[ip+ng]*2*angArray[12]+rad[ip+2*ng]*5*angArray[42]-rad[ip+3*ng]*angArray[100];
         bas[4] = -rad[ip+ng]*3*angArray[11]+rad[ip+2*ng]*(3*angArray[43]+angArray[41])-rad[ip+3*ng]*angArray[101];
         bas[5] = -rad[ip+ng]*2*angArray[12]+rad[ip+2*ng]*(angArray[44]+2*angArray[42])-rad[ip+3*ng]*angArray[102];
         bas[6] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*7*angArray[46]-rad[ip+3*ng]*angArray[106];
         bas[7] = rad[ip]*2*angArray[1]-rad[ip+ng]*(5*angArray[13]+2*angArray[15])+rad[ip+2*ng]*(5*angArray[47]+angArray[45])-rad[ip+3*ng]*angArray[107];
         bas[8] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*(3*angArray[48]+2*angArray[46])-rad[ip+3*ng]*angArray[108];
         bas[9] = -rad[ip+ng]*3*angArray[15]+rad[ip+2*ng]*(angArray[49]+3*angArray[47])-rad[ip+3*ng]*angArray[109];
         bas[10] = -rad[ip+ng]*12*angArray[17]+rad[ip+2*ng]*9*angArray[51]-rad[ip+3*ng]*angArray[113];
         bas[11] = rad[ip]*6*angArray[2]-rad[ip+ng]*(7*angArray[16]+6*angArray[18])+rad[ip+2*ng]*(7*angArray[52]+angArray[50])-rad[ip+3*ng]*angArray[114];
         bas[12] = rad[ip]*4*angArray[3]-rad[ip+ng]*(10*angArray[17]+2*angArray[19])+rad[ip+2*ng]*(5*angArray[53]+2*angArray[51])-rad[ip+3*ng]*angArray[115];
         bas[13] = -rad[ip+ng]*9*angArray[18]+rad[ip+2*ng]*(3*angArray[54]+3*angArray[52])-rad[ip+3*ng]*angArray[116];
         bas[14] = -rad[ip+ng]*4*angArray[19]+rad[ip+2*ng]*(4*angArray[53]+angArray[55])-rad[ip+3*ng]*angArray[117];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[58]-rad[ip+3*ng]*angArray[127];
         bas[1] = rad[ip+2*ng]*3*angArray[60]-rad[ip+3*ng]*angArray[131];
         bas[2] = -rad[ip+ng]*angArray[20]+rad[ip+2*ng]*(angArray[61]+angArray[59])-rad[ip+3*ng]*angArray[132];
         bas[3] = -rad[ip+ng]*2*angArray[22]+rad[ip+2*ng]*5*angArray[63]-rad[ip+3*ng]*angArray[136];
         bas[4] = -rad[ip+ng]*3*angArray[21]+rad[ip+2*ng]*(3*angArray[64]+angArray[62])-rad[ip+3*ng]*angArray[137];
         bas[5] = -rad[ip+ng]*2*angArray[22]+rad[ip+2*ng]*(angArray[65]+2*angArray[63])-rad[ip+3*ng]*angArray[138];
         bas[6] = -rad[ip+ng]*6*angArray[24]+rad[ip+2*ng]*7*angArray[67]-rad[ip+3*ng]*angArray[142];
         bas[7] = rad[ip]*2*angArray[4]-rad[ip+ng]*(5*angArray[23]+2*angArray[25])+rad[ip+2*ng]*(5*angArray[68]+angArray[66])-rad[ip+3*ng]*angArray[143];
         bas[8] = -rad[ip+ng]*6*angArray[24]+rad[ip+2*ng]*(3*angArray[69]+2*angArray[67])-rad[ip+3*ng]*angArray[144];
         bas[9] = -rad[ip+ng]*3*angArray[25]+rad[ip+2*ng]*(angArray[70]+3*angArray[68])-rad[ip+3*ng]*angArray[145];
         bas[10] = -rad[ip+ng]*12*angArray[27]+rad[ip+2*ng]*9*angArray[72]-rad[ip+3*ng]*angArray[149];
         bas[11] = rad[ip]*6*angArray[5]-rad[ip+ng]*(7*angArray[26]+6*angArray[28])+rad[ip+2*ng]*(7*angArray[73]+angArray[71])-rad[ip+3*ng]*angArray[150];
         bas[12] = rad[ip]*4*angArray[6]-rad[ip+ng]*(10*angArray[27]+2*angArray[29])+rad[ip+2*ng]*(5*angArray[74]+2*angArray[72])-rad[ip+3*ng]*angArray[151];
         bas[13] = -rad[ip+ng]*9*angArray[28]+rad[ip+2*ng]*(3*angArray[75]+3*angArray[73])-rad[ip+3*ng]*angArray[152];
         bas[14] = -rad[ip+ng]*4*angArray[29]+rad[ip+2*ng]*(4*angArray[74]+angArray[76])-rad[ip+3*ng]*angArray[153];
         bas[15] = -rad[ip+ng]*20*angArray[31]+rad[ip+2*ng]*11*angArray[78]-rad[ip+3*ng]*angArray[157];
         bas[16] = rad[ip]*12*angArray[7]-rad[ip+ng]*(9*angArray[30]+12*angArray[32])+rad[ip+2*ng]*(9*angArray[79]+angArray[77])-rad[ip+3*ng]*angArray[158];
         bas[17] = rad[ip]*12*angArray[8]-rad[ip+ng]*(14*angArray[31]+6*angArray[33])+rad[ip+2*ng]*(7*angArray[80]+2*angArray[78])-rad[ip+3*ng]*angArray[159];
         bas[18] = rad[ip]*6*angArray[9]-rad[ip+ng]*(15*angArray[32]+2*angArray[34])+rad[ip+2*ng]*(5*angArray[81]+3*angArray[79])-rad[ip+3*ng]*angArray[160];
         bas[19] = -rad[ip+ng]*12*angArray[33]+rad[ip+2*ng]*(3*angArray[82]+4*angArray[80])-rad[ip+3*ng]*angArray[161];
         bas[20] = -rad[ip+ng]*5*angArray[34]+rad[ip+2*ng]*(5*angArray[81]+angArray[83])-rad[ip+3*ng]*angArray[162];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[86]-rad[ip+3*ng]*angArray[172];
         bas[1] = rad[ip+2*ng]*3*angArray[88]-rad[ip+3*ng]*angArray[176];
         bas[2] = -rad[ip+ng]*angArray[35]+rad[ip+2*ng]*(angArray[89]+angArray[87])-rad[ip+3*ng]*angArray[177];
         bas[3] = -rad[ip+ng]*2*angArray[37]+rad[ip+2*ng]*5*angArray[91]-rad[ip+3*ng]*angArray[181];
         bas[4] = -rad[ip+ng]*3*angArray[36]+rad[ip+2*ng]*(3*angArray[92]+angArray[90])-rad[ip+3*ng]*angArray[182];
         bas[5] = -rad[ip+ng]*2*angArray[37]+rad[ip+2*ng]*(angArray[93]+2*angArray[91])-rad[ip+3*ng]*angArray[183];
         bas[6] = -rad[ip+ng]*6*angArray[39]+rad[ip+2*ng]*7*angArray[95]-rad[ip+3*ng]*angArray[187];
         bas[7] = rad[ip]*2*angArray[10]-rad[ip+ng]*(5*angArray[38]+2*angArray[40])+rad[ip+2*ng]*(5*angArray[96]+angArray[94])-rad[ip+3*ng]*angArray[188];
         bas[8] = -rad[ip+ng]*6*angArray[39]+rad[ip+2*ng]*(3*angArray[97]+2*angArray[95])-rad[ip+3*ng]*angArray[189];
         bas[9] = -rad[ip+ng]*3*angArray[40]+rad[ip+2*ng]*(angArray[98]+3*angArray[96])-rad[ip+3*ng]*angArray[190];
         bas[10] = -rad[ip+ng]*12*angArray[42]+rad[ip+2*ng]*9*angArray[100]-rad[ip+3*ng]*angArray[194];
         bas[11] = rad[ip]*6*angArray[11]-rad[ip+ng]*(7*angArray[41]+6*angArray[43])+rad[ip+2*ng]*(7*angArray[101]+angArray[99])-rad[ip+3*ng]*angArray[195];
         bas[12] = rad[ip]*4*angArray[12]-rad[ip+ng]*(10*angArray[42]+2*angArray[44])+rad[ip+2*ng]*(5*angArray[102]+2*angArray[100])-rad[ip+3*ng]*angArray[196];
         bas[13] = -rad[ip+ng]*9*angArray[43]+rad[ip+2*ng]*(3*angArray[103]+3*angArray[101])-rad[ip+3*ng]*angArray[197];
         bas[14] = -rad[ip+ng]*4*angArray[44]+rad[ip+2*ng]*(4*angArray[102]+angArray[104])-rad[ip+3*ng]*angArray[198];
         bas[15] = -rad[ip+ng]*20*angArray[46]+rad[ip+2*ng]*11*angArray[106]-rad[ip+3*ng]*angArray[202];
         bas[16] = rad[ip]*12*angArray[13]-rad[ip+ng]*(9*angArray[45]+12*angArray[47])+rad[ip+2*ng]*(9*angArray[107]+angArray[105])-rad[ip+3*ng]*angArray[203];
         bas[17] = rad[ip]*12*angArray[14]-rad[ip+ng]*(14*angArray[46]+6*angArray[48])+rad[ip+2*ng]*(7*angArray[108]+2*angArray[106])-rad[ip+3*ng]*angArray[204];
         bas[18] = rad[ip]*6*angArray[15]-rad[ip+ng]*(15*angArray[47]+2*angArray[49])+rad[ip+2*ng]*(5*angArray[109]+3*angArray[107])-rad[ip+3*ng]*angArray[205];
         bas[19] = -rad[ip+ng]*12*angArray[48]+rad[ip+2*ng]*(3*angArray[110]+4*angArray[108])-rad[ip+3*ng]*angArray[206];
         bas[20] = -rad[ip+ng]*5*angArray[49]+rad[ip+2*ng]*(5*angArray[109]+angArray[111])-rad[ip+3*ng]*angArray[207];
         bas[21] = -rad[ip+ng]*30*angArray[51]+rad[ip+2*ng]*13*angArray[113]-rad[ip+3*ng]*angArray[211];
         bas[22] = rad[ip]*20*angArray[16]-rad[ip+ng]*(11*angArray[50]+20*angArray[52])+rad[ip+2*ng]*(11*angArray[114]+angArray[112])-rad[ip+3*ng]*angArray[212];
         bas[23] = rad[ip]*24*angArray[17]-rad[ip+ng]*(12*angArray[53]+18*angArray[51])+rad[ip+2*ng]*(2*angArray[113]+9*angArray[115])-rad[ip+3*ng]*angArray[213];
         bas[24] = rad[ip]*18*angArray[18]-rad[ip+ng]*(21*angArray[52]+6*angArray[54])+rad[ip+2*ng]*(7*angArray[116]+3*angArray[114])-rad[ip+3*ng]*angArray[214];
         bas[25] = rad[ip]*8*angArray[19]-rad[ip+ng]*(20*angArray[53]+2*angArray[55])+rad[ip+2*ng]*(5*angArray[117]+4*angArray[115])-rad[ip+3*ng]*angArray[215];
         bas[26] = -rad[ip+ng]*15*angArray[54]+rad[ip+2*ng]*(3*angArray[118]+5*angArray[116])-rad[ip+3*ng]*angArray[216];
         bas[27] = -rad[ip+ng]*6*angArray[55]+rad[ip+2*ng]*(6*angArray[117]+angArray[119])-rad[ip+3*ng]*angArray[217];
      }

   }


   // now we do derivatives for the given basis set to XZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[1]-rad[ip+3*ng]*angArray[15];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[9]+angArray[4])-rad[ip+3*ng]*angArray[25];
         bas[1] = rad[ip+2*ng]*angArray[5]-rad[ip+3*ng]*angArray[28];
         bas[2] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*angArray[29];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(2*angArray[15]+angArray[10])-rad[ip+3*ng]*angArray[40];
         bas[1] = -rad[ip+ng]*angArray[2]+rad[ip+2*ng]*(angArray[18]+angArray[11])-rad[ip+3*ng]*angArray[43];
         bas[2] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*(angArray[19]+3*angArray[12])-rad[ip+3*ng]*angArray[44];
         bas[3] = rad[ip+2*ng]*angArray[13]-rad[ip+3*ng]*angArray[47];
         bas[4] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[48];
         bas[5] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*5*angArray[15]-rad[ip+3*ng]*angArray[49];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*(3*angArray[25]+angArray[20])-rad[ip+3*ng]*angArray[61];
         bas[1] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*(2*angArray[28]+angArray[21])-rad[ip+3*ng]*angArray[64];
         bas[2] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(2*angArray[29]+3*angArray[22])-rad[ip+3*ng]*angArray[65];
         bas[3] = -rad[ip+ng]*angArray[7]+rad[ip+2*ng]*(angArray[32]+angArray[23])-rad[ip+3*ng]*angArray[68];
         bas[4] = -rad[ip+ng]*3*angArray[8]+rad[ip+2*ng]*(angArray[33]+3*angArray[24])-rad[ip+3*ng]*angArray[69];
         bas[5] = rad[ip]*2*angArray[0]-rad[ip+ng]*(5*angArray[9]+2*angArray[4])+rad[ip+2*ng]*(5*angArray[25]+angArray[34])-rad[ip+3*ng]*angArray[70];
         bas[6] = rad[ip+2*ng]*angArray[26]-rad[ip+3*ng]*angArray[73];
         bas[7] = rad[ip+2*ng]*3*angArray[27]-rad[ip+3*ng]*angArray[74];
         bas[8] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*5*angArray[28]-rad[ip+3*ng]*angArray[75];
         bas[9] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*7*angArray[29]-rad[ip+3*ng]*angArray[76];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*4*angArray[10]+rad[ip+2*ng]*(4*angArray[40]+angArray[35])-rad[ip+3*ng]*angArray[89];
         bas[1] = -rad[ip+ng]*3*angArray[11]+rad[ip+2*ng]*(3*angArray[43]+angArray[36])-rad[ip+3*ng]*angArray[92];
         bas[2] = -rad[ip+ng]*9*angArray[12]+rad[ip+2*ng]*(3*angArray[44]+3*angArray[37])-rad[ip+3*ng]*angArray[93];
         bas[3] = -rad[ip+ng]*2*angArray[13]+rad[ip+2*ng]*(2*angArray[47]+angArray[38])-rad[ip+3*ng]*angArray[96];
         bas[4] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*(2*angArray[48]+3*angArray[39])-rad[ip+3*ng]*angArray[97];
         bas[5] = rad[ip]*4*angArray[1]-rad[ip+ng]*(10*angArray[15]+2*angArray[10])+rad[ip+2*ng]*(2*angArray[49]+5*angArray[40])-rad[ip+3*ng]*angArray[98];
         bas[6] = -rad[ip+ng]*angArray[16]+rad[ip+2*ng]*(angArray[52]+angArray[41])-rad[ip+3*ng]*angArray[101];
         bas[7] = -rad[ip+ng]*3*angArray[17]+rad[ip+2*ng]*(angArray[53]+3*angArray[42])-rad[ip+3*ng]*angArray[102];
         bas[8] = rad[ip]*2*angArray[2]-rad[ip+ng]*(5*angArray[18]+2*angArray[11])+rad[ip+2*ng]*(5*angArray[43]+angArray[54])-rad[ip+3*ng]*angArray[103];
         bas[9] = rad[ip]*6*angArray[3]-rad[ip+ng]*(6*angArray[12]+7*angArray[19])+rad[ip+2*ng]*(angArray[55]+7*angArray[44])-rad[ip+3*ng]*angArray[104];
         bas[10] = rad[ip+2*ng]*angArray[45]-rad[ip+3*ng]*angArray[107];
         bas[11] = rad[ip+2*ng]*3*angArray[46]-rad[ip+3*ng]*angArray[108];
         bas[12] = -rad[ip+ng]*2*angArray[13]+rad[ip+2*ng]*5*angArray[47]-rad[ip+3*ng]*angArray[109];
         bas[13] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*7*angArray[48]-rad[ip+3*ng]*angArray[110];
         bas[14] = -rad[ip+ng]*12*angArray[15]+rad[ip+2*ng]*9*angArray[49]-rad[ip+3*ng]*angArray[111];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*5*angArray[20]+rad[ip+2*ng]*(angArray[56]+5*angArray[61])-rad[ip+3*ng]*angArray[125];
         bas[1] = -rad[ip+ng]*4*angArray[21]+rad[ip+2*ng]*(4*angArray[64]+angArray[57])-rad[ip+3*ng]*angArray[128];
         bas[2] = -rad[ip+ng]*12*angArray[22]+rad[ip+2*ng]*(4*angArray[65]+3*angArray[58])-rad[ip+3*ng]*angArray[129];
         bas[3] = -rad[ip+ng]*3*angArray[23]+rad[ip+2*ng]*(3*angArray[68]+angArray[59])-rad[ip+3*ng]*angArray[132];
         bas[4] = -rad[ip+ng]*9*angArray[24]+rad[ip+2*ng]*(3*angArray[69]+3*angArray[60])-rad[ip+3*ng]*angArray[133];
         bas[5] = rad[ip]*6*angArray[4]-rad[ip+ng]*(15*angArray[25]+2*angArray[20])+rad[ip+2*ng]*(5*angArray[61]+3*angArray[70])-rad[ip+3*ng]*angArray[134];
         bas[6] = -rad[ip+ng]*2*angArray[26]+rad[ip+2*ng]*(2*angArray[73]+angArray[62])-rad[ip+3*ng]*angArray[137];
         bas[7] = -rad[ip+ng]*6*angArray[27]+rad[ip+2*ng]*(2*angArray[74]+3*angArray[63])-rad[ip+3*ng]*angArray[138];
         bas[8] = rad[ip]*4*angArray[5]-rad[ip+ng]*(10*angArray[28]+2*angArray[21])+rad[ip+2*ng]*(2*angArray[75]+5*angArray[64])-rad[ip+3*ng]*angArray[139];
         bas[9] = rad[ip]*12*angArray[6]-rad[ip+ng]*(6*angArray[22]+14*angArray[29])+rad[ip+2*ng]*(2*angArray[76]+7*angArray[65])-rad[ip+3*ng]*angArray[140];
         bas[10] = -rad[ip+ng]*angArray[30]+rad[ip+2*ng]*(angArray[79]+angArray[66])-rad[ip+3*ng]*angArray[143];
         bas[11] = -rad[ip+ng]*3*angArray[31]+rad[ip+2*ng]*(angArray[80]+3*angArray[67])-rad[ip+3*ng]*angArray[144];
         bas[12] = rad[ip]*2*angArray[7]-rad[ip+ng]*(5*angArray[32]+2*angArray[23])+rad[ip+2*ng]*(5*angArray[68]+angArray[81])-rad[ip+3*ng]*angArray[145];
         bas[13] = rad[ip]*6*angArray[8]-rad[ip+ng]*(6*angArray[24]+7*angArray[33])+rad[ip+2*ng]*(angArray[82]+7*angArray[69])-rad[ip+3*ng]*angArray[146];
         bas[14] = rad[ip]*12*angArray[9]-rad[ip+ng]*(12*angArray[25]+9*angArray[34])+rad[ip+2*ng]*(9*angArray[70]+angArray[83])-rad[ip+3*ng]*angArray[147];
         bas[15] = rad[ip+2*ng]*angArray[71]-rad[ip+3*ng]*angArray[150];
         bas[16] = rad[ip+2*ng]*3*angArray[72]-rad[ip+3*ng]*angArray[151];
         bas[17] = -rad[ip+ng]*2*angArray[26]+rad[ip+2*ng]*5*angArray[73]-rad[ip+3*ng]*angArray[152];
         bas[18] = -rad[ip+ng]*6*angArray[27]+rad[ip+2*ng]*7*angArray[74]-rad[ip+3*ng]*angArray[153];
         bas[19] = -rad[ip+ng]*12*angArray[28]+rad[ip+2*ng]*9*angArray[75]-rad[ip+3*ng]*angArray[154];
         bas[20] = -rad[ip+ng]*20*angArray[29]+rad[ip+2*ng]*11*angArray[76]-rad[ip+3*ng]*angArray[155];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[35]+rad[ip+2*ng]*(6*angArray[89]+angArray[84])-rad[ip+3*ng]*angArray[170];
         bas[1] = -rad[ip+ng]*5*angArray[36]+rad[ip+2*ng]*(angArray[85]+5*angArray[92])-rad[ip+3*ng]*angArray[173];
         bas[2] = -rad[ip+ng]*15*angArray[37]+rad[ip+2*ng]*(5*angArray[93]+3*angArray[86])-rad[ip+3*ng]*angArray[174];
         bas[3] = -rad[ip+ng]*4*angArray[38]+rad[ip+2*ng]*(4*angArray[96]+angArray[87])-rad[ip+3*ng]*angArray[177];
         bas[4] = -rad[ip+ng]*12*angArray[39]+rad[ip+2*ng]*(4*angArray[97]+3*angArray[88])-rad[ip+3*ng]*angArray[178];
         bas[5] = rad[ip]*8*angArray[10]-rad[ip+ng]*(20*angArray[40]+2*angArray[35])+rad[ip+2*ng]*(4*angArray[98]+5*angArray[89])-rad[ip+3*ng]*angArray[179];
         bas[6] = -rad[ip+ng]*3*angArray[41]+rad[ip+2*ng]*(3*angArray[101]+angArray[90])-rad[ip+3*ng]*angArray[182];
         bas[7] = -rad[ip+ng]*9*angArray[42]+rad[ip+2*ng]*(3*angArray[102]+3*angArray[91])-rad[ip+3*ng]*angArray[183];
         bas[8] = rad[ip]*6*angArray[11]-rad[ip+ng]*(15*angArray[43]+2*angArray[36])+rad[ip+2*ng]*(5*angArray[92]+3*angArray[103])-rad[ip+3*ng]*angArray[184];
         bas[9] = rad[ip]*18*angArray[12]-rad[ip+ng]*(6*angArray[37]+21*angArray[44])+rad[ip+2*ng]*(3*angArray[104]+7*angArray[93])-rad[ip+3*ng]*angArray[185];
         bas[10] = -rad[ip+ng]*2*angArray[45]+rad[ip+2*ng]*(2*angArray[107]+angArray[94])-rad[ip+3*ng]*angArray[188];
         bas[11] = -rad[ip+ng]*6*angArray[46]+rad[ip+2*ng]*(2*angArray[108]+3*angArray[95])-rad[ip+3*ng]*angArray[189];
         bas[12] = rad[ip]*4*angArray[13]-rad[ip+ng]*(10*angArray[47]+2*angArray[38])+rad[ip+2*ng]*(2*angArray[109]+5*angArray[96])-rad[ip+3*ng]*angArray[190];
         bas[13] = rad[ip]*12*angArray[14]-rad[ip+ng]*(6*angArray[39]+14*angArray[48])+rad[ip+2*ng]*(2*angArray[110]+7*angArray[97])-rad[ip+3*ng]*angArray[191];
         bas[14] = rad[ip]*24*angArray[15]-rad[ip+ng]*(12*angArray[40]+18*angArray[49])+rad[ip+2*ng]*(9*angArray[98]+2*angArray[111])-rad[ip+3*ng]*angArray[192];
         bas[15] = -rad[ip+ng]*angArray[50]+rad[ip+2*ng]*(angArray[114]+angArray[99])-rad[ip+3*ng]*angArray[195];
         bas[16] = -rad[ip+ng]*3*angArray[51]+rad[ip+2*ng]*(angArray[115]+3*angArray[100])-rad[ip+3*ng]*angArray[196];
         bas[17] = rad[ip]*2*angArray[16]-rad[ip+ng]*(5*angArray[52]+2*angArray[41])+rad[ip+2*ng]*(5*angArray[101]+angArray[116])-rad[ip+3*ng]*angArray[197];
         bas[18] = rad[ip]*6*angArray[17]-rad[ip+ng]*(6*angArray[42]+7*angArray[53])+rad[ip+2*ng]*(angArray[117]+7*angArray[102])-rad[ip+3*ng]*angArray[198];
         bas[19] = rad[ip]*12*angArray[18]-rad[ip+ng]*(12*angArray[43]+9*angArray[54])+rad[ip+2*ng]*(9*angArray[103]+angArray[118])-rad[ip+3*ng]*angArray[199];
         bas[20] = rad[ip]*20*angArray[19]-rad[ip+ng]*(11*angArray[55]+20*angArray[44])+rad[ip+2*ng]*(angArray[119]+11*angArray[104])-rad[ip+3*ng]*angArray[200];
         bas[21] = rad[ip+2*ng]*angArray[105]-rad[ip+3*ng]*angArray[203];
         bas[22] = rad[ip+2*ng]*3*angArray[106]-rad[ip+3*ng]*angArray[204];
         bas[23] = -rad[ip+ng]*2*angArray[45]+rad[ip+2*ng]*5*angArray[107]-rad[ip+3*ng]*angArray[205];
         bas[24] = -rad[ip+ng]*6*angArray[46]+rad[ip+2*ng]*7*angArray[108]-rad[ip+3*ng]*angArray[206];
         bas[25] = -rad[ip+ng]*12*angArray[47]+rad[ip+2*ng]*9*angArray[109]-rad[ip+3*ng]*angArray[207];
         bas[26] = -rad[ip+ng]*20*angArray[48]+rad[ip+2*ng]*11*angArray[110]-rad[ip+3*ng]*angArray[208];
         bas[27] = -rad[ip+ng]*30*angArray[49]+rad[ip+2*ng]*13*angArray[111]-rad[ip+3*ng]*angArray[209];
      }

   }


   // now we do derivatives for the given basis set to YZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[2]-rad[ip+3*ng]*angArray[18];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[5]-rad[ip+3*ng]*angArray[28];
         bas[1] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[9]+angArray[7])-rad[ip+3*ng]*angArray[32];
         bas[2] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*angArray[33];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[11]-rad[ip+3*ng]*angArray[43];
         bas[1] = -rad[ip+ng]*angArray[1]+rad[ip+2*ng]*(angArray[15]+angArray[13])-rad[ip+3*ng]*angArray[47];
         bas[2] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[48];
         bas[3] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(2*angArray[18]+angArray[16])-rad[ip+3*ng]*angArray[52];
         bas[4] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*(angArray[19]+3*angArray[17])-rad[ip+3*ng]*angArray[53];
         bas[5] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*5*angArray[18]-rad[ip+3*ng]*angArray[54];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[21]-rad[ip+3*ng]*angArray[64];
         bas[1] = -rad[ip+ng]*angArray[4]+rad[ip+2*ng]*(angArray[25]+angArray[23])-rad[ip+3*ng]*angArray[68];
         bas[2] = rad[ip+2*ng]*3*angArray[24]-rad[ip+3*ng]*angArray[69];
         bas[3] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*(2*angArray[28]+angArray[26])-rad[ip+3*ng]*angArray[73];
         bas[4] = -rad[ip+ng]*3*angArray[6]+rad[ip+2*ng]*(angArray[29]+3*angArray[27])-rad[ip+3*ng]*angArray[74];
         bas[5] = -rad[ip+ng]*2*angArray[5]+rad[ip+2*ng]*5*angArray[28]-rad[ip+3*ng]*angArray[75];
         bas[6] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*(3*angArray[32]+angArray[30])-rad[ip+3*ng]*angArray[79];
         bas[7] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(2*angArray[33]+3*angArray[31])-rad[ip+3*ng]*angArray[80];
         bas[8] = rad[ip]*2*angArray[0]-rad[ip+ng]*(5*angArray[9]+2*angArray[7])+rad[ip+2*ng]*(5*angArray[32]+angArray[34])-rad[ip+3*ng]*angArray[81];
         bas[9] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*7*angArray[33]-rad[ip+3*ng]*angArray[82];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[36]-rad[ip+3*ng]*angArray[92];
         bas[1] = -rad[ip+ng]*angArray[10]+rad[ip+2*ng]*(angArray[40]+angArray[38])-rad[ip+3*ng]*angArray[96];
         bas[2] = rad[ip+2*ng]*3*angArray[39]-rad[ip+3*ng]*angArray[97];
         bas[3] = -rad[ip+ng]*2*angArray[11]+rad[ip+2*ng]*(2*angArray[43]+angArray[41])-rad[ip+3*ng]*angArray[101];
         bas[4] = -rad[ip+ng]*3*angArray[12]+rad[ip+2*ng]*(angArray[44]+3*angArray[42])-rad[ip+3*ng]*angArray[102];
         bas[5] = -rad[ip+ng]*2*angArray[11]+rad[ip+2*ng]*5*angArray[43]-rad[ip+3*ng]*angArray[103];
         bas[6] = -rad[ip+ng]*3*angArray[13]+rad[ip+2*ng]*(3*angArray[47]+angArray[45])-rad[ip+3*ng]*angArray[107];
         bas[7] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*(2*angArray[48]+3*angArray[46])-rad[ip+3*ng]*angArray[108];
         bas[8] = rad[ip]*2*angArray[1]-rad[ip+ng]*(5*angArray[15]+2*angArray[13])+rad[ip+2*ng]*(5*angArray[47]+angArray[49])-rad[ip+3*ng]*angArray[109];
         bas[9] = -rad[ip+ng]*6*angArray[14]+rad[ip+2*ng]*7*angArray[48]-rad[ip+3*ng]*angArray[110];
         bas[10] = -rad[ip+ng]*4*angArray[16]+rad[ip+2*ng]*(4*angArray[52]+angArray[50])-rad[ip+3*ng]*angArray[114];
         bas[11] = -rad[ip+ng]*9*angArray[17]+rad[ip+2*ng]*(3*angArray[53]+3*angArray[51])-rad[ip+3*ng]*angArray[115];
         bas[12] = rad[ip]*4*angArray[2]-rad[ip+ng]*(10*angArray[18]+2*angArray[16])+rad[ip+2*ng]*(2*angArray[54]+5*angArray[52])-rad[ip+3*ng]*angArray[116];
         bas[13] = rad[ip]*6*angArray[3]-rad[ip+ng]*(6*angArray[17]+7*angArray[19])+rad[ip+2*ng]*(angArray[55]+7*angArray[53])-rad[ip+3*ng]*angArray[117];
         bas[14] = -rad[ip+ng]*12*angArray[18]+rad[ip+2*ng]*9*angArray[54]-rad[ip+3*ng]*angArray[118];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[57]-rad[ip+3*ng]*angArray[128];
         bas[1] = -rad[ip+ng]*angArray[20]+rad[ip+2*ng]*(angArray[61]+angArray[59])-rad[ip+3*ng]*angArray[132];
         bas[2] = rad[ip+2*ng]*3*angArray[60]-rad[ip+3*ng]*angArray[133];
         bas[3] = -rad[ip+ng]*2*angArray[21]+rad[ip+2*ng]*(2*angArray[64]+angArray[62])-rad[ip+3*ng]*angArray[137];
         bas[4] = -rad[ip+ng]*3*angArray[22]+rad[ip+2*ng]*(angArray[65]+3*angArray[63])-rad[ip+3*ng]*angArray[138];
         bas[5] = -rad[ip+ng]*2*angArray[21]+rad[ip+2*ng]*5*angArray[64]-rad[ip+3*ng]*angArray[139];
         bas[6] = -rad[ip+ng]*3*angArray[23]+rad[ip+2*ng]*(3*angArray[68]+angArray[66])-rad[ip+3*ng]*angArray[143];
         bas[7] = -rad[ip+ng]*6*angArray[24]+rad[ip+2*ng]*(2*angArray[69]+3*angArray[67])-rad[ip+3*ng]*angArray[144];
         bas[8] = rad[ip]*2*angArray[4]-rad[ip+ng]*(5*angArray[25]+2*angArray[23])+rad[ip+2*ng]*(5*angArray[68]+angArray[70])-rad[ip+3*ng]*angArray[145];
         bas[9] = -rad[ip+ng]*6*angArray[24]+rad[ip+2*ng]*7*angArray[69]-rad[ip+3*ng]*angArray[146];
         bas[10] = -rad[ip+ng]*4*angArray[26]+rad[ip+2*ng]*(4*angArray[73]+angArray[71])-rad[ip+3*ng]*angArray[150];
         bas[11] = -rad[ip+ng]*9*angArray[27]+rad[ip+2*ng]*(3*angArray[74]+3*angArray[72])-rad[ip+3*ng]*angArray[151];
         bas[12] = rad[ip]*4*angArray[5]-rad[ip+ng]*(10*angArray[28]+2*angArray[26])+rad[ip+2*ng]*(2*angArray[75]+5*angArray[73])-rad[ip+3*ng]*angArray[152];
         bas[13] = rad[ip]*6*angArray[6]-rad[ip+ng]*(6*angArray[27]+7*angArray[29])+rad[ip+2*ng]*(angArray[76]+7*angArray[74])-rad[ip+3*ng]*angArray[153];
         bas[14] = -rad[ip+ng]*12*angArray[28]+rad[ip+2*ng]*9*angArray[75]-rad[ip+3*ng]*angArray[154];
         bas[15] = -rad[ip+ng]*5*angArray[30]+rad[ip+2*ng]*(angArray[77]+5*angArray[79])-rad[ip+3*ng]*angArray[158];
         bas[16] = -rad[ip+ng]*12*angArray[31]+rad[ip+2*ng]*(4*angArray[80]+3*angArray[78])-rad[ip+3*ng]*angArray[159];
         bas[17] = rad[ip]*6*angArray[7]-rad[ip+ng]*(15*angArray[32]+2*angArray[30])+rad[ip+2*ng]*(5*angArray[79]+3*angArray[81])-rad[ip+3*ng]*angArray[160];
         bas[18] = rad[ip]*12*angArray[8]-rad[ip+ng]*(6*angArray[31]+14*angArray[33])+rad[ip+2*ng]*(2*angArray[82]+7*angArray[80])-rad[ip+3*ng]*angArray[161];
         bas[19] = rad[ip]*12*angArray[9]-rad[ip+ng]*(12*angArray[32]+9*angArray[34])+rad[ip+2*ng]*(9*angArray[81]+angArray[83])-rad[ip+3*ng]*angArray[162];
         bas[20] = -rad[ip+ng]*20*angArray[33]+rad[ip+2*ng]*11*angArray[82]-rad[ip+3*ng]*angArray[163];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[85]-rad[ip+3*ng]*angArray[173];
         bas[1] = -rad[ip+ng]*angArray[35]+rad[ip+2*ng]*(angArray[89]+angArray[87])-rad[ip+3*ng]*angArray[177];
         bas[2] = rad[ip+2*ng]*3*angArray[88]-rad[ip+3*ng]*angArray[178];
         bas[3] = -rad[ip+ng]*2*angArray[36]+rad[ip+2*ng]*(2*angArray[92]+angArray[90])-rad[ip+3*ng]*angArray[182];
         bas[4] = -rad[ip+ng]*3*angArray[37]+rad[ip+2*ng]*(angArray[93]+3*angArray[91])-rad[ip+3*ng]*angArray[183];
         bas[5] = -rad[ip+ng]*2*angArray[36]+rad[ip+2*ng]*5*angArray[92]-rad[ip+3*ng]*angArray[184];
         bas[6] = -rad[ip+ng]*3*angArray[38]+rad[ip+2*ng]*(3*angArray[96]+angArray[94])-rad[ip+3*ng]*angArray[188];
         bas[7] = -rad[ip+ng]*6*angArray[39]+rad[ip+2*ng]*(2*angArray[97]+3*angArray[95])-rad[ip+3*ng]*angArray[189];
         bas[8] = rad[ip]*2*angArray[10]-rad[ip+ng]*(5*angArray[40]+2*angArray[38])+rad[ip+2*ng]*(5*angArray[96]+angArray[98])-rad[ip+3*ng]*angArray[190];
         bas[9] = -rad[ip+ng]*6*angArray[39]+rad[ip+2*ng]*7*angArray[97]-rad[ip+3*ng]*angArray[191];
         bas[10] = -rad[ip+ng]*4*angArray[41]+rad[ip+2*ng]*(4*angArray[101]+angArray[99])-rad[ip+3*ng]*angArray[195];
         bas[11] = -rad[ip+ng]*9*angArray[42]+rad[ip+2*ng]*(3*angArray[102]+3*angArray[100])-rad[ip+3*ng]*angArray[196];
         bas[12] = rad[ip]*4*angArray[11]-rad[ip+ng]*(10*angArray[43]+2*angArray[41])+rad[ip+2*ng]*(2*angArray[103]+5*angArray[101])-rad[ip+3*ng]*angArray[197];
         bas[13] = rad[ip]*6*angArray[12]-rad[ip+ng]*(6*angArray[42]+7*angArray[44])+rad[ip+2*ng]*(angArray[104]+7*angArray[102])-rad[ip+3*ng]*angArray[198];
         bas[14] = -rad[ip+ng]*12*angArray[43]+rad[ip+2*ng]*9*angArray[103]-rad[ip+3*ng]*angArray[199];
         bas[15] = -rad[ip+ng]*5*angArray[45]+rad[ip+2*ng]*(angArray[105]+5*angArray[107])-rad[ip+3*ng]*angArray[203];
         bas[16] = -rad[ip+ng]*12*angArray[46]+rad[ip+2*ng]*(4*angArray[108]+3*angArray[106])-rad[ip+3*ng]*angArray[204];
         bas[17] = rad[ip]*6*angArray[13]-rad[ip+ng]*(15*angArray[47]+2*angArray[45])+rad[ip+2*ng]*(5*angArray[107]+3*angArray[109])-rad[ip+3*ng]*angArray[205];
         bas[18] = rad[ip]*12*angArray[14]-rad[ip+ng]*(6*angArray[46]+14*angArray[48])+rad[ip+2*ng]*(2*angArray[110]+7*angArray[108])-rad[ip+3*ng]*angArray[206];
         bas[19] = rad[ip]*12*angArray[15]-rad[ip+ng]*(12*angArray[47]+9*angArray[49])+rad[ip+2*ng]*(9*angArray[109]+angArray[111])-rad[ip+3*ng]*angArray[207];
         bas[20] = -rad[ip+ng]*20*angArray[48]+rad[ip+2*ng]*11*angArray[110]-rad[ip+3*ng]*angArray[208];
         bas[21] = -rad[ip+ng]*6*angArray[50]+rad[ip+2*ng]*(6*angArray[114]+angArray[112])-rad[ip+3*ng]*angArray[212];
         bas[22] = -rad[ip+ng]*15*angArray[51]+rad[ip+2*ng]*(5*angArray[115]+3*angArray[113])-rad[ip+3*ng]*angArray[213];
         bas[23] = rad[ip]*8*angArray[16]-rad[ip+ng]*(20*angArray[52]+2*angArray[50])+rad[ip+2*ng]*(4*angArray[116]+5*angArray[114])-rad[ip+3*ng]*angArray[214];
         bas[24] = rad[ip]*18*angArray[17]-rad[ip+ng]*(6*angArray[51]+21*angArray[53])+rad[ip+2*ng]*(3*angArray[117]+7*angArray[115])-rad[ip+3*ng]*angArray[215];
         bas[25] = rad[ip]*24*angArray[18]-rad[ip+ng]*(12*angArray[52]+18*angArray[54])+rad[ip+2*ng]*(9*angArray[116]+2*angArray[118])-rad[ip+3*ng]*angArray[216];
         bas[26] = rad[ip]*20*angArray[19]-rad[ip+ng]*(11*angArray[55]+20*angArray[53])+rad[ip+2*ng]*(angArray[119]+11*angArray[117])-rad[ip+3*ng]*angArray[217];
         bas[27] = -rad[ip+ng]*30*angArray[54]+rad[ip+2*ng]*13*angArray[118]-rad[ip+3*ng]*angArray[218];
      }

   }


   // now we do derivatives for the given basis set to ZZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*angArray[19];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*angArray[29];
         bas[1] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*angArray[33];
         bas[2] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*6*angArray[9]-rad[ip+3*ng]*angArray[34];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[12]-rad[ip+3*ng]*angArray[44];
         bas[1] = rad[ip+2*ng]*3*angArray[14]-rad[ip+3*ng]*angArray[48];
         bas[2] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*6*angArray[15]-rad[ip+3*ng]*angArray[49];
         bas[3] = rad[ip+2*ng]*3*angArray[17]-rad[ip+3*ng]*angArray[53];
         bas[4] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*6*angArray[18]-rad[ip+3*ng]*angArray[54];
         bas[5] = -rad[ip+ng]*12*angArray[3]+rad[ip+2*ng]*9*angArray[19]-rad[ip+3*ng]*angArray[55];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[22]-rad[ip+3*ng]*angArray[65];
         bas[1] = rad[ip+2*ng]*3*angArray[24]-rad[ip+3*ng]*angArray[69];
         bas[2] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*6*angArray[25]-rad[ip+3*ng]*angArray[70];
         bas[3] = rad[ip+2*ng]*3*angArray[27]-rad[ip+3*ng]*angArray[74];
         bas[4] = -rad[ip+ng]*3*angArray[5]+rad[ip+2*ng]*6*angArray[28]-rad[ip+3*ng]*angArray[75];
         bas[5] = -rad[ip+ng]*12*angArray[6]+rad[ip+2*ng]*9*angArray[29]-rad[ip+3*ng]*angArray[76];
         bas[6] = rad[ip+2*ng]*3*angArray[31]-rad[ip+3*ng]*angArray[80];
         bas[7] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*6*angArray[32]-rad[ip+3*ng]*angArray[81];
         bas[8] = -rad[ip+ng]*12*angArray[8]+rad[ip+2*ng]*9*angArray[33]-rad[ip+3*ng]*angArray[82];
         bas[9] = rad[ip]*6*angArray[0]-rad[ip+ng]*27*angArray[9]+rad[ip+2*ng]*12*angArray[34]-rad[ip+3*ng]*angArray[83];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[37]-rad[ip+3*ng]*angArray[93];
         bas[1] = rad[ip+2*ng]*3*angArray[39]-rad[ip+3*ng]*angArray[97];
         bas[2] = -rad[ip+ng]*3*angArray[10]+rad[ip+2*ng]*6*angArray[40]-rad[ip+3*ng]*angArray[98];
         bas[3] = rad[ip+2*ng]*3*angArray[42]-rad[ip+3*ng]*angArray[102];
         bas[4] = -rad[ip+ng]*3*angArray[11]+rad[ip+2*ng]*6*angArray[43]-rad[ip+3*ng]*angArray[103];
         bas[5] = -rad[ip+ng]*12*angArray[12]+rad[ip+2*ng]*9*angArray[44]-rad[ip+3*ng]*angArray[104];
         bas[6] = rad[ip+2*ng]*3*angArray[46]-rad[ip+3*ng]*angArray[108];
         bas[7] = -rad[ip+ng]*3*angArray[13]+rad[ip+2*ng]*6*angArray[47]-rad[ip+3*ng]*angArray[109];
         bas[8] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*9*angArray[48]-rad[ip+3*ng]*angArray[110];
         bas[9] = rad[ip]*6*angArray[1]-rad[ip+ng]*27*angArray[15]+rad[ip+2*ng]*12*angArray[49]-rad[ip+3*ng]*angArray[111];
         bas[10] = rad[ip+2*ng]*3*angArray[51]-rad[ip+3*ng]*angArray[115];
         bas[11] = -rad[ip+ng]*3*angArray[16]+rad[ip+2*ng]*6*angArray[52]-rad[ip+3*ng]*angArray[116];
         bas[12] = -rad[ip+ng]*12*angArray[17]+rad[ip+2*ng]*9*angArray[53]-rad[ip+3*ng]*angArray[117];
         bas[13] = rad[ip]*6*angArray[2]-rad[ip+ng]*27*angArray[18]+rad[ip+2*ng]*12*angArray[54]-rad[ip+3*ng]*angArray[118];
         bas[14] = rad[ip]*24*angArray[3]-rad[ip+ng]*48*angArray[19]+rad[ip+2*ng]*15*angArray[55]-rad[ip+3*ng]*angArray[119];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[58]-rad[ip+3*ng]*angArray[129];
         bas[1] = rad[ip+2*ng]*3*angArray[60]-rad[ip+3*ng]*angArray[133];
         bas[2] = -rad[ip+ng]*3*angArray[20]+rad[ip+2*ng]*6*angArray[61]-rad[ip+3*ng]*angArray[134];
         bas[3] = rad[ip+2*ng]*3*angArray[63]-rad[ip+3*ng]*angArray[138];
         bas[4] = -rad[ip+ng]*3*angArray[21]+rad[ip+2*ng]*6*angArray[64]-rad[ip+3*ng]*angArray[139];
         bas[5] = -rad[ip+ng]*12*angArray[22]+rad[ip+2*ng]*9*angArray[65]-rad[ip+3*ng]*angArray[140];
         bas[6] = rad[ip+2*ng]*3*angArray[67]-rad[ip+3*ng]*angArray[144];
         bas[7] = -rad[ip+ng]*3*angArray[23]+rad[ip+2*ng]*6*angArray[68]-rad[ip+3*ng]*angArray[145];
         bas[8] = -rad[ip+ng]*12*angArray[24]+rad[ip+2*ng]*9*angArray[69]-rad[ip+3*ng]*angArray[146];
         bas[9] = rad[ip]*6*angArray[4]-rad[ip+ng]*27*angArray[25]+rad[ip+2*ng]*12*angArray[70]-rad[ip+3*ng]*angArray[147];
         bas[10] = rad[ip+2*ng]*3*angArray[72]-rad[ip+3*ng]*angArray[151];
         bas[11] = -rad[ip+ng]*3*angArray[26]+rad[ip+2*ng]*6*angArray[73]-rad[ip+3*ng]*angArray[152];
         bas[12] = -rad[ip+ng]*12*angArray[27]+rad[ip+2*ng]*9*angArray[74]-rad[ip+3*ng]*angArray[153];
         bas[13] = rad[ip]*6*angArray[5]-rad[ip+ng]*27*angArray[28]+rad[ip+2*ng]*12*angArray[75]-rad[ip+3*ng]*angArray[154];
         bas[14] = rad[ip]*24*angArray[6]-rad[ip+ng]*48*angArray[29]+rad[ip+2*ng]*15*angArray[76]-rad[ip+3*ng]*angArray[155];
         bas[15] = rad[ip+2*ng]*3*angArray[78]-rad[ip+3*ng]*angArray[159];
         bas[16] = -rad[ip+ng]*3*angArray[30]+rad[ip+2*ng]*6*angArray[79]-rad[ip+3*ng]*angArray[160];
         bas[17] = -rad[ip+ng]*12*angArray[31]+rad[ip+2*ng]*9*angArray[80]-rad[ip+3*ng]*angArray[161];
         bas[18] = rad[ip]*6*angArray[7]-rad[ip+ng]*27*angArray[32]+rad[ip+2*ng]*12*angArray[81]-rad[ip+3*ng]*angArray[162];
         bas[19] = rad[ip]*24*angArray[8]-rad[ip+ng]*48*angArray[33]+rad[ip+2*ng]*15*angArray[82]-rad[ip+3*ng]*angArray[163];
         bas[20] = rad[ip]*60*angArray[9]-rad[ip+ng]*75*angArray[34]+rad[ip+2*ng]*18*angArray[83]-rad[ip+3*ng]*angArray[164];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[86]-rad[ip+3*ng]*angArray[174];
         bas[1] = rad[ip+2*ng]*3*angArray[88]-rad[ip+3*ng]*angArray[178];
         bas[2] = -rad[ip+ng]*3*angArray[35]+rad[ip+2*ng]*6*angArray[89]-rad[ip+3*ng]*angArray[179];
         bas[3] = rad[ip+2*ng]*3*angArray[91]-rad[ip+3*ng]*angArray[183];
         bas[4] = -rad[ip+ng]*3*angArray[36]+rad[ip+2*ng]*6*angArray[92]-rad[ip+3*ng]*angArray[184];
         bas[5] = -rad[ip+ng]*12*angArray[37]+rad[ip+2*ng]*9*angArray[93]-rad[ip+3*ng]*angArray[185];
         bas[6] = rad[ip+2*ng]*3*angArray[95]-rad[ip+3*ng]*angArray[189];
         bas[7] = -rad[ip+ng]*3*angArray[38]+rad[ip+2*ng]*6*angArray[96]-rad[ip+3*ng]*angArray[190];
         bas[8] = -rad[ip+ng]*12*angArray[39]+rad[ip+2*ng]*9*angArray[97]-rad[ip+3*ng]*angArray[191];
         bas[9] = rad[ip]*6*angArray[10]-rad[ip+ng]*27*angArray[40]+rad[ip+2*ng]*12*angArray[98]-rad[ip+3*ng]*angArray[192];
         bas[10] = rad[ip+2*ng]*3*angArray[100]-rad[ip+3*ng]*angArray[196];
         bas[11] = -rad[ip+ng]*3*angArray[41]+rad[ip+2*ng]*6*angArray[101]-rad[ip+3*ng]*angArray[197];
         bas[12] = -rad[ip+ng]*12*angArray[42]+rad[ip+2*ng]*9*angArray[102]-rad[ip+3*ng]*angArray[198];
         bas[13] = rad[ip]*6*angArray[11]-rad[ip+ng]*27*angArray[43]+rad[ip+2*ng]*12*angArray[103]-rad[ip+3*ng]*angArray[199];
         bas[14] = rad[ip]*24*angArray[12]-rad[ip+ng]*48*angArray[44]+rad[ip+2*ng]*15*angArray[104]-rad[ip+3*ng]*angArray[200];
         bas[15] = rad[ip+2*ng]*3*angArray[106]-rad[ip+3*ng]*angArray[204];
         bas[16] = -rad[ip+ng]*3*angArray[45]+rad[ip+2*ng]*6*angArray[107]-rad[ip+3*ng]*angArray[205];
         bas[17] = -rad[ip+ng]*12*angArray[46]+rad[ip+2*ng]*9*angArray[108]-rad[ip+3*ng]*angArray[206];
         bas[18] = rad[ip]*6*angArray[13]-rad[ip+ng]*27*angArray[47]+rad[ip+2*ng]*12*angArray[109]-rad[ip+3*ng]*angArray[207];
         bas[19] = rad[ip]*24*angArray[14]-rad[ip+ng]*48*angArray[48]+rad[ip+2*ng]*15*angArray[110]-rad[ip+3*ng]*angArray[208];
         bas[20] = rad[ip]*60*angArray[15]-rad[ip+ng]*75*angArray[49]+rad[ip+2*ng]*18*angArray[111]-rad[ip+3*ng]*angArray[209];
         bas[21] = rad[ip+2*ng]*3*angArray[113]-rad[ip+3*ng]*angArray[213];
         bas[22] = -rad[ip+ng]*3*angArray[50]+rad[ip+2*ng]*6*angArray[114]-rad[ip+3*ng]*angArray[214];
         bas[23] = -rad[ip+ng]*12*angArray[51]+rad[ip+2*ng]*9*angArray[115]-rad[ip+3*ng]*angArray[215];
         bas[24] = rad[ip]*6*angArray[16]-rad[ip+ng]*27*angArray[52]+rad[ip+2*ng]*12*angArray[116]-rad[ip+3*ng]*angArray[216];
         bas[25] = rad[ip]*24*angArray[17]-rad[ip+ng]*48*angArray[53]+rad[ip+2*ng]*15*angArray[117]-rad[ip+3*ng]*angArray[217];
         bas[26] = rad[ip]*60*angArray[18]-rad[ip+ng]*75*angArray[54]+rad[ip+2*ng]*18*angArray[118]-rad[ip+3*ng]*angArray[218];
         bas[27] = rad[ip]*120*angArray[19]-rad[ip+ng]*108*angArray[55]+rad[ip+2*ng]*21*angArray[119]-rad[ip+3*ng]*angArray[219];
      }

   }


}


