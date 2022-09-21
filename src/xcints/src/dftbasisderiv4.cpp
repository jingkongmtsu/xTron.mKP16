/**
 * This function is used to generate 4 derivatives for basis set 
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

void BatchBasis::dftbasisderiv4(const UInt& ng, const UInt& L, const UInt& nTolCarBas, const Double* ang, const Double* rad, Double* basis) const 
{

   // now we set up the nBas for the computation
   UInt nBas = (L+1)*(L+2)/2;

   // now we do derivatives for the given basis set to XXXX
   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[0]-rad[ip+3*ng]*6*angArray[4]+rad[ip+4*ng]*angArray[20];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*15*angArray[1]-rad[ip+3*ng]*10*angArray[10]+rad[ip+4*ng]*angArray[35];
         bas[1] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*6*angArray[11]+rad[ip+4*ng]*angArray[36];
         bas[2] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*6*angArray[12]+rad[ip+4*ng]*angArray[37];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*12*angArray[0]+rad[ip+2*ng]*39*angArray[4]-rad[ip+3*ng]*14*angArray[20]+rad[ip+4*ng]*angArray[56];
         bas[1] = rad[ip+2*ng]*15*angArray[5]-rad[ip+3*ng]*10*angArray[21]+rad[ip+4*ng]*angArray[57];
         bas[2] = rad[ip+2*ng]*15*angArray[6]-rad[ip+3*ng]*10*angArray[22]+rad[ip+4*ng]*angArray[58];
         bas[3] = rad[ip+2*ng]*3*angArray[7]-rad[ip+3*ng]*6*angArray[23]+rad[ip+4*ng]*angArray[59];
         bas[4] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*6*angArray[24]+rad[ip+4*ng]*angArray[60];
         bas[5] = rad[ip+2*ng]*3*angArray[9]-rad[ip+3*ng]*6*angArray[25]+rad[ip+4*ng]*angArray[61];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*60*angArray[1]+rad[ip+2*ng]*75*angArray[10]-rad[ip+3*ng]*18*angArray[35]+rad[ip+4*ng]*angArray[84];
         bas[1] = -rad[ip+ng]*12*angArray[2]+rad[ip+2*ng]*39*angArray[11]-rad[ip+3*ng]*14*angArray[36]+rad[ip+4*ng]*angArray[85];
         bas[2] = -rad[ip+ng]*12*angArray[3]+rad[ip+2*ng]*39*angArray[12]-rad[ip+3*ng]*14*angArray[37]+rad[ip+4*ng]*angArray[86];
         bas[3] = rad[ip+2*ng]*15*angArray[13]-rad[ip+3*ng]*10*angArray[38]+rad[ip+4*ng]*angArray[87];
         bas[4] = rad[ip+2*ng]*15*angArray[14]-rad[ip+3*ng]*10*angArray[39]+rad[ip+4*ng]*angArray[88];
         bas[5] = rad[ip+2*ng]*15*angArray[15]-rad[ip+3*ng]*10*angArray[40]+rad[ip+4*ng]*angArray[89];
         bas[6] = rad[ip+2*ng]*3*angArray[16]-rad[ip+3*ng]*6*angArray[41]+rad[ip+4*ng]*angArray[90];
         bas[7] = rad[ip+2*ng]*3*angArray[17]-rad[ip+3*ng]*6*angArray[42]+rad[ip+4*ng]*angArray[91];
         bas[8] = rad[ip+2*ng]*3*angArray[18]-rad[ip+3*ng]*6*angArray[43]+rad[ip+4*ng]*angArray[92];
         bas[9] = rad[ip+2*ng]*3*angArray[19]-rad[ip+3*ng]*6*angArray[44]+rad[ip+4*ng]*angArray[93];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*24*angArray[0]-rad[ip+ng]*168*angArray[4]+rad[ip+2*ng]*123*angArray[20]-rad[ip+3*ng]*22*angArray[56]+rad[ip+4*ng]*angArray[120];
         bas[1] = -rad[ip+ng]*60*angArray[5]+rad[ip+2*ng]*75*angArray[21]-rad[ip+3*ng]*18*angArray[57]+rad[ip+4*ng]*angArray[121];
         bas[2] = -rad[ip+ng]*60*angArray[6]+rad[ip+2*ng]*75*angArray[22]-rad[ip+3*ng]*18*angArray[58]+rad[ip+4*ng]*angArray[122];
         bas[3] = -rad[ip+ng]*12*angArray[7]+rad[ip+2*ng]*39*angArray[23]-rad[ip+3*ng]*14*angArray[59]+rad[ip+4*ng]*angArray[123];
         bas[4] = -rad[ip+ng]*12*angArray[8]+rad[ip+2*ng]*39*angArray[24]-rad[ip+3*ng]*14*angArray[60]+rad[ip+4*ng]*angArray[124];
         bas[5] = -rad[ip+ng]*12*angArray[9]+rad[ip+2*ng]*39*angArray[25]-rad[ip+3*ng]*14*angArray[61]+rad[ip+4*ng]*angArray[125];
         bas[6] = rad[ip+2*ng]*15*angArray[26]-rad[ip+3*ng]*10*angArray[62]+rad[ip+4*ng]*angArray[126];
         bas[7] = rad[ip+2*ng]*15*angArray[27]-rad[ip+3*ng]*10*angArray[63]+rad[ip+4*ng]*angArray[127];
         bas[8] = rad[ip+2*ng]*15*angArray[28]-rad[ip+3*ng]*10*angArray[64]+rad[ip+4*ng]*angArray[128];
         bas[9] = rad[ip+2*ng]*15*angArray[29]-rad[ip+3*ng]*10*angArray[65]+rad[ip+4*ng]*angArray[129];
         bas[10] = rad[ip+2*ng]*3*angArray[30]-rad[ip+3*ng]*6*angArray[66]+rad[ip+4*ng]*angArray[130];
         bas[11] = rad[ip+2*ng]*3*angArray[31]-rad[ip+3*ng]*6*angArray[67]+rad[ip+4*ng]*angArray[131];
         bas[12] = rad[ip+2*ng]*3*angArray[32]-rad[ip+3*ng]*6*angArray[68]+rad[ip+4*ng]*angArray[132];
         bas[13] = rad[ip+2*ng]*3*angArray[33]-rad[ip+3*ng]*6*angArray[69]+rad[ip+4*ng]*angArray[133];
         bas[14] = rad[ip+2*ng]*3*angArray[34]-rad[ip+3*ng]*6*angArray[70]+rad[ip+4*ng]*angArray[134];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*120*angArray[1]-rad[ip+ng]*360*angArray[10]+rad[ip+2*ng]*183*angArray[35]-rad[ip+3*ng]*26*angArray[84]+rad[ip+4*ng]*angArray[165];
         bas[1] = rad[ip]*24*angArray[2]-rad[ip+ng]*168*angArray[11]+rad[ip+2*ng]*123*angArray[36]-rad[ip+3*ng]*22*angArray[85]+rad[ip+4*ng]*angArray[166];
         bas[2] = rad[ip]*24*angArray[3]-rad[ip+ng]*168*angArray[12]+rad[ip+2*ng]*123*angArray[37]-rad[ip+3*ng]*22*angArray[86]+rad[ip+4*ng]*angArray[167];
         bas[3] = -rad[ip+ng]*60*angArray[13]+rad[ip+2*ng]*75*angArray[38]-rad[ip+3*ng]*18*angArray[87]+rad[ip+4*ng]*angArray[168];
         bas[4] = -rad[ip+ng]*60*angArray[14]+rad[ip+2*ng]*75*angArray[39]-rad[ip+3*ng]*18*angArray[88]+rad[ip+4*ng]*angArray[169];
         bas[5] = -rad[ip+ng]*60*angArray[15]+rad[ip+2*ng]*75*angArray[40]-rad[ip+3*ng]*18*angArray[89]+rad[ip+4*ng]*angArray[170];
         bas[6] = -rad[ip+ng]*12*angArray[16]+rad[ip+2*ng]*39*angArray[41]-rad[ip+3*ng]*14*angArray[90]+rad[ip+4*ng]*angArray[171];
         bas[7] = -rad[ip+ng]*12*angArray[17]+rad[ip+2*ng]*39*angArray[42]-rad[ip+3*ng]*14*angArray[91]+rad[ip+4*ng]*angArray[172];
         bas[8] = -rad[ip+ng]*12*angArray[18]+rad[ip+2*ng]*39*angArray[43]-rad[ip+3*ng]*14*angArray[92]+rad[ip+4*ng]*angArray[173];
         bas[9] = -rad[ip+ng]*12*angArray[19]+rad[ip+2*ng]*39*angArray[44]-rad[ip+3*ng]*14*angArray[93]+rad[ip+4*ng]*angArray[174];
         bas[10] = rad[ip+2*ng]*15*angArray[45]-rad[ip+3*ng]*10*angArray[94]+rad[ip+4*ng]*angArray[175];
         bas[11] = rad[ip+2*ng]*15*angArray[46]-rad[ip+3*ng]*10*angArray[95]+rad[ip+4*ng]*angArray[176];
         bas[12] = rad[ip+2*ng]*15*angArray[47]-rad[ip+3*ng]*10*angArray[96]+rad[ip+4*ng]*angArray[177];
         bas[13] = rad[ip+2*ng]*15*angArray[48]-rad[ip+3*ng]*10*angArray[97]+rad[ip+4*ng]*angArray[178];
         bas[14] = rad[ip+2*ng]*15*angArray[49]-rad[ip+3*ng]*10*angArray[98]+rad[ip+4*ng]*angArray[179];
         bas[15] = rad[ip+2*ng]*3*angArray[50]-rad[ip+3*ng]*6*angArray[99]+rad[ip+4*ng]*angArray[180];
         bas[16] = rad[ip+2*ng]*3*angArray[51]-rad[ip+3*ng]*6*angArray[100]+rad[ip+4*ng]*angArray[181];
         bas[17] = rad[ip+2*ng]*3*angArray[52]-rad[ip+3*ng]*6*angArray[101]+rad[ip+4*ng]*angArray[182];
         bas[18] = rad[ip+2*ng]*3*angArray[53]-rad[ip+3*ng]*6*angArray[102]+rad[ip+4*ng]*angArray[183];
         bas[19] = rad[ip+2*ng]*3*angArray[54]-rad[ip+3*ng]*6*angArray[103]+rad[ip+4*ng]*angArray[184];
         bas[20] = rad[ip+2*ng]*3*angArray[55]-rad[ip+3*ng]*6*angArray[104]+rad[ip+4*ng]*angArray[185];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip]*360*angArray[4]-rad[ip+ng]*660*angArray[20]+rad[ip+2*ng]*255*angArray[56]-rad[ip+3*ng]*30*angArray[120]+rad[ip+4*ng]*angArray[220];
         bas[1] = rad[ip]*120*angArray[5]-rad[ip+ng]*360*angArray[21]+rad[ip+2*ng]*183*angArray[57]-rad[ip+3*ng]*26*angArray[121]+rad[ip+4*ng]*angArray[221];
         bas[2] = rad[ip]*120*angArray[6]-rad[ip+ng]*360*angArray[22]+rad[ip+2*ng]*183*angArray[58]-rad[ip+3*ng]*26*angArray[122]+rad[ip+4*ng]*angArray[222];
         bas[3] = rad[ip]*24*angArray[7]-rad[ip+ng]*168*angArray[23]+rad[ip+2*ng]*123*angArray[59]-rad[ip+3*ng]*22*angArray[123]+rad[ip+4*ng]*angArray[223];
         bas[4] = rad[ip]*24*angArray[8]-rad[ip+ng]*168*angArray[24]+rad[ip+2*ng]*123*angArray[60]-rad[ip+3*ng]*22*angArray[124]+rad[ip+4*ng]*angArray[224];
         bas[5] = rad[ip]*24*angArray[9]-rad[ip+ng]*168*angArray[25]+rad[ip+2*ng]*123*angArray[61]-rad[ip+3*ng]*22*angArray[125]+rad[ip+4*ng]*angArray[225];
         bas[6] = -rad[ip+ng]*60*angArray[26]+rad[ip+2*ng]*75*angArray[62]-rad[ip+3*ng]*18*angArray[126]+rad[ip+4*ng]*angArray[226];
         bas[7] = -rad[ip+ng]*60*angArray[27]+rad[ip+2*ng]*75*angArray[63]-rad[ip+3*ng]*18*angArray[127]+rad[ip+4*ng]*angArray[227];
         bas[8] = -rad[ip+ng]*60*angArray[28]+rad[ip+2*ng]*75*angArray[64]-rad[ip+3*ng]*18*angArray[128]+rad[ip+4*ng]*angArray[228];
         bas[9] = -rad[ip+ng]*60*angArray[29]+rad[ip+2*ng]*75*angArray[65]-rad[ip+3*ng]*18*angArray[129]+rad[ip+4*ng]*angArray[229];
         bas[10] = -rad[ip+ng]*12*angArray[30]+rad[ip+2*ng]*39*angArray[66]-rad[ip+3*ng]*14*angArray[130]+rad[ip+4*ng]*angArray[230];
         bas[11] = -rad[ip+ng]*12*angArray[31]+rad[ip+2*ng]*39*angArray[67]-rad[ip+3*ng]*14*angArray[131]+rad[ip+4*ng]*angArray[231];
         bas[12] = -rad[ip+ng]*12*angArray[32]+rad[ip+2*ng]*39*angArray[68]-rad[ip+3*ng]*14*angArray[132]+rad[ip+4*ng]*angArray[232];
         bas[13] = -rad[ip+ng]*12*angArray[33]+rad[ip+2*ng]*39*angArray[69]-rad[ip+3*ng]*14*angArray[133]+rad[ip+4*ng]*angArray[233];
         bas[14] = -rad[ip+ng]*12*angArray[34]+rad[ip+2*ng]*39*angArray[70]-rad[ip+3*ng]*14*angArray[134]+rad[ip+4*ng]*angArray[234];
         bas[15] = rad[ip+2*ng]*15*angArray[71]-rad[ip+3*ng]*10*angArray[135]+rad[ip+4*ng]*angArray[235];
         bas[16] = rad[ip+2*ng]*15*angArray[72]-rad[ip+3*ng]*10*angArray[136]+rad[ip+4*ng]*angArray[236];
         bas[17] = rad[ip+2*ng]*15*angArray[73]-rad[ip+3*ng]*10*angArray[137]+rad[ip+4*ng]*angArray[237];
         bas[18] = rad[ip+2*ng]*15*angArray[74]-rad[ip+3*ng]*10*angArray[138]+rad[ip+4*ng]*angArray[238];
         bas[19] = rad[ip+2*ng]*15*angArray[75]-rad[ip+3*ng]*10*angArray[139]+rad[ip+4*ng]*angArray[239];
         bas[20] = rad[ip+2*ng]*15*angArray[76]-rad[ip+3*ng]*10*angArray[140]+rad[ip+4*ng]*angArray[240];
         bas[21] = rad[ip+2*ng]*3*angArray[77]-rad[ip+3*ng]*6*angArray[141]+rad[ip+4*ng]*angArray[241];
         bas[22] = rad[ip+2*ng]*3*angArray[78]-rad[ip+3*ng]*6*angArray[142]+rad[ip+4*ng]*angArray[242];
         bas[23] = rad[ip+2*ng]*3*angArray[79]-rad[ip+3*ng]*6*angArray[143]+rad[ip+4*ng]*angArray[243];
         bas[24] = rad[ip+2*ng]*3*angArray[80]-rad[ip+3*ng]*6*angArray[144]+rad[ip+4*ng]*angArray[244];
         bas[25] = rad[ip+2*ng]*3*angArray[81]-rad[ip+3*ng]*6*angArray[145]+rad[ip+4*ng]*angArray[245];
         bas[26] = rad[ip+2*ng]*3*angArray[82]-rad[ip+3*ng]*6*angArray[146]+rad[ip+4*ng]*angArray[246];
         bas[27] = rad[ip+2*ng]*3*angArray[83]-rad[ip+3*ng]*6*angArray[147]+rad[ip+4*ng]*angArray[247];
      }

   }


   // now we do derivatives for the given basis set to XXXY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[5]+rad[ip+4*ng]*angArray[21];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*6*angArray[11]+rad[ip+4*ng]*angArray[36];
         bas[1] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*(3*angArray[13]+angArray[10])+rad[ip+4*ng]*angArray[38];
         bas[2] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[39];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*12*angArray[5]-rad[ip+3*ng]*9*angArray[21]+rad[ip+4*ng]*angArray[57];
         bas[1] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*(6*angArray[4]+3*angArray[7])-rad[ip+3*ng]*(angArray[20]+6*angArray[23])+rad[ip+4*ng]*angArray[59];
         bas[2] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*6*angArray[24]+rad[ip+4*ng]*angArray[60];
         bas[3] = rad[ip+2*ng]*6*angArray[5]-rad[ip+3*ng]*(2*angArray[21]+3*angArray[26])+rad[ip+4*ng]*angArray[62];
         bas[4] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*(3*angArray[27]+angArray[22])+rad[ip+4*ng]*angArray[63];
         bas[5] = -rad[ip+3*ng]*3*angArray[28]+rad[ip+4*ng]*angArray[64];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*27*angArray[11]-rad[ip+3*ng]*12*angArray[36]+rad[ip+4*ng]*angArray[85];
         bas[1] = -rad[ip+ng]*12*angArray[1]+rad[ip+2*ng]*(12*angArray[13]+9*angArray[10])-rad[ip+3*ng]*(9*angArray[38]+angArray[35])+rad[ip+4*ng]*angArray[87];
         bas[2] = rad[ip+2*ng]*12*angArray[14]-rad[ip+3*ng]*9*angArray[39]+rad[ip+4*ng]*angArray[88];
         bas[3] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*(12*angArray[11]+3*angArray[16])-rad[ip+3*ng]*(2*angArray[36]+6*angArray[41])+rad[ip+4*ng]*angArray[90];
         bas[4] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*(6*angArray[12]+3*angArray[17])-rad[ip+3*ng]*(angArray[37]+6*angArray[42])+rad[ip+4*ng]*angArray[91];
         bas[5] = rad[ip+2*ng]*3*angArray[18]-rad[ip+3*ng]*6*angArray[43]+rad[ip+4*ng]*angArray[92];
         bas[6] = rad[ip+2*ng]*9*angArray[13]-rad[ip+3*ng]*(3*angArray[38]+3*angArray[45])+rad[ip+4*ng]*angArray[94];
         bas[7] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[39]+3*angArray[46])+rad[ip+4*ng]*angArray[95];
         bas[8] = rad[ip+2*ng]*3*angArray[15]-rad[ip+3*ng]*(3*angArray[47]+angArray[40])+rad[ip+4*ng]*angArray[96];
         bas[9] = -rad[ip+3*ng]*3*angArray[48]+rad[ip+4*ng]*angArray[97];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*24*angArray[5]+rad[ip+2*ng]*48*angArray[21]-rad[ip+3*ng]*15*angArray[57]+rad[ip+4*ng]*angArray[121];
         bas[1] = rad[ip]*6*angArray[0]-rad[ip+ng]*(27*angArray[4]+6*angArray[7])+rad[ip+2*ng]*(12*angArray[20]+27*angArray[23])-rad[ip+3*ng]*(12*angArray[59]+angArray[56])+rad[ip+4*ng]*angArray[123];
         bas[2] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*27*angArray[24]-rad[ip+3*ng]*12*angArray[60]+rad[ip+4*ng]*angArray[124];
         bas[3] = -rad[ip+ng]*24*angArray[5]+rad[ip+2*ng]*(12*angArray[26]+18*angArray[21])-rad[ip+3*ng]*(9*angArray[62]+2*angArray[57])+rad[ip+4*ng]*angArray[126];
         bas[4] = -rad[ip+ng]*12*angArray[6]+rad[ip+2*ng]*(12*angArray[27]+9*angArray[22])-rad[ip+3*ng]*(9*angArray[63]+angArray[58])+rad[ip+4*ng]*angArray[127];
         bas[5] = rad[ip+2*ng]*12*angArray[28]-rad[ip+3*ng]*9*angArray[64]+rad[ip+4*ng]*angArray[128];
         bas[6] = -rad[ip+ng]*9*angArray[7]+rad[ip+2*ng]*(18*angArray[23]+3*angArray[30])-rad[ip+3*ng]*(3*angArray[59]+6*angArray[66])+rad[ip+4*ng]*angArray[130];
         bas[7] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(12*angArray[24]+3*angArray[31])-rad[ip+3*ng]*(2*angArray[60]+6*angArray[67])+rad[ip+4*ng]*angArray[131];
         bas[8] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*(6*angArray[25]+3*angArray[32])-rad[ip+3*ng]*(angArray[61]+6*angArray[68])+rad[ip+4*ng]*angArray[132];
         bas[9] = rad[ip+2*ng]*3*angArray[33]-rad[ip+3*ng]*6*angArray[69]+rad[ip+4*ng]*angArray[133];
         bas[10] = rad[ip+2*ng]*12*angArray[26]-rad[ip+3*ng]*(3*angArray[71]+4*angArray[62])+rad[ip+4*ng]*angArray[135];
         bas[11] = rad[ip+2*ng]*9*angArray[27]-rad[ip+3*ng]*(3*angArray[63]+3*angArray[72])+rad[ip+4*ng]*angArray[136];
         bas[12] = rad[ip+2*ng]*6*angArray[28]-rad[ip+3*ng]*(2*angArray[64]+3*angArray[73])+rad[ip+4*ng]*angArray[137];
         bas[13] = rad[ip+2*ng]*3*angArray[29]-rad[ip+3*ng]*(3*angArray[74]+angArray[65])+rad[ip+4*ng]*angArray[138];
         bas[14] = -rad[ip+3*ng]*3*angArray[75]+rad[ip+4*ng]*angArray[139];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*60*angArray[11]+rad[ip+2*ng]*75*angArray[36]-rad[ip+3*ng]*18*angArray[85]+rad[ip+4*ng]*angArray[166];
         bas[1] = rad[ip]*24*angArray[1]-rad[ip+ng]*(48*angArray[10]+24*angArray[13])+rad[ip+2*ng]*(15*angArray[35]+48*angArray[38])-rad[ip+3*ng]*(angArray[84]+15*angArray[87])+rad[ip+4*ng]*angArray[168];
         bas[2] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*48*angArray[39]-rad[ip+3*ng]*15*angArray[88]+rad[ip+4*ng]*angArray[169];
         bas[3] = rad[ip]*12*angArray[2]-rad[ip+ng]*(54*angArray[11]+6*angArray[16])+rad[ip+2*ng]*(24*angArray[36]+27*angArray[41])-rad[ip+3*ng]*(12*angArray[90]+2*angArray[85])+rad[ip+4*ng]*angArray[171];
         bas[4] = rad[ip]*6*angArray[3]-rad[ip+ng]*(27*angArray[12]+6*angArray[17])+rad[ip+2*ng]*(12*angArray[37]+27*angArray[42])-rad[ip+3*ng]*(12*angArray[91]+angArray[86])+rad[ip+4*ng]*angArray[172];
         bas[5] = -rad[ip+ng]*6*angArray[18]+rad[ip+2*ng]*27*angArray[43]-rad[ip+3*ng]*12*angArray[92]+rad[ip+4*ng]*angArray[173];
         bas[6] = -rad[ip+ng]*36*angArray[13]+rad[ip+2*ng]*(12*angArray[45]+27*angArray[38])-rad[ip+3*ng]*(9*angArray[94]+3*angArray[87])+rad[ip+4*ng]*angArray[175];
         bas[7] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*(12*angArray[46]+18*angArray[39])-rad[ip+3*ng]*(9*angArray[95]+2*angArray[88])+rad[ip+4*ng]*angArray[176];
         bas[8] = -rad[ip+ng]*12*angArray[15]+rad[ip+2*ng]*(12*angArray[47]+9*angArray[40])-rad[ip+3*ng]*(9*angArray[96]+angArray[89])+rad[ip+4*ng]*angArray[177];
         bas[9] = rad[ip+2*ng]*12*angArray[48]-rad[ip+3*ng]*9*angArray[97]+rad[ip+4*ng]*angArray[178];
         bas[10] = -rad[ip+ng]*12*angArray[16]+rad[ip+2*ng]*(3*angArray[50]+24*angArray[41])-rad[ip+3*ng]*(4*angArray[90]+6*angArray[99])+rad[ip+4*ng]*angArray[180];
         bas[11] = -rad[ip+ng]*9*angArray[17]+rad[ip+2*ng]*(18*angArray[42]+3*angArray[51])-rad[ip+3*ng]*(3*angArray[91]+6*angArray[100])+rad[ip+4*ng]*angArray[181];
         bas[12] = -rad[ip+ng]*6*angArray[18]+rad[ip+2*ng]*(12*angArray[43]+3*angArray[52])-rad[ip+3*ng]*(2*angArray[92]+6*angArray[101])+rad[ip+4*ng]*angArray[182];
         bas[13] = -rad[ip+ng]*3*angArray[19]+rad[ip+2*ng]*(6*angArray[44]+3*angArray[53])-rad[ip+3*ng]*(angArray[93]+6*angArray[102])+rad[ip+4*ng]*angArray[183];
         bas[14] = rad[ip+2*ng]*3*angArray[54]-rad[ip+3*ng]*6*angArray[103]+rad[ip+4*ng]*angArray[184];
         bas[15] = rad[ip+2*ng]*15*angArray[45]-rad[ip+3*ng]*(3*angArray[105]+5*angArray[94])+rad[ip+4*ng]*angArray[186];
         bas[16] = rad[ip+2*ng]*12*angArray[46]-rad[ip+3*ng]*(3*angArray[106]+4*angArray[95])+rad[ip+4*ng]*angArray[187];
         bas[17] = rad[ip+2*ng]*9*angArray[47]-rad[ip+3*ng]*(3*angArray[96]+3*angArray[107])+rad[ip+4*ng]*angArray[188];
         bas[18] = rad[ip+2*ng]*6*angArray[48]-rad[ip+3*ng]*(2*angArray[97]+3*angArray[108])+rad[ip+4*ng]*angArray[189];
         bas[19] = rad[ip+2*ng]*3*angArray[49]-rad[ip+3*ng]*(3*angArray[109]+angArray[98])+rad[ip+4*ng]*angArray[190];
         bas[20] = -rad[ip+3*ng]*3*angArray[110]+rad[ip+4*ng]*angArray[191];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*120*angArray[21]+rad[ip+2*ng]*108*angArray[57]-rad[ip+3*ng]*21*angArray[121]+rad[ip+4*ng]*angArray[221];
         bas[1] = rad[ip]*60*angArray[4]-rad[ip+ng]*(75*angArray[20]+60*angArray[23])+rad[ip+2*ng]*(75*angArray[59]+18*angArray[56])-rad[ip+3*ng]*(18*angArray[123]+angArray[120])+rad[ip+4*ng]*angArray[223];
         bas[2] = -rad[ip+ng]*60*angArray[24]+rad[ip+2*ng]*75*angArray[60]-rad[ip+3*ng]*18*angArray[124]+rad[ip+4*ng]*angArray[224];
         bas[3] = rad[ip]*48*angArray[5]-rad[ip+ng]*(96*angArray[21]+24*angArray[26])+rad[ip+2*ng]*(48*angArray[62]+30*angArray[57])-rad[ip+3*ng]*(2*angArray[121]+15*angArray[126])+rad[ip+4*ng]*angArray[226];
         bas[4] = rad[ip]*24*angArray[6]-rad[ip+ng]*(48*angArray[22]+24*angArray[27])+rad[ip+2*ng]*(15*angArray[58]+48*angArray[63])-rad[ip+3*ng]*(angArray[122]+15*angArray[127])+rad[ip+4*ng]*angArray[227];
         bas[5] = -rad[ip+ng]*24*angArray[28]+rad[ip+2*ng]*48*angArray[64]-rad[ip+3*ng]*15*angArray[128]+rad[ip+4*ng]*angArray[228];
         bas[6] = rad[ip]*18*angArray[7]-rad[ip+ng]*(81*angArray[23]+6*angArray[30])+rad[ip+2*ng]*(36*angArray[59]+27*angArray[66])-rad[ip+3*ng]*(12*angArray[130]+3*angArray[123])+rad[ip+4*ng]*angArray[230];
         bas[7] = rad[ip]*12*angArray[8]-rad[ip+ng]*(54*angArray[24]+6*angArray[31])+rad[ip+2*ng]*(24*angArray[60]+27*angArray[67])-rad[ip+3*ng]*(12*angArray[131]+2*angArray[124])+rad[ip+4*ng]*angArray[231];
         bas[8] = rad[ip]*6*angArray[9]-rad[ip+ng]*(27*angArray[25]+6*angArray[32])+rad[ip+2*ng]*(12*angArray[61]+27*angArray[68])-rad[ip+3*ng]*(12*angArray[132]+angArray[125])+rad[ip+4*ng]*angArray[232];
         bas[9] = -rad[ip+ng]*6*angArray[33]+rad[ip+2*ng]*27*angArray[69]-rad[ip+3*ng]*12*angArray[133]+rad[ip+4*ng]*angArray[233];
         bas[10] = -rad[ip+ng]*48*angArray[26]+rad[ip+2*ng]*(36*angArray[62]+12*angArray[71])-rad[ip+3*ng]*(9*angArray[135]+4*angArray[126])+rad[ip+4*ng]*angArray[235];
         bas[11] = -rad[ip+ng]*36*angArray[27]+rad[ip+2*ng]*(12*angArray[72]+27*angArray[63])-rad[ip+3*ng]*(9*angArray[136]+3*angArray[127])+rad[ip+4*ng]*angArray[236];
         bas[12] = -rad[ip+ng]*24*angArray[28]+rad[ip+2*ng]*(12*angArray[73]+18*angArray[64])-rad[ip+3*ng]*(9*angArray[137]+2*angArray[128])+rad[ip+4*ng]*angArray[237];
         bas[13] = -rad[ip+ng]*12*angArray[29]+rad[ip+2*ng]*(12*angArray[74]+9*angArray[65])-rad[ip+3*ng]*(9*angArray[138]+angArray[129])+rad[ip+4*ng]*angArray[238];
         bas[14] = rad[ip+2*ng]*12*angArray[75]-rad[ip+3*ng]*9*angArray[139]+rad[ip+4*ng]*angArray[239];
         bas[15] = -rad[ip+ng]*15*angArray[30]+rad[ip+2*ng]*(30*angArray[66]+3*angArray[77])-rad[ip+3*ng]*(5*angArray[130]+6*angArray[141])+rad[ip+4*ng]*angArray[241];
         bas[16] = -rad[ip+ng]*12*angArray[31]+rad[ip+2*ng]*(3*angArray[78]+24*angArray[67])-rad[ip+3*ng]*(4*angArray[131]+6*angArray[142])+rad[ip+4*ng]*angArray[242];
         bas[17] = -rad[ip+ng]*9*angArray[32]+rad[ip+2*ng]*(18*angArray[68]+3*angArray[79])-rad[ip+3*ng]*(3*angArray[132]+6*angArray[143])+rad[ip+4*ng]*angArray[243];
         bas[18] = -rad[ip+ng]*6*angArray[33]+rad[ip+2*ng]*(12*angArray[69]+3*angArray[80])-rad[ip+3*ng]*(2*angArray[133]+6*angArray[144])+rad[ip+4*ng]*angArray[244];
         bas[19] = -rad[ip+ng]*3*angArray[34]+rad[ip+2*ng]*(6*angArray[70]+3*angArray[81])-rad[ip+3*ng]*(angArray[134]+6*angArray[145])+rad[ip+4*ng]*angArray[245];
         bas[20] = rad[ip+2*ng]*3*angArray[82]-rad[ip+3*ng]*6*angArray[146]+rad[ip+4*ng]*angArray[246];
         bas[21] = rad[ip+2*ng]*18*angArray[71]-rad[ip+3*ng]*(3*angArray[148]+6*angArray[135])+rad[ip+4*ng]*angArray[248];
         bas[22] = rad[ip+2*ng]*15*angArray[72]-rad[ip+3*ng]*(3*angArray[149]+5*angArray[136])+rad[ip+4*ng]*angArray[249];
         bas[23] = rad[ip+2*ng]*12*angArray[73]-rad[ip+3*ng]*(3*angArray[150]+4*angArray[137])+rad[ip+4*ng]*angArray[250];
         bas[24] = rad[ip+2*ng]*9*angArray[74]-rad[ip+3*ng]*(3*angArray[138]+3*angArray[151])+rad[ip+4*ng]*angArray[251];
         bas[25] = rad[ip+2*ng]*6*angArray[75]-rad[ip+3*ng]*(2*angArray[139]+3*angArray[152])+rad[ip+4*ng]*angArray[252];
         bas[26] = rad[ip+2*ng]*3*angArray[76]-rad[ip+3*ng]*(3*angArray[153]+angArray[140])+rad[ip+4*ng]*angArray[253];
         bas[27] = -rad[ip+3*ng]*3*angArray[154]+rad[ip+4*ng]*angArray[254];
      }

   }


   // now we do derivatives for the given basis set to XXYY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[0]-rad[ip+3*ng]*(angArray[4]+angArray[7])+rad[ip+4*ng]*angArray[23];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*(angArray[10]+3*angArray[13])+rad[ip+4*ng]*angArray[38];
         bas[1] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*(3*angArray[11]+angArray[16])+rad[ip+4*ng]*angArray[41];
         bas[2] = rad[ip+2*ng]*angArray[3]-rad[ip+3*ng]*(angArray[12]+angArray[17])+rad[ip+4*ng]*angArray[42];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[0]+rad[ip+2*ng]*(5*angArray[4]+2*angArray[7])-rad[ip+3*ng]*(5*angArray[23]+angArray[20])+rad[ip+4*ng]*angArray[59];
         bas[1] = rad[ip+2*ng]*9*angArray[5]-rad[ip+3*ng]*(3*angArray[21]+3*angArray[26])+rad[ip+4*ng]*angArray[62];
         bas[2] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*(angArray[22]+3*angArray[27])+rad[ip+4*ng]*angArray[63];
         bas[3] = -rad[ip+ng]*2*angArray[0]+rad[ip+2*ng]*(5*angArray[7]+2*angArray[4])-rad[ip+3*ng]*(5*angArray[23]+angArray[30])+rad[ip+4*ng]*angArray[66];
         bas[4] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*(3*angArray[24]+angArray[31])+rad[ip+4*ng]*angArray[67];
         bas[5] = rad[ip+2*ng]*angArray[9]-rad[ip+3*ng]*(angArray[25]+angArray[32])+rad[ip+4*ng]*angArray[68];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*(7*angArray[10]+6*angArray[13])-rad[ip+3*ng]*(7*angArray[38]+angArray[35])+rad[ip+4*ng]*angArray[87];
         bas[1] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*(15*angArray[11]+2*angArray[16])-rad[ip+3*ng]*(5*angArray[41]+3*angArray[36])+rad[ip+4*ng]*angArray[90];
         bas[2] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(5*angArray[12]+2*angArray[17])-rad[ip+3*ng]*(5*angArray[42]+angArray[37])+rad[ip+4*ng]*angArray[91];
         bas[3] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*(2*angArray[10]+15*angArray[13])-rad[ip+3*ng]*(5*angArray[38]+3*angArray[45])+rad[ip+4*ng]*angArray[94];
         bas[4] = rad[ip+2*ng]*9*angArray[14]-rad[ip+3*ng]*(3*angArray[39]+3*angArray[46])+rad[ip+4*ng]*angArray[95];
         bas[5] = rad[ip+2*ng]*3*angArray[15]-rad[ip+3*ng]*(angArray[40]+3*angArray[47])+rad[ip+4*ng]*angArray[96];
         bas[6] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*(7*angArray[16]+6*angArray[11])-rad[ip+3*ng]*(7*angArray[41]+angArray[50])+rad[ip+4*ng]*angArray[99];
         bas[7] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(5*angArray[17]+2*angArray[12])-rad[ip+3*ng]*(5*angArray[42]+angArray[51])+rad[ip+4*ng]*angArray[100];
         bas[8] = rad[ip+2*ng]*3*angArray[18]-rad[ip+3*ng]*(3*angArray[43]+angArray[52])+rad[ip+4*ng]*angArray[101];
         bas[9] = rad[ip+2*ng]*angArray[19]-rad[ip+3*ng]*(angArray[44]+angArray[53])+rad[ip+4*ng]*angArray[102];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*12*angArray[4]+rad[ip+2*ng]*(9*angArray[20]+12*angArray[23])-rad[ip+3*ng]*(angArray[56]+9*angArray[59])+rad[ip+4*ng]*angArray[123];
         bas[1] = -rad[ip+ng]*18*angArray[5]+rad[ip+2*ng]*(21*angArray[21]+6*angArray[26])-rad[ip+3*ng]*(7*angArray[62]+3*angArray[57])+rad[ip+4*ng]*angArray[126];
         bas[2] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(7*angArray[22]+6*angArray[27])-rad[ip+3*ng]*(7*angArray[63]+angArray[58])+rad[ip+4*ng]*angArray[127];
         bas[3] = rad[ip]*4*angArray[0]-rad[ip+ng]*(10*angArray[4]+10*angArray[7])+rad[ip+2*ng]*(25*angArray[23]+2*angArray[20]+2*angArray[30])-rad[ip+3*ng]*(5*angArray[66]+5*angArray[59])+rad[ip+4*ng]*angArray[130];
         bas[4] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(15*angArray[24]+2*angArray[31])-rad[ip+3*ng]*(5*angArray[67]+3*angArray[60])+rad[ip+4*ng]*angArray[131];
         bas[5] = -rad[ip+ng]*2*angArray[9]+rad[ip+2*ng]*(5*angArray[25]+2*angArray[32])-rad[ip+3*ng]*(5*angArray[68]+angArray[61])+rad[ip+4*ng]*angArray[132];
         bas[6] = -rad[ip+ng]*18*angArray[5]+rad[ip+2*ng]*(21*angArray[26]+6*angArray[21])-rad[ip+3*ng]*(7*angArray[62]+3*angArray[71])+rad[ip+4*ng]*angArray[135];
         bas[7] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(2*angArray[22]+15*angArray[27])-rad[ip+3*ng]*(5*angArray[63]+3*angArray[72])+rad[ip+4*ng]*angArray[136];
         bas[8] = rad[ip+2*ng]*9*angArray[28]-rad[ip+3*ng]*(3*angArray[64]+3*angArray[73])+rad[ip+4*ng]*angArray[137];
         bas[9] = rad[ip+2*ng]*3*angArray[29]-rad[ip+3*ng]*(angArray[65]+3*angArray[74])+rad[ip+4*ng]*angArray[138];
         bas[10] = -rad[ip+ng]*12*angArray[7]+rad[ip+2*ng]*(9*angArray[30]+12*angArray[23])-rad[ip+3*ng]*(9*angArray[66]+angArray[77])+rad[ip+4*ng]*angArray[141];
         bas[11] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(7*angArray[31]+6*angArray[24])-rad[ip+3*ng]*(7*angArray[67]+angArray[78])+rad[ip+4*ng]*angArray[142];
         bas[12] = -rad[ip+ng]*2*angArray[9]+rad[ip+2*ng]*(5*angArray[32]+2*angArray[25])-rad[ip+3*ng]*(5*angArray[68]+angArray[79])+rad[ip+4*ng]*angArray[143];
         bas[13] = rad[ip+2*ng]*3*angArray[33]-rad[ip+3*ng]*(3*angArray[69]+angArray[80])+rad[ip+4*ng]*angArray[144];
         bas[14] = rad[ip+2*ng]*angArray[34]-rad[ip+3*ng]*(angArray[70]+angArray[81])+rad[ip+4*ng]*angArray[145];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*20*angArray[10]+rad[ip+2*ng]*(20*angArray[38]+11*angArray[35])-rad[ip+3*ng]*(11*angArray[87]+angArray[84])+rad[ip+4*ng]*angArray[168];
         bas[1] = -rad[ip+ng]*36*angArray[11]+rad[ip+2*ng]*(27*angArray[36]+12*angArray[41])-rad[ip+3*ng]*(3*angArray[85]+9*angArray[90])+rad[ip+4*ng]*angArray[171];
         bas[2] = -rad[ip+ng]*12*angArray[12]+rad[ip+2*ng]*(9*angArray[37]+12*angArray[42])-rad[ip+3*ng]*(angArray[86]+9*angArray[91])+rad[ip+4*ng]*angArray[172];
         bas[3] = rad[ip]*12*angArray[1]-rad[ip+ng]*(14*angArray[10]+30*angArray[13])+rad[ip+2*ng]*(35*angArray[38]+2*angArray[35]+6*angArray[45])-rad[ip+3*ng]*(7*angArray[94]+5*angArray[87])+rad[ip+4*ng]*angArray[175];
         bas[4] = -rad[ip+ng]*18*angArray[14]+rad[ip+2*ng]*(21*angArray[39]+6*angArray[46])-rad[ip+3*ng]*(7*angArray[95]+3*angArray[88])+rad[ip+4*ng]*angArray[176];
         bas[5] = -rad[ip+ng]*6*angArray[15]+rad[ip+2*ng]*(7*angArray[40]+6*angArray[47])-rad[ip+3*ng]*(7*angArray[96]+angArray[89])+rad[ip+4*ng]*angArray[177];
         bas[6] = rad[ip]*12*angArray[2]-rad[ip+ng]*(30*angArray[11]+14*angArray[16])+rad[ip+2*ng]*(35*angArray[41]+6*angArray[36]+2*angArray[50])-rad[ip+3*ng]*(5*angArray[99]+7*angArray[90])+rad[ip+4*ng]*angArray[180];
         bas[7] = rad[ip]*4*angArray[3]-rad[ip+ng]*(10*angArray[12]+10*angArray[17])+rad[ip+2*ng]*(25*angArray[42]+2*angArray[37]+2*angArray[51])-rad[ip+3*ng]*(5*angArray[100]+5*angArray[91])+rad[ip+4*ng]*angArray[181];
         bas[8] = -rad[ip+ng]*6*angArray[18]+rad[ip+2*ng]*(15*angArray[43]+2*angArray[52])-rad[ip+3*ng]*(5*angArray[101]+3*angArray[92])+rad[ip+4*ng]*angArray[182];
         bas[9] = -rad[ip+ng]*2*angArray[19]+rad[ip+2*ng]*(5*angArray[44]+2*angArray[53])-rad[ip+3*ng]*(5*angArray[102]+angArray[93])+rad[ip+4*ng]*angArray[183];
         bas[10] = -rad[ip+ng]*36*angArray[13]+rad[ip+2*ng]*(27*angArray[45]+12*angArray[38])-rad[ip+3*ng]*(9*angArray[94]+3*angArray[105])+rad[ip+4*ng]*angArray[186];
         bas[11] = -rad[ip+ng]*18*angArray[14]+rad[ip+2*ng]*(21*angArray[46]+6*angArray[39])-rad[ip+3*ng]*(7*angArray[95]+3*angArray[106])+rad[ip+4*ng]*angArray[187];
         bas[12] = -rad[ip+ng]*6*angArray[15]+rad[ip+2*ng]*(2*angArray[40]+15*angArray[47])-rad[ip+3*ng]*(5*angArray[96]+3*angArray[107])+rad[ip+4*ng]*angArray[188];
         bas[13] = rad[ip+2*ng]*9*angArray[48]-rad[ip+3*ng]*(3*angArray[97]+3*angArray[108])+rad[ip+4*ng]*angArray[189];
         bas[14] = rad[ip+2*ng]*3*angArray[49]-rad[ip+3*ng]*(angArray[98]+3*angArray[109])+rad[ip+4*ng]*angArray[190];
         bas[15] = -rad[ip+ng]*20*angArray[16]+rad[ip+2*ng]*(20*angArray[41]+11*angArray[50])-rad[ip+3*ng]*(11*angArray[99]+angArray[112])+rad[ip+4*ng]*angArray[193];
         bas[16] = -rad[ip+ng]*12*angArray[17]+rad[ip+2*ng]*(9*angArray[51]+12*angArray[42])-rad[ip+3*ng]*(9*angArray[100]+angArray[113])+rad[ip+4*ng]*angArray[194];
         bas[17] = -rad[ip+ng]*6*angArray[18]+rad[ip+2*ng]*(7*angArray[52]+6*angArray[43])-rad[ip+3*ng]*(7*angArray[101]+angArray[114])+rad[ip+4*ng]*angArray[195];
         bas[18] = -rad[ip+ng]*2*angArray[19]+rad[ip+2*ng]*(5*angArray[53]+2*angArray[44])-rad[ip+3*ng]*(5*angArray[102]+angArray[115])+rad[ip+4*ng]*angArray[196];
         bas[19] = rad[ip+2*ng]*3*angArray[54]-rad[ip+3*ng]*(3*angArray[103]+angArray[116])+rad[ip+4*ng]*angArray[197];
         bas[20] = rad[ip+2*ng]*angArray[55]-rad[ip+3*ng]*(angArray[104]+angArray[117])+rad[ip+4*ng]*angArray[198];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*30*angArray[20]+rad[ip+2*ng]*(30*angArray[59]+13*angArray[56])-rad[ip+3*ng]*(angArray[120]+13*angArray[123])+rad[ip+4*ng]*angArray[223];
         bas[1] = -rad[ip+ng]*60*angArray[21]+rad[ip+2*ng]*(20*angArray[62]+33*angArray[57])-rad[ip+3*ng]*(3*angArray[121]+11*angArray[126])+rad[ip+4*ng]*angArray[226];
         bas[2] = -rad[ip+ng]*20*angArray[22]+rad[ip+2*ng]*(20*angArray[63]+11*angArray[58])-rad[ip+3*ng]*(11*angArray[127]+angArray[122])+rad[ip+4*ng]*angArray[227];
         bas[3] = rad[ip]*24*angArray[4]-rad[ip+ng]*(60*angArray[23]+18*angArray[20])+rad[ip+2*ng]*(45*angArray[59]+12*angArray[66]+2*angArray[56])-rad[ip+3*ng]*(5*angArray[123]+9*angArray[130])+rad[ip+4*ng]*angArray[230];
         bas[4] = -rad[ip+ng]*36*angArray[24]+rad[ip+2*ng]*(27*angArray[60]+12*angArray[67])-rad[ip+3*ng]*(3*angArray[124]+9*angArray[131])+rad[ip+4*ng]*angArray[231];
         bas[5] = -rad[ip+ng]*12*angArray[25]+rad[ip+2*ng]*(9*angArray[61]+12*angArray[68])-rad[ip+3*ng]*(angArray[125]+9*angArray[132])+rad[ip+4*ng]*angArray[232];
         bas[6] = rad[ip]*36*angArray[5]-rad[ip+ng]*(42*angArray[21]+42*angArray[26])+rad[ip+2*ng]*(49*angArray[62]+6*angArray[71]+6*angArray[57])-rad[ip+3*ng]*(7*angArray[135]+7*angArray[126])+rad[ip+4*ng]*angArray[235];
         bas[7] = rad[ip]*12*angArray[6]-rad[ip+ng]*(14*angArray[22]+30*angArray[27])+rad[ip+2*ng]*(35*angArray[63]+2*angArray[58]+6*angArray[72])-rad[ip+3*ng]*(7*angArray[136]+5*angArray[127])+rad[ip+4*ng]*angArray[236];
         bas[8] = -rad[ip+ng]*18*angArray[28]+rad[ip+2*ng]*(21*angArray[64]+6*angArray[73])-rad[ip+3*ng]*(7*angArray[137]+3*angArray[128])+rad[ip+4*ng]*angArray[237];
         bas[9] = -rad[ip+ng]*6*angArray[29]+rad[ip+2*ng]*(7*angArray[65]+6*angArray[74])-rad[ip+3*ng]*(7*angArray[138]+angArray[129])+rad[ip+4*ng]*angArray[238];
         bas[10] = rad[ip]*24*angArray[7]-rad[ip+ng]*(60*angArray[23]+18*angArray[30])+rad[ip+2*ng]*(45*angArray[66]+2*angArray[77]+12*angArray[59])-rad[ip+3*ng]*(5*angArray[141]+9*angArray[130])+rad[ip+4*ng]*angArray[241];
         bas[11] = rad[ip]*12*angArray[8]-rad[ip+ng]*(30*angArray[24]+14*angArray[31])+rad[ip+2*ng]*(35*angArray[67]+6*angArray[60]+2*angArray[78])-rad[ip+3*ng]*(5*angArray[142]+7*angArray[131])+rad[ip+4*ng]*angArray[242];
         bas[12] = rad[ip]*4*angArray[9]-rad[ip+ng]*(10*angArray[25]+10*angArray[32])+rad[ip+2*ng]*(25*angArray[68]+2*angArray[61]+2*angArray[79])-rad[ip+3*ng]*(5*angArray[143]+5*angArray[132])+rad[ip+4*ng]*angArray[243];
         bas[13] = -rad[ip+ng]*6*angArray[33]+rad[ip+2*ng]*(15*angArray[69]+2*angArray[80])-rad[ip+3*ng]*(5*angArray[144]+3*angArray[133])+rad[ip+4*ng]*angArray[244];
         bas[14] = -rad[ip+ng]*2*angArray[34]+rad[ip+2*ng]*(5*angArray[70]+2*angArray[81])-rad[ip+3*ng]*(5*angArray[145]+angArray[134])+rad[ip+4*ng]*angArray[245];
         bas[15] = -rad[ip+ng]*60*angArray[26]+rad[ip+2*ng]*(20*angArray[62]+33*angArray[71])-rad[ip+3*ng]*(11*angArray[135]+3*angArray[148])+rad[ip+4*ng]*angArray[248];
         bas[16] = -rad[ip+ng]*36*angArray[27]+rad[ip+2*ng]*(27*angArray[72]+12*angArray[63])-rad[ip+3*ng]*(9*angArray[136]+3*angArray[149])+rad[ip+4*ng]*angArray[249];
         bas[17] = -rad[ip+ng]*18*angArray[28]+rad[ip+2*ng]*(21*angArray[73]+6*angArray[64])-rad[ip+3*ng]*(7*angArray[137]+3*angArray[150])+rad[ip+4*ng]*angArray[250];
         bas[18] = -rad[ip+ng]*6*angArray[29]+rad[ip+2*ng]*(2*angArray[65]+15*angArray[74])-rad[ip+3*ng]*(5*angArray[138]+3*angArray[151])+rad[ip+4*ng]*angArray[251];
         bas[19] = rad[ip+2*ng]*9*angArray[75]-rad[ip+3*ng]*(3*angArray[139]+3*angArray[152])+rad[ip+4*ng]*angArray[252];
         bas[20] = rad[ip+2*ng]*3*angArray[76]-rad[ip+3*ng]*(angArray[140]+3*angArray[153])+rad[ip+4*ng]*angArray[253];
         bas[21] = -rad[ip+ng]*30*angArray[30]+rad[ip+2*ng]*(13*angArray[77]+30*angArray[66])-rad[ip+3*ng]*(13*angArray[141]+angArray[156])+rad[ip+4*ng]*angArray[256];
         bas[22] = -rad[ip+ng]*20*angArray[31]+rad[ip+2*ng]*(20*angArray[67]+11*angArray[78])-rad[ip+3*ng]*(11*angArray[142]+angArray[157])+rad[ip+4*ng]*angArray[257];
         bas[23] = -rad[ip+ng]*12*angArray[32]+rad[ip+2*ng]*(9*angArray[79]+12*angArray[68])-rad[ip+3*ng]*(9*angArray[143]+angArray[158])+rad[ip+4*ng]*angArray[258];
         bas[24] = -rad[ip+ng]*6*angArray[33]+rad[ip+2*ng]*(7*angArray[80]+6*angArray[69])-rad[ip+3*ng]*(7*angArray[144]+angArray[159])+rad[ip+4*ng]*angArray[259];
         bas[25] = -rad[ip+ng]*2*angArray[34]+rad[ip+2*ng]*(5*angArray[81]+2*angArray[70])-rad[ip+3*ng]*(5*angArray[145]+angArray[160])+rad[ip+4*ng]*angArray[260];
         bas[26] = rad[ip+2*ng]*3*angArray[82]-rad[ip+3*ng]*(3*angArray[146]+angArray[161])+rad[ip+4*ng]*angArray[261];
         bas[27] = rad[ip+2*ng]*angArray[83]-rad[ip+3*ng]*(angArray[147]+angArray[162])+rad[ip+4*ng]*angArray[262];
      }

   }


   // now we do derivatives for the given basis set to XYYY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[5]+rad[ip+4*ng]*angArray[26];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*(3*angArray[11]+angArray[16])+rad[ip+4*ng]*angArray[41];
         bas[1] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*6*angArray[13]+rad[ip+4*ng]*angArray[45];
         bas[2] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[46];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*6*angArray[5]-rad[ip+3*ng]*(2*angArray[26]+3*angArray[21])+rad[ip+4*ng]*angArray[62];
         bas[1] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*(3*angArray[4]+6*angArray[7])-rad[ip+3*ng]*(6*angArray[23]+angArray[30])+rad[ip+4*ng]*angArray[66];
         bas[2] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*(3*angArray[24]+angArray[31])+rad[ip+4*ng]*angArray[67];
         bas[3] = rad[ip+2*ng]*12*angArray[5]-rad[ip+3*ng]*9*angArray[26]+rad[ip+4*ng]*angArray[71];
         bas[4] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*6*angArray[27]+rad[ip+4*ng]*angArray[72];
         bas[5] = -rad[ip+3*ng]*3*angArray[28]+rad[ip+4*ng]*angArray[73];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*9*angArray[11]-rad[ip+3*ng]*(3*angArray[36]+3*angArray[41])+rad[ip+4*ng]*angArray[90];
         bas[1] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*(3*angArray[10]+12*angArray[13])-rad[ip+3*ng]*(2*angArray[45]+6*angArray[38])+rad[ip+4*ng]*angArray[94];
         bas[2] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[46]+3*angArray[39])+rad[ip+4*ng]*angArray[95];
         bas[3] = -rad[ip+ng]*12*angArray[2]+rad[ip+2*ng]*(9*angArray[16]+12*angArray[11])-rad[ip+3*ng]*(9*angArray[41]+angArray[50])+rad[ip+4*ng]*angArray[99];
         bas[4] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*(3*angArray[12]+6*angArray[17])-rad[ip+3*ng]*(6*angArray[42]+angArray[51])+rad[ip+4*ng]*angArray[100];
         bas[5] = rad[ip+2*ng]*3*angArray[18]-rad[ip+3*ng]*(3*angArray[43]+angArray[52])+rad[ip+4*ng]*angArray[101];
         bas[6] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*27*angArray[13]-rad[ip+3*ng]*12*angArray[45]+rad[ip+4*ng]*angArray[105];
         bas[7] = rad[ip+2*ng]*12*angArray[14]-rad[ip+3*ng]*9*angArray[46]+rad[ip+4*ng]*angArray[106];
         bas[8] = rad[ip+2*ng]*3*angArray[15]-rad[ip+3*ng]*6*angArray[47]+rad[ip+4*ng]*angArray[107];
         bas[9] = -rad[ip+3*ng]*3*angArray[48]+rad[ip+4*ng]*angArray[108];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*12*angArray[21]-rad[ip+3*ng]*(3*angArray[57]+4*angArray[62])+rad[ip+4*ng]*angArray[126];
         bas[1] = -rad[ip+ng]*9*angArray[4]+rad[ip+2*ng]*(18*angArray[23]+3*angArray[20])-rad[ip+3*ng]*(6*angArray[59]+3*angArray[66])+rad[ip+4*ng]*angArray[130];
         bas[2] = rad[ip+2*ng]*9*angArray[24]-rad[ip+3*ng]*(3*angArray[60]+3*angArray[67])+rad[ip+4*ng]*angArray[131];
         bas[3] = -rad[ip+ng]*24*angArray[5]+rad[ip+2*ng]*(18*angArray[26]+12*angArray[21])-rad[ip+3*ng]*(2*angArray[71]+9*angArray[62])+rad[ip+4*ng]*angArray[135];
         bas[4] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(3*angArray[22]+12*angArray[27])-rad[ip+3*ng]*(2*angArray[72]+6*angArray[63])+rad[ip+4*ng]*angArray[136];
         bas[5] = rad[ip+2*ng]*6*angArray[28]-rad[ip+3*ng]*(2*angArray[73]+3*angArray[64])+rad[ip+4*ng]*angArray[137];
         bas[6] = rad[ip]*6*angArray[0]-rad[ip+ng]*(6*angArray[4]+27*angArray[7])+rad[ip+2*ng]*(12*angArray[30]+27*angArray[23])-rad[ip+3*ng]*(12*angArray[66]+angArray[77])+rad[ip+4*ng]*angArray[141];
         bas[7] = -rad[ip+ng]*12*angArray[8]+rad[ip+2*ng]*(9*angArray[31]+12*angArray[24])-rad[ip+3*ng]*(9*angArray[67]+angArray[78])+rad[ip+4*ng]*angArray[142];
         bas[8] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*(3*angArray[25]+6*angArray[32])-rad[ip+3*ng]*(6*angArray[68]+angArray[79])+rad[ip+4*ng]*angArray[143];
         bas[9] = rad[ip+2*ng]*3*angArray[33]-rad[ip+3*ng]*(3*angArray[69]+angArray[80])+rad[ip+4*ng]*angArray[144];
         bas[10] = -rad[ip+ng]*24*angArray[5]+rad[ip+2*ng]*48*angArray[26]-rad[ip+3*ng]*15*angArray[71]+rad[ip+4*ng]*angArray[148];
         bas[11] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*27*angArray[27]-rad[ip+3*ng]*12*angArray[72]+rad[ip+4*ng]*angArray[149];
         bas[12] = rad[ip+2*ng]*12*angArray[28]-rad[ip+3*ng]*9*angArray[73]+rad[ip+4*ng]*angArray[150];
         bas[13] = rad[ip+2*ng]*3*angArray[29]-rad[ip+3*ng]*6*angArray[74]+rad[ip+4*ng]*angArray[151];
         bas[14] = -rad[ip+3*ng]*3*angArray[75]+rad[ip+4*ng]*angArray[152];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*15*angArray[36]-rad[ip+3*ng]*(3*angArray[85]+5*angArray[90])+rad[ip+4*ng]*angArray[171];
         bas[1] = -rad[ip+ng]*12*angArray[10]+rad[ip+2*ng]*(3*angArray[35]+24*angArray[38])-rad[ip+3*ng]*(6*angArray[87]+4*angArray[94])+rad[ip+4*ng]*angArray[175];
         bas[2] = rad[ip+2*ng]*12*angArray[39]-rad[ip+3*ng]*(3*angArray[88]+4*angArray[95])+rad[ip+4*ng]*angArray[176];
         bas[3] = -rad[ip+ng]*36*angArray[11]+rad[ip+2*ng]*(27*angArray[41]+12*angArray[36])-rad[ip+3*ng]*(3*angArray[99]+9*angArray[90])+rad[ip+4*ng]*angArray[180];
         bas[4] = -rad[ip+ng]*9*angArray[12]+rad[ip+2*ng]*(18*angArray[42]+3*angArray[37])-rad[ip+3*ng]*(6*angArray[91]+3*angArray[100])+rad[ip+4*ng]*angArray[181];
         bas[5] = rad[ip+2*ng]*9*angArray[43]-rad[ip+3*ng]*(3*angArray[92]+3*angArray[101])+rad[ip+4*ng]*angArray[182];
         bas[6] = rad[ip]*12*angArray[1]-rad[ip+ng]*(54*angArray[13]+6*angArray[10])+rad[ip+2*ng]*(24*angArray[45]+27*angArray[38])-rad[ip+3*ng]*(12*angArray[94]+2*angArray[105])+rad[ip+4*ng]*angArray[186];
         bas[7] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*(18*angArray[46]+12*angArray[39])-rad[ip+3*ng]*(2*angArray[106]+9*angArray[95])+rad[ip+4*ng]*angArray[187];
         bas[8] = -rad[ip+ng]*6*angArray[15]+rad[ip+2*ng]*(3*angArray[40]+12*angArray[47])-rad[ip+3*ng]*(2*angArray[107]+6*angArray[96])+rad[ip+4*ng]*angArray[188];
         bas[9] = rad[ip+2*ng]*6*angArray[48]-rad[ip+3*ng]*(2*angArray[108]+3*angArray[97])+rad[ip+4*ng]*angArray[189];
         bas[10] = rad[ip]*24*angArray[2]-rad[ip+ng]*(24*angArray[11]+48*angArray[16])+rad[ip+2*ng]*(48*angArray[41]+15*angArray[50])-rad[ip+3*ng]*(15*angArray[99]+angArray[112])+rad[ip+4*ng]*angArray[193];
         bas[11] = rad[ip]*6*angArray[3]-rad[ip+ng]*(6*angArray[12]+27*angArray[17])+rad[ip+2*ng]*(12*angArray[51]+27*angArray[42])-rad[ip+3*ng]*(12*angArray[100]+angArray[113])+rad[ip+4*ng]*angArray[194];
         bas[12] = -rad[ip+ng]*12*angArray[18]+rad[ip+2*ng]*(9*angArray[52]+12*angArray[43])-rad[ip+3*ng]*(9*angArray[101]+angArray[114])+rad[ip+4*ng]*angArray[195];
         bas[13] = -rad[ip+ng]*3*angArray[19]+rad[ip+2*ng]*(3*angArray[44]+6*angArray[53])-rad[ip+3*ng]*(6*angArray[102]+angArray[115])+rad[ip+4*ng]*angArray[196];
         bas[14] = rad[ip+2*ng]*3*angArray[54]-rad[ip+3*ng]*(3*angArray[103]+angArray[116])+rad[ip+4*ng]*angArray[197];
         bas[15] = -rad[ip+ng]*60*angArray[13]+rad[ip+2*ng]*75*angArray[45]-rad[ip+3*ng]*18*angArray[105]+rad[ip+4*ng]*angArray[201];
         bas[16] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*48*angArray[46]-rad[ip+3*ng]*15*angArray[106]+rad[ip+4*ng]*angArray[202];
         bas[17] = -rad[ip+ng]*6*angArray[15]+rad[ip+2*ng]*27*angArray[47]-rad[ip+3*ng]*12*angArray[107]+rad[ip+4*ng]*angArray[203];
         bas[18] = rad[ip+2*ng]*12*angArray[48]-rad[ip+3*ng]*9*angArray[108]+rad[ip+4*ng]*angArray[204];
         bas[19] = rad[ip+2*ng]*3*angArray[49]-rad[ip+3*ng]*6*angArray[109]+rad[ip+4*ng]*angArray[205];
         bas[20] = -rad[ip+3*ng]*3*angArray[110]+rad[ip+4*ng]*angArray[206];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*18*angArray[57]-rad[ip+3*ng]*(3*angArray[121]+6*angArray[126])+rad[ip+4*ng]*angArray[226];
         bas[1] = -rad[ip+ng]*15*angArray[20]+rad[ip+2*ng]*(30*angArray[59]+3*angArray[56])-rad[ip+3*ng]*(6*angArray[123]+5*angArray[130])+rad[ip+4*ng]*angArray[230];
         bas[2] = rad[ip+2*ng]*15*angArray[60]-rad[ip+3*ng]*(3*angArray[124]+5*angArray[131])+rad[ip+4*ng]*angArray[231];
         bas[3] = -rad[ip+ng]*48*angArray[21]+rad[ip+2*ng]*(12*angArray[57]+36*angArray[62])-rad[ip+3*ng]*(9*angArray[126]+4*angArray[135])+rad[ip+4*ng]*angArray[235];
         bas[4] = -rad[ip+ng]*12*angArray[22]+rad[ip+2*ng]*(3*angArray[58]+24*angArray[63])-rad[ip+3*ng]*(6*angArray[127]+4*angArray[136])+rad[ip+4*ng]*angArray[236];
         bas[5] = rad[ip+2*ng]*12*angArray[64]-rad[ip+3*ng]*(3*angArray[128]+4*angArray[137])+rad[ip+4*ng]*angArray[237];
         bas[6] = rad[ip]*18*angArray[4]-rad[ip+ng]*(81*angArray[23]+6*angArray[20])+rad[ip+2*ng]*(36*angArray[66]+27*angArray[59])-rad[ip+3*ng]*(12*angArray[130]+3*angArray[141])+rad[ip+4*ng]*angArray[241];
         bas[7] = -rad[ip+ng]*36*angArray[24]+rad[ip+2*ng]*(27*angArray[67]+12*angArray[60])-rad[ip+3*ng]*(3*angArray[142]+9*angArray[131])+rad[ip+4*ng]*angArray[242];
         bas[8] = -rad[ip+ng]*9*angArray[25]+rad[ip+2*ng]*(18*angArray[68]+3*angArray[61])-rad[ip+3*ng]*(6*angArray[132]+3*angArray[143])+rad[ip+4*ng]*angArray[243];
         bas[9] = rad[ip+2*ng]*9*angArray[69]-rad[ip+3*ng]*(3*angArray[133]+3*angArray[144])+rad[ip+4*ng]*angArray[244];
         bas[10] = rad[ip]*48*angArray[5]-rad[ip+ng]*(24*angArray[21]+96*angArray[26])+rad[ip+2*ng]*(30*angArray[71]+48*angArray[62])-rad[ip+3*ng]*(15*angArray[135]+2*angArray[148])+rad[ip+4*ng]*angArray[248];
         bas[11] = rad[ip]*12*angArray[6]-rad[ip+ng]*(54*angArray[27]+6*angArray[22])+rad[ip+2*ng]*(24*angArray[72]+27*angArray[63])-rad[ip+3*ng]*(12*angArray[136]+2*angArray[149])+rad[ip+4*ng]*angArray[249];
         bas[12] = -rad[ip+ng]*24*angArray[28]+rad[ip+2*ng]*(18*angArray[73]+12*angArray[64])-rad[ip+3*ng]*(2*angArray[150]+9*angArray[137])+rad[ip+4*ng]*angArray[250];
         bas[13] = -rad[ip+ng]*6*angArray[29]+rad[ip+2*ng]*(3*angArray[65]+12*angArray[74])-rad[ip+3*ng]*(2*angArray[151]+6*angArray[138])+rad[ip+4*ng]*angArray[251];
         bas[14] = rad[ip+2*ng]*6*angArray[75]-rad[ip+3*ng]*(2*angArray[152]+3*angArray[139])+rad[ip+4*ng]*angArray[252];
         bas[15] = rad[ip]*60*angArray[7]-rad[ip+ng]*(60*angArray[23]+75*angArray[30])+rad[ip+2*ng]*(18*angArray[77]+75*angArray[66])-rad[ip+3*ng]*(18*angArray[141]+angArray[156])+rad[ip+4*ng]*angArray[256];
         bas[16] = rad[ip]*24*angArray[8]-rad[ip+ng]*(24*angArray[24]+48*angArray[31])+rad[ip+2*ng]*(48*angArray[67]+15*angArray[78])-rad[ip+3*ng]*(15*angArray[142]+angArray[157])+rad[ip+4*ng]*angArray[257];
         bas[17] = rad[ip]*6*angArray[9]-rad[ip+ng]*(6*angArray[25]+27*angArray[32])+rad[ip+2*ng]*(12*angArray[79]+27*angArray[68])-rad[ip+3*ng]*(12*angArray[143]+angArray[158])+rad[ip+4*ng]*angArray[258];
         bas[18] = -rad[ip+ng]*12*angArray[33]+rad[ip+2*ng]*(9*angArray[80]+12*angArray[69])-rad[ip+3*ng]*(9*angArray[144]+angArray[159])+rad[ip+4*ng]*angArray[259];
         bas[19] = -rad[ip+ng]*3*angArray[34]+rad[ip+2*ng]*(3*angArray[70]+6*angArray[81])-rad[ip+3*ng]*(6*angArray[145]+angArray[160])+rad[ip+4*ng]*angArray[260];
         bas[20] = rad[ip+2*ng]*3*angArray[82]-rad[ip+3*ng]*(3*angArray[146]+angArray[161])+rad[ip+4*ng]*angArray[261];
         bas[21] = -rad[ip+ng]*120*angArray[26]+rad[ip+2*ng]*108*angArray[71]-rad[ip+3*ng]*21*angArray[148]+rad[ip+4*ng]*angArray[265];
         bas[22] = -rad[ip+ng]*60*angArray[27]+rad[ip+2*ng]*75*angArray[72]-rad[ip+3*ng]*18*angArray[149]+rad[ip+4*ng]*angArray[266];
         bas[23] = -rad[ip+ng]*24*angArray[28]+rad[ip+2*ng]*48*angArray[73]-rad[ip+3*ng]*15*angArray[150]+rad[ip+4*ng]*angArray[267];
         bas[24] = -rad[ip+ng]*6*angArray[29]+rad[ip+2*ng]*27*angArray[74]-rad[ip+3*ng]*12*angArray[151]+rad[ip+4*ng]*angArray[268];
         bas[25] = rad[ip+2*ng]*12*angArray[75]-rad[ip+3*ng]*9*angArray[152]+rad[ip+4*ng]*angArray[269];
         bas[26] = rad[ip+2*ng]*3*angArray[76]-rad[ip+3*ng]*6*angArray[153]+rad[ip+4*ng]*angArray[270];
         bas[27] = -rad[ip+3*ng]*3*angArray[154]+rad[ip+4*ng]*angArray[271];
      }

   }


   // now we do derivatives for the given basis set to YYYY
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[0]-rad[ip+3*ng]*6*angArray[7]+rad[ip+4*ng]*angArray[30];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*6*angArray[13]+rad[ip+4*ng]*angArray[45];
         bas[1] = rad[ip+2*ng]*15*angArray[2]-rad[ip+3*ng]*10*angArray[16]+rad[ip+4*ng]*angArray[50];
         bas[2] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*6*angArray[17]+rad[ip+4*ng]*angArray[51];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[4]-rad[ip+3*ng]*6*angArray[23]+rad[ip+4*ng]*angArray[66];
         bas[1] = rad[ip+2*ng]*15*angArray[5]-rad[ip+3*ng]*10*angArray[26]+rad[ip+4*ng]*angArray[71];
         bas[2] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*6*angArray[27]+rad[ip+4*ng]*angArray[72];
         bas[3] = -rad[ip+ng]*12*angArray[0]+rad[ip+2*ng]*39*angArray[7]-rad[ip+3*ng]*14*angArray[30]+rad[ip+4*ng]*angArray[77];
         bas[4] = rad[ip+2*ng]*15*angArray[8]-rad[ip+3*ng]*10*angArray[31]+rad[ip+4*ng]*angArray[78];
         bas[5] = rad[ip+2*ng]*3*angArray[9]-rad[ip+3*ng]*6*angArray[32]+rad[ip+4*ng]*angArray[79];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[10]-rad[ip+3*ng]*6*angArray[38]+rad[ip+4*ng]*angArray[94];
         bas[1] = rad[ip+2*ng]*15*angArray[11]-rad[ip+3*ng]*10*angArray[41]+rad[ip+4*ng]*angArray[99];
         bas[2] = rad[ip+2*ng]*3*angArray[12]-rad[ip+3*ng]*6*angArray[42]+rad[ip+4*ng]*angArray[100];
         bas[3] = -rad[ip+ng]*12*angArray[1]+rad[ip+2*ng]*39*angArray[13]-rad[ip+3*ng]*14*angArray[45]+rad[ip+4*ng]*angArray[105];
         bas[4] = rad[ip+2*ng]*15*angArray[14]-rad[ip+3*ng]*10*angArray[46]+rad[ip+4*ng]*angArray[106];
         bas[5] = rad[ip+2*ng]*3*angArray[15]-rad[ip+3*ng]*6*angArray[47]+rad[ip+4*ng]*angArray[107];
         bas[6] = -rad[ip+ng]*60*angArray[2]+rad[ip+2*ng]*75*angArray[16]-rad[ip+3*ng]*18*angArray[50]+rad[ip+4*ng]*angArray[112];
         bas[7] = -rad[ip+ng]*12*angArray[3]+rad[ip+2*ng]*39*angArray[17]-rad[ip+3*ng]*14*angArray[51]+rad[ip+4*ng]*angArray[113];
         bas[8] = rad[ip+2*ng]*15*angArray[18]-rad[ip+3*ng]*10*angArray[52]+rad[ip+4*ng]*angArray[114];
         bas[9] = rad[ip+2*ng]*3*angArray[19]-rad[ip+3*ng]*6*angArray[53]+rad[ip+4*ng]*angArray[115];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[20]-rad[ip+3*ng]*6*angArray[59]+rad[ip+4*ng]*angArray[130];
         bas[1] = rad[ip+2*ng]*15*angArray[21]-rad[ip+3*ng]*10*angArray[62]+rad[ip+4*ng]*angArray[135];
         bas[2] = rad[ip+2*ng]*3*angArray[22]-rad[ip+3*ng]*6*angArray[63]+rad[ip+4*ng]*angArray[136];
         bas[3] = -rad[ip+ng]*12*angArray[4]+rad[ip+2*ng]*39*angArray[23]-rad[ip+3*ng]*14*angArray[66]+rad[ip+4*ng]*angArray[141];
         bas[4] = rad[ip+2*ng]*15*angArray[24]-rad[ip+3*ng]*10*angArray[67]+rad[ip+4*ng]*angArray[142];
         bas[5] = rad[ip+2*ng]*3*angArray[25]-rad[ip+3*ng]*6*angArray[68]+rad[ip+4*ng]*angArray[143];
         bas[6] = -rad[ip+ng]*60*angArray[5]+rad[ip+2*ng]*75*angArray[26]-rad[ip+3*ng]*18*angArray[71]+rad[ip+4*ng]*angArray[148];
         bas[7] = -rad[ip+ng]*12*angArray[6]+rad[ip+2*ng]*39*angArray[27]-rad[ip+3*ng]*14*angArray[72]+rad[ip+4*ng]*angArray[149];
         bas[8] = rad[ip+2*ng]*15*angArray[28]-rad[ip+3*ng]*10*angArray[73]+rad[ip+4*ng]*angArray[150];
         bas[9] = rad[ip+2*ng]*3*angArray[29]-rad[ip+3*ng]*6*angArray[74]+rad[ip+4*ng]*angArray[151];
         bas[10] = rad[ip]*24*angArray[0]-rad[ip+ng]*168*angArray[7]+rad[ip+2*ng]*123*angArray[30]-rad[ip+3*ng]*22*angArray[77]+rad[ip+4*ng]*angArray[156];
         bas[11] = -rad[ip+ng]*60*angArray[8]+rad[ip+2*ng]*75*angArray[31]-rad[ip+3*ng]*18*angArray[78]+rad[ip+4*ng]*angArray[157];
         bas[12] = -rad[ip+ng]*12*angArray[9]+rad[ip+2*ng]*39*angArray[32]-rad[ip+3*ng]*14*angArray[79]+rad[ip+4*ng]*angArray[158];
         bas[13] = rad[ip+2*ng]*15*angArray[33]-rad[ip+3*ng]*10*angArray[80]+rad[ip+4*ng]*angArray[159];
         bas[14] = rad[ip+2*ng]*3*angArray[34]-rad[ip+3*ng]*6*angArray[81]+rad[ip+4*ng]*angArray[160];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[35]-rad[ip+3*ng]*6*angArray[87]+rad[ip+4*ng]*angArray[175];
         bas[1] = rad[ip+2*ng]*15*angArray[36]-rad[ip+3*ng]*10*angArray[90]+rad[ip+4*ng]*angArray[180];
         bas[2] = rad[ip+2*ng]*3*angArray[37]-rad[ip+3*ng]*6*angArray[91]+rad[ip+4*ng]*angArray[181];
         bas[3] = -rad[ip+ng]*12*angArray[10]+rad[ip+2*ng]*39*angArray[38]-rad[ip+3*ng]*14*angArray[94]+rad[ip+4*ng]*angArray[186];
         bas[4] = rad[ip+2*ng]*15*angArray[39]-rad[ip+3*ng]*10*angArray[95]+rad[ip+4*ng]*angArray[187];
         bas[5] = rad[ip+2*ng]*3*angArray[40]-rad[ip+3*ng]*6*angArray[96]+rad[ip+4*ng]*angArray[188];
         bas[6] = -rad[ip+ng]*60*angArray[11]+rad[ip+2*ng]*75*angArray[41]-rad[ip+3*ng]*18*angArray[99]+rad[ip+4*ng]*angArray[193];
         bas[7] = -rad[ip+ng]*12*angArray[12]+rad[ip+2*ng]*39*angArray[42]-rad[ip+3*ng]*14*angArray[100]+rad[ip+4*ng]*angArray[194];
         bas[8] = rad[ip+2*ng]*15*angArray[43]-rad[ip+3*ng]*10*angArray[101]+rad[ip+4*ng]*angArray[195];
         bas[9] = rad[ip+2*ng]*3*angArray[44]-rad[ip+3*ng]*6*angArray[102]+rad[ip+4*ng]*angArray[196];
         bas[10] = rad[ip]*24*angArray[1]-rad[ip+ng]*168*angArray[13]+rad[ip+2*ng]*123*angArray[45]-rad[ip+3*ng]*22*angArray[105]+rad[ip+4*ng]*angArray[201];
         bas[11] = -rad[ip+ng]*60*angArray[14]+rad[ip+2*ng]*75*angArray[46]-rad[ip+3*ng]*18*angArray[106]+rad[ip+4*ng]*angArray[202];
         bas[12] = -rad[ip+ng]*12*angArray[15]+rad[ip+2*ng]*39*angArray[47]-rad[ip+3*ng]*14*angArray[107]+rad[ip+4*ng]*angArray[203];
         bas[13] = rad[ip+2*ng]*15*angArray[48]-rad[ip+3*ng]*10*angArray[108]+rad[ip+4*ng]*angArray[204];
         bas[14] = rad[ip+2*ng]*3*angArray[49]-rad[ip+3*ng]*6*angArray[109]+rad[ip+4*ng]*angArray[205];
         bas[15] = rad[ip]*120*angArray[2]-rad[ip+ng]*360*angArray[16]+rad[ip+2*ng]*183*angArray[50]-rad[ip+3*ng]*26*angArray[112]+rad[ip+4*ng]*angArray[210];
         bas[16] = rad[ip]*24*angArray[3]-rad[ip+ng]*168*angArray[17]+rad[ip+2*ng]*123*angArray[51]-rad[ip+3*ng]*22*angArray[113]+rad[ip+4*ng]*angArray[211];
         bas[17] = -rad[ip+ng]*60*angArray[18]+rad[ip+2*ng]*75*angArray[52]-rad[ip+3*ng]*18*angArray[114]+rad[ip+4*ng]*angArray[212];
         bas[18] = -rad[ip+ng]*12*angArray[19]+rad[ip+2*ng]*39*angArray[53]-rad[ip+3*ng]*14*angArray[115]+rad[ip+4*ng]*angArray[213];
         bas[19] = rad[ip+2*ng]*15*angArray[54]-rad[ip+3*ng]*10*angArray[116]+rad[ip+4*ng]*angArray[214];
         bas[20] = rad[ip+2*ng]*3*angArray[55]-rad[ip+3*ng]*6*angArray[117]+rad[ip+4*ng]*angArray[215];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[56]-rad[ip+3*ng]*6*angArray[123]+rad[ip+4*ng]*angArray[230];
         bas[1] = rad[ip+2*ng]*15*angArray[57]-rad[ip+3*ng]*10*angArray[126]+rad[ip+4*ng]*angArray[235];
         bas[2] = rad[ip+2*ng]*3*angArray[58]-rad[ip+3*ng]*6*angArray[127]+rad[ip+4*ng]*angArray[236];
         bas[3] = -rad[ip+ng]*12*angArray[20]+rad[ip+2*ng]*39*angArray[59]-rad[ip+3*ng]*14*angArray[130]+rad[ip+4*ng]*angArray[241];
         bas[4] = rad[ip+2*ng]*15*angArray[60]-rad[ip+3*ng]*10*angArray[131]+rad[ip+4*ng]*angArray[242];
         bas[5] = rad[ip+2*ng]*3*angArray[61]-rad[ip+3*ng]*6*angArray[132]+rad[ip+4*ng]*angArray[243];
         bas[6] = -rad[ip+ng]*60*angArray[21]+rad[ip+2*ng]*75*angArray[62]-rad[ip+3*ng]*18*angArray[135]+rad[ip+4*ng]*angArray[248];
         bas[7] = -rad[ip+ng]*12*angArray[22]+rad[ip+2*ng]*39*angArray[63]-rad[ip+3*ng]*14*angArray[136]+rad[ip+4*ng]*angArray[249];
         bas[8] = rad[ip+2*ng]*15*angArray[64]-rad[ip+3*ng]*10*angArray[137]+rad[ip+4*ng]*angArray[250];
         bas[9] = rad[ip+2*ng]*3*angArray[65]-rad[ip+3*ng]*6*angArray[138]+rad[ip+4*ng]*angArray[251];
         bas[10] = rad[ip]*24*angArray[4]-rad[ip+ng]*168*angArray[23]+rad[ip+2*ng]*123*angArray[66]-rad[ip+3*ng]*22*angArray[141]+rad[ip+4*ng]*angArray[256];
         bas[11] = -rad[ip+ng]*60*angArray[24]+rad[ip+2*ng]*75*angArray[67]-rad[ip+3*ng]*18*angArray[142]+rad[ip+4*ng]*angArray[257];
         bas[12] = -rad[ip+ng]*12*angArray[25]+rad[ip+2*ng]*39*angArray[68]-rad[ip+3*ng]*14*angArray[143]+rad[ip+4*ng]*angArray[258];
         bas[13] = rad[ip+2*ng]*15*angArray[69]-rad[ip+3*ng]*10*angArray[144]+rad[ip+4*ng]*angArray[259];
         bas[14] = rad[ip+2*ng]*3*angArray[70]-rad[ip+3*ng]*6*angArray[145]+rad[ip+4*ng]*angArray[260];
         bas[15] = rad[ip]*120*angArray[5]-rad[ip+ng]*360*angArray[26]+rad[ip+2*ng]*183*angArray[71]-rad[ip+3*ng]*26*angArray[148]+rad[ip+4*ng]*angArray[265];
         bas[16] = rad[ip]*24*angArray[6]-rad[ip+ng]*168*angArray[27]+rad[ip+2*ng]*123*angArray[72]-rad[ip+3*ng]*22*angArray[149]+rad[ip+4*ng]*angArray[266];
         bas[17] = -rad[ip+ng]*60*angArray[28]+rad[ip+2*ng]*75*angArray[73]-rad[ip+3*ng]*18*angArray[150]+rad[ip+4*ng]*angArray[267];
         bas[18] = -rad[ip+ng]*12*angArray[29]+rad[ip+2*ng]*39*angArray[74]-rad[ip+3*ng]*14*angArray[151]+rad[ip+4*ng]*angArray[268];
         bas[19] = rad[ip+2*ng]*15*angArray[75]-rad[ip+3*ng]*10*angArray[152]+rad[ip+4*ng]*angArray[269];
         bas[20] = rad[ip+2*ng]*3*angArray[76]-rad[ip+3*ng]*6*angArray[153]+rad[ip+4*ng]*angArray[270];
         bas[21] = rad[ip]*360*angArray[7]-rad[ip+ng]*660*angArray[30]+rad[ip+2*ng]*255*angArray[77]-rad[ip+3*ng]*30*angArray[156]+rad[ip+4*ng]*angArray[275];
         bas[22] = rad[ip]*120*angArray[8]-rad[ip+ng]*360*angArray[31]+rad[ip+2*ng]*183*angArray[78]-rad[ip+3*ng]*26*angArray[157]+rad[ip+4*ng]*angArray[276];
         bas[23] = rad[ip]*24*angArray[9]-rad[ip+ng]*168*angArray[32]+rad[ip+2*ng]*123*angArray[79]-rad[ip+3*ng]*22*angArray[158]+rad[ip+4*ng]*angArray[277];
         bas[24] = -rad[ip+ng]*60*angArray[33]+rad[ip+2*ng]*75*angArray[80]-rad[ip+3*ng]*18*angArray[159]+rad[ip+4*ng]*angArray[278];
         bas[25] = -rad[ip+ng]*12*angArray[34]+rad[ip+2*ng]*39*angArray[81]-rad[ip+3*ng]*14*angArray[160]+rad[ip+4*ng]*angArray[279];
         bas[26] = rad[ip+2*ng]*15*angArray[82]-rad[ip+3*ng]*10*angArray[161]+rad[ip+4*ng]*angArray[280];
         bas[27] = rad[ip+2*ng]*3*angArray[83]-rad[ip+3*ng]*6*angArray[162]+rad[ip+4*ng]*angArray[281];
      }

   }


   // now we do derivatives for the given basis set to XXXZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[6]+rad[ip+4*ng]*angArray[22];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*6*angArray[12]+rad[ip+4*ng]*angArray[37];
         bas[1] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[39];
         bas[2] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*(3*angArray[15]+angArray[10])+rad[ip+4*ng]*angArray[40];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*12*angArray[6]-rad[ip+3*ng]*9*angArray[22]+rad[ip+4*ng]*angArray[58];
         bas[1] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*6*angArray[24]+rad[ip+4*ng]*angArray[60];
         bas[2] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*(6*angArray[4]+3*angArray[9])-rad[ip+3*ng]*(angArray[20]+6*angArray[25])+rad[ip+4*ng]*angArray[61];
         bas[3] = -rad[ip+3*ng]*3*angArray[27]+rad[ip+4*ng]*angArray[63];
         bas[4] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*(3*angArray[28]+angArray[21])+rad[ip+4*ng]*angArray[64];
         bas[5] = rad[ip+2*ng]*6*angArray[6]-rad[ip+3*ng]*(2*angArray[22]+3*angArray[29])+rad[ip+4*ng]*angArray[65];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*27*angArray[12]-rad[ip+3*ng]*12*angArray[37]+rad[ip+4*ng]*angArray[86];
         bas[1] = rad[ip+2*ng]*12*angArray[14]-rad[ip+3*ng]*9*angArray[39]+rad[ip+4*ng]*angArray[88];
         bas[2] = -rad[ip+ng]*12*angArray[1]+rad[ip+2*ng]*(12*angArray[15]+9*angArray[10])-rad[ip+3*ng]*(9*angArray[40]+angArray[35])+rad[ip+4*ng]*angArray[89];
         bas[3] = rad[ip+2*ng]*3*angArray[17]-rad[ip+3*ng]*6*angArray[42]+rad[ip+4*ng]*angArray[91];
         bas[4] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*(6*angArray[11]+3*angArray[18])-rad[ip+3*ng]*(angArray[36]+6*angArray[43])+rad[ip+4*ng]*angArray[92];
         bas[5] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*(12*angArray[12]+3*angArray[19])-rad[ip+3*ng]*(2*angArray[37]+6*angArray[44])+rad[ip+4*ng]*angArray[93];
         bas[6] = -rad[ip+3*ng]*3*angArray[46]+rad[ip+4*ng]*angArray[95];
         bas[7] = rad[ip+2*ng]*3*angArray[13]-rad[ip+3*ng]*(3*angArray[47]+angArray[38])+rad[ip+4*ng]*angArray[96];
         bas[8] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[39]+3*angArray[48])+rad[ip+4*ng]*angArray[97];
         bas[9] = rad[ip+2*ng]*9*angArray[15]-rad[ip+3*ng]*(3*angArray[40]+3*angArray[49])+rad[ip+4*ng]*angArray[98];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*24*angArray[6]+rad[ip+2*ng]*48*angArray[22]-rad[ip+3*ng]*15*angArray[58]+rad[ip+4*ng]*angArray[122];
         bas[1] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*27*angArray[24]-rad[ip+3*ng]*12*angArray[60]+rad[ip+4*ng]*angArray[124];
         bas[2] = rad[ip]*6*angArray[0]-rad[ip+ng]*(27*angArray[4]+6*angArray[9])+rad[ip+2*ng]*(12*angArray[20]+27*angArray[25])-rad[ip+3*ng]*(12*angArray[61]+angArray[56])+rad[ip+4*ng]*angArray[125];
         bas[3] = rad[ip+2*ng]*12*angArray[27]-rad[ip+3*ng]*9*angArray[63]+rad[ip+4*ng]*angArray[127];
         bas[4] = -rad[ip+ng]*12*angArray[5]+rad[ip+2*ng]*(12*angArray[28]+9*angArray[21])-rad[ip+3*ng]*(9*angArray[64]+angArray[57])+rad[ip+4*ng]*angArray[128];
         bas[5] = -rad[ip+ng]*24*angArray[6]+rad[ip+2*ng]*(12*angArray[29]+18*angArray[22])-rad[ip+3*ng]*(9*angArray[65]+2*angArray[58])+rad[ip+4*ng]*angArray[129];
         bas[6] = rad[ip+2*ng]*3*angArray[31]-rad[ip+3*ng]*6*angArray[67]+rad[ip+4*ng]*angArray[131];
         bas[7] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*(6*angArray[23]+3*angArray[32])-rad[ip+3*ng]*(angArray[59]+6*angArray[68])+rad[ip+4*ng]*angArray[132];
         bas[8] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(12*angArray[24]+3*angArray[33])-rad[ip+3*ng]*(2*angArray[60]+6*angArray[69])+rad[ip+4*ng]*angArray[133];
         bas[9] = -rad[ip+ng]*9*angArray[9]+rad[ip+2*ng]*(18*angArray[25]+3*angArray[34])-rad[ip+3*ng]*(3*angArray[61]+6*angArray[70])+rad[ip+4*ng]*angArray[134];
         bas[10] = -rad[ip+3*ng]*3*angArray[72]+rad[ip+4*ng]*angArray[136];
         bas[11] = rad[ip+2*ng]*3*angArray[26]-rad[ip+3*ng]*(3*angArray[73]+angArray[62])+rad[ip+4*ng]*angArray[137];
         bas[12] = rad[ip+2*ng]*6*angArray[27]-rad[ip+3*ng]*(2*angArray[63]+3*angArray[74])+rad[ip+4*ng]*angArray[138];
         bas[13] = rad[ip+2*ng]*9*angArray[28]-rad[ip+3*ng]*(3*angArray[64]+3*angArray[75])+rad[ip+4*ng]*angArray[139];
         bas[14] = rad[ip+2*ng]*12*angArray[29]-rad[ip+3*ng]*(3*angArray[76]+4*angArray[65])+rad[ip+4*ng]*angArray[140];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*60*angArray[12]+rad[ip+2*ng]*75*angArray[37]-rad[ip+3*ng]*18*angArray[86]+rad[ip+4*ng]*angArray[167];
         bas[1] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*48*angArray[39]-rad[ip+3*ng]*15*angArray[88]+rad[ip+4*ng]*angArray[169];
         bas[2] = rad[ip]*24*angArray[1]-rad[ip+ng]*(48*angArray[10]+24*angArray[15])+rad[ip+2*ng]*(15*angArray[35]+48*angArray[40])-rad[ip+3*ng]*(angArray[84]+15*angArray[89])+rad[ip+4*ng]*angArray[170];
         bas[3] = -rad[ip+ng]*6*angArray[17]+rad[ip+2*ng]*27*angArray[42]-rad[ip+3*ng]*12*angArray[91]+rad[ip+4*ng]*angArray[172];
         bas[4] = rad[ip]*6*angArray[2]-rad[ip+ng]*(27*angArray[11]+6*angArray[18])+rad[ip+2*ng]*(12*angArray[36]+27*angArray[43])-rad[ip+3*ng]*(12*angArray[92]+angArray[85])+rad[ip+4*ng]*angArray[173];
         bas[5] = rad[ip]*12*angArray[3]-rad[ip+ng]*(54*angArray[12]+6*angArray[19])+rad[ip+2*ng]*(24*angArray[37]+27*angArray[44])-rad[ip+3*ng]*(12*angArray[93]+2*angArray[86])+rad[ip+4*ng]*angArray[174];
         bas[6] = rad[ip+2*ng]*12*angArray[46]-rad[ip+3*ng]*9*angArray[95]+rad[ip+4*ng]*angArray[176];
         bas[7] = -rad[ip+ng]*12*angArray[13]+rad[ip+2*ng]*(12*angArray[47]+9*angArray[38])-rad[ip+3*ng]*(9*angArray[96]+angArray[87])+rad[ip+4*ng]*angArray[177];
         bas[8] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*(12*angArray[48]+18*angArray[39])-rad[ip+3*ng]*(9*angArray[97]+2*angArray[88])+rad[ip+4*ng]*angArray[178];
         bas[9] = -rad[ip+ng]*36*angArray[15]+rad[ip+2*ng]*(12*angArray[49]+27*angArray[40])-rad[ip+3*ng]*(9*angArray[98]+3*angArray[89])+rad[ip+4*ng]*angArray[179];
         bas[10] = rad[ip+2*ng]*3*angArray[51]-rad[ip+3*ng]*6*angArray[100]+rad[ip+4*ng]*angArray[181];
         bas[11] = -rad[ip+ng]*3*angArray[16]+rad[ip+2*ng]*(6*angArray[41]+3*angArray[52])-rad[ip+3*ng]*(angArray[90]+6*angArray[101])+rad[ip+4*ng]*angArray[182];
         bas[12] = -rad[ip+ng]*6*angArray[17]+rad[ip+2*ng]*(12*angArray[42]+3*angArray[53])-rad[ip+3*ng]*(2*angArray[91]+6*angArray[102])+rad[ip+4*ng]*angArray[183];
         bas[13] = -rad[ip+ng]*9*angArray[18]+rad[ip+2*ng]*(18*angArray[43]+3*angArray[54])-rad[ip+3*ng]*(3*angArray[92]+6*angArray[103])+rad[ip+4*ng]*angArray[184];
         bas[14] = -rad[ip+ng]*12*angArray[19]+rad[ip+2*ng]*(3*angArray[55]+24*angArray[44])-rad[ip+3*ng]*(4*angArray[93]+6*angArray[104])+rad[ip+4*ng]*angArray[185];
         bas[15] = -rad[ip+3*ng]*3*angArray[106]+rad[ip+4*ng]*angArray[187];
         bas[16] = rad[ip+2*ng]*3*angArray[45]-rad[ip+3*ng]*(3*angArray[107]+angArray[94])+rad[ip+4*ng]*angArray[188];
         bas[17] = rad[ip+2*ng]*6*angArray[46]-rad[ip+3*ng]*(2*angArray[95]+3*angArray[108])+rad[ip+4*ng]*angArray[189];
         bas[18] = rad[ip+2*ng]*9*angArray[47]-rad[ip+3*ng]*(3*angArray[96]+3*angArray[109])+rad[ip+4*ng]*angArray[190];
         bas[19] = rad[ip+2*ng]*12*angArray[48]-rad[ip+3*ng]*(3*angArray[110]+4*angArray[97])+rad[ip+4*ng]*angArray[191];
         bas[20] = rad[ip+2*ng]*15*angArray[49]-rad[ip+3*ng]*(3*angArray[111]+5*angArray[98])+rad[ip+4*ng]*angArray[192];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*120*angArray[22]+rad[ip+2*ng]*108*angArray[58]-rad[ip+3*ng]*21*angArray[122]+rad[ip+4*ng]*angArray[222];
         bas[1] = -rad[ip+ng]*60*angArray[24]+rad[ip+2*ng]*75*angArray[60]-rad[ip+3*ng]*18*angArray[124]+rad[ip+4*ng]*angArray[224];
         bas[2] = rad[ip]*60*angArray[4]-rad[ip+ng]*(75*angArray[20]+60*angArray[25])+rad[ip+2*ng]*(75*angArray[61]+18*angArray[56])-rad[ip+3*ng]*(18*angArray[125]+angArray[120])+rad[ip+4*ng]*angArray[225];
         bas[3] = -rad[ip+ng]*24*angArray[27]+rad[ip+2*ng]*48*angArray[63]-rad[ip+3*ng]*15*angArray[127]+rad[ip+4*ng]*angArray[227];
         bas[4] = rad[ip]*24*angArray[5]-rad[ip+ng]*(48*angArray[21]+24*angArray[28])+rad[ip+2*ng]*(15*angArray[57]+48*angArray[64])-rad[ip+3*ng]*(angArray[121]+15*angArray[128])+rad[ip+4*ng]*angArray[228];
         bas[5] = rad[ip]*48*angArray[6]-rad[ip+ng]*(96*angArray[22]+24*angArray[29])+rad[ip+2*ng]*(48*angArray[65]+30*angArray[58])-rad[ip+3*ng]*(2*angArray[122]+15*angArray[129])+rad[ip+4*ng]*angArray[229];
         bas[6] = -rad[ip+ng]*6*angArray[31]+rad[ip+2*ng]*27*angArray[67]-rad[ip+3*ng]*12*angArray[131]+rad[ip+4*ng]*angArray[231];
         bas[7] = rad[ip]*6*angArray[7]-rad[ip+ng]*(27*angArray[23]+6*angArray[32])+rad[ip+2*ng]*(12*angArray[59]+27*angArray[68])-rad[ip+3*ng]*(12*angArray[132]+angArray[123])+rad[ip+4*ng]*angArray[232];
         bas[8] = rad[ip]*12*angArray[8]-rad[ip+ng]*(54*angArray[24]+6*angArray[33])+rad[ip+2*ng]*(24*angArray[60]+27*angArray[69])-rad[ip+3*ng]*(12*angArray[133]+2*angArray[124])+rad[ip+4*ng]*angArray[233];
         bas[9] = rad[ip]*18*angArray[9]-rad[ip+ng]*(81*angArray[25]+6*angArray[34])+rad[ip+2*ng]*(36*angArray[61]+27*angArray[70])-rad[ip+3*ng]*(12*angArray[134]+3*angArray[125])+rad[ip+4*ng]*angArray[234];
         bas[10] = rad[ip+2*ng]*12*angArray[72]-rad[ip+3*ng]*9*angArray[136]+rad[ip+4*ng]*angArray[236];
         bas[11] = -rad[ip+ng]*12*angArray[26]+rad[ip+2*ng]*(12*angArray[73]+9*angArray[62])-rad[ip+3*ng]*(9*angArray[137]+angArray[126])+rad[ip+4*ng]*angArray[237];
         bas[12] = -rad[ip+ng]*24*angArray[27]+rad[ip+2*ng]*(12*angArray[74]+18*angArray[63])-rad[ip+3*ng]*(9*angArray[138]+2*angArray[127])+rad[ip+4*ng]*angArray[238];
         bas[13] = -rad[ip+ng]*36*angArray[28]+rad[ip+2*ng]*(12*angArray[75]+27*angArray[64])-rad[ip+3*ng]*(9*angArray[139]+3*angArray[128])+rad[ip+4*ng]*angArray[239];
         bas[14] = -rad[ip+ng]*48*angArray[29]+rad[ip+2*ng]*(36*angArray[65]+12*angArray[76])-rad[ip+3*ng]*(9*angArray[140]+4*angArray[129])+rad[ip+4*ng]*angArray[240];
         bas[15] = rad[ip+2*ng]*3*angArray[78]-rad[ip+3*ng]*6*angArray[142]+rad[ip+4*ng]*angArray[242];
         bas[16] = -rad[ip+ng]*3*angArray[30]+rad[ip+2*ng]*(6*angArray[66]+3*angArray[79])-rad[ip+3*ng]*(angArray[130]+6*angArray[143])+rad[ip+4*ng]*angArray[243];
         bas[17] = -rad[ip+ng]*6*angArray[31]+rad[ip+2*ng]*(12*angArray[67]+3*angArray[80])-rad[ip+3*ng]*(2*angArray[131]+6*angArray[144])+rad[ip+4*ng]*angArray[244];
         bas[18] = -rad[ip+ng]*9*angArray[32]+rad[ip+2*ng]*(18*angArray[68]+3*angArray[81])-rad[ip+3*ng]*(3*angArray[132]+6*angArray[145])+rad[ip+4*ng]*angArray[245];
         bas[19] = -rad[ip+ng]*12*angArray[33]+rad[ip+2*ng]*(3*angArray[82]+24*angArray[69])-rad[ip+3*ng]*(4*angArray[133]+6*angArray[146])+rad[ip+4*ng]*angArray[246];
         bas[20] = -rad[ip+ng]*15*angArray[34]+rad[ip+2*ng]*(30*angArray[70]+3*angArray[83])-rad[ip+3*ng]*(5*angArray[134]+6*angArray[147])+rad[ip+4*ng]*angArray[247];
         bas[21] = -rad[ip+3*ng]*3*angArray[149]+rad[ip+4*ng]*angArray[249];
         bas[22] = rad[ip+2*ng]*3*angArray[71]-rad[ip+3*ng]*(3*angArray[150]+angArray[135])+rad[ip+4*ng]*angArray[250];
         bas[23] = rad[ip+2*ng]*6*angArray[72]-rad[ip+3*ng]*(2*angArray[136]+3*angArray[151])+rad[ip+4*ng]*angArray[251];
         bas[24] = rad[ip+2*ng]*9*angArray[73]-rad[ip+3*ng]*(3*angArray[137]+3*angArray[152])+rad[ip+4*ng]*angArray[252];
         bas[25] = rad[ip+2*ng]*12*angArray[74]-rad[ip+3*ng]*(3*angArray[153]+4*angArray[138])+rad[ip+4*ng]*angArray[253];
         bas[26] = rad[ip+2*ng]*15*angArray[75]-rad[ip+3*ng]*(3*angArray[154]+5*angArray[139])+rad[ip+4*ng]*angArray[254];
         bas[27] = rad[ip+2*ng]*18*angArray[76]-rad[ip+3*ng]*(3*angArray[155]+6*angArray[140])+rad[ip+4*ng]*angArray[255];
      }

   }


   // now we do derivatives for the given basis set to XXYZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*angArray[8]+rad[ip+4*ng]*angArray[24];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[39];
         bas[1] = rad[ip+2*ng]*angArray[3]-rad[ip+3*ng]*(angArray[12]+angArray[17])+rad[ip+4*ng]*angArray[42];
         bas[2] = rad[ip+2*ng]*angArray[2]-rad[ip+3*ng]*(angArray[11]+angArray[18])+rad[ip+4*ng]*angArray[43];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*2*angArray[8]-rad[ip+3*ng]*5*angArray[24]+rad[ip+4*ng]*angArray[60];
         bas[1] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*(angArray[22]+3*angArray[27])+rad[ip+4*ng]*angArray[63];
         bas[2] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*(angArray[21]+3*angArray[28])+rad[ip+4*ng]*angArray[64];
         bas[3] = rad[ip+2*ng]*2*angArray[8]-rad[ip+3*ng]*(2*angArray[24]+angArray[31])+rad[ip+4*ng]*angArray[67];
         bas[4] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[7]+angArray[9]+angArray[4])-rad[ip+3*ng]*(angArray[23]+angArray[32]+angArray[25])+rad[ip+4*ng]*angArray[68];
         bas[5] = rad[ip+2*ng]*2*angArray[8]-rad[ip+3*ng]*(2*angArray[24]+angArray[33])+rad[ip+4*ng]*angArray[69];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*7*angArray[39]+rad[ip+4*ng]*angArray[88];
         bas[1] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(5*angArray[12]+2*angArray[17])-rad[ip+3*ng]*(5*angArray[42]+angArray[37])+rad[ip+4*ng]*angArray[91];
         bas[2] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(5*angArray[11]+2*angArray[18])-rad[ip+3*ng]*(5*angArray[43]+angArray[36])+rad[ip+4*ng]*angArray[92];
         bas[3] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[39]+3*angArray[46])+rad[ip+4*ng]*angArray[95];
         bas[4] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*(3*angArray[13]+3*angArray[15]+angArray[10])-rad[ip+3*ng]*(angArray[38]+angArray[40]+3*angArray[47])+rad[ip+4*ng]*angArray[96];
         bas[5] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[39]+3*angArray[48])+rad[ip+4*ng]*angArray[97];
         bas[6] = rad[ip+2*ng]*3*angArray[17]-rad[ip+3*ng]*(angArray[51]+3*angArray[42])+rad[ip+4*ng]*angArray[100];
         bas[7] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(angArray[16]+2*angArray[11]+2*angArray[18])-rad[ip+3*ng]*(angArray[41]+angArray[52]+2*angArray[43])+rad[ip+4*ng]*angArray[101];
         bas[8] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(2*angArray[17]+angArray[19]+2*angArray[12])-rad[ip+3*ng]*(2*angArray[42]+angArray[53]+angArray[44])+rad[ip+4*ng]*angArray[102];
         bas[9] = rad[ip+2*ng]*3*angArray[18]-rad[ip+3*ng]*(3*angArray[43]+angArray[54])+rad[ip+4*ng]*angArray[103];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*12*angArray[24]-rad[ip+3*ng]*9*angArray[60]+rad[ip+4*ng]*angArray[124];
         bas[1] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(7*angArray[22]+6*angArray[27])-rad[ip+3*ng]*(7*angArray[63]+angArray[58])+rad[ip+4*ng]*angArray[127];
         bas[2] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(7*angArray[21]+6*angArray[28])-rad[ip+3*ng]*(7*angArray[64]+angArray[57])+rad[ip+4*ng]*angArray[128];
         bas[3] = -rad[ip+ng]*4*angArray[8]+rad[ip+2*ng]*(10*angArray[24]+2*angArray[31])-rad[ip+3*ng]*(5*angArray[67]+2*angArray[60])+rad[ip+4*ng]*angArray[131];
         bas[4] = rad[ip]*2*angArray[0]-rad[ip+ng]*(5*angArray[4]+2*angArray[9]+2*angArray[7])+rad[ip+2*ng]*(5*angArray[25]+5*angArray[23]+2*angArray[32]+angArray[20])-rad[ip+3*ng]*(5*angArray[68]+angArray[59]+angArray[61])+rad[ip+4*ng]*angArray[132];
         bas[5] = -rad[ip+ng]*4*angArray[8]+rad[ip+2*ng]*(10*angArray[24]+2*angArray[33])-rad[ip+3*ng]*(5*angArray[69]+2*angArray[60])+rad[ip+4*ng]*angArray[133];
         bas[6] = rad[ip+2*ng]*9*angArray[27]-rad[ip+3*ng]*(3*angArray[72]+3*angArray[63])+rad[ip+4*ng]*angArray[136];
         bas[7] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(3*angArray[26]+6*angArray[28]+2*angArray[21])-rad[ip+3*ng]*(angArray[62]+3*angArray[73]+2*angArray[64])+rad[ip+4*ng]*angArray[137];
         bas[8] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(6*angArray[27]+3*angArray[29]+2*angArray[22])-rad[ip+3*ng]*(2*angArray[63]+angArray[65]+3*angArray[74])+rad[ip+4*ng]*angArray[138];
         bas[9] = rad[ip+2*ng]*9*angArray[28]-rad[ip+3*ng]*(3*angArray[64]+3*angArray[75])+rad[ip+4*ng]*angArray[139];
         bas[10] = rad[ip+2*ng]*4*angArray[31]-rad[ip+3*ng]*(4*angArray[67]+angArray[78])+rad[ip+4*ng]*angArray[142];
         bas[11] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*(angArray[30]+3*angArray[23]+3*angArray[32])-rad[ip+3*ng]*(angArray[66]+angArray[79]+3*angArray[68])+rad[ip+4*ng]*angArray[143];
         bas[12] = -rad[ip+ng]*4*angArray[8]+rad[ip+2*ng]*(2*angArray[31]+4*angArray[24]+2*angArray[33])-rad[ip+3*ng]*(2*angArray[67]+2*angArray[69]+angArray[80])+rad[ip+4*ng]*angArray[144];
         bas[13] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*(3*angArray[32]+angArray[34]+3*angArray[25])-rad[ip+3*ng]*(3*angArray[68]+angArray[81]+angArray[70])+rad[ip+4*ng]*angArray[145];
         bas[14] = rad[ip+2*ng]*4*angArray[33]-rad[ip+3*ng]*(4*angArray[69]+angArray[82])+rad[ip+4*ng]*angArray[146];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*20*angArray[39]-rad[ip+3*ng]*11*angArray[88]+rad[ip+4*ng]*angArray[169];
         bas[1] = -rad[ip+ng]*12*angArray[12]+rad[ip+2*ng]*(12*angArray[42]+9*angArray[37])-rad[ip+3*ng]*(angArray[86]+9*angArray[91])+rad[ip+4*ng]*angArray[172];
         bas[2] = -rad[ip+ng]*12*angArray[11]+rad[ip+2*ng]*(9*angArray[36]+12*angArray[43])-rad[ip+3*ng]*(angArray[85]+9*angArray[92])+rad[ip+4*ng]*angArray[173];
         bas[3] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(14*angArray[39]+6*angArray[46])-rad[ip+3*ng]*(7*angArray[95]+2*angArray[88])+rad[ip+4*ng]*angArray[176];
         bas[4] = rad[ip]*6*angArray[1]-rad[ip+ng]*(7*angArray[10]+6*angArray[15]+6*angArray[13])+rad[ip+2*ng]*(7*angArray[38]+7*angArray[40]+6*angArray[47]+angArray[35])-rad[ip+3*ng]*(7*angArray[96]+angArray[87]+angArray[89])+rad[ip+4*ng]*angArray[177];
         bas[5] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(14*angArray[39]+6*angArray[48])-rad[ip+3*ng]*(7*angArray[97]+2*angArray[88])+rad[ip+4*ng]*angArray[178];
         bas[6] = -rad[ip+ng]*6*angArray[17]+rad[ip+2*ng]*(15*angArray[42]+2*angArray[51])-rad[ip+3*ng]*(5*angArray[100]+3*angArray[91])+rad[ip+4*ng]*angArray[181];
         bas[7] = rad[ip]*4*angArray[2]-rad[ip+ng]*(10*angArray[11]+2*angArray[16]+4*angArray[18])+rad[ip+2*ng]*(10*angArray[43]+5*angArray[41]+2*angArray[52]+2*angArray[36])-rad[ip+3*ng]*(5*angArray[101]+2*angArray[92]+angArray[90])+rad[ip+4*ng]*angArray[182];
         bas[8] = rad[ip]*4*angArray[3]-rad[ip+ng]*(10*angArray[12]+2*angArray[19]+4*angArray[17])+rad[ip+2*ng]*(5*angArray[44]+10*angArray[42]+2*angArray[53]+2*angArray[37])-rad[ip+3*ng]*(5*angArray[102]+2*angArray[91]+angArray[93])+rad[ip+4*ng]*angArray[183];
         bas[9] = -rad[ip+ng]*6*angArray[18]+rad[ip+2*ng]*(15*angArray[43]+2*angArray[54])-rad[ip+3*ng]*(5*angArray[103]+3*angArray[92])+rad[ip+4*ng]*angArray[184];
         bas[10] = rad[ip+2*ng]*12*angArray[46]-rad[ip+3*ng]*(4*angArray[95]+3*angArray[106])+rad[ip+4*ng]*angArray[187];
         bas[11] = -rad[ip+ng]*9*angArray[13]+rad[ip+2*ng]*(3*angArray[45]+9*angArray[47]+3*angArray[38])-rad[ip+3*ng]*(angArray[94]+3*angArray[107]+3*angArray[96])+rad[ip+4*ng]*angArray[188];
         bas[12] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(6*angArray[48]+6*angArray[46]+4*angArray[39])-rad[ip+3*ng]*(2*angArray[95]+2*angArray[97]+3*angArray[108])+rad[ip+4*ng]*angArray[189];
         bas[13] = -rad[ip+ng]*9*angArray[15]+rad[ip+2*ng]*(9*angArray[47]+3*angArray[49]+3*angArray[40])-rad[ip+3*ng]*(3*angArray[109]+3*angArray[96]+angArray[98])+rad[ip+4*ng]*angArray[190];
         bas[14] = rad[ip+2*ng]*12*angArray[48]-rad[ip+3*ng]*(4*angArray[97]+3*angArray[110])+rad[ip+4*ng]*angArray[191];
         bas[15] = rad[ip+2*ng]*5*angArray[51]-rad[ip+3*ng]*(5*angArray[100]+angArray[113])+rad[ip+4*ng]*angArray[194];
         bas[16] = -rad[ip+ng]*4*angArray[16]+rad[ip+2*ng]*(angArray[50]+4*angArray[52]+4*angArray[41])-rad[ip+3*ng]*(angArray[99]+4*angArray[101]+angArray[114])+rad[ip+4*ng]*angArray[195];
         bas[17] = -rad[ip+ng]*6*angArray[17]+rad[ip+2*ng]*(2*angArray[51]+6*angArray[42]+3*angArray[53])-rad[ip+3*ng]*(2*angArray[100]+angArray[115]+3*angArray[102])+rad[ip+4*ng]*angArray[196];
         bas[18] = -rad[ip+ng]*6*angArray[18]+rad[ip+2*ng]*(3*angArray[52]+6*angArray[43]+2*angArray[54])-rad[ip+3*ng]*(3*angArray[101]+2*angArray[103]+angArray[116])+rad[ip+4*ng]*angArray[197];
         bas[19] = -rad[ip+ng]*4*angArray[19]+rad[ip+2*ng]*(4*angArray[53]+4*angArray[44]+angArray[55])-rad[ip+3*ng]*(4*angArray[102]+angArray[117]+angArray[104])+rad[ip+4*ng]*angArray[198];
         bas[20] = rad[ip+2*ng]*5*angArray[54]-rad[ip+3*ng]*(angArray[118]+5*angArray[103])+rad[ip+4*ng]*angArray[199];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*30*angArray[60]-rad[ip+3*ng]*13*angArray[124]+rad[ip+4*ng]*angArray[224];
         bas[1] = -rad[ip+ng]*20*angArray[22]+rad[ip+2*ng]*(20*angArray[63]+11*angArray[58])-rad[ip+3*ng]*(11*angArray[127]+angArray[122])+rad[ip+4*ng]*angArray[227];
         bas[2] = -rad[ip+ng]*20*angArray[21]+rad[ip+2*ng]*(20*angArray[64]+11*angArray[57])-rad[ip+3*ng]*(11*angArray[128]+angArray[121])+rad[ip+4*ng]*angArray[228];
         bas[3] = -rad[ip+ng]*24*angArray[24]+rad[ip+2*ng]*(18*angArray[60]+12*angArray[67])-rad[ip+3*ng]*(2*angArray[124]+9*angArray[131])+rad[ip+4*ng]*angArray[231];
         bas[4] = rad[ip]*12*angArray[4]-rad[ip+ng]*(12*angArray[25]+12*angArray[23]+9*angArray[20])+rad[ip+2*ng]*(angArray[56]+9*angArray[59]+12*angArray[68]+9*angArray[61])-rad[ip+3*ng]*(angArray[125]+angArray[123]+9*angArray[132])+rad[ip+4*ng]*angArray[232];
         bas[5] = -rad[ip+ng]*24*angArray[24]+rad[ip+2*ng]*(18*angArray[60]+12*angArray[69])-rad[ip+3*ng]*(2*angArray[124]+9*angArray[133])+rad[ip+4*ng]*angArray[233];
         bas[6] = -rad[ip+ng]*18*angArray[27]+rad[ip+2*ng]*(21*angArray[63]+6*angArray[72])-rad[ip+3*ng]*(7*angArray[136]+3*angArray[127])+rad[ip+4*ng]*angArray[236];
         bas[7] = rad[ip]*12*angArray[5]-rad[ip+ng]*(14*angArray[21]+6*angArray[26]+12*angArray[28])+rad[ip+2*ng]*(14*angArray[64]+7*angArray[62]+2*angArray[57]+6*angArray[73])-rad[ip+3*ng]*(7*angArray[137]+angArray[126]+2*angArray[128])+rad[ip+4*ng]*angArray[237];
         bas[8] = rad[ip]*12*angArray[6]-rad[ip+ng]*(14*angArray[22]+12*angArray[27]+6*angArray[29])+rad[ip+2*ng]*(14*angArray[63]+7*angArray[65]+6*angArray[74]+2*angArray[58])-rad[ip+3*ng]*(7*angArray[138]+2*angArray[127]+angArray[129])+rad[ip+4*ng]*angArray[238];
         bas[9] = -rad[ip+ng]*18*angArray[28]+rad[ip+2*ng]*(21*angArray[64]+6*angArray[75])-rad[ip+3*ng]*(7*angArray[139]+3*angArray[128])+rad[ip+4*ng]*angArray[239];
         bas[10] = -rad[ip+ng]*8*angArray[31]+rad[ip+2*ng]*(20*angArray[67]+2*angArray[78])-rad[ip+3*ng]*(5*angArray[142]+4*angArray[131])+rad[ip+4*ng]*angArray[242];
         bas[11] = rad[ip]*6*angArray[7]-rad[ip+ng]*(15*angArray[23]+6*angArray[32]+2*angArray[30])+rad[ip+2*ng]*(15*angArray[68]+5*angArray[66]+2*angArray[79]+3*angArray[59])-rad[ip+3*ng]*(5*angArray[143]+angArray[130]+3*angArray[132])+rad[ip+4*ng]*angArray[243];
         bas[12] = rad[ip]*8*angArray[8]-rad[ip+ng]*(20*angArray[24]+4*angArray[31]+4*angArray[33])+rad[ip+2*ng]*(10*angArray[69]+10*angArray[67]+4*angArray[60]+2*angArray[80])-rad[ip+3*ng]*(5*angArray[144]+2*angArray[131]+2*angArray[133])+rad[ip+4*ng]*angArray[244];
         bas[13] = rad[ip]*6*angArray[9]-rad[ip+ng]*(15*angArray[25]+2*angArray[34]+6*angArray[32])+rad[ip+2*ng]*(5*angArray[70]+15*angArray[68]+2*angArray[81]+3*angArray[61])-rad[ip+3*ng]*(5*angArray[145]+3*angArray[132]+angArray[134])+rad[ip+4*ng]*angArray[245];
         bas[14] = -rad[ip+ng]*8*angArray[33]+rad[ip+2*ng]*(20*angArray[69]+2*angArray[82])-rad[ip+3*ng]*(5*angArray[146]+4*angArray[133])+rad[ip+4*ng]*angArray[246];
         bas[15] = rad[ip+2*ng]*15*angArray[72]-rad[ip+3*ng]*(5*angArray[136]+3*angArray[149])+rad[ip+4*ng]*angArray[249];
         bas[16] = -rad[ip+ng]*12*angArray[26]+rad[ip+2*ng]*(3*angArray[71]+12*angArray[73]+4*angArray[62])-rad[ip+3*ng]*(angArray[135]+4*angArray[137]+3*angArray[150])+rad[ip+4*ng]*angArray[250];
         bas[17] = -rad[ip+ng]*18*angArray[27]+rad[ip+2*ng]*(6*angArray[72]+9*angArray[74]+6*angArray[63])-rad[ip+3*ng]*(2*angArray[136]+3*angArray[151]+3*angArray[138])+rad[ip+4*ng]*angArray[251];
         bas[18] = -rad[ip+ng]*18*angArray[28]+rad[ip+2*ng]*(6*angArray[64]+6*angArray[75]+9*angArray[73])-rad[ip+3*ng]*(3*angArray[137]+2*angArray[139]+3*angArray[152])+rad[ip+4*ng]*angArray[252];
         bas[19] = -rad[ip+ng]*12*angArray[29]+rad[ip+2*ng]*(3*angArray[76]+12*angArray[74]+4*angArray[65])-rad[ip+3*ng]*(4*angArray[138]+angArray[140]+3*angArray[153])+rad[ip+4*ng]*angArray[253];
         bas[20] = rad[ip+2*ng]*15*angArray[75]-rad[ip+3*ng]*(5*angArray[139]+3*angArray[154])+rad[ip+4*ng]*angArray[254];
         bas[21] = rad[ip+2*ng]*6*angArray[78]-rad[ip+3*ng]*(6*angArray[142]+angArray[157])+rad[ip+4*ng]*angArray[257];
         bas[22] = -rad[ip+ng]*5*angArray[30]+rad[ip+2*ng]*(5*angArray[66]+angArray[77]+5*angArray[79])-rad[ip+3*ng]*(angArray[141]+5*angArray[143]+angArray[158])+rad[ip+4*ng]*angArray[258];
         bas[23] = -rad[ip+ng]*8*angArray[31]+rad[ip+2*ng]*(2*angArray[78]+4*angArray[80]+8*angArray[67])-rad[ip+3*ng]*(2*angArray[142]+4*angArray[144]+angArray[159])+rad[ip+4*ng]*angArray[259];
         bas[24] = -rad[ip+ng]*9*angArray[32]+rad[ip+2*ng]*(3*angArray[79]+9*angArray[68]+3*angArray[81])-rad[ip+3*ng]*(3*angArray[143]+angArray[160]+3*angArray[145])+rad[ip+4*ng]*angArray[260];
         bas[25] = -rad[ip+ng]*8*angArray[33]+rad[ip+2*ng]*(4*angArray[80]+8*angArray[69]+2*angArray[82])-rad[ip+3*ng]*(2*angArray[146]+4*angArray[144]+angArray[161])+rad[ip+4*ng]*angArray[261];
         bas[26] = -rad[ip+ng]*5*angArray[34]+rad[ip+2*ng]*(5*angArray[81]+5*angArray[70]+angArray[83])-rad[ip+3*ng]*(5*angArray[145]+angArray[147]+angArray[162])+rad[ip+4*ng]*angArray[262];
         bas[27] = rad[ip+2*ng]*6*angArray[82]-rad[ip+3*ng]*(6*angArray[146]+angArray[163])+rad[ip+4*ng]*angArray[263];
      }

   }


   // now we do derivatives for the given basis set to XYYZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*angArray[6]+rad[ip+4*ng]*angArray[27];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[3]-rad[ip+3*ng]*(angArray[12]+angArray[17])+rad[ip+4*ng]*angArray[42];
         bas[1] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[46];
         bas[2] = rad[ip+2*ng]*angArray[1]-rad[ip+3*ng]*(angArray[13]+angArray[15])+rad[ip+4*ng]*angArray[47];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*2*angArray[6]-rad[ip+3*ng]*(angArray[22]+2*angArray[27])+rad[ip+4*ng]*angArray[63];
         bas[1] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*(3*angArray[24]+angArray[31])+rad[ip+4*ng]*angArray[67];
         bas[2] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[7]+angArray[9]+angArray[4])-rad[ip+3*ng]*(angArray[23]+angArray[25]+angArray[32])+rad[ip+4*ng]*angArray[68];
         bas[3] = rad[ip+2*ng]*2*angArray[6]-rad[ip+3*ng]*5*angArray[27]+rad[ip+4*ng]*angArray[72];
         bas[4] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*(3*angArray[28]+angArray[26])+rad[ip+4*ng]*angArray[73];
         bas[5] = rad[ip+2*ng]*2*angArray[6]-rad[ip+3*ng]*(2*angArray[27]+angArray[29])+rad[ip+4*ng]*angArray[74];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[12]-rad[ip+3*ng]*(angArray[37]+3*angArray[42])+rad[ip+4*ng]*angArray[91];
         bas[1] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[46]+3*angArray[39])+rad[ip+4*ng]*angArray[95];
         bas[2] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(2*angArray[13]+2*angArray[15]+angArray[10])-rad[ip+3*ng]*(2*angArray[47]+angArray[38]+angArray[40])+rad[ip+4*ng]*angArray[96];
         bas[3] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(5*angArray[17]+2*angArray[12])-rad[ip+3*ng]*(5*angArray[42]+angArray[51])+rad[ip+4*ng]*angArray[100];
         bas[4] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*(3*angArray[11]+angArray[16]+3*angArray[18])-rad[ip+3*ng]*(angArray[41]+3*angArray[43]+angArray[52])+rad[ip+4*ng]*angArray[101];
         bas[5] = -rad[ip+ng]*2*angArray[3]+rad[ip+2*ng]*(angArray[19]+2*angArray[17]+2*angArray[12])-rad[ip+3*ng]*(2*angArray[42]+angArray[44]+angArray[53])+rad[ip+4*ng]*angArray[102];
         bas[6] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*7*angArray[46]+rad[ip+4*ng]*angArray[106];
         bas[7] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(5*angArray[13]+2*angArray[15])-rad[ip+3*ng]*(angArray[45]+5*angArray[47])+rad[ip+4*ng]*angArray[107];
         bas[8] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[46]+3*angArray[48])+rad[ip+4*ng]*angArray[108];
         bas[9] = rad[ip+2*ng]*3*angArray[15]-rad[ip+3*ng]*(3*angArray[47]+angArray[49])+rad[ip+4*ng]*angArray[109];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*4*angArray[22]-rad[ip+3*ng]*(angArray[58]+4*angArray[63])+rad[ip+4*ng]*angArray[127];
         bas[1] = rad[ip+2*ng]*9*angArray[24]-rad[ip+3*ng]*(3*angArray[60]+3*angArray[67])+rad[ip+4*ng]*angArray[131];
         bas[2] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*(3*angArray[23]+angArray[20]+3*angArray[25])-rad[ip+3*ng]*(angArray[59]+3*angArray[68]+angArray[61])+rad[ip+4*ng]*angArray[132];
         bas[3] = -rad[ip+ng]*4*angArray[6]+rad[ip+2*ng]*(10*angArray[27]+2*angArray[22])-rad[ip+3*ng]*(2*angArray[72]+5*angArray[63])+rad[ip+4*ng]*angArray[136];
         bas[4] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(3*angArray[21]+6*angArray[28]+2*angArray[26])-rad[ip+3*ng]*(2*angArray[73]+3*angArray[64]+angArray[62])+rad[ip+4*ng]*angArray[137];
         bas[5] = -rad[ip+ng]*4*angArray[6]+rad[ip+2*ng]*(4*angArray[27]+2*angArray[29]+2*angArray[22])-rad[ip+3*ng]*(2*angArray[74]+2*angArray[63]+angArray[65])+rad[ip+4*ng]*angArray[138];
         bas[6] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(7*angArray[31]+6*angArray[24])-rad[ip+3*ng]*(7*angArray[67]+angArray[78])+rad[ip+4*ng]*angArray[142];
         bas[7] = rad[ip]*2*angArray[0]-rad[ip+ng]*(2*angArray[9]+5*angArray[7]+2*angArray[4])+rad[ip+2*ng]*(5*angArray[32]+angArray[30]+5*angArray[23]+2*angArray[25])-rad[ip+3*ng]*(angArray[66]+5*angArray[68]+angArray[79])+rad[ip+4*ng]*angArray[143];
         bas[8] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(6*angArray[24]+3*angArray[33]+2*angArray[31])-rad[ip+3*ng]*(2*angArray[67]+3*angArray[69]+angArray[80])+rad[ip+4*ng]*angArray[144];
         bas[9] = -rad[ip+ng]*3*angArray[9]+rad[ip+2*ng]*(angArray[34]+3*angArray[32]+3*angArray[25])-rad[ip+3*ng]*(angArray[81]+3*angArray[68]+angArray[70])+rad[ip+4*ng]*angArray[145];
         bas[10] = rad[ip+2*ng]*12*angArray[27]-rad[ip+3*ng]*9*angArray[72]+rad[ip+4*ng]*angArray[149];
         bas[11] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(6*angArray[28]+7*angArray[26])-rad[ip+3*ng]*(7*angArray[73]+angArray[71])+rad[ip+4*ng]*angArray[150];
         bas[12] = -rad[ip+ng]*4*angArray[6]+rad[ip+2*ng]*(10*angArray[27]+2*angArray[29])-rad[ip+3*ng]*(2*angArray[72]+5*angArray[74])+rad[ip+4*ng]*angArray[151];
         bas[13] = rad[ip+2*ng]*9*angArray[28]-rad[ip+3*ng]*(3*angArray[73]+3*angArray[75])+rad[ip+4*ng]*angArray[152];
         bas[14] = rad[ip+2*ng]*4*angArray[29]-rad[ip+3*ng]*(angArray[76]+4*angArray[74])+rad[ip+4*ng]*angArray[153];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*5*angArray[37]-rad[ip+3*ng]*(5*angArray[91]+angArray[86])+rad[ip+4*ng]*angArray[172];
         bas[1] = rad[ip+2*ng]*12*angArray[39]-rad[ip+3*ng]*(3*angArray[88]+4*angArray[95])+rad[ip+4*ng]*angArray[176];
         bas[2] = -rad[ip+ng]*4*angArray[10]+rad[ip+2*ng]*(4*angArray[38]+4*angArray[40]+angArray[35])-rad[ip+3*ng]*(angArray[87]+angArray[89]+4*angArray[96])+rad[ip+4*ng]*angArray[177];
         bas[3] = -rad[ip+ng]*6*angArray[12]+rad[ip+2*ng]*(15*angArray[42]+2*angArray[37])-rad[ip+3*ng]*(3*angArray[100]+5*angArray[91])+rad[ip+4*ng]*angArray[181];
         bas[4] = -rad[ip+ng]*9*angArray[11]+rad[ip+2*ng]*(3*angArray[36]+3*angArray[41]+9*angArray[43])-rad[ip+3*ng]*(angArray[90]+3*angArray[101]+3*angArray[92])+rad[ip+4*ng]*angArray[182];
         bas[5] = -rad[ip+ng]*6*angArray[12]+rad[ip+2*ng]*(6*angArray[42]+2*angArray[37]+3*angArray[44])-rad[ip+3*ng]*(2*angArray[91]+3*angArray[102]+angArray[93])+rad[ip+4*ng]*angArray[183];
         bas[6] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(14*angArray[46]+6*angArray[39])-rad[ip+3*ng]*(7*angArray[95]+2*angArray[106])+rad[ip+4*ng]*angArray[187];
         bas[7] = rad[ip]*4*angArray[1]-rad[ip+ng]*(10*angArray[13]+4*angArray[15]+2*angArray[10])+rad[ip+2*ng]*(10*angArray[47]+2*angArray[45]+5*angArray[38]+2*angArray[40])-rad[ip+3*ng]*(2*angArray[107]+angArray[94]+5*angArray[96])+rad[ip+4*ng]*angArray[188];
         bas[8] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(6*angArray[39]+6*angArray[48]+4*angArray[46])-rad[ip+3*ng]*(2*angArray[108]+2*angArray[95]+3*angArray[97])+rad[ip+4*ng]*angArray[189];
         bas[9] = -rad[ip+ng]*6*angArray[15]+rad[ip+2*ng]*(6*angArray[47]+2*angArray[49]+3*angArray[40])-rad[ip+3*ng]*(2*angArray[109]+3*angArray[96]+angArray[98])+rad[ip+4*ng]*angArray[190];
         bas[10] = -rad[ip+ng]*12*angArray[17]+rad[ip+2*ng]*(9*angArray[51]+12*angArray[42])-rad[ip+3*ng]*(9*angArray[100]+angArray[113])+rad[ip+4*ng]*angArray[194];
         bas[11] = rad[ip]*6*angArray[2]-rad[ip+ng]*(6*angArray[11]+7*angArray[16]+6*angArray[18])+rad[ip+2*ng]*(6*angArray[43]+angArray[50]+7*angArray[52]+7*angArray[41])-rad[ip+3*ng]*(7*angArray[101]+angArray[99]+angArray[114])+rad[ip+4*ng]*angArray[195];
         bas[12] = rad[ip]*4*angArray[3]-rad[ip+ng]*(2*angArray[19]+10*angArray[17]+4*angArray[12])+rad[ip+2*ng]*(10*angArray[42]+5*angArray[53]+2*angArray[51]+2*angArray[44])-rad[ip+3*ng]*(2*angArray[100]+5*angArray[102]+angArray[115])+rad[ip+4*ng]*angArray[196];
         bas[13] = -rad[ip+ng]*9*angArray[18]+rad[ip+2*ng]*(9*angArray[43]+3*angArray[54]+3*angArray[52])-rad[ip+3*ng]*(3*angArray[103]+3*angArray[101]+angArray[116])+rad[ip+4*ng]*angArray[197];
         bas[14] = -rad[ip+ng]*4*angArray[19]+rad[ip+2*ng]*(angArray[55]+4*angArray[53]+4*angArray[44])-rad[ip+3*ng]*(4*angArray[102]+angArray[104]+angArray[117])+rad[ip+4*ng]*angArray[198];
         bas[15] = rad[ip+2*ng]*20*angArray[46]-rad[ip+3*ng]*11*angArray[106]+rad[ip+4*ng]*angArray[202];
         bas[16] = -rad[ip+ng]*12*angArray[13]+rad[ip+2*ng]*(9*angArray[45]+12*angArray[47])-rad[ip+3*ng]*(9*angArray[107]+angArray[105])+rad[ip+4*ng]*angArray[203];
         bas[17] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(14*angArray[46]+6*angArray[48])-rad[ip+3*ng]*(7*angArray[108]+2*angArray[106])+rad[ip+4*ng]*angArray[204];
         bas[18] = -rad[ip+ng]*6*angArray[15]+rad[ip+2*ng]*(15*angArray[47]+2*angArray[49])-rad[ip+3*ng]*(3*angArray[107]+5*angArray[109])+rad[ip+4*ng]*angArray[205];
         bas[19] = rad[ip+2*ng]*12*angArray[48]-rad[ip+3*ng]*(3*angArray[110]+4*angArray[108])+rad[ip+4*ng]*angArray[206];
         bas[20] = rad[ip+2*ng]*5*angArray[49]-rad[ip+3*ng]*(5*angArray[109]+angArray[111])+rad[ip+4*ng]*angArray[207];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*6*angArray[58]-rad[ip+3*ng]*(6*angArray[127]+angArray[122])+rad[ip+4*ng]*angArray[227];
         bas[1] = rad[ip+2*ng]*15*angArray[60]-rad[ip+3*ng]*(3*angArray[124]+5*angArray[131])+rad[ip+4*ng]*angArray[231];
         bas[2] = -rad[ip+ng]*5*angArray[20]+rad[ip+2*ng]*(5*angArray[61]+5*angArray[59]+angArray[56])-rad[ip+3*ng]*(angArray[123]+5*angArray[132]+angArray[125])+rad[ip+4*ng]*angArray[232];
         bas[3] = -rad[ip+ng]*8*angArray[22]+rad[ip+2*ng]*(2*angArray[58]+20*angArray[63])-rad[ip+3*ng]*(4*angArray[136]+5*angArray[127])+rad[ip+4*ng]*angArray[236];
         bas[4] = -rad[ip+ng]*12*angArray[21]+rad[ip+2*ng]*(3*angArray[57]+4*angArray[62]+12*angArray[64])-rad[ip+3*ng]*(angArray[126]+3*angArray[128]+4*angArray[137])+rad[ip+4*ng]*angArray[237];
         bas[5] = -rad[ip+ng]*8*angArray[22]+rad[ip+2*ng]*(8*angArray[63]+4*angArray[65]+2*angArray[58])-rad[ip+3*ng]*(2*angArray[127]+angArray[129]+4*angArray[138])+rad[ip+4*ng]*angArray[238];
         bas[6] = -rad[ip+ng]*18*angArray[24]+rad[ip+2*ng]*(21*angArray[67]+6*angArray[60])-rad[ip+3*ng]*(7*angArray[131]+3*angArray[142])+rad[ip+4*ng]*angArray[242];
         bas[7] = rad[ip]*6*angArray[4]-rad[ip+ng]*(15*angArray[23]+2*angArray[20]+6*angArray[25])+rad[ip+2*ng]*(5*angArray[59]+15*angArray[68]+3*angArray[66]+2*angArray[61])-rad[ip+3*ng]*(angArray[130]+3*angArray[143]+5*angArray[132])+rad[ip+4*ng]*angArray[243];
         bas[8] = -rad[ip+ng]*18*angArray[24]+rad[ip+2*ng]*(6*angArray[67]+9*angArray[69]+6*angArray[60])-rad[ip+3*ng]*(2*angArray[131]+3*angArray[133]+3*angArray[144])+rad[ip+4*ng]*angArray[244];
         bas[9] = -rad[ip+ng]*9*angArray[25]+rad[ip+2*ng]*(9*angArray[68]+3*angArray[70]+3*angArray[61])-rad[ip+3*ng]*(3*angArray[132]+3*angArray[145]+angArray[134])+rad[ip+4*ng]*angArray[245];
         bas[10] = -rad[ip+ng]*24*angArray[27]+rad[ip+2*ng]*(18*angArray[72]+12*angArray[63])-rad[ip+3*ng]*(9*angArray[136]+2*angArray[149])+rad[ip+4*ng]*angArray[249];
         bas[11] = rad[ip]*12*angArray[5]-rad[ip+ng]*(12*angArray[28]+14*angArray[26]+6*angArray[21])+rad[ip+2*ng]*(2*angArray[71]+14*angArray[73]+6*angArray[64]+7*angArray[62])-rad[ip+3*ng]*(7*angArray[137]+2*angArray[150]+angArray[135])+rad[ip+4*ng]*angArray[250];
         bas[12] = rad[ip]*8*angArray[6]-rad[ip+ng]*(20*angArray[27]+4*angArray[29]+4*angArray[22])+rad[ip+2*ng]*(10*angArray[74]+4*angArray[72]+10*angArray[63]+2*angArray[65])-rad[ip+3*ng]*(2*angArray[151]+2*angArray[136]+5*angArray[138])+rad[ip+4*ng]*angArray[251];
         bas[13] = -rad[ip+ng]*18*angArray[28]+rad[ip+2*ng]*(6*angArray[75]+6*angArray[73]+9*angArray[64])-rad[ip+3*ng]*(2*angArray[152]+3*angArray[137]+3*angArray[139])+rad[ip+4*ng]*angArray[252];
         bas[14] = -rad[ip+ng]*8*angArray[29]+rad[ip+2*ng]*(8*angArray[74]+2*angArray[76]+4*angArray[65])-rad[ip+3*ng]*(2*angArray[153]+4*angArray[138]+angArray[140])+rad[ip+4*ng]*angArray[253];
         bas[15] = -rad[ip+ng]*20*angArray[31]+rad[ip+2*ng]*(11*angArray[78]+20*angArray[67])-rad[ip+3*ng]*(11*angArray[142]+angArray[157])+rad[ip+4*ng]*angArray[257];
         bas[16] = rad[ip]*12*angArray[7]-rad[ip+ng]*(12*angArray[23]+9*angArray[30]+12*angArray[32])+rad[ip+2*ng]*(9*angArray[66]+angArray[77]+9*angArray[79]+12*angArray[68])-rad[ip+3*ng]*(9*angArray[143]+angArray[141]+angArray[158])+rad[ip+4*ng]*angArray[258];
         bas[17] = rad[ip]*12*angArray[8]-rad[ip+ng]*(12*angArray[24]+14*angArray[31]+6*angArray[33])+rad[ip+2*ng]*(2*angArray[78]+7*angArray[80]+14*angArray[67]+6*angArray[69])-rad[ip+3*ng]*(7*angArray[144]+2*angArray[142]+angArray[159])+rad[ip+4*ng]*angArray[259];
         bas[18] = rad[ip]*6*angArray[9]-rad[ip+ng]*(2*angArray[34]+6*angArray[25]+15*angArray[32])+rad[ip+2*ng]*(5*angArray[81]+15*angArray[68]+3*angArray[79]+2*angArray[70])-rad[ip+3*ng]*(3*angArray[143]+5*angArray[145]+angArray[160])+rad[ip+4*ng]*angArray[260];
         bas[19] = -rad[ip+ng]*12*angArray[33]+rad[ip+2*ng]*(12*angArray[69]+3*angArray[82]+4*angArray[80])-rad[ip+3*ng]*(3*angArray[146]+4*angArray[144]+angArray[161])+rad[ip+4*ng]*angArray[261];
         bas[20] = -rad[ip+ng]*5*angArray[34]+rad[ip+2*ng]*(angArray[83]+5*angArray[81]+5*angArray[70])-rad[ip+3*ng]*(5*angArray[145]+angArray[147]+angArray[162])+rad[ip+4*ng]*angArray[262];
         bas[21] = rad[ip+2*ng]*30*angArray[72]-rad[ip+3*ng]*13*angArray[149]+rad[ip+4*ng]*angArray[266];
         bas[22] = -rad[ip+ng]*20*angArray[26]+rad[ip+2*ng]*(20*angArray[73]+11*angArray[71])-rad[ip+3*ng]*(angArray[148]+11*angArray[150])+rad[ip+4*ng]*angArray[267];
         bas[23] = -rad[ip+ng]*24*angArray[27]+rad[ip+2*ng]*(18*angArray[72]+12*angArray[74])-rad[ip+3*ng]*(9*angArray[151]+2*angArray[149])+rad[ip+4*ng]*angArray[268];
         bas[24] = -rad[ip+ng]*18*angArray[28]+rad[ip+2*ng]*(21*angArray[73]+6*angArray[75])-rad[ip+3*ng]*(7*angArray[152]+3*angArray[150])+rad[ip+4*ng]*angArray[269];
         bas[25] = -rad[ip+ng]*8*angArray[29]+rad[ip+2*ng]*(20*angArray[74]+2*angArray[76])-rad[ip+3*ng]*(4*angArray[151]+5*angArray[153])+rad[ip+4*ng]*angArray[270];
         bas[26] = rad[ip+2*ng]*15*angArray[75]-rad[ip+3*ng]*(3*angArray[154]+5*angArray[152])+rad[ip+4*ng]*angArray[271];
         bas[27] = rad[ip+2*ng]*6*angArray[76]-rad[ip+3*ng]*(angArray[155]+6*angArray[153])+rad[ip+4*ng]*angArray[272];
      }

   }


   // now we do derivatives for the given basis set to YYYZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[8]+rad[ip+4*ng]*angArray[31];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[46];
         bas[1] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*6*angArray[17]+rad[ip+4*ng]*angArray[51];
         bas[2] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*(3*angArray[18]+angArray[16])+rad[ip+4*ng]*angArray[52];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[24]+rad[ip+4*ng]*angArray[67];
         bas[1] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*6*angArray[27]+rad[ip+4*ng]*angArray[72];
         bas[2] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*(3*angArray[28]+angArray[26])+rad[ip+4*ng]*angArray[73];
         bas[3] = rad[ip+2*ng]*12*angArray[8]-rad[ip+3*ng]*9*angArray[31]+rad[ip+4*ng]*angArray[78];
         bas[4] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*(6*angArray[7]+3*angArray[9])-rad[ip+3*ng]*(angArray[30]+6*angArray[32])+rad[ip+4*ng]*angArray[79];
         bas[5] = rad[ip+2*ng]*6*angArray[8]-rad[ip+3*ng]*(2*angArray[31]+3*angArray[33])+rad[ip+4*ng]*angArray[80];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[39]+rad[ip+4*ng]*angArray[95];
         bas[1] = rad[ip+2*ng]*3*angArray[12]-rad[ip+3*ng]*6*angArray[42]+rad[ip+4*ng]*angArray[100];
         bas[2] = rad[ip+2*ng]*3*angArray[11]-rad[ip+3*ng]*(3*angArray[43]+angArray[41])+rad[ip+4*ng]*angArray[101];
         bas[3] = rad[ip+2*ng]*12*angArray[14]-rad[ip+3*ng]*9*angArray[46]+rad[ip+4*ng]*angArray[106];
         bas[4] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*(6*angArray[13]+3*angArray[15])-rad[ip+3*ng]*(angArray[45]+6*angArray[47])+rad[ip+4*ng]*angArray[107];
         bas[5] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[46]+3*angArray[48])+rad[ip+4*ng]*angArray[108];
         bas[6] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*27*angArray[17]-rad[ip+3*ng]*12*angArray[51]+rad[ip+4*ng]*angArray[113];
         bas[7] = -rad[ip+ng]*12*angArray[2]+rad[ip+2*ng]*(12*angArray[18]+9*angArray[16])-rad[ip+3*ng]*(9*angArray[52]+angArray[50])+rad[ip+4*ng]*angArray[114];
         bas[8] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*(12*angArray[17]+3*angArray[19])-rad[ip+3*ng]*(2*angArray[51]+6*angArray[53])+rad[ip+4*ng]*angArray[115];
         bas[9] = rad[ip+2*ng]*9*angArray[18]-rad[ip+3*ng]*(3*angArray[52]+3*angArray[54])+rad[ip+4*ng]*angArray[116];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[60]+rad[ip+4*ng]*angArray[131];
         bas[1] = rad[ip+2*ng]*3*angArray[22]-rad[ip+3*ng]*6*angArray[63]+rad[ip+4*ng]*angArray[136];
         bas[2] = rad[ip+2*ng]*3*angArray[21]-rad[ip+3*ng]*(3*angArray[64]+angArray[62])+rad[ip+4*ng]*angArray[137];
         bas[3] = rad[ip+2*ng]*12*angArray[24]-rad[ip+3*ng]*9*angArray[67]+rad[ip+4*ng]*angArray[142];
         bas[4] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*(6*angArray[23]+3*angArray[25])-rad[ip+3*ng]*(angArray[66]+6*angArray[68])+rad[ip+4*ng]*angArray[143];
         bas[5] = rad[ip+2*ng]*6*angArray[24]-rad[ip+3*ng]*(2*angArray[67]+3*angArray[69])+rad[ip+4*ng]*angArray[144];
         bas[6] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*27*angArray[27]-rad[ip+3*ng]*12*angArray[72]+rad[ip+4*ng]*angArray[149];
         bas[7] = -rad[ip+ng]*12*angArray[5]+rad[ip+2*ng]*(12*angArray[28]+9*angArray[26])-rad[ip+3*ng]*(9*angArray[73]+angArray[71])+rad[ip+4*ng]*angArray[150];
         bas[8] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(12*angArray[27]+3*angArray[29])-rad[ip+3*ng]*(2*angArray[72]+6*angArray[74])+rad[ip+4*ng]*angArray[151];
         bas[9] = rad[ip+2*ng]*9*angArray[28]-rad[ip+3*ng]*(3*angArray[73]+3*angArray[75])+rad[ip+4*ng]*angArray[152];
         bas[10] = -rad[ip+ng]*24*angArray[8]+rad[ip+2*ng]*48*angArray[31]-rad[ip+3*ng]*15*angArray[78]+rad[ip+4*ng]*angArray[157];
         bas[11] = rad[ip]*6*angArray[0]-rad[ip+ng]*(27*angArray[7]+6*angArray[9])+rad[ip+2*ng]*(12*angArray[30]+27*angArray[32])-rad[ip+3*ng]*(12*angArray[79]+angArray[77])+rad[ip+4*ng]*angArray[158];
         bas[12] = -rad[ip+ng]*24*angArray[8]+rad[ip+2*ng]*(12*angArray[33]+18*angArray[31])-rad[ip+3*ng]*(9*angArray[80]+2*angArray[78])+rad[ip+4*ng]*angArray[159];
         bas[13] = -rad[ip+ng]*9*angArray[9]+rad[ip+2*ng]*(18*angArray[32]+3*angArray[34])-rad[ip+3*ng]*(3*angArray[79]+6*angArray[81])+rad[ip+4*ng]*angArray[160];
         bas[14] = rad[ip+2*ng]*12*angArray[33]-rad[ip+3*ng]*(3*angArray[82]+4*angArray[80])+rad[ip+4*ng]*angArray[161];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[88]+rad[ip+4*ng]*angArray[176];
         bas[1] = rad[ip+2*ng]*3*angArray[37]-rad[ip+3*ng]*6*angArray[91]+rad[ip+4*ng]*angArray[181];
         bas[2] = rad[ip+2*ng]*3*angArray[36]-rad[ip+3*ng]*(3*angArray[92]+angArray[90])+rad[ip+4*ng]*angArray[182];
         bas[3] = rad[ip+2*ng]*12*angArray[39]-rad[ip+3*ng]*9*angArray[95]+rad[ip+4*ng]*angArray[187];
         bas[4] = -rad[ip+ng]*3*angArray[10]+rad[ip+2*ng]*(6*angArray[38]+3*angArray[40])-rad[ip+3*ng]*(angArray[94]+6*angArray[96])+rad[ip+4*ng]*angArray[188];
         bas[5] = rad[ip+2*ng]*6*angArray[39]-rad[ip+3*ng]*(2*angArray[95]+3*angArray[97])+rad[ip+4*ng]*angArray[189];
         bas[6] = -rad[ip+ng]*6*angArray[12]+rad[ip+2*ng]*27*angArray[42]-rad[ip+3*ng]*12*angArray[100]+rad[ip+4*ng]*angArray[194];
         bas[7] = -rad[ip+ng]*12*angArray[11]+rad[ip+2*ng]*(12*angArray[43]+9*angArray[41])-rad[ip+3*ng]*(9*angArray[101]+angArray[99])+rad[ip+4*ng]*angArray[195];
         bas[8] = -rad[ip+ng]*6*angArray[12]+rad[ip+2*ng]*(12*angArray[42]+3*angArray[44])-rad[ip+3*ng]*(2*angArray[100]+6*angArray[102])+rad[ip+4*ng]*angArray[196];
         bas[9] = rad[ip+2*ng]*9*angArray[43]-rad[ip+3*ng]*(3*angArray[101]+3*angArray[103])+rad[ip+4*ng]*angArray[197];
         bas[10] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*48*angArray[46]-rad[ip+3*ng]*15*angArray[106]+rad[ip+4*ng]*angArray[202];
         bas[11] = rad[ip]*6*angArray[1]-rad[ip+ng]*(27*angArray[13]+6*angArray[15])+rad[ip+2*ng]*(12*angArray[45]+27*angArray[47])-rad[ip+3*ng]*(12*angArray[107]+angArray[105])+rad[ip+4*ng]*angArray[203];
         bas[12] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*(12*angArray[48]+18*angArray[46])-rad[ip+3*ng]*(9*angArray[108]+2*angArray[106])+rad[ip+4*ng]*angArray[204];
         bas[13] = -rad[ip+ng]*9*angArray[15]+rad[ip+2*ng]*(18*angArray[47]+3*angArray[49])-rad[ip+3*ng]*(3*angArray[107]+6*angArray[109])+rad[ip+4*ng]*angArray[205];
         bas[14] = rad[ip+2*ng]*12*angArray[48]-rad[ip+3*ng]*(3*angArray[110]+4*angArray[108])+rad[ip+4*ng]*angArray[206];
         bas[15] = -rad[ip+ng]*60*angArray[17]+rad[ip+2*ng]*75*angArray[51]-rad[ip+3*ng]*18*angArray[113]+rad[ip+4*ng]*angArray[211];
         bas[16] = rad[ip]*24*angArray[2]-rad[ip+ng]*(48*angArray[16]+24*angArray[18])+rad[ip+2*ng]*(15*angArray[50]+48*angArray[52])-rad[ip+3*ng]*(angArray[112]+15*angArray[114])+rad[ip+4*ng]*angArray[212];
         bas[17] = rad[ip]*12*angArray[3]-rad[ip+ng]*(54*angArray[17]+6*angArray[19])+rad[ip+2*ng]*(24*angArray[51]+27*angArray[53])-rad[ip+3*ng]*(12*angArray[115]+2*angArray[113])+rad[ip+4*ng]*angArray[213];
         bas[18] = -rad[ip+ng]*36*angArray[18]+rad[ip+2*ng]*(12*angArray[54]+27*angArray[52])-rad[ip+3*ng]*(9*angArray[116]+3*angArray[114])+rad[ip+4*ng]*angArray[214];
         bas[19] = -rad[ip+ng]*12*angArray[19]+rad[ip+2*ng]*(3*angArray[55]+24*angArray[53])-rad[ip+3*ng]*(4*angArray[115]+6*angArray[117])+rad[ip+4*ng]*angArray[215];
         bas[20] = rad[ip+2*ng]*15*angArray[54]-rad[ip+3*ng]*(3*angArray[118]+5*angArray[116])+rad[ip+4*ng]*angArray[216];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[124]+rad[ip+4*ng]*angArray[231];
         bas[1] = rad[ip+2*ng]*3*angArray[58]-rad[ip+3*ng]*6*angArray[127]+rad[ip+4*ng]*angArray[236];
         bas[2] = rad[ip+2*ng]*3*angArray[57]-rad[ip+3*ng]*(3*angArray[128]+angArray[126])+rad[ip+4*ng]*angArray[237];
         bas[3] = rad[ip+2*ng]*12*angArray[60]-rad[ip+3*ng]*9*angArray[131]+rad[ip+4*ng]*angArray[242];
         bas[4] = -rad[ip+ng]*3*angArray[20]+rad[ip+2*ng]*(6*angArray[59]+3*angArray[61])-rad[ip+3*ng]*(angArray[130]+6*angArray[132])+rad[ip+4*ng]*angArray[243];
         bas[5] = rad[ip+2*ng]*6*angArray[60]-rad[ip+3*ng]*(2*angArray[131]+3*angArray[133])+rad[ip+4*ng]*angArray[244];
         bas[6] = -rad[ip+ng]*6*angArray[22]+rad[ip+2*ng]*27*angArray[63]-rad[ip+3*ng]*12*angArray[136]+rad[ip+4*ng]*angArray[249];
         bas[7] = -rad[ip+ng]*12*angArray[21]+rad[ip+2*ng]*(12*angArray[64]+9*angArray[62])-rad[ip+3*ng]*(9*angArray[137]+angArray[135])+rad[ip+4*ng]*angArray[250];
         bas[8] = -rad[ip+ng]*6*angArray[22]+rad[ip+2*ng]*(12*angArray[63]+3*angArray[65])-rad[ip+3*ng]*(2*angArray[136]+6*angArray[138])+rad[ip+4*ng]*angArray[251];
         bas[9] = rad[ip+2*ng]*9*angArray[64]-rad[ip+3*ng]*(3*angArray[137]+3*angArray[139])+rad[ip+4*ng]*angArray[252];
         bas[10] = -rad[ip+ng]*24*angArray[24]+rad[ip+2*ng]*48*angArray[67]-rad[ip+3*ng]*15*angArray[142]+rad[ip+4*ng]*angArray[257];
         bas[11] = rad[ip]*6*angArray[4]-rad[ip+ng]*(27*angArray[23]+6*angArray[25])+rad[ip+2*ng]*(12*angArray[66]+27*angArray[68])-rad[ip+3*ng]*(12*angArray[143]+angArray[141])+rad[ip+4*ng]*angArray[258];
         bas[12] = -rad[ip+ng]*24*angArray[24]+rad[ip+2*ng]*(12*angArray[69]+18*angArray[67])-rad[ip+3*ng]*(9*angArray[144]+2*angArray[142])+rad[ip+4*ng]*angArray[259];
         bas[13] = -rad[ip+ng]*9*angArray[25]+rad[ip+2*ng]*(18*angArray[68]+3*angArray[70])-rad[ip+3*ng]*(3*angArray[143]+6*angArray[145])+rad[ip+4*ng]*angArray[260];
         bas[14] = rad[ip+2*ng]*12*angArray[69]-rad[ip+3*ng]*(3*angArray[146]+4*angArray[144])+rad[ip+4*ng]*angArray[261];
         bas[15] = -rad[ip+ng]*60*angArray[27]+rad[ip+2*ng]*75*angArray[72]-rad[ip+3*ng]*18*angArray[149]+rad[ip+4*ng]*angArray[266];
         bas[16] = rad[ip]*24*angArray[5]-rad[ip+ng]*(48*angArray[26]+24*angArray[28])+rad[ip+2*ng]*(15*angArray[71]+48*angArray[73])-rad[ip+3*ng]*(angArray[148]+15*angArray[150])+rad[ip+4*ng]*angArray[267];
         bas[17] = rad[ip]*12*angArray[6]-rad[ip+ng]*(54*angArray[27]+6*angArray[29])+rad[ip+2*ng]*(24*angArray[72]+27*angArray[74])-rad[ip+3*ng]*(12*angArray[151]+2*angArray[149])+rad[ip+4*ng]*angArray[268];
         bas[18] = -rad[ip+ng]*36*angArray[28]+rad[ip+2*ng]*(12*angArray[75]+27*angArray[73])-rad[ip+3*ng]*(9*angArray[152]+3*angArray[150])+rad[ip+4*ng]*angArray[269];
         bas[19] = -rad[ip+ng]*12*angArray[29]+rad[ip+2*ng]*(3*angArray[76]+24*angArray[74])-rad[ip+3*ng]*(4*angArray[151]+6*angArray[153])+rad[ip+4*ng]*angArray[270];
         bas[20] = rad[ip+2*ng]*15*angArray[75]-rad[ip+3*ng]*(3*angArray[154]+5*angArray[152])+rad[ip+4*ng]*angArray[271];
         bas[21] = -rad[ip+ng]*120*angArray[31]+rad[ip+2*ng]*108*angArray[78]-rad[ip+3*ng]*21*angArray[157]+rad[ip+4*ng]*angArray[276];
         bas[22] = rad[ip]*60*angArray[7]-rad[ip+ng]*(75*angArray[30]+60*angArray[32])+rad[ip+2*ng]*(75*angArray[79]+18*angArray[77])-rad[ip+3*ng]*(18*angArray[158]+angArray[156])+rad[ip+4*ng]*angArray[277];
         bas[23] = rad[ip]*48*angArray[8]-rad[ip+ng]*(96*angArray[31]+24*angArray[33])+rad[ip+2*ng]*(48*angArray[80]+30*angArray[78])-rad[ip+3*ng]*(2*angArray[157]+15*angArray[159])+rad[ip+4*ng]*angArray[278];
         bas[24] = rad[ip]*18*angArray[9]-rad[ip+ng]*(81*angArray[32]+6*angArray[34])+rad[ip+2*ng]*(36*angArray[79]+27*angArray[81])-rad[ip+3*ng]*(12*angArray[160]+3*angArray[158])+rad[ip+4*ng]*angArray[279];
         bas[25] = -rad[ip+ng]*48*angArray[33]+rad[ip+2*ng]*(36*angArray[80]+12*angArray[82])-rad[ip+3*ng]*(9*angArray[161]+4*angArray[159])+rad[ip+4*ng]*angArray[280];
         bas[26] = -rad[ip+ng]*15*angArray[34]+rad[ip+2*ng]*(30*angArray[81]+3*angArray[83])-rad[ip+3*ng]*(5*angArray[160]+6*angArray[162])+rad[ip+4*ng]*angArray[281];
         bas[27] = rad[ip+2*ng]*18*angArray[82]-rad[ip+3*ng]*(3*angArray[163]+6*angArray[161])+rad[ip+4*ng]*angArray[282];
      }

   }


   // now we do derivatives for the given basis set to XXZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[0]-rad[ip+3*ng]*(angArray[4]+angArray[9])+rad[ip+4*ng]*angArray[25];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*(angArray[10]+3*angArray[15])+rad[ip+4*ng]*angArray[40];
         bas[1] = rad[ip+2*ng]*angArray[2]-rad[ip+3*ng]*(angArray[11]+angArray[18])+rad[ip+4*ng]*angArray[43];
         bas[2] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*(3*angArray[12]+angArray[19])+rad[ip+4*ng]*angArray[44];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*2*angArray[0]+rad[ip+2*ng]*(5*angArray[4]+2*angArray[9])-rad[ip+3*ng]*(5*angArray[25]+angArray[20])+rad[ip+4*ng]*angArray[61];
         bas[1] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*(angArray[21]+3*angArray[28])+rad[ip+4*ng]*angArray[64];
         bas[2] = rad[ip+2*ng]*9*angArray[6]-rad[ip+3*ng]*(3*angArray[22]+3*angArray[29])+rad[ip+4*ng]*angArray[65];
         bas[3] = rad[ip+2*ng]*angArray[7]-rad[ip+3*ng]*(angArray[23]+angArray[32])+rad[ip+4*ng]*angArray[68];
         bas[4] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*(3*angArray[24]+angArray[33])+rad[ip+4*ng]*angArray[69];
         bas[5] = -rad[ip+ng]*2*angArray[0]+rad[ip+2*ng]*(5*angArray[9]+2*angArray[4])-rad[ip+3*ng]*(5*angArray[25]+angArray[34])+rad[ip+4*ng]*angArray[70];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*(7*angArray[10]+6*angArray[15])-rad[ip+3*ng]*(7*angArray[40]+angArray[35])+rad[ip+4*ng]*angArray[89];
         bas[1] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(5*angArray[11]+2*angArray[18])-rad[ip+3*ng]*(5*angArray[43]+angArray[36])+rad[ip+4*ng]*angArray[92];
         bas[2] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*(15*angArray[12]+2*angArray[19])-rad[ip+3*ng]*(5*angArray[44]+3*angArray[37])+rad[ip+4*ng]*angArray[93];
         bas[3] = rad[ip+2*ng]*3*angArray[13]-rad[ip+3*ng]*(angArray[38]+3*angArray[47])+rad[ip+4*ng]*angArray[96];
         bas[4] = rad[ip+2*ng]*9*angArray[14]-rad[ip+3*ng]*(3*angArray[39]+3*angArray[48])+rad[ip+4*ng]*angArray[97];
         bas[5] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*(2*angArray[10]+15*angArray[15])-rad[ip+3*ng]*(5*angArray[40]+3*angArray[49])+rad[ip+4*ng]*angArray[98];
         bas[6] = rad[ip+2*ng]*angArray[16]-rad[ip+3*ng]*(angArray[41]+angArray[52])+rad[ip+4*ng]*angArray[101];
         bas[7] = rad[ip+2*ng]*3*angArray[17]-rad[ip+3*ng]*(3*angArray[42]+angArray[53])+rad[ip+4*ng]*angArray[102];
         bas[8] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(5*angArray[18]+2*angArray[11])-rad[ip+3*ng]*(5*angArray[43]+angArray[54])+rad[ip+4*ng]*angArray[103];
         bas[9] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*(7*angArray[19]+6*angArray[12])-rad[ip+3*ng]*(7*angArray[44]+angArray[55])+rad[ip+4*ng]*angArray[104];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*12*angArray[4]+rad[ip+2*ng]*(9*angArray[20]+12*angArray[25])-rad[ip+3*ng]*(angArray[56]+9*angArray[61])+rad[ip+4*ng]*angArray[125];
         bas[1] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(7*angArray[21]+6*angArray[28])-rad[ip+3*ng]*(7*angArray[64]+angArray[57])+rad[ip+4*ng]*angArray[128];
         bas[2] = -rad[ip+ng]*18*angArray[6]+rad[ip+2*ng]*(21*angArray[22]+6*angArray[29])-rad[ip+3*ng]*(7*angArray[65]+3*angArray[58])+rad[ip+4*ng]*angArray[129];
         bas[3] = -rad[ip+ng]*2*angArray[7]+rad[ip+2*ng]*(5*angArray[23]+2*angArray[32])-rad[ip+3*ng]*(5*angArray[68]+angArray[59])+rad[ip+4*ng]*angArray[132];
         bas[4] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(15*angArray[24]+2*angArray[33])-rad[ip+3*ng]*(5*angArray[69]+3*angArray[60])+rad[ip+4*ng]*angArray[133];
         bas[5] = rad[ip]*4*angArray[0]-rad[ip+ng]*(10*angArray[4]+10*angArray[9])+rad[ip+2*ng]*(25*angArray[25]+2*angArray[20]+2*angArray[34])-rad[ip+3*ng]*(5*angArray[70]+5*angArray[61])+rad[ip+4*ng]*angArray[134];
         bas[6] = rad[ip+2*ng]*3*angArray[26]-rad[ip+3*ng]*(angArray[62]+3*angArray[73])+rad[ip+4*ng]*angArray[137];
         bas[7] = rad[ip+2*ng]*9*angArray[27]-rad[ip+3*ng]*(3*angArray[63]+3*angArray[74])+rad[ip+4*ng]*angArray[138];
         bas[8] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(2*angArray[21]+15*angArray[28])-rad[ip+3*ng]*(5*angArray[64]+3*angArray[75])+rad[ip+4*ng]*angArray[139];
         bas[9] = -rad[ip+ng]*18*angArray[6]+rad[ip+2*ng]*(21*angArray[29]+6*angArray[22])-rad[ip+3*ng]*(7*angArray[65]+3*angArray[76])+rad[ip+4*ng]*angArray[140];
         bas[10] = rad[ip+2*ng]*angArray[30]-rad[ip+3*ng]*(angArray[66]+angArray[79])+rad[ip+4*ng]*angArray[143];
         bas[11] = rad[ip+2*ng]*3*angArray[31]-rad[ip+3*ng]*(3*angArray[67]+angArray[80])+rad[ip+4*ng]*angArray[144];
         bas[12] = -rad[ip+ng]*2*angArray[7]+rad[ip+2*ng]*(5*angArray[32]+2*angArray[23])-rad[ip+3*ng]*(5*angArray[68]+angArray[81])+rad[ip+4*ng]*angArray[145];
         bas[13] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(7*angArray[33]+6*angArray[24])-rad[ip+3*ng]*(7*angArray[69]+angArray[82])+rad[ip+4*ng]*angArray[146];
         bas[14] = -rad[ip+ng]*12*angArray[9]+rad[ip+2*ng]*(9*angArray[34]+12*angArray[25])-rad[ip+3*ng]*(9*angArray[70]+angArray[83])+rad[ip+4*ng]*angArray[147];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*20*angArray[10]+rad[ip+2*ng]*(20*angArray[40]+11*angArray[35])-rad[ip+3*ng]*(11*angArray[89]+angArray[84])+rad[ip+4*ng]*angArray[170];
         bas[1] = -rad[ip+ng]*12*angArray[11]+rad[ip+2*ng]*(9*angArray[36]+12*angArray[43])-rad[ip+3*ng]*(angArray[85]+9*angArray[92])+rad[ip+4*ng]*angArray[173];
         bas[2] = -rad[ip+ng]*36*angArray[12]+rad[ip+2*ng]*(27*angArray[37]+12*angArray[44])-rad[ip+3*ng]*(3*angArray[86]+9*angArray[93])+rad[ip+4*ng]*angArray[174];
         bas[3] = -rad[ip+ng]*6*angArray[13]+rad[ip+2*ng]*(7*angArray[38]+6*angArray[47])-rad[ip+3*ng]*(7*angArray[96]+angArray[87])+rad[ip+4*ng]*angArray[177];
         bas[4] = -rad[ip+ng]*18*angArray[14]+rad[ip+2*ng]*(21*angArray[39]+6*angArray[48])-rad[ip+3*ng]*(7*angArray[97]+3*angArray[88])+rad[ip+4*ng]*angArray[178];
         bas[5] = rad[ip]*12*angArray[1]-rad[ip+ng]*(14*angArray[10]+30*angArray[15])+rad[ip+2*ng]*(35*angArray[40]+2*angArray[35]+6*angArray[49])-rad[ip+3*ng]*(7*angArray[98]+5*angArray[89])+rad[ip+4*ng]*angArray[179];
         bas[6] = -rad[ip+ng]*2*angArray[16]+rad[ip+2*ng]*(5*angArray[41]+2*angArray[52])-rad[ip+3*ng]*(5*angArray[101]+angArray[90])+rad[ip+4*ng]*angArray[182];
         bas[7] = -rad[ip+ng]*6*angArray[17]+rad[ip+2*ng]*(15*angArray[42]+2*angArray[53])-rad[ip+3*ng]*(5*angArray[102]+3*angArray[91])+rad[ip+4*ng]*angArray[183];
         bas[8] = rad[ip]*4*angArray[2]-rad[ip+ng]*(10*angArray[11]+10*angArray[18])+rad[ip+2*ng]*(25*angArray[43]+2*angArray[36]+2*angArray[54])-rad[ip+3*ng]*(5*angArray[103]+5*angArray[92])+rad[ip+4*ng]*angArray[184];
         bas[9] = rad[ip]*12*angArray[3]-rad[ip+ng]*(30*angArray[12]+14*angArray[19])+rad[ip+2*ng]*(35*angArray[44]+6*angArray[37]+2*angArray[55])-rad[ip+3*ng]*(5*angArray[104]+7*angArray[93])+rad[ip+4*ng]*angArray[185];
         bas[10] = rad[ip+2*ng]*3*angArray[45]-rad[ip+3*ng]*(angArray[94]+3*angArray[107])+rad[ip+4*ng]*angArray[188];
         bas[11] = rad[ip+2*ng]*9*angArray[46]-rad[ip+3*ng]*(3*angArray[95]+3*angArray[108])+rad[ip+4*ng]*angArray[189];
         bas[12] = -rad[ip+ng]*6*angArray[13]+rad[ip+2*ng]*(2*angArray[38]+15*angArray[47])-rad[ip+3*ng]*(5*angArray[96]+3*angArray[109])+rad[ip+4*ng]*angArray[190];
         bas[13] = -rad[ip+ng]*18*angArray[14]+rad[ip+2*ng]*(21*angArray[48]+6*angArray[39])-rad[ip+3*ng]*(7*angArray[97]+3*angArray[110])+rad[ip+4*ng]*angArray[191];
         bas[14] = -rad[ip+ng]*36*angArray[15]+rad[ip+2*ng]*(27*angArray[49]+12*angArray[40])-rad[ip+3*ng]*(9*angArray[98]+3*angArray[111])+rad[ip+4*ng]*angArray[192];
         bas[15] = rad[ip+2*ng]*angArray[50]-rad[ip+3*ng]*(angArray[99]+angArray[114])+rad[ip+4*ng]*angArray[195];
         bas[16] = rad[ip+2*ng]*3*angArray[51]-rad[ip+3*ng]*(3*angArray[100]+angArray[115])+rad[ip+4*ng]*angArray[196];
         bas[17] = -rad[ip+ng]*2*angArray[16]+rad[ip+2*ng]*(5*angArray[52]+2*angArray[41])-rad[ip+3*ng]*(5*angArray[101]+angArray[116])+rad[ip+4*ng]*angArray[197];
         bas[18] = -rad[ip+ng]*6*angArray[17]+rad[ip+2*ng]*(7*angArray[53]+6*angArray[42])-rad[ip+3*ng]*(7*angArray[102]+angArray[117])+rad[ip+4*ng]*angArray[198];
         bas[19] = -rad[ip+ng]*12*angArray[18]+rad[ip+2*ng]*(9*angArray[54]+12*angArray[43])-rad[ip+3*ng]*(9*angArray[103]+angArray[118])+rad[ip+4*ng]*angArray[199];
         bas[20] = -rad[ip+ng]*20*angArray[19]+rad[ip+2*ng]*(20*angArray[44]+11*angArray[55])-rad[ip+3*ng]*(11*angArray[104]+angArray[119])+rad[ip+4*ng]*angArray[200];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+ng]*30*angArray[20]+rad[ip+2*ng]*(30*angArray[61]+13*angArray[56])-rad[ip+3*ng]*(angArray[120]+13*angArray[125])+rad[ip+4*ng]*angArray[225];
         bas[1] = -rad[ip+ng]*20*angArray[21]+rad[ip+2*ng]*(20*angArray[64]+11*angArray[57])-rad[ip+3*ng]*(11*angArray[128]+angArray[121])+rad[ip+4*ng]*angArray[228];
         bas[2] = -rad[ip+ng]*60*angArray[22]+rad[ip+2*ng]*(20*angArray[65]+33*angArray[58])-rad[ip+3*ng]*(3*angArray[122]+11*angArray[129])+rad[ip+4*ng]*angArray[229];
         bas[3] = -rad[ip+ng]*12*angArray[23]+rad[ip+2*ng]*(9*angArray[59]+12*angArray[68])-rad[ip+3*ng]*(angArray[123]+9*angArray[132])+rad[ip+4*ng]*angArray[232];
         bas[4] = -rad[ip+ng]*36*angArray[24]+rad[ip+2*ng]*(27*angArray[60]+12*angArray[69])-rad[ip+3*ng]*(3*angArray[124]+9*angArray[133])+rad[ip+4*ng]*angArray[233];
         bas[5] = rad[ip]*24*angArray[4]-rad[ip+ng]*(60*angArray[25]+18*angArray[20])+rad[ip+2*ng]*(45*angArray[61]+12*angArray[70]+2*angArray[56])-rad[ip+3*ng]*(5*angArray[125]+9*angArray[134])+rad[ip+4*ng]*angArray[234];
         bas[6] = -rad[ip+ng]*6*angArray[26]+rad[ip+2*ng]*(7*angArray[62]+6*angArray[73])-rad[ip+3*ng]*(7*angArray[137]+angArray[126])+rad[ip+4*ng]*angArray[237];
         bas[7] = -rad[ip+ng]*18*angArray[27]+rad[ip+2*ng]*(21*angArray[63]+6*angArray[74])-rad[ip+3*ng]*(7*angArray[138]+3*angArray[127])+rad[ip+4*ng]*angArray[238];
         bas[8] = rad[ip]*12*angArray[5]-rad[ip+ng]*(14*angArray[21]+30*angArray[28])+rad[ip+2*ng]*(35*angArray[64]+2*angArray[57]+6*angArray[75])-rad[ip+3*ng]*(7*angArray[139]+5*angArray[128])+rad[ip+4*ng]*angArray[239];
         bas[9] = rad[ip]*36*angArray[6]-rad[ip+ng]*(42*angArray[22]+42*angArray[29])+rad[ip+2*ng]*(49*angArray[65]+6*angArray[76]+6*angArray[58])-rad[ip+3*ng]*(7*angArray[140]+7*angArray[129])+rad[ip+4*ng]*angArray[240];
         bas[10] = -rad[ip+ng]*2*angArray[30]+rad[ip+2*ng]*(5*angArray[66]+2*angArray[79])-rad[ip+3*ng]*(5*angArray[143]+angArray[130])+rad[ip+4*ng]*angArray[243];
         bas[11] = -rad[ip+ng]*6*angArray[31]+rad[ip+2*ng]*(15*angArray[67]+2*angArray[80])-rad[ip+3*ng]*(5*angArray[144]+3*angArray[131])+rad[ip+4*ng]*angArray[244];
         bas[12] = rad[ip]*4*angArray[7]-rad[ip+ng]*(10*angArray[23]+10*angArray[32])+rad[ip+2*ng]*(25*angArray[68]+2*angArray[59]+2*angArray[81])-rad[ip+3*ng]*(5*angArray[145]+5*angArray[132])+rad[ip+4*ng]*angArray[245];
         bas[13] = rad[ip]*12*angArray[8]-rad[ip+ng]*(30*angArray[24]+14*angArray[33])+rad[ip+2*ng]*(35*angArray[69]+6*angArray[60]+2*angArray[82])-rad[ip+3*ng]*(5*angArray[146]+7*angArray[133])+rad[ip+4*ng]*angArray[246];
         bas[14] = rad[ip]*24*angArray[9]-rad[ip+ng]*(60*angArray[25]+18*angArray[34])+rad[ip+2*ng]*(45*angArray[70]+2*angArray[83]+12*angArray[61])-rad[ip+3*ng]*(5*angArray[147]+9*angArray[134])+rad[ip+4*ng]*angArray[247];
         bas[15] = rad[ip+2*ng]*3*angArray[71]-rad[ip+3*ng]*(angArray[135]+3*angArray[150])+rad[ip+4*ng]*angArray[250];
         bas[16] = rad[ip+2*ng]*9*angArray[72]-rad[ip+3*ng]*(3*angArray[136]+3*angArray[151])+rad[ip+4*ng]*angArray[251];
         bas[17] = -rad[ip+ng]*6*angArray[26]+rad[ip+2*ng]*(2*angArray[62]+15*angArray[73])-rad[ip+3*ng]*(5*angArray[137]+3*angArray[152])+rad[ip+4*ng]*angArray[252];
         bas[18] = -rad[ip+ng]*18*angArray[27]+rad[ip+2*ng]*(21*angArray[74]+6*angArray[63])-rad[ip+3*ng]*(7*angArray[138]+3*angArray[153])+rad[ip+4*ng]*angArray[253];
         bas[19] = -rad[ip+ng]*36*angArray[28]+rad[ip+2*ng]*(27*angArray[75]+12*angArray[64])-rad[ip+3*ng]*(9*angArray[139]+3*angArray[154])+rad[ip+4*ng]*angArray[254];
         bas[20] = -rad[ip+ng]*60*angArray[29]+rad[ip+2*ng]*(20*angArray[65]+33*angArray[76])-rad[ip+3*ng]*(11*angArray[140]+3*angArray[155])+rad[ip+4*ng]*angArray[255];
         bas[21] = rad[ip+2*ng]*angArray[77]-rad[ip+3*ng]*(angArray[141]+angArray[158])+rad[ip+4*ng]*angArray[258];
         bas[22] = rad[ip+2*ng]*3*angArray[78]-rad[ip+3*ng]*(3*angArray[142]+angArray[159])+rad[ip+4*ng]*angArray[259];
         bas[23] = -rad[ip+ng]*2*angArray[30]+rad[ip+2*ng]*(5*angArray[79]+2*angArray[66])-rad[ip+3*ng]*(5*angArray[143]+angArray[160])+rad[ip+4*ng]*angArray[260];
         bas[24] = -rad[ip+ng]*6*angArray[31]+rad[ip+2*ng]*(7*angArray[80]+6*angArray[67])-rad[ip+3*ng]*(7*angArray[144]+angArray[161])+rad[ip+4*ng]*angArray[261];
         bas[25] = -rad[ip+ng]*12*angArray[32]+rad[ip+2*ng]*(9*angArray[81]+12*angArray[68])-rad[ip+3*ng]*(9*angArray[145]+angArray[162])+rad[ip+4*ng]*angArray[262];
         bas[26] = -rad[ip+ng]*20*angArray[33]+rad[ip+2*ng]*(20*angArray[69]+11*angArray[82])-rad[ip+3*ng]*(11*angArray[146]+angArray[163])+rad[ip+4*ng]*angArray[263];
         bas[27] = -rad[ip+ng]*30*angArray[34]+rad[ip+2*ng]*(13*angArray[83]+30*angArray[70])-rad[ip+3*ng]*(13*angArray[147]+angArray[164])+rad[ip+4*ng]*angArray[264];
      }

   }


   // now we do derivatives for the given basis set to XYZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*angArray[5]+rad[ip+4*ng]*angArray[28];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[2]-rad[ip+3*ng]*(angArray[11]+angArray[18])+rad[ip+4*ng]*angArray[43];
         bas[1] = rad[ip+2*ng]*angArray[1]-rad[ip+3*ng]*(angArray[13]+angArray[15])+rad[ip+4*ng]*angArray[47];
         bas[2] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[48];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*2*angArray[5]-rad[ip+3*ng]*(angArray[21]+2*angArray[28])+rad[ip+4*ng]*angArray[64];
         bas[1] = -rad[ip+ng]*angArray[0]+rad[ip+2*ng]*(angArray[4]+angArray[7]+angArray[9])-rad[ip+3*ng]*(angArray[23]+angArray[25]+angArray[32])+rad[ip+4*ng]*angArray[68];
         bas[2] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*(3*angArray[24]+angArray[33])+rad[ip+4*ng]*angArray[69];
         bas[3] = rad[ip+2*ng]*2*angArray[5]-rad[ip+3*ng]*(angArray[26]+2*angArray[28])+rad[ip+4*ng]*angArray[73];
         bas[4] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*(3*angArray[27]+angArray[29])+rad[ip+4*ng]*angArray[74];
         bas[5] = rad[ip+2*ng]*2*angArray[5]-rad[ip+3*ng]*5*angArray[28]+rad[ip+4*ng]*angArray[75];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[11]-rad[ip+3*ng]*(angArray[36]+3*angArray[43])+rad[ip+4*ng]*angArray[92];
         bas[1] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(angArray[10]+2*angArray[13]+2*angArray[15])-rad[ip+3*ng]*(2*angArray[47]+angArray[38]+angArray[40])+rad[ip+4*ng]*angArray[96];
         bas[2] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[48]+3*angArray[39])+rad[ip+4*ng]*angArray[97];
         bas[3] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(2*angArray[18]+angArray[16]+2*angArray[11])-rad[ip+3*ng]*(angArray[41]+2*angArray[43]+angArray[52])+rad[ip+4*ng]*angArray[101];
         bas[4] = -rad[ip+ng]*3*angArray[3]+rad[ip+2*ng]*(3*angArray[12]+3*angArray[17]+angArray[19])-rad[ip+3*ng]*(3*angArray[42]+angArray[44]+angArray[53])+rad[ip+4*ng]*angArray[102];
         bas[5] = -rad[ip+ng]*2*angArray[2]+rad[ip+2*ng]*(2*angArray[11]+5*angArray[18])-rad[ip+3*ng]*(5*angArray[43]+angArray[54])+rad[ip+4*ng]*angArray[103];
         bas[6] = rad[ip+2*ng]*3*angArray[13]-rad[ip+3*ng]*(3*angArray[47]+angArray[45])+rad[ip+4*ng]*angArray[107];
         bas[7] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(3*angArray[46]+2*angArray[48])+rad[ip+4*ng]*angArray[108];
         bas[8] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(5*angArray[15]+2*angArray[13])-rad[ip+3*ng]*(5*angArray[47]+angArray[49])+rad[ip+4*ng]*angArray[109];
         bas[9] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*7*angArray[48]+rad[ip+4*ng]*angArray[110];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*4*angArray[21]-rad[ip+3*ng]*(angArray[57]+4*angArray[64])+rad[ip+4*ng]*angArray[128];
         bas[1] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*(angArray[20]+3*angArray[23]+3*angArray[25])-rad[ip+3*ng]*(angArray[59]+3*angArray[68]+angArray[61])+rad[ip+4*ng]*angArray[132];
         bas[2] = rad[ip+2*ng]*9*angArray[24]-rad[ip+3*ng]*(3*angArray[60]+3*angArray[69])+rad[ip+4*ng]*angArray[133];
         bas[3] = -rad[ip+ng]*4*angArray[5]+rad[ip+2*ng]*(4*angArray[28]+2*angArray[26]+2*angArray[21])-rad[ip+3*ng]*(2*angArray[73]+angArray[62]+2*angArray[64])+rad[ip+4*ng]*angArray[137];
         bas[4] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(3*angArray[22]+6*angArray[27]+2*angArray[29])-rad[ip+3*ng]*(2*angArray[74]+3*angArray[63]+angArray[65])+rad[ip+4*ng]*angArray[138];
         bas[5] = -rad[ip+ng]*4*angArray[5]+rad[ip+2*ng]*(10*angArray[28]+2*angArray[21])-rad[ip+3*ng]*(2*angArray[75]+5*angArray[64])+rad[ip+4*ng]*angArray[139];
         bas[6] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*(angArray[30]+3*angArray[32]+3*angArray[23])-rad[ip+3*ng]*(3*angArray[68]+angArray[66]+angArray[79])+rad[ip+4*ng]*angArray[143];
         bas[7] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(2*angArray[33]+6*angArray[24]+3*angArray[31])-rad[ip+3*ng]*(3*angArray[67]+2*angArray[69]+angArray[80])+rad[ip+4*ng]*angArray[144];
         bas[8] = rad[ip]*2*angArray[0]-rad[ip+ng]*(2*angArray[7]+5*angArray[9]+2*angArray[4])+rad[ip+2*ng]*(5*angArray[25]+angArray[34]+5*angArray[32]+2*angArray[23])-rad[ip+3*ng]*(5*angArray[68]+angArray[70]+angArray[81])+rad[ip+4*ng]*angArray[145];
         bas[9] = -rad[ip+ng]*6*angArray[8]+rad[ip+2*ng]*(7*angArray[33]+6*angArray[24])-rad[ip+3*ng]*(7*angArray[69]+angArray[82])+rad[ip+4*ng]*angArray[146];
         bas[10] = rad[ip+2*ng]*4*angArray[26]-rad[ip+3*ng]*(angArray[71]+4*angArray[73])+rad[ip+4*ng]*angArray[150];
         bas[11] = rad[ip+2*ng]*9*angArray[27]-rad[ip+3*ng]*(3*angArray[74]+3*angArray[72])+rad[ip+4*ng]*angArray[151];
         bas[12] = -rad[ip+ng]*4*angArray[5]+rad[ip+2*ng]*(10*angArray[28]+2*angArray[26])-rad[ip+3*ng]*(5*angArray[73]+2*angArray[75])+rad[ip+4*ng]*angArray[152];
         bas[13] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(7*angArray[29]+6*angArray[27])-rad[ip+3*ng]*(7*angArray[74]+angArray[76])+rad[ip+4*ng]*angArray[153];
         bas[14] = rad[ip+2*ng]*12*angArray[28]-rad[ip+3*ng]*9*angArray[75]+rad[ip+4*ng]*angArray[154];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*5*angArray[36]-rad[ip+3*ng]*(angArray[85]+5*angArray[92])+rad[ip+4*ng]*angArray[173];
         bas[1] = -rad[ip+ng]*4*angArray[10]+rad[ip+2*ng]*(angArray[35]+4*angArray[38]+4*angArray[40])-rad[ip+3*ng]*(angArray[87]+angArray[89]+4*angArray[96])+rad[ip+4*ng]*angArray[177];
         bas[2] = rad[ip+2*ng]*12*angArray[39]-rad[ip+3*ng]*(3*angArray[88]+4*angArray[97])+rad[ip+4*ng]*angArray[178];
         bas[3] = -rad[ip+ng]*6*angArray[11]+rad[ip+2*ng]*(3*angArray[41]+2*angArray[36]+6*angArray[43])-rad[ip+3*ng]*(angArray[90]+3*angArray[101]+2*angArray[92])+rad[ip+4*ng]*angArray[182];
         bas[4] = -rad[ip+ng]*9*angArray[12]+rad[ip+2*ng]*(9*angArray[42]+3*angArray[37]+3*angArray[44])-rad[ip+3*ng]*(3*angArray[91]+3*angArray[102]+angArray[93])+rad[ip+4*ng]*angArray[183];
         bas[5] = -rad[ip+ng]*6*angArray[11]+rad[ip+2*ng]*(15*angArray[43]+2*angArray[36])-rad[ip+3*ng]*(5*angArray[92]+3*angArray[103])+rad[ip+4*ng]*angArray[184];
         bas[6] = -rad[ip+ng]*6*angArray[13]+rad[ip+2*ng]*(2*angArray[45]+6*angArray[47]+3*angArray[38])-rad[ip+3*ng]*(3*angArray[96]+2*angArray[107]+angArray[94])+rad[ip+4*ng]*angArray[188];
         bas[7] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(4*angArray[48]+6*angArray[46]+6*angArray[39])-rad[ip+3*ng]*(2*angArray[108]+3*angArray[95]+2*angArray[97])+rad[ip+4*ng]*angArray[189];
         bas[8] = rad[ip]*4*angArray[1]-rad[ip+ng]*(4*angArray[13]+2*angArray[10]+10*angArray[15])+rad[ip+2*ng]*(10*angArray[47]+2*angArray[38]+5*angArray[40]+2*angArray[49])-rad[ip+3*ng]*(2*angArray[109]+angArray[98]+5*angArray[96])+rad[ip+4*ng]*angArray[190];
         bas[9] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(14*angArray[48]+6*angArray[39])-rad[ip+3*ng]*(2*angArray[110]+7*angArray[97])+rad[ip+4*ng]*angArray[191];
         bas[10] = -rad[ip+ng]*4*angArray[16]+rad[ip+2*ng]*(4*angArray[41]+angArray[50]+4*angArray[52])-rad[ip+3*ng]*(4*angArray[101]+angArray[99]+angArray[114])+rad[ip+4*ng]*angArray[195];
         bas[11] = -rad[ip+ng]*9*angArray[17]+rad[ip+2*ng]*(9*angArray[42]+3*angArray[51]+3*angArray[53])-rad[ip+3*ng]*(3*angArray[102]+3*angArray[100]+angArray[115])+rad[ip+4*ng]*angArray[196];
         bas[12] = rad[ip]*4*angArray[2]-rad[ip+ng]*(2*angArray[16]+10*angArray[18]+4*angArray[11])+rad[ip+2*ng]*(2*angArray[41]+2*angArray[54]+10*angArray[43]+5*angArray[52])-rad[ip+3*ng]*(2*angArray[103]+angArray[116]+5*angArray[101])+rad[ip+4*ng]*angArray[197];
         bas[13] = rad[ip]*6*angArray[3]-rad[ip+ng]*(7*angArray[19]+6*angArray[17]+6*angArray[12])+rad[ip+2*ng]*(7*angArray[44]+angArray[55]+7*angArray[53]+6*angArray[42])-rad[ip+3*ng]*(7*angArray[102]+angArray[104]+angArray[117])+rad[ip+4*ng]*angArray[198];
         bas[14] = -rad[ip+ng]*12*angArray[18]+rad[ip+2*ng]*(9*angArray[54]+12*angArray[43])-rad[ip+3*ng]*(9*angArray[103]+angArray[118])+rad[ip+4*ng]*angArray[199];
         bas[15] = rad[ip+2*ng]*5*angArray[45]-rad[ip+3*ng]*(angArray[105]+5*angArray[107])+rad[ip+4*ng]*angArray[203];
         bas[16] = rad[ip+2*ng]*12*angArray[46]-rad[ip+3*ng]*(4*angArray[108]+3*angArray[106])+rad[ip+4*ng]*angArray[204];
         bas[17] = -rad[ip+ng]*6*angArray[13]+rad[ip+2*ng]*(2*angArray[45]+15*angArray[47])-rad[ip+3*ng]*(3*angArray[109]+5*angArray[107])+rad[ip+4*ng]*angArray[205];
         bas[18] = -rad[ip+ng]*12*angArray[14]+rad[ip+2*ng]*(14*angArray[48]+6*angArray[46])-rad[ip+3*ng]*(7*angArray[108]+2*angArray[110])+rad[ip+4*ng]*angArray[206];
         bas[19] = -rad[ip+ng]*12*angArray[15]+rad[ip+2*ng]*(9*angArray[49]+12*angArray[47])-rad[ip+3*ng]*(9*angArray[109]+angArray[111])+rad[ip+4*ng]*angArray[207];
         bas[20] = rad[ip+2*ng]*20*angArray[48]-rad[ip+3*ng]*11*angArray[110]+rad[ip+4*ng]*angArray[208];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*6*angArray[57]-rad[ip+3*ng]*(angArray[121]+6*angArray[128])+rad[ip+4*ng]*angArray[228];
         bas[1] = -rad[ip+ng]*5*angArray[20]+rad[ip+2*ng]*(5*angArray[61]+angArray[56]+5*angArray[59])-rad[ip+3*ng]*(angArray[123]+angArray[125]+5*angArray[132])+rad[ip+4*ng]*angArray[232];
         bas[2] = rad[ip+2*ng]*15*angArray[60]-rad[ip+3*ng]*(3*angArray[124]+5*angArray[133])+rad[ip+4*ng]*angArray[233];
         bas[3] = -rad[ip+ng]*8*angArray[21]+rad[ip+2*ng]*(4*angArray[62]+2*angArray[57]+8*angArray[64])-rad[ip+3*ng]*(2*angArray[128]+angArray[126]+4*angArray[137])+rad[ip+4*ng]*angArray[237];
         bas[4] = -rad[ip+ng]*12*angArray[22]+rad[ip+2*ng]*(3*angArray[58]+12*angArray[63]+4*angArray[65])-rad[ip+3*ng]*(3*angArray[127]+angArray[129]+4*angArray[138])+rad[ip+4*ng]*angArray[238];
         bas[5] = -rad[ip+ng]*8*angArray[21]+rad[ip+2*ng]*(20*angArray[64]+2*angArray[57])-rad[ip+3*ng]*(5*angArray[128]+4*angArray[139])+rad[ip+4*ng]*angArray[239];
         bas[6] = -rad[ip+ng]*9*angArray[23]+rad[ip+2*ng]*(3*angArray[66]+9*angArray[68]+3*angArray[59])-rad[ip+3*ng]*(3*angArray[132]+angArray[130]+3*angArray[143])+rad[ip+4*ng]*angArray[243];
         bas[7] = -rad[ip+ng]*18*angArray[24]+rad[ip+2*ng]*(9*angArray[67]+6*angArray[69]+6*angArray[60])-rad[ip+3*ng]*(3*angArray[131]+3*angArray[144]+2*angArray[133])+rad[ip+4*ng]*angArray[244];
         bas[8] = rad[ip]*6*angArray[4]-rad[ip+ng]*(6*angArray[23]+15*angArray[25]+2*angArray[20])+rad[ip+2*ng]*(15*angArray[68]+5*angArray[61]+2*angArray[59]+3*angArray[70])-rad[ip+3*ng]*(5*angArray[132]+3*angArray[145]+angArray[134])+rad[ip+4*ng]*angArray[245];
         bas[9] = -rad[ip+ng]*18*angArray[24]+rad[ip+2*ng]*(21*angArray[69]+6*angArray[60])-rad[ip+3*ng]*(3*angArray[146]+7*angArray[133])+rad[ip+4*ng]*angArray[246];
         bas[10] = -rad[ip+ng]*8*angArray[26]+rad[ip+2*ng]*(4*angArray[62]+8*angArray[73]+2*angArray[71])-rad[ip+3*ng]*(2*angArray[150]+angArray[135]+4*angArray[137])+rad[ip+4*ng]*angArray[250];
         bas[11] = -rad[ip+ng]*18*angArray[27]+rad[ip+2*ng]*(9*angArray[63]+6*angArray[72]+6*angArray[74])-rad[ip+3*ng]*(3*angArray[138]+2*angArray[151]+3*angArray[136])+rad[ip+4*ng]*angArray[251];
         bas[12] = rad[ip]*8*angArray[5]-rad[ip+ng]*(4*angArray[26]+20*angArray[28]+4*angArray[21])+rad[ip+2*ng]*(4*angArray[75]+10*angArray[73]+10*angArray[64]+2*angArray[62])-rad[ip+3*ng]*(2*angArray[152]+5*angArray[137]+2*angArray[139])+rad[ip+4*ng]*angArray[252];
         bas[13] = rad[ip]*12*angArray[6]-rad[ip+ng]*(12*angArray[27]+14*angArray[29]+6*angArray[22])+rad[ip+2*ng]*(14*angArray[74]+7*angArray[65]+6*angArray[63]+2*angArray[76])-rad[ip+3*ng]*(2*angArray[153]+7*angArray[138]+angArray[140])+rad[ip+4*ng]*angArray[253];
         bas[14] = -rad[ip+ng]*24*angArray[28]+rad[ip+2*ng]*(18*angArray[75]+12*angArray[64])-rad[ip+3*ng]*(2*angArray[154]+9*angArray[139])+rad[ip+4*ng]*angArray[254];
         bas[15] = -rad[ip+ng]*5*angArray[30]+rad[ip+2*ng]*(5*angArray[79]+angArray[77]+5*angArray[66])-rad[ip+3*ng]*(angArray[141]+5*angArray[143]+angArray[158])+rad[ip+4*ng]*angArray[258];
         bas[16] = -rad[ip+ng]*12*angArray[31]+rad[ip+2*ng]*(12*angArray[67]+3*angArray[78]+4*angArray[80])-rad[ip+3*ng]*(4*angArray[144]+3*angArray[142]+angArray[159])+rad[ip+4*ng]*angArray[259];
         bas[17] = rad[ip]*6*angArray[7]-rad[ip+ng]*(2*angArray[30]+6*angArray[23]+15*angArray[32])+rad[ip+2*ng]*(5*angArray[79]+2*angArray[66]+3*angArray[81]+15*angArray[68])-rad[ip+3*ng]*(3*angArray[145]+5*angArray[143]+angArray[160])+rad[ip+4*ng]*angArray[260];
         bas[18] = rad[ip]*12*angArray[8]-rad[ip+ng]*(14*angArray[33]+6*angArray[31]+12*angArray[24])+rad[ip+2*ng]*(2*angArray[82]+14*angArray[69]+6*angArray[67]+7*angArray[80])-rad[ip+3*ng]*(angArray[161]+7*angArray[144]+2*angArray[146])+rad[ip+4*ng]*angArray[261];
         bas[19] = rad[ip]*12*angArray[9]-rad[ip+ng]*(12*angArray[25]+9*angArray[34]+12*angArray[32])+rad[ip+2*ng]*(9*angArray[70]+9*angArray[81]+12*angArray[68]+angArray[83])-rad[ip+3*ng]*(9*angArray[145]+angArray[147]+angArray[162])+rad[ip+4*ng]*angArray[262];
         bas[20] = -rad[ip+ng]*20*angArray[33]+rad[ip+2*ng]*(20*angArray[69]+11*angArray[82])-rad[ip+3*ng]*(11*angArray[146]+angArray[163])+rad[ip+4*ng]*angArray[263];
         bas[21] = rad[ip+2*ng]*6*angArray[71]-rad[ip+3*ng]*(6*angArray[150]+angArray[148])+rad[ip+4*ng]*angArray[267];
         bas[22] = rad[ip+2*ng]*15*angArray[72]-rad[ip+3*ng]*(3*angArray[149]+5*angArray[151])+rad[ip+4*ng]*angArray[268];
         bas[23] = -rad[ip+ng]*8*angArray[26]+rad[ip+2*ng]*(20*angArray[73]+2*angArray[71])-rad[ip+3*ng]*(4*angArray[152]+5*angArray[150])+rad[ip+4*ng]*angArray[269];
         bas[24] = -rad[ip+ng]*18*angArray[27]+rad[ip+2*ng]*(21*angArray[74]+6*angArray[72])-rad[ip+3*ng]*(3*angArray[153]+7*angArray[151])+rad[ip+4*ng]*angArray[270];
         bas[25] = -rad[ip+ng]*24*angArray[28]+rad[ip+2*ng]*(18*angArray[75]+12*angArray[73])-rad[ip+3*ng]*(9*angArray[152]+2*angArray[154])+rad[ip+4*ng]*angArray[271];
         bas[26] = -rad[ip+ng]*20*angArray[29]+rad[ip+2*ng]*(20*angArray[74]+11*angArray[76])-rad[ip+3*ng]*(11*angArray[153]+angArray[155])+rad[ip+4*ng]*angArray[272];
         bas[27] = rad[ip+2*ng]*30*angArray[75]-rad[ip+3*ng]*13*angArray[154]+rad[ip+4*ng]*angArray[273];
      }

   }


   // now we do derivatives for the given basis set to YYZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[0]-rad[ip+3*ng]*(angArray[7]+angArray[9])+rad[ip+4*ng]*angArray[32];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[1]-rad[ip+3*ng]*(angArray[13]+angArray[15])+rad[ip+4*ng]*angArray[47];
         bas[1] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*(angArray[16]+3*angArray[18])+rad[ip+4*ng]*angArray[52];
         bas[2] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*(3*angArray[17]+angArray[19])+rad[ip+4*ng]*angArray[53];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[4]-rad[ip+3*ng]*(angArray[23]+angArray[25])+rad[ip+4*ng]*angArray[68];
         bas[1] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*(angArray[26]+3*angArray[28])+rad[ip+4*ng]*angArray[73];
         bas[2] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*(3*angArray[27]+angArray[29])+rad[ip+4*ng]*angArray[74];
         bas[3] = -rad[ip+ng]*2*angArray[0]+rad[ip+2*ng]*(5*angArray[7]+2*angArray[9])-rad[ip+3*ng]*(5*angArray[32]+angArray[30])+rad[ip+4*ng]*angArray[79];
         bas[4] = rad[ip+2*ng]*9*angArray[8]-rad[ip+3*ng]*(3*angArray[31]+3*angArray[33])+rad[ip+4*ng]*angArray[80];
         bas[5] = -rad[ip+ng]*2*angArray[0]+rad[ip+2*ng]*(5*angArray[9]+2*angArray[7])-rad[ip+3*ng]*(5*angArray[32]+angArray[34])+rad[ip+4*ng]*angArray[81];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[10]-rad[ip+3*ng]*(angArray[38]+angArray[40])+rad[ip+4*ng]*angArray[96];
         bas[1] = rad[ip+2*ng]*3*angArray[11]-rad[ip+3*ng]*(angArray[41]+3*angArray[43])+rad[ip+4*ng]*angArray[101];
         bas[2] = rad[ip+2*ng]*3*angArray[12]-rad[ip+3*ng]*(3*angArray[42]+angArray[44])+rad[ip+4*ng]*angArray[102];
         bas[3] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(5*angArray[13]+2*angArray[15])-rad[ip+3*ng]*(5*angArray[47]+angArray[45])+rad[ip+4*ng]*angArray[107];
         bas[4] = rad[ip+2*ng]*9*angArray[14]-rad[ip+3*ng]*(3*angArray[46]+3*angArray[48])+rad[ip+4*ng]*angArray[108];
         bas[5] = -rad[ip+ng]*2*angArray[1]+rad[ip+2*ng]*(5*angArray[15]+2*angArray[13])-rad[ip+3*ng]*(5*angArray[47]+angArray[49])+rad[ip+4*ng]*angArray[109];
         bas[6] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*(7*angArray[16]+6*angArray[18])-rad[ip+3*ng]*(7*angArray[52]+angArray[50])+rad[ip+4*ng]*angArray[114];
         bas[7] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*(15*angArray[17]+2*angArray[19])-rad[ip+3*ng]*(5*angArray[53]+3*angArray[51])+rad[ip+4*ng]*angArray[115];
         bas[8] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*(2*angArray[16]+15*angArray[18])-rad[ip+3*ng]*(5*angArray[52]+3*angArray[54])+rad[ip+4*ng]*angArray[116];
         bas[9] = -rad[ip+ng]*6*angArray[3]+rad[ip+2*ng]*(7*angArray[19]+6*angArray[17])-rad[ip+3*ng]*(7*angArray[53]+angArray[55])+rad[ip+4*ng]*angArray[117];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[20]-rad[ip+3*ng]*(angArray[59]+angArray[61])+rad[ip+4*ng]*angArray[132];
         bas[1] = rad[ip+2*ng]*3*angArray[21]-rad[ip+3*ng]*(angArray[62]+3*angArray[64])+rad[ip+4*ng]*angArray[137];
         bas[2] = rad[ip+2*ng]*3*angArray[22]-rad[ip+3*ng]*(3*angArray[63]+angArray[65])+rad[ip+4*ng]*angArray[138];
         bas[3] = -rad[ip+ng]*2*angArray[4]+rad[ip+2*ng]*(5*angArray[23]+2*angArray[25])-rad[ip+3*ng]*(5*angArray[68]+angArray[66])+rad[ip+4*ng]*angArray[143];
         bas[4] = rad[ip+2*ng]*9*angArray[24]-rad[ip+3*ng]*(3*angArray[67]+3*angArray[69])+rad[ip+4*ng]*angArray[144];
         bas[5] = -rad[ip+ng]*2*angArray[4]+rad[ip+2*ng]*(5*angArray[25]+2*angArray[23])-rad[ip+3*ng]*(5*angArray[68]+angArray[70])+rad[ip+4*ng]*angArray[145];
         bas[6] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(7*angArray[26]+6*angArray[28])-rad[ip+3*ng]*(7*angArray[73]+angArray[71])+rad[ip+4*ng]*angArray[150];
         bas[7] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(15*angArray[27]+2*angArray[29])-rad[ip+3*ng]*(5*angArray[74]+3*angArray[72])+rad[ip+4*ng]*angArray[151];
         bas[8] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(2*angArray[26]+15*angArray[28])-rad[ip+3*ng]*(5*angArray[73]+3*angArray[75])+rad[ip+4*ng]*angArray[152];
         bas[9] = -rad[ip+ng]*6*angArray[6]+rad[ip+2*ng]*(7*angArray[29]+6*angArray[27])-rad[ip+3*ng]*(7*angArray[74]+angArray[76])+rad[ip+4*ng]*angArray[153];
         bas[10] = -rad[ip+ng]*12*angArray[7]+rad[ip+2*ng]*(9*angArray[30]+12*angArray[32])-rad[ip+3*ng]*(angArray[77]+9*angArray[79])+rad[ip+4*ng]*angArray[158];
         bas[11] = -rad[ip+ng]*18*angArray[8]+rad[ip+2*ng]*(21*angArray[31]+6*angArray[33])-rad[ip+3*ng]*(7*angArray[80]+3*angArray[78])+rad[ip+4*ng]*angArray[159];
         bas[12] = rad[ip]*4*angArray[0]-rad[ip+ng]*(10*angArray[7]+10*angArray[9])+rad[ip+2*ng]*(25*angArray[32]+2*angArray[30]+2*angArray[34])-rad[ip+3*ng]*(5*angArray[81]+5*angArray[79])+rad[ip+4*ng]*angArray[160];
         bas[13] = -rad[ip+ng]*18*angArray[8]+rad[ip+2*ng]*(21*angArray[33]+6*angArray[31])-rad[ip+3*ng]*(7*angArray[80]+3*angArray[82])+rad[ip+4*ng]*angArray[161];
         bas[14] = -rad[ip+ng]*12*angArray[9]+rad[ip+2*ng]*(9*angArray[34]+12*angArray[32])-rad[ip+3*ng]*(9*angArray[81]+angArray[83])+rad[ip+4*ng]*angArray[162];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[35]-rad[ip+3*ng]*(angArray[87]+angArray[89])+rad[ip+4*ng]*angArray[177];
         bas[1] = rad[ip+2*ng]*3*angArray[36]-rad[ip+3*ng]*(angArray[90]+3*angArray[92])+rad[ip+4*ng]*angArray[182];
         bas[2] = rad[ip+2*ng]*3*angArray[37]-rad[ip+3*ng]*(3*angArray[91]+angArray[93])+rad[ip+4*ng]*angArray[183];
         bas[3] = -rad[ip+ng]*2*angArray[10]+rad[ip+2*ng]*(5*angArray[38]+2*angArray[40])-rad[ip+3*ng]*(5*angArray[96]+angArray[94])+rad[ip+4*ng]*angArray[188];
         bas[4] = rad[ip+2*ng]*9*angArray[39]-rad[ip+3*ng]*(3*angArray[95]+3*angArray[97])+rad[ip+4*ng]*angArray[189];
         bas[5] = -rad[ip+ng]*2*angArray[10]+rad[ip+2*ng]*(5*angArray[40]+2*angArray[38])-rad[ip+3*ng]*(5*angArray[96]+angArray[98])+rad[ip+4*ng]*angArray[190];
         bas[6] = -rad[ip+ng]*6*angArray[11]+rad[ip+2*ng]*(7*angArray[41]+6*angArray[43])-rad[ip+3*ng]*(7*angArray[101]+angArray[99])+rad[ip+4*ng]*angArray[195];
         bas[7] = -rad[ip+ng]*6*angArray[12]+rad[ip+2*ng]*(15*angArray[42]+2*angArray[44])-rad[ip+3*ng]*(5*angArray[102]+3*angArray[100])+rad[ip+4*ng]*angArray[196];
         bas[8] = -rad[ip+ng]*6*angArray[11]+rad[ip+2*ng]*(2*angArray[41]+15*angArray[43])-rad[ip+3*ng]*(5*angArray[101]+3*angArray[103])+rad[ip+4*ng]*angArray[197];
         bas[9] = -rad[ip+ng]*6*angArray[12]+rad[ip+2*ng]*(7*angArray[44]+6*angArray[42])-rad[ip+3*ng]*(7*angArray[102]+angArray[104])+rad[ip+4*ng]*angArray[198];
         bas[10] = -rad[ip+ng]*12*angArray[13]+rad[ip+2*ng]*(9*angArray[45]+12*angArray[47])-rad[ip+3*ng]*(angArray[105]+9*angArray[107])+rad[ip+4*ng]*angArray[203];
         bas[11] = -rad[ip+ng]*18*angArray[14]+rad[ip+2*ng]*(21*angArray[46]+6*angArray[48])-rad[ip+3*ng]*(7*angArray[108]+3*angArray[106])+rad[ip+4*ng]*angArray[204];
         bas[12] = rad[ip]*4*angArray[1]-rad[ip+ng]*(10*angArray[13]+10*angArray[15])+rad[ip+2*ng]*(25*angArray[47]+2*angArray[45]+2*angArray[49])-rad[ip+3*ng]*(5*angArray[109]+5*angArray[107])+rad[ip+4*ng]*angArray[205];
         bas[13] = -rad[ip+ng]*18*angArray[14]+rad[ip+2*ng]*(21*angArray[48]+6*angArray[46])-rad[ip+3*ng]*(7*angArray[108]+3*angArray[110])+rad[ip+4*ng]*angArray[206];
         bas[14] = -rad[ip+ng]*12*angArray[15]+rad[ip+2*ng]*(9*angArray[49]+12*angArray[47])-rad[ip+3*ng]*(9*angArray[109]+angArray[111])+rad[ip+4*ng]*angArray[207];
         bas[15] = -rad[ip+ng]*20*angArray[16]+rad[ip+2*ng]*(20*angArray[52]+11*angArray[50])-rad[ip+3*ng]*(11*angArray[114]+angArray[112])+rad[ip+4*ng]*angArray[212];
         bas[16] = -rad[ip+ng]*36*angArray[17]+rad[ip+2*ng]*(27*angArray[51]+12*angArray[53])-rad[ip+3*ng]*(3*angArray[113]+9*angArray[115])+rad[ip+4*ng]*angArray[213];
         bas[17] = rad[ip]*12*angArray[2]-rad[ip+ng]*(14*angArray[16]+30*angArray[18])+rad[ip+2*ng]*(35*angArray[52]+2*angArray[50]+6*angArray[54])-rad[ip+3*ng]*(7*angArray[116]+5*angArray[114])+rad[ip+4*ng]*angArray[214];
         bas[18] = rad[ip]*12*angArray[3]-rad[ip+ng]*(30*angArray[17]+14*angArray[19])+rad[ip+2*ng]*(35*angArray[53]+6*angArray[51]+2*angArray[55])-rad[ip+3*ng]*(5*angArray[117]+7*angArray[115])+rad[ip+4*ng]*angArray[215];
         bas[19] = -rad[ip+ng]*36*angArray[18]+rad[ip+2*ng]*(27*angArray[54]+12*angArray[52])-rad[ip+3*ng]*(9*angArray[116]+3*angArray[118])+rad[ip+4*ng]*angArray[216];
         bas[20] = -rad[ip+ng]*20*angArray[19]+rad[ip+2*ng]*(20*angArray[53]+11*angArray[55])-rad[ip+3*ng]*(11*angArray[117]+angArray[119])+rad[ip+4*ng]*angArray[217];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*angArray[56]-rad[ip+3*ng]*(angArray[123]+angArray[125])+rad[ip+4*ng]*angArray[232];
         bas[1] = rad[ip+2*ng]*3*angArray[57]-rad[ip+3*ng]*(angArray[126]+3*angArray[128])+rad[ip+4*ng]*angArray[237];
         bas[2] = rad[ip+2*ng]*3*angArray[58]-rad[ip+3*ng]*(3*angArray[127]+angArray[129])+rad[ip+4*ng]*angArray[238];
         bas[3] = -rad[ip+ng]*2*angArray[20]+rad[ip+2*ng]*(5*angArray[59]+2*angArray[61])-rad[ip+3*ng]*(5*angArray[132]+angArray[130])+rad[ip+4*ng]*angArray[243];
         bas[4] = rad[ip+2*ng]*9*angArray[60]-rad[ip+3*ng]*(3*angArray[131]+3*angArray[133])+rad[ip+4*ng]*angArray[244];
         bas[5] = -rad[ip+ng]*2*angArray[20]+rad[ip+2*ng]*(5*angArray[61]+2*angArray[59])-rad[ip+3*ng]*(5*angArray[132]+angArray[134])+rad[ip+4*ng]*angArray[245];
         bas[6] = -rad[ip+ng]*6*angArray[21]+rad[ip+2*ng]*(7*angArray[62]+6*angArray[64])-rad[ip+3*ng]*(7*angArray[137]+angArray[135])+rad[ip+4*ng]*angArray[250];
         bas[7] = -rad[ip+ng]*6*angArray[22]+rad[ip+2*ng]*(15*angArray[63]+2*angArray[65])-rad[ip+3*ng]*(5*angArray[138]+3*angArray[136])+rad[ip+4*ng]*angArray[251];
         bas[8] = -rad[ip+ng]*6*angArray[21]+rad[ip+2*ng]*(2*angArray[62]+15*angArray[64])-rad[ip+3*ng]*(5*angArray[137]+3*angArray[139])+rad[ip+4*ng]*angArray[252];
         bas[9] = -rad[ip+ng]*6*angArray[22]+rad[ip+2*ng]*(7*angArray[65]+6*angArray[63])-rad[ip+3*ng]*(7*angArray[138]+angArray[140])+rad[ip+4*ng]*angArray[253];
         bas[10] = -rad[ip+ng]*12*angArray[23]+rad[ip+2*ng]*(9*angArray[66]+12*angArray[68])-rad[ip+3*ng]*(angArray[141]+9*angArray[143])+rad[ip+4*ng]*angArray[258];
         bas[11] = -rad[ip+ng]*18*angArray[24]+rad[ip+2*ng]*(21*angArray[67]+6*angArray[69])-rad[ip+3*ng]*(7*angArray[144]+3*angArray[142])+rad[ip+4*ng]*angArray[259];
         bas[12] = rad[ip]*4*angArray[4]-rad[ip+ng]*(10*angArray[23]+10*angArray[25])+rad[ip+2*ng]*(25*angArray[68]+2*angArray[66]+2*angArray[70])-rad[ip+3*ng]*(5*angArray[145]+5*angArray[143])+rad[ip+4*ng]*angArray[260];
         bas[13] = -rad[ip+ng]*18*angArray[24]+rad[ip+2*ng]*(21*angArray[69]+6*angArray[67])-rad[ip+3*ng]*(7*angArray[144]+3*angArray[146])+rad[ip+4*ng]*angArray[261];
         bas[14] = -rad[ip+ng]*12*angArray[25]+rad[ip+2*ng]*(9*angArray[70]+12*angArray[68])-rad[ip+3*ng]*(9*angArray[145]+angArray[147])+rad[ip+4*ng]*angArray[262];
         bas[15] = -rad[ip+ng]*20*angArray[26]+rad[ip+2*ng]*(20*angArray[73]+11*angArray[71])-rad[ip+3*ng]*(11*angArray[150]+angArray[148])+rad[ip+4*ng]*angArray[267];
         bas[16] = -rad[ip+ng]*36*angArray[27]+rad[ip+2*ng]*(27*angArray[72]+12*angArray[74])-rad[ip+3*ng]*(3*angArray[149]+9*angArray[151])+rad[ip+4*ng]*angArray[268];
         bas[17] = rad[ip]*12*angArray[5]-rad[ip+ng]*(14*angArray[26]+30*angArray[28])+rad[ip+2*ng]*(35*angArray[73]+2*angArray[71]+6*angArray[75])-rad[ip+3*ng]*(7*angArray[152]+5*angArray[150])+rad[ip+4*ng]*angArray[269];
         bas[18] = rad[ip]*12*angArray[6]-rad[ip+ng]*(30*angArray[27]+14*angArray[29])+rad[ip+2*ng]*(35*angArray[74]+6*angArray[72]+2*angArray[76])-rad[ip+3*ng]*(5*angArray[153]+7*angArray[151])+rad[ip+4*ng]*angArray[270];
         bas[19] = -rad[ip+ng]*36*angArray[28]+rad[ip+2*ng]*(27*angArray[75]+12*angArray[73])-rad[ip+3*ng]*(9*angArray[152]+3*angArray[154])+rad[ip+4*ng]*angArray[271];
         bas[20] = -rad[ip+ng]*20*angArray[29]+rad[ip+2*ng]*(20*angArray[74]+11*angArray[76])-rad[ip+3*ng]*(11*angArray[153]+angArray[155])+rad[ip+4*ng]*angArray[272];
         bas[21] = -rad[ip+ng]*30*angArray[30]+rad[ip+2*ng]*(30*angArray[79]+13*angArray[77])-rad[ip+3*ng]*(angArray[156]+13*angArray[158])+rad[ip+4*ng]*angArray[277];
         bas[22] = -rad[ip+ng]*60*angArray[31]+rad[ip+2*ng]*(20*angArray[80]+33*angArray[78])-rad[ip+3*ng]*(3*angArray[157]+11*angArray[159])+rad[ip+4*ng]*angArray[278];
         bas[23] = rad[ip]*24*angArray[7]-rad[ip+ng]*(60*angArray[32]+18*angArray[30])+rad[ip+2*ng]*(45*angArray[79]+12*angArray[81]+2*angArray[77])-rad[ip+3*ng]*(5*angArray[158]+9*angArray[160])+rad[ip+4*ng]*angArray[279];
         bas[24] = rad[ip]*36*angArray[8]-rad[ip+ng]*(42*angArray[31]+42*angArray[33])+rad[ip+2*ng]*(49*angArray[80]+6*angArray[82]+6*angArray[78])-rad[ip+3*ng]*(7*angArray[161]+7*angArray[159])+rad[ip+4*ng]*angArray[280];
         bas[25] = rad[ip]*24*angArray[9]-rad[ip+ng]*(60*angArray[32]+18*angArray[34])+rad[ip+2*ng]*(45*angArray[81]+2*angArray[83]+12*angArray[79])-rad[ip+3*ng]*(5*angArray[162]+9*angArray[160])+rad[ip+4*ng]*angArray[281];
         bas[26] = -rad[ip+ng]*60*angArray[33]+rad[ip+2*ng]*(20*angArray[80]+33*angArray[82])-rad[ip+3*ng]*(11*angArray[161]+3*angArray[163])+rad[ip+4*ng]*angArray[282];
         bas[27] = -rad[ip+ng]*30*angArray[34]+rad[ip+2*ng]*(13*angArray[83]+30*angArray[81])-rad[ip+3*ng]*(13*angArray[162]+angArray[164])+rad[ip+4*ng]*angArray[283];
      }

   }


   // now we do derivatives for the given basis set to XZZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[6]+rad[ip+4*ng]*angArray[29];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*(3*angArray[12]+angArray[19])+rad[ip+4*ng]*angArray[44];
         bas[1] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[48];
         bas[2] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*6*angArray[15]+rad[ip+4*ng]*angArray[49];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*6*angArray[6]-rad[ip+3*ng]*(2*angArray[29]+3*angArray[22])+rad[ip+4*ng]*angArray[65];
         bas[1] = rad[ip+2*ng]*3*angArray[8]-rad[ip+3*ng]*(3*angArray[24]+angArray[33])+rad[ip+4*ng]*angArray[69];
         bas[2] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*(3*angArray[4]+6*angArray[9])-rad[ip+3*ng]*(6*angArray[25]+angArray[34])+rad[ip+4*ng]*angArray[70];
         bas[3] = -rad[ip+3*ng]*3*angArray[27]+rad[ip+4*ng]*angArray[74];
         bas[4] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*6*angArray[28]+rad[ip+4*ng]*angArray[75];
         bas[5] = rad[ip+2*ng]*12*angArray[6]-rad[ip+3*ng]*9*angArray[29]+rad[ip+4*ng]*angArray[76];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*9*angArray[12]-rad[ip+3*ng]*(3*angArray[37]+3*angArray[44])+rad[ip+4*ng]*angArray[93];
         bas[1] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[48]+3*angArray[39])+rad[ip+4*ng]*angArray[97];
         bas[2] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*(3*angArray[10]+12*angArray[15])-rad[ip+3*ng]*(2*angArray[49]+6*angArray[40])+rad[ip+4*ng]*angArray[98];
         bas[3] = rad[ip+2*ng]*3*angArray[17]-rad[ip+3*ng]*(3*angArray[42]+angArray[53])+rad[ip+4*ng]*angArray[102];
         bas[4] = -rad[ip+ng]*3*angArray[2]+rad[ip+2*ng]*(3*angArray[11]+6*angArray[18])-rad[ip+3*ng]*(6*angArray[43]+angArray[54])+rad[ip+4*ng]*angArray[103];
         bas[5] = -rad[ip+ng]*12*angArray[3]+rad[ip+2*ng]*(9*angArray[19]+12*angArray[12])-rad[ip+3*ng]*(9*angArray[44]+angArray[55])+rad[ip+4*ng]*angArray[104];
         bas[6] = -rad[ip+3*ng]*3*angArray[46]+rad[ip+4*ng]*angArray[108];
         bas[7] = rad[ip+2*ng]*3*angArray[13]-rad[ip+3*ng]*6*angArray[47]+rad[ip+4*ng]*angArray[109];
         bas[8] = rad[ip+2*ng]*12*angArray[14]-rad[ip+3*ng]*9*angArray[48]+rad[ip+4*ng]*angArray[110];
         bas[9] = -rad[ip+ng]*6*angArray[1]+rad[ip+2*ng]*27*angArray[15]-rad[ip+3*ng]*12*angArray[49]+rad[ip+4*ng]*angArray[111];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*12*angArray[22]-rad[ip+3*ng]*(3*angArray[58]+4*angArray[65])+rad[ip+4*ng]*angArray[129];
         bas[1] = rad[ip+2*ng]*9*angArray[24]-rad[ip+3*ng]*(3*angArray[60]+3*angArray[69])+rad[ip+4*ng]*angArray[133];
         bas[2] = -rad[ip+ng]*9*angArray[4]+rad[ip+2*ng]*(18*angArray[25]+3*angArray[20])-rad[ip+3*ng]*(6*angArray[61]+3*angArray[70])+rad[ip+4*ng]*angArray[134];
         bas[3] = rad[ip+2*ng]*6*angArray[27]-rad[ip+3*ng]*(2*angArray[74]+3*angArray[63])+rad[ip+4*ng]*angArray[138];
         bas[4] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(3*angArray[21]+12*angArray[28])-rad[ip+3*ng]*(2*angArray[75]+6*angArray[64])+rad[ip+4*ng]*angArray[139];
         bas[5] = -rad[ip+ng]*24*angArray[6]+rad[ip+2*ng]*(18*angArray[29]+12*angArray[22])-rad[ip+3*ng]*(2*angArray[76]+9*angArray[65])+rad[ip+4*ng]*angArray[140];
         bas[6] = rad[ip+2*ng]*3*angArray[31]-rad[ip+3*ng]*(3*angArray[67]+angArray[80])+rad[ip+4*ng]*angArray[144];
         bas[7] = -rad[ip+ng]*3*angArray[7]+rad[ip+2*ng]*(3*angArray[23]+6*angArray[32])-rad[ip+3*ng]*(6*angArray[68]+angArray[81])+rad[ip+4*ng]*angArray[145];
         bas[8] = -rad[ip+ng]*12*angArray[8]+rad[ip+2*ng]*(9*angArray[33]+12*angArray[24])-rad[ip+3*ng]*(9*angArray[69]+angArray[82])+rad[ip+4*ng]*angArray[146];
         bas[9] = rad[ip]*6*angArray[0]-rad[ip+ng]*(6*angArray[4]+27*angArray[9])+rad[ip+2*ng]*(12*angArray[34]+27*angArray[25])-rad[ip+3*ng]*(12*angArray[70]+angArray[83])+rad[ip+4*ng]*angArray[147];
         bas[10] = -rad[ip+3*ng]*3*angArray[72]+rad[ip+4*ng]*angArray[151];
         bas[11] = rad[ip+2*ng]*3*angArray[26]-rad[ip+3*ng]*6*angArray[73]+rad[ip+4*ng]*angArray[152];
         bas[12] = rad[ip+2*ng]*12*angArray[27]-rad[ip+3*ng]*9*angArray[74]+rad[ip+4*ng]*angArray[153];
         bas[13] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*27*angArray[28]-rad[ip+3*ng]*12*angArray[75]+rad[ip+4*ng]*angArray[154];
         bas[14] = -rad[ip+ng]*24*angArray[6]+rad[ip+2*ng]*48*angArray[29]-rad[ip+3*ng]*15*angArray[76]+rad[ip+4*ng]*angArray[155];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*15*angArray[37]-rad[ip+3*ng]*(3*angArray[86]+5*angArray[93])+rad[ip+4*ng]*angArray[174];
         bas[1] = rad[ip+2*ng]*12*angArray[39]-rad[ip+3*ng]*(3*angArray[88]+4*angArray[97])+rad[ip+4*ng]*angArray[178];
         bas[2] = -rad[ip+ng]*12*angArray[10]+rad[ip+2*ng]*(3*angArray[35]+24*angArray[40])-rad[ip+3*ng]*(6*angArray[89]+4*angArray[98])+rad[ip+4*ng]*angArray[179];
         bas[3] = rad[ip+2*ng]*9*angArray[42]-rad[ip+3*ng]*(3*angArray[91]+3*angArray[102])+rad[ip+4*ng]*angArray[183];
         bas[4] = -rad[ip+ng]*9*angArray[11]+rad[ip+2*ng]*(18*angArray[43]+3*angArray[36])-rad[ip+3*ng]*(6*angArray[92]+3*angArray[103])+rad[ip+4*ng]*angArray[184];
         bas[5] = -rad[ip+ng]*36*angArray[12]+rad[ip+2*ng]*(27*angArray[44]+12*angArray[37])-rad[ip+3*ng]*(3*angArray[104]+9*angArray[93])+rad[ip+4*ng]*angArray[185];
         bas[6] = rad[ip+2*ng]*6*angArray[46]-rad[ip+3*ng]*(2*angArray[108]+3*angArray[95])+rad[ip+4*ng]*angArray[189];
         bas[7] = -rad[ip+ng]*6*angArray[13]+rad[ip+2*ng]*(3*angArray[38]+12*angArray[47])-rad[ip+3*ng]*(2*angArray[109]+6*angArray[96])+rad[ip+4*ng]*angArray[190];
         bas[8] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*(18*angArray[48]+12*angArray[39])-rad[ip+3*ng]*(2*angArray[110]+9*angArray[97])+rad[ip+4*ng]*angArray[191];
         bas[9] = rad[ip]*12*angArray[1]-rad[ip+ng]*(54*angArray[15]+6*angArray[10])+rad[ip+2*ng]*(24*angArray[49]+27*angArray[40])-rad[ip+3*ng]*(12*angArray[98]+2*angArray[111])+rad[ip+4*ng]*angArray[192];
         bas[10] = rad[ip+2*ng]*3*angArray[51]-rad[ip+3*ng]*(3*angArray[100]+angArray[115])+rad[ip+4*ng]*angArray[196];
         bas[11] = -rad[ip+ng]*3*angArray[16]+rad[ip+2*ng]*(3*angArray[41]+6*angArray[52])-rad[ip+3*ng]*(6*angArray[101]+angArray[116])+rad[ip+4*ng]*angArray[197];
         bas[12] = -rad[ip+ng]*12*angArray[17]+rad[ip+2*ng]*(9*angArray[53]+12*angArray[42])-rad[ip+3*ng]*(9*angArray[102]+angArray[117])+rad[ip+4*ng]*angArray[198];
         bas[13] = rad[ip]*6*angArray[2]-rad[ip+ng]*(6*angArray[11]+27*angArray[18])+rad[ip+2*ng]*(12*angArray[54]+27*angArray[43])-rad[ip+3*ng]*(12*angArray[103]+angArray[118])+rad[ip+4*ng]*angArray[199];
         bas[14] = rad[ip]*24*angArray[3]-rad[ip+ng]*(24*angArray[12]+48*angArray[19])+rad[ip+2*ng]*(48*angArray[44]+15*angArray[55])-rad[ip+3*ng]*(15*angArray[104]+angArray[119])+rad[ip+4*ng]*angArray[200];
         bas[15] = -rad[ip+3*ng]*3*angArray[106]+rad[ip+4*ng]*angArray[204];
         bas[16] = rad[ip+2*ng]*3*angArray[45]-rad[ip+3*ng]*6*angArray[107]+rad[ip+4*ng]*angArray[205];
         bas[17] = rad[ip+2*ng]*12*angArray[46]-rad[ip+3*ng]*9*angArray[108]+rad[ip+4*ng]*angArray[206];
         bas[18] = -rad[ip+ng]*6*angArray[13]+rad[ip+2*ng]*27*angArray[47]-rad[ip+3*ng]*12*angArray[109]+rad[ip+4*ng]*angArray[207];
         bas[19] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*48*angArray[48]-rad[ip+3*ng]*15*angArray[110]+rad[ip+4*ng]*angArray[208];
         bas[20] = -rad[ip+ng]*60*angArray[15]+rad[ip+2*ng]*75*angArray[49]-rad[ip+3*ng]*18*angArray[111]+rad[ip+4*ng]*angArray[209];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*18*angArray[58]-rad[ip+3*ng]*(3*angArray[122]+6*angArray[129])+rad[ip+4*ng]*angArray[229];
         bas[1] = rad[ip+2*ng]*15*angArray[60]-rad[ip+3*ng]*(3*angArray[124]+5*angArray[133])+rad[ip+4*ng]*angArray[233];
         bas[2] = -rad[ip+ng]*15*angArray[20]+rad[ip+2*ng]*(30*angArray[61]+3*angArray[56])-rad[ip+3*ng]*(6*angArray[125]+5*angArray[134])+rad[ip+4*ng]*angArray[234];
         bas[3] = rad[ip+2*ng]*12*angArray[63]-rad[ip+3*ng]*(3*angArray[127]+4*angArray[138])+rad[ip+4*ng]*angArray[238];
         bas[4] = -rad[ip+ng]*12*angArray[21]+rad[ip+2*ng]*(3*angArray[57]+24*angArray[64])-rad[ip+3*ng]*(6*angArray[128]+4*angArray[139])+rad[ip+4*ng]*angArray[239];
         bas[5] = -rad[ip+ng]*48*angArray[22]+rad[ip+2*ng]*(12*angArray[58]+36*angArray[65])-rad[ip+3*ng]*(9*angArray[129]+4*angArray[140])+rad[ip+4*ng]*angArray[240];
         bas[6] = rad[ip+2*ng]*9*angArray[67]-rad[ip+3*ng]*(3*angArray[131]+3*angArray[144])+rad[ip+4*ng]*angArray[244];
         bas[7] = -rad[ip+ng]*9*angArray[23]+rad[ip+2*ng]*(18*angArray[68]+3*angArray[59])-rad[ip+3*ng]*(6*angArray[132]+3*angArray[145])+rad[ip+4*ng]*angArray[245];
         bas[8] = -rad[ip+ng]*36*angArray[24]+rad[ip+2*ng]*(27*angArray[69]+12*angArray[60])-rad[ip+3*ng]*(3*angArray[146]+9*angArray[133])+rad[ip+4*ng]*angArray[246];
         bas[9] = rad[ip]*18*angArray[4]-rad[ip+ng]*(81*angArray[25]+6*angArray[20])+rad[ip+2*ng]*(36*angArray[70]+27*angArray[61])-rad[ip+3*ng]*(12*angArray[134]+3*angArray[147])+rad[ip+4*ng]*angArray[247];
         bas[10] = rad[ip+2*ng]*6*angArray[72]-rad[ip+3*ng]*(2*angArray[151]+3*angArray[136])+rad[ip+4*ng]*angArray[251];
         bas[11] = -rad[ip+ng]*6*angArray[26]+rad[ip+2*ng]*(3*angArray[62]+12*angArray[73])-rad[ip+3*ng]*(2*angArray[152]+6*angArray[137])+rad[ip+4*ng]*angArray[252];
         bas[12] = -rad[ip+ng]*24*angArray[27]+rad[ip+2*ng]*(18*angArray[74]+12*angArray[63])-rad[ip+3*ng]*(2*angArray[153]+9*angArray[138])+rad[ip+4*ng]*angArray[253];
         bas[13] = rad[ip]*12*angArray[5]-rad[ip+ng]*(54*angArray[28]+6*angArray[21])+rad[ip+2*ng]*(24*angArray[75]+27*angArray[64])-rad[ip+3*ng]*(12*angArray[139]+2*angArray[154])+rad[ip+4*ng]*angArray[254];
         bas[14] = rad[ip]*48*angArray[6]-rad[ip+ng]*(24*angArray[22]+96*angArray[29])+rad[ip+2*ng]*(30*angArray[76]+48*angArray[65])-rad[ip+3*ng]*(15*angArray[140]+2*angArray[155])+rad[ip+4*ng]*angArray[255];
         bas[15] = rad[ip+2*ng]*3*angArray[78]-rad[ip+3*ng]*(3*angArray[142]+angArray[159])+rad[ip+4*ng]*angArray[259];
         bas[16] = -rad[ip+ng]*3*angArray[30]+rad[ip+2*ng]*(3*angArray[66]+6*angArray[79])-rad[ip+3*ng]*(6*angArray[143]+angArray[160])+rad[ip+4*ng]*angArray[260];
         bas[17] = -rad[ip+ng]*12*angArray[31]+rad[ip+2*ng]*(9*angArray[80]+12*angArray[67])-rad[ip+3*ng]*(9*angArray[144]+angArray[161])+rad[ip+4*ng]*angArray[261];
         bas[18] = rad[ip]*6*angArray[7]-rad[ip+ng]*(6*angArray[23]+27*angArray[32])+rad[ip+2*ng]*(12*angArray[81]+27*angArray[68])-rad[ip+3*ng]*(12*angArray[145]+angArray[162])+rad[ip+4*ng]*angArray[262];
         bas[19] = rad[ip]*24*angArray[8]-rad[ip+ng]*(24*angArray[24]+48*angArray[33])+rad[ip+2*ng]*(48*angArray[69]+15*angArray[82])-rad[ip+3*ng]*(15*angArray[146]+angArray[163])+rad[ip+4*ng]*angArray[263];
         bas[20] = rad[ip]*60*angArray[9]-rad[ip+ng]*(60*angArray[25]+75*angArray[34])+rad[ip+2*ng]*(18*angArray[83]+75*angArray[70])-rad[ip+3*ng]*(18*angArray[147]+angArray[164])+rad[ip+4*ng]*angArray[264];
         bas[21] = -rad[ip+3*ng]*3*angArray[149]+rad[ip+4*ng]*angArray[268];
         bas[22] = rad[ip+2*ng]*3*angArray[71]-rad[ip+3*ng]*6*angArray[150]+rad[ip+4*ng]*angArray[269];
         bas[23] = rad[ip+2*ng]*12*angArray[72]-rad[ip+3*ng]*9*angArray[151]+rad[ip+4*ng]*angArray[270];
         bas[24] = -rad[ip+ng]*6*angArray[26]+rad[ip+2*ng]*27*angArray[73]-rad[ip+3*ng]*12*angArray[152]+rad[ip+4*ng]*angArray[271];
         bas[25] = -rad[ip+ng]*24*angArray[27]+rad[ip+2*ng]*48*angArray[74]-rad[ip+3*ng]*15*angArray[153]+rad[ip+4*ng]*angArray[272];
         bas[26] = -rad[ip+ng]*60*angArray[28]+rad[ip+2*ng]*75*angArray[75]-rad[ip+3*ng]*18*angArray[154]+rad[ip+4*ng]*angArray[273];
         bas[27] = -rad[ip+ng]*120*angArray[29]+rad[ip+2*ng]*108*angArray[76]-rad[ip+3*ng]*21*angArray[155]+rad[ip+4*ng]*angArray[274];
      }

   }


   // now we do derivatives for the given basis set to YZZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[8]+rad[ip+4*ng]*angArray[33];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[14]+rad[ip+4*ng]*angArray[48];
         bas[1] = rad[ip+2*ng]*3*angArray[3]-rad[ip+3*ng]*(3*angArray[17]+angArray[19])+rad[ip+4*ng]*angArray[53];
         bas[2] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*6*angArray[18]+rad[ip+4*ng]*angArray[54];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[24]+rad[ip+4*ng]*angArray[69];
         bas[1] = rad[ip+2*ng]*3*angArray[6]-rad[ip+3*ng]*(3*angArray[27]+angArray[29])+rad[ip+4*ng]*angArray[74];
         bas[2] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*6*angArray[28]+rad[ip+4*ng]*angArray[75];
         bas[3] = rad[ip+2*ng]*6*angArray[8]-rad[ip+3*ng]*(2*angArray[33]+3*angArray[31])+rad[ip+4*ng]*angArray[80];
         bas[4] = -rad[ip+ng]*3*angArray[0]+rad[ip+2*ng]*(3*angArray[7]+6*angArray[9])-rad[ip+3*ng]*(6*angArray[32]+angArray[34])+rad[ip+4*ng]*angArray[81];
         bas[5] = rad[ip+2*ng]*12*angArray[8]-rad[ip+3*ng]*9*angArray[33]+rad[ip+4*ng]*angArray[82];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[39]+rad[ip+4*ng]*angArray[97];
         bas[1] = rad[ip+2*ng]*3*angArray[12]-rad[ip+3*ng]*(3*angArray[42]+angArray[44])+rad[ip+4*ng]*angArray[102];
         bas[2] = rad[ip+2*ng]*3*angArray[11]-rad[ip+3*ng]*6*angArray[43]+rad[ip+4*ng]*angArray[103];
         bas[3] = rad[ip+2*ng]*6*angArray[14]-rad[ip+3*ng]*(2*angArray[48]+3*angArray[46])+rad[ip+4*ng]*angArray[108];
         bas[4] = -rad[ip+ng]*3*angArray[1]+rad[ip+2*ng]*(3*angArray[13]+6*angArray[15])-rad[ip+3*ng]*(6*angArray[47]+angArray[49])+rad[ip+4*ng]*angArray[109];
         bas[5] = rad[ip+2*ng]*12*angArray[14]-rad[ip+3*ng]*9*angArray[48]+rad[ip+4*ng]*angArray[110];
         bas[6] = rad[ip+2*ng]*9*angArray[17]-rad[ip+3*ng]*(3*angArray[51]+3*angArray[53])+rad[ip+4*ng]*angArray[115];
         bas[7] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*(3*angArray[16]+12*angArray[18])-rad[ip+3*ng]*(2*angArray[54]+6*angArray[52])+rad[ip+4*ng]*angArray[116];
         bas[8] = -rad[ip+ng]*12*angArray[3]+rad[ip+2*ng]*(9*angArray[19]+12*angArray[17])-rad[ip+3*ng]*(9*angArray[53]+angArray[55])+rad[ip+4*ng]*angArray[117];
         bas[9] = -rad[ip+ng]*6*angArray[2]+rad[ip+2*ng]*27*angArray[18]-rad[ip+3*ng]*12*angArray[54]+rad[ip+4*ng]*angArray[118];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[60]+rad[ip+4*ng]*angArray[133];
         bas[1] = rad[ip+2*ng]*3*angArray[22]-rad[ip+3*ng]*(3*angArray[63]+angArray[65])+rad[ip+4*ng]*angArray[138];
         bas[2] = rad[ip+2*ng]*3*angArray[21]-rad[ip+3*ng]*6*angArray[64]+rad[ip+4*ng]*angArray[139];
         bas[3] = rad[ip+2*ng]*6*angArray[24]-rad[ip+3*ng]*(2*angArray[69]+3*angArray[67])+rad[ip+4*ng]*angArray[144];
         bas[4] = -rad[ip+ng]*3*angArray[4]+rad[ip+2*ng]*(3*angArray[23]+6*angArray[25])-rad[ip+3*ng]*(6*angArray[68]+angArray[70])+rad[ip+4*ng]*angArray[145];
         bas[5] = rad[ip+2*ng]*12*angArray[24]-rad[ip+3*ng]*9*angArray[69]+rad[ip+4*ng]*angArray[146];
         bas[6] = rad[ip+2*ng]*9*angArray[27]-rad[ip+3*ng]*(3*angArray[72]+3*angArray[74])+rad[ip+4*ng]*angArray[151];
         bas[7] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*(3*angArray[26]+12*angArray[28])-rad[ip+3*ng]*(2*angArray[75]+6*angArray[73])+rad[ip+4*ng]*angArray[152];
         bas[8] = -rad[ip+ng]*12*angArray[6]+rad[ip+2*ng]*(9*angArray[29]+12*angArray[27])-rad[ip+3*ng]*(9*angArray[74]+angArray[76])+rad[ip+4*ng]*angArray[153];
         bas[9] = -rad[ip+ng]*6*angArray[5]+rad[ip+2*ng]*27*angArray[28]-rad[ip+3*ng]*12*angArray[75]+rad[ip+4*ng]*angArray[154];
         bas[10] = rad[ip+2*ng]*12*angArray[31]-rad[ip+3*ng]*(3*angArray[78]+4*angArray[80])+rad[ip+4*ng]*angArray[159];
         bas[11] = -rad[ip+ng]*9*angArray[7]+rad[ip+2*ng]*(18*angArray[32]+3*angArray[30])-rad[ip+3*ng]*(6*angArray[79]+3*angArray[81])+rad[ip+4*ng]*angArray[160];
         bas[12] = -rad[ip+ng]*24*angArray[8]+rad[ip+2*ng]*(18*angArray[33]+12*angArray[31])-rad[ip+3*ng]*(2*angArray[82]+9*angArray[80])+rad[ip+4*ng]*angArray[161];
         bas[13] = rad[ip]*6*angArray[0]-rad[ip+ng]*(6*angArray[7]+27*angArray[9])+rad[ip+2*ng]*(12*angArray[34]+27*angArray[32])-rad[ip+3*ng]*(12*angArray[81]+angArray[83])+rad[ip+4*ng]*angArray[162];
         bas[14] = -rad[ip+ng]*24*angArray[8]+rad[ip+2*ng]*48*angArray[33]-rad[ip+3*ng]*15*angArray[82]+rad[ip+4*ng]*angArray[163];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[88]+rad[ip+4*ng]*angArray[178];
         bas[1] = rad[ip+2*ng]*3*angArray[37]-rad[ip+3*ng]*(3*angArray[91]+angArray[93])+rad[ip+4*ng]*angArray[183];
         bas[2] = rad[ip+2*ng]*3*angArray[36]-rad[ip+3*ng]*6*angArray[92]+rad[ip+4*ng]*angArray[184];
         bas[3] = rad[ip+2*ng]*6*angArray[39]-rad[ip+3*ng]*(2*angArray[97]+3*angArray[95])+rad[ip+4*ng]*angArray[189];
         bas[4] = -rad[ip+ng]*3*angArray[10]+rad[ip+2*ng]*(3*angArray[38]+6*angArray[40])-rad[ip+3*ng]*(6*angArray[96]+angArray[98])+rad[ip+4*ng]*angArray[190];
         bas[5] = rad[ip+2*ng]*12*angArray[39]-rad[ip+3*ng]*9*angArray[97]+rad[ip+4*ng]*angArray[191];
         bas[6] = rad[ip+2*ng]*9*angArray[42]-rad[ip+3*ng]*(3*angArray[100]+3*angArray[102])+rad[ip+4*ng]*angArray[196];
         bas[7] = -rad[ip+ng]*6*angArray[11]+rad[ip+2*ng]*(3*angArray[41]+12*angArray[43])-rad[ip+3*ng]*(2*angArray[103]+6*angArray[101])+rad[ip+4*ng]*angArray[197];
         bas[8] = -rad[ip+ng]*12*angArray[12]+rad[ip+2*ng]*(9*angArray[44]+12*angArray[42])-rad[ip+3*ng]*(9*angArray[102]+angArray[104])+rad[ip+4*ng]*angArray[198];
         bas[9] = -rad[ip+ng]*6*angArray[11]+rad[ip+2*ng]*27*angArray[43]-rad[ip+3*ng]*12*angArray[103]+rad[ip+4*ng]*angArray[199];
         bas[10] = rad[ip+2*ng]*12*angArray[46]-rad[ip+3*ng]*(3*angArray[106]+4*angArray[108])+rad[ip+4*ng]*angArray[204];
         bas[11] = -rad[ip+ng]*9*angArray[13]+rad[ip+2*ng]*(18*angArray[47]+3*angArray[45])-rad[ip+3*ng]*(6*angArray[107]+3*angArray[109])+rad[ip+4*ng]*angArray[205];
         bas[12] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*(18*angArray[48]+12*angArray[46])-rad[ip+3*ng]*(2*angArray[110]+9*angArray[108])+rad[ip+4*ng]*angArray[206];
         bas[13] = rad[ip]*6*angArray[1]-rad[ip+ng]*(6*angArray[13]+27*angArray[15])+rad[ip+2*ng]*(12*angArray[49]+27*angArray[47])-rad[ip+3*ng]*(12*angArray[109]+angArray[111])+rad[ip+4*ng]*angArray[207];
         bas[14] = -rad[ip+ng]*24*angArray[14]+rad[ip+2*ng]*48*angArray[48]-rad[ip+3*ng]*15*angArray[110]+rad[ip+4*ng]*angArray[208];
         bas[15] = rad[ip+2*ng]*15*angArray[51]-rad[ip+3*ng]*(3*angArray[113]+5*angArray[115])+rad[ip+4*ng]*angArray[213];
         bas[16] = -rad[ip+ng]*12*angArray[16]+rad[ip+2*ng]*(3*angArray[50]+24*angArray[52])-rad[ip+3*ng]*(6*angArray[114]+4*angArray[116])+rad[ip+4*ng]*angArray[214];
         bas[17] = -rad[ip+ng]*36*angArray[17]+rad[ip+2*ng]*(27*angArray[53]+12*angArray[51])-rad[ip+3*ng]*(3*angArray[117]+9*angArray[115])+rad[ip+4*ng]*angArray[215];
         bas[18] = rad[ip]*12*angArray[2]-rad[ip+ng]*(54*angArray[18]+6*angArray[16])+rad[ip+2*ng]*(24*angArray[54]+27*angArray[52])-rad[ip+3*ng]*(12*angArray[116]+2*angArray[118])+rad[ip+4*ng]*angArray[216];
         bas[19] = rad[ip]*24*angArray[3]-rad[ip+ng]*(24*angArray[17]+48*angArray[19])+rad[ip+2*ng]*(48*angArray[53]+15*angArray[55])-rad[ip+3*ng]*(15*angArray[117]+angArray[119])+rad[ip+4*ng]*angArray[217];
         bas[20] = -rad[ip+ng]*60*angArray[18]+rad[ip+2*ng]*75*angArray[54]-rad[ip+3*ng]*18*angArray[118]+rad[ip+4*ng]*angArray[218];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = -rad[ip+3*ng]*3*angArray[124]+rad[ip+4*ng]*angArray[233];
         bas[1] = rad[ip+2*ng]*3*angArray[58]-rad[ip+3*ng]*(3*angArray[127]+angArray[129])+rad[ip+4*ng]*angArray[238];
         bas[2] = rad[ip+2*ng]*3*angArray[57]-rad[ip+3*ng]*6*angArray[128]+rad[ip+4*ng]*angArray[239];
         bas[3] = rad[ip+2*ng]*6*angArray[60]-rad[ip+3*ng]*(2*angArray[133]+3*angArray[131])+rad[ip+4*ng]*angArray[244];
         bas[4] = -rad[ip+ng]*3*angArray[20]+rad[ip+2*ng]*(3*angArray[59]+6*angArray[61])-rad[ip+3*ng]*(6*angArray[132]+angArray[134])+rad[ip+4*ng]*angArray[245];
         bas[5] = rad[ip+2*ng]*12*angArray[60]-rad[ip+3*ng]*9*angArray[133]+rad[ip+4*ng]*angArray[246];
         bas[6] = rad[ip+2*ng]*9*angArray[63]-rad[ip+3*ng]*(3*angArray[136]+3*angArray[138])+rad[ip+4*ng]*angArray[251];
         bas[7] = -rad[ip+ng]*6*angArray[21]+rad[ip+2*ng]*(3*angArray[62]+12*angArray[64])-rad[ip+3*ng]*(2*angArray[139]+6*angArray[137])+rad[ip+4*ng]*angArray[252];
         bas[8] = -rad[ip+ng]*12*angArray[22]+rad[ip+2*ng]*(9*angArray[65]+12*angArray[63])-rad[ip+3*ng]*(9*angArray[138]+angArray[140])+rad[ip+4*ng]*angArray[253];
         bas[9] = -rad[ip+ng]*6*angArray[21]+rad[ip+2*ng]*27*angArray[64]-rad[ip+3*ng]*12*angArray[139]+rad[ip+4*ng]*angArray[254];
         bas[10] = rad[ip+2*ng]*12*angArray[67]-rad[ip+3*ng]*(3*angArray[142]+4*angArray[144])+rad[ip+4*ng]*angArray[259];
         bas[11] = -rad[ip+ng]*9*angArray[23]+rad[ip+2*ng]*(18*angArray[68]+3*angArray[66])-rad[ip+3*ng]*(6*angArray[143]+3*angArray[145])+rad[ip+4*ng]*angArray[260];
         bas[12] = -rad[ip+ng]*24*angArray[24]+rad[ip+2*ng]*(18*angArray[69]+12*angArray[67])-rad[ip+3*ng]*(2*angArray[146]+9*angArray[144])+rad[ip+4*ng]*angArray[261];
         bas[13] = rad[ip]*6*angArray[4]-rad[ip+ng]*(6*angArray[23]+27*angArray[25])+rad[ip+2*ng]*(12*angArray[70]+27*angArray[68])-rad[ip+3*ng]*(12*angArray[145]+angArray[147])+rad[ip+4*ng]*angArray[262];
         bas[14] = -rad[ip+ng]*24*angArray[24]+rad[ip+2*ng]*48*angArray[69]-rad[ip+3*ng]*15*angArray[146]+rad[ip+4*ng]*angArray[263];
         bas[15] = rad[ip+2*ng]*15*angArray[72]-rad[ip+3*ng]*(3*angArray[149]+5*angArray[151])+rad[ip+4*ng]*angArray[268];
         bas[16] = -rad[ip+ng]*12*angArray[26]+rad[ip+2*ng]*(3*angArray[71]+24*angArray[73])-rad[ip+3*ng]*(6*angArray[150]+4*angArray[152])+rad[ip+4*ng]*angArray[269];
         bas[17] = -rad[ip+ng]*36*angArray[27]+rad[ip+2*ng]*(27*angArray[74]+12*angArray[72])-rad[ip+3*ng]*(3*angArray[153]+9*angArray[151])+rad[ip+4*ng]*angArray[270];
         bas[18] = rad[ip]*12*angArray[5]-rad[ip+ng]*(54*angArray[28]+6*angArray[26])+rad[ip+2*ng]*(24*angArray[75]+27*angArray[73])-rad[ip+3*ng]*(12*angArray[152]+2*angArray[154])+rad[ip+4*ng]*angArray[271];
         bas[19] = rad[ip]*24*angArray[6]-rad[ip+ng]*(24*angArray[27]+48*angArray[29])+rad[ip+2*ng]*(48*angArray[74]+15*angArray[76])-rad[ip+3*ng]*(15*angArray[153]+angArray[155])+rad[ip+4*ng]*angArray[272];
         bas[20] = -rad[ip+ng]*60*angArray[28]+rad[ip+2*ng]*75*angArray[75]-rad[ip+3*ng]*18*angArray[154]+rad[ip+4*ng]*angArray[273];
         bas[21] = rad[ip+2*ng]*18*angArray[78]-rad[ip+3*ng]*(3*angArray[157]+6*angArray[159])+rad[ip+4*ng]*angArray[278];
         bas[22] = -rad[ip+ng]*15*angArray[30]+rad[ip+2*ng]*(30*angArray[79]+3*angArray[77])-rad[ip+3*ng]*(6*angArray[158]+5*angArray[160])+rad[ip+4*ng]*angArray[279];
         bas[23] = -rad[ip+ng]*48*angArray[31]+rad[ip+2*ng]*(12*angArray[78]+36*angArray[80])-rad[ip+3*ng]*(9*angArray[159]+4*angArray[161])+rad[ip+4*ng]*angArray[280];
         bas[24] = rad[ip]*18*angArray[7]-rad[ip+ng]*(81*angArray[32]+6*angArray[30])+rad[ip+2*ng]*(36*angArray[81]+27*angArray[79])-rad[ip+3*ng]*(12*angArray[160]+3*angArray[162])+rad[ip+4*ng]*angArray[281];
         bas[25] = rad[ip]*48*angArray[8]-rad[ip+ng]*(24*angArray[31]+96*angArray[33])+rad[ip+2*ng]*(30*angArray[82]+48*angArray[80])-rad[ip+3*ng]*(15*angArray[161]+2*angArray[163])+rad[ip+4*ng]*angArray[282];
         bas[26] = rad[ip]*60*angArray[9]-rad[ip+ng]*(60*angArray[32]+75*angArray[34])+rad[ip+2*ng]*(18*angArray[83]+75*angArray[81])-rad[ip+3*ng]*(18*angArray[162]+angArray[164])+rad[ip+4*ng]*angArray[283];
         bas[27] = -rad[ip+ng]*120*angArray[33]+rad[ip+2*ng]*108*angArray[82]-rad[ip+3*ng]*21*angArray[163]+rad[ip+4*ng]*angArray[284];
      }

   }


   // now we do derivatives for the given basis set to ZZZZ
   basis = basis + ng*nBas; 

   if(L == 0) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[0]-rad[ip+3*ng]*6*angArray[9]+rad[ip+4*ng]*angArray[34];
      }

   } else if(L == 1) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[1]-rad[ip+3*ng]*6*angArray[15]+rad[ip+4*ng]*angArray[49];
         bas[1] = rad[ip+2*ng]*3*angArray[2]-rad[ip+3*ng]*6*angArray[18]+rad[ip+4*ng]*angArray[54];
         bas[2] = rad[ip+2*ng]*15*angArray[3]-rad[ip+3*ng]*10*angArray[19]+rad[ip+4*ng]*angArray[55];
      }

   } else if(L == 2) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[4]-rad[ip+3*ng]*6*angArray[25]+rad[ip+4*ng]*angArray[70];
         bas[1] = rad[ip+2*ng]*3*angArray[5]-rad[ip+3*ng]*6*angArray[28]+rad[ip+4*ng]*angArray[75];
         bas[2] = rad[ip+2*ng]*15*angArray[6]-rad[ip+3*ng]*10*angArray[29]+rad[ip+4*ng]*angArray[76];
         bas[3] = rad[ip+2*ng]*3*angArray[7]-rad[ip+3*ng]*6*angArray[32]+rad[ip+4*ng]*angArray[81];
         bas[4] = rad[ip+2*ng]*15*angArray[8]-rad[ip+3*ng]*10*angArray[33]+rad[ip+4*ng]*angArray[82];
         bas[5] = -rad[ip+ng]*12*angArray[0]+rad[ip+2*ng]*39*angArray[9]-rad[ip+3*ng]*14*angArray[34]+rad[ip+4*ng]*angArray[83];
      }

   } else if(L == 3) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[10]-rad[ip+3*ng]*6*angArray[40]+rad[ip+4*ng]*angArray[98];
         bas[1] = rad[ip+2*ng]*3*angArray[11]-rad[ip+3*ng]*6*angArray[43]+rad[ip+4*ng]*angArray[103];
         bas[2] = rad[ip+2*ng]*15*angArray[12]-rad[ip+3*ng]*10*angArray[44]+rad[ip+4*ng]*angArray[104];
         bas[3] = rad[ip+2*ng]*3*angArray[13]-rad[ip+3*ng]*6*angArray[47]+rad[ip+4*ng]*angArray[109];
         bas[4] = rad[ip+2*ng]*15*angArray[14]-rad[ip+3*ng]*10*angArray[48]+rad[ip+4*ng]*angArray[110];
         bas[5] = -rad[ip+ng]*12*angArray[1]+rad[ip+2*ng]*39*angArray[15]-rad[ip+3*ng]*14*angArray[49]+rad[ip+4*ng]*angArray[111];
         bas[6] = rad[ip+2*ng]*3*angArray[16]-rad[ip+3*ng]*6*angArray[52]+rad[ip+4*ng]*angArray[116];
         bas[7] = rad[ip+2*ng]*15*angArray[17]-rad[ip+3*ng]*10*angArray[53]+rad[ip+4*ng]*angArray[117];
         bas[8] = -rad[ip+ng]*12*angArray[2]+rad[ip+2*ng]*39*angArray[18]-rad[ip+3*ng]*14*angArray[54]+rad[ip+4*ng]*angArray[118];
         bas[9] = -rad[ip+ng]*60*angArray[3]+rad[ip+2*ng]*75*angArray[19]-rad[ip+3*ng]*18*angArray[55]+rad[ip+4*ng]*angArray[119];
      }

   } else if(L == 4) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[20]-rad[ip+3*ng]*6*angArray[61]+rad[ip+4*ng]*angArray[134];
         bas[1] = rad[ip+2*ng]*3*angArray[21]-rad[ip+3*ng]*6*angArray[64]+rad[ip+4*ng]*angArray[139];
         bas[2] = rad[ip+2*ng]*15*angArray[22]-rad[ip+3*ng]*10*angArray[65]+rad[ip+4*ng]*angArray[140];
         bas[3] = rad[ip+2*ng]*3*angArray[23]-rad[ip+3*ng]*6*angArray[68]+rad[ip+4*ng]*angArray[145];
         bas[4] = rad[ip+2*ng]*15*angArray[24]-rad[ip+3*ng]*10*angArray[69]+rad[ip+4*ng]*angArray[146];
         bas[5] = -rad[ip+ng]*12*angArray[4]+rad[ip+2*ng]*39*angArray[25]-rad[ip+3*ng]*14*angArray[70]+rad[ip+4*ng]*angArray[147];
         bas[6] = rad[ip+2*ng]*3*angArray[26]-rad[ip+3*ng]*6*angArray[73]+rad[ip+4*ng]*angArray[152];
         bas[7] = rad[ip+2*ng]*15*angArray[27]-rad[ip+3*ng]*10*angArray[74]+rad[ip+4*ng]*angArray[153];
         bas[8] = -rad[ip+ng]*12*angArray[5]+rad[ip+2*ng]*39*angArray[28]-rad[ip+3*ng]*14*angArray[75]+rad[ip+4*ng]*angArray[154];
         bas[9] = -rad[ip+ng]*60*angArray[6]+rad[ip+2*ng]*75*angArray[29]-rad[ip+3*ng]*18*angArray[76]+rad[ip+4*ng]*angArray[155];
         bas[10] = rad[ip+2*ng]*3*angArray[30]-rad[ip+3*ng]*6*angArray[79]+rad[ip+4*ng]*angArray[160];
         bas[11] = rad[ip+2*ng]*15*angArray[31]-rad[ip+3*ng]*10*angArray[80]+rad[ip+4*ng]*angArray[161];
         bas[12] = -rad[ip+ng]*12*angArray[7]+rad[ip+2*ng]*39*angArray[32]-rad[ip+3*ng]*14*angArray[81]+rad[ip+4*ng]*angArray[162];
         bas[13] = -rad[ip+ng]*60*angArray[8]+rad[ip+2*ng]*75*angArray[33]-rad[ip+3*ng]*18*angArray[82]+rad[ip+4*ng]*angArray[163];
         bas[14] = rad[ip]*24*angArray[0]-rad[ip+ng]*168*angArray[9]+rad[ip+2*ng]*123*angArray[34]-rad[ip+3*ng]*22*angArray[83]+rad[ip+4*ng]*angArray[164];
      }

   } else if(L == 5) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[35]-rad[ip+3*ng]*6*angArray[89]+rad[ip+4*ng]*angArray[179];
         bas[1] = rad[ip+2*ng]*3*angArray[36]-rad[ip+3*ng]*6*angArray[92]+rad[ip+4*ng]*angArray[184];
         bas[2] = rad[ip+2*ng]*15*angArray[37]-rad[ip+3*ng]*10*angArray[93]+rad[ip+4*ng]*angArray[185];
         bas[3] = rad[ip+2*ng]*3*angArray[38]-rad[ip+3*ng]*6*angArray[96]+rad[ip+4*ng]*angArray[190];
         bas[4] = rad[ip+2*ng]*15*angArray[39]-rad[ip+3*ng]*10*angArray[97]+rad[ip+4*ng]*angArray[191];
         bas[5] = -rad[ip+ng]*12*angArray[10]+rad[ip+2*ng]*39*angArray[40]-rad[ip+3*ng]*14*angArray[98]+rad[ip+4*ng]*angArray[192];
         bas[6] = rad[ip+2*ng]*3*angArray[41]-rad[ip+3*ng]*6*angArray[101]+rad[ip+4*ng]*angArray[197];
         bas[7] = rad[ip+2*ng]*15*angArray[42]-rad[ip+3*ng]*10*angArray[102]+rad[ip+4*ng]*angArray[198];
         bas[8] = -rad[ip+ng]*12*angArray[11]+rad[ip+2*ng]*39*angArray[43]-rad[ip+3*ng]*14*angArray[103]+rad[ip+4*ng]*angArray[199];
         bas[9] = -rad[ip+ng]*60*angArray[12]+rad[ip+2*ng]*75*angArray[44]-rad[ip+3*ng]*18*angArray[104]+rad[ip+4*ng]*angArray[200];
         bas[10] = rad[ip+2*ng]*3*angArray[45]-rad[ip+3*ng]*6*angArray[107]+rad[ip+4*ng]*angArray[205];
         bas[11] = rad[ip+2*ng]*15*angArray[46]-rad[ip+3*ng]*10*angArray[108]+rad[ip+4*ng]*angArray[206];
         bas[12] = -rad[ip+ng]*12*angArray[13]+rad[ip+2*ng]*39*angArray[47]-rad[ip+3*ng]*14*angArray[109]+rad[ip+4*ng]*angArray[207];
         bas[13] = -rad[ip+ng]*60*angArray[14]+rad[ip+2*ng]*75*angArray[48]-rad[ip+3*ng]*18*angArray[110]+rad[ip+4*ng]*angArray[208];
         bas[14] = rad[ip]*24*angArray[1]-rad[ip+ng]*168*angArray[15]+rad[ip+2*ng]*123*angArray[49]-rad[ip+3*ng]*22*angArray[111]+rad[ip+4*ng]*angArray[209];
         bas[15] = rad[ip+2*ng]*3*angArray[50]-rad[ip+3*ng]*6*angArray[114]+rad[ip+4*ng]*angArray[214];
         bas[16] = rad[ip+2*ng]*15*angArray[51]-rad[ip+3*ng]*10*angArray[115]+rad[ip+4*ng]*angArray[215];
         bas[17] = -rad[ip+ng]*12*angArray[16]+rad[ip+2*ng]*39*angArray[52]-rad[ip+3*ng]*14*angArray[116]+rad[ip+4*ng]*angArray[216];
         bas[18] = -rad[ip+ng]*60*angArray[17]+rad[ip+2*ng]*75*angArray[53]-rad[ip+3*ng]*18*angArray[117]+rad[ip+4*ng]*angArray[217];
         bas[19] = rad[ip]*24*angArray[2]-rad[ip+ng]*168*angArray[18]+rad[ip+2*ng]*123*angArray[54]-rad[ip+3*ng]*22*angArray[118]+rad[ip+4*ng]*angArray[218];
         bas[20] = rad[ip]*120*angArray[3]-rad[ip+ng]*360*angArray[19]+rad[ip+2*ng]*183*angArray[55]-rad[ip+3*ng]*26*angArray[119]+rad[ip+4*ng]*angArray[219];
      }

   } else if(L == 6) {

      for(UInt ip = 0; ip<ng; ip++) {
         Double* bas = &basis[ip*nBas];
         const Double* angArray = &ang[ip*nTolCarBas];
         bas[0] = rad[ip+2*ng]*3*angArray[56]-rad[ip+3*ng]*6*angArray[125]+rad[ip+4*ng]*angArray[234];
         bas[1] = rad[ip+2*ng]*3*angArray[57]-rad[ip+3*ng]*6*angArray[128]+rad[ip+4*ng]*angArray[239];
         bas[2] = rad[ip+2*ng]*15*angArray[58]-rad[ip+3*ng]*10*angArray[129]+rad[ip+4*ng]*angArray[240];
         bas[3] = rad[ip+2*ng]*3*angArray[59]-rad[ip+3*ng]*6*angArray[132]+rad[ip+4*ng]*angArray[245];
         bas[4] = rad[ip+2*ng]*15*angArray[60]-rad[ip+3*ng]*10*angArray[133]+rad[ip+4*ng]*angArray[246];
         bas[5] = -rad[ip+ng]*12*angArray[20]+rad[ip+2*ng]*39*angArray[61]-rad[ip+3*ng]*14*angArray[134]+rad[ip+4*ng]*angArray[247];
         bas[6] = rad[ip+2*ng]*3*angArray[62]-rad[ip+3*ng]*6*angArray[137]+rad[ip+4*ng]*angArray[252];
         bas[7] = rad[ip+2*ng]*15*angArray[63]-rad[ip+3*ng]*10*angArray[138]+rad[ip+4*ng]*angArray[253];
         bas[8] = -rad[ip+ng]*12*angArray[21]+rad[ip+2*ng]*39*angArray[64]-rad[ip+3*ng]*14*angArray[139]+rad[ip+4*ng]*angArray[254];
         bas[9] = -rad[ip+ng]*60*angArray[22]+rad[ip+2*ng]*75*angArray[65]-rad[ip+3*ng]*18*angArray[140]+rad[ip+4*ng]*angArray[255];
         bas[10] = rad[ip+2*ng]*3*angArray[66]-rad[ip+3*ng]*6*angArray[143]+rad[ip+4*ng]*angArray[260];
         bas[11] = rad[ip+2*ng]*15*angArray[67]-rad[ip+3*ng]*10*angArray[144]+rad[ip+4*ng]*angArray[261];
         bas[12] = -rad[ip+ng]*12*angArray[23]+rad[ip+2*ng]*39*angArray[68]-rad[ip+3*ng]*14*angArray[145]+rad[ip+4*ng]*angArray[262];
         bas[13] = -rad[ip+ng]*60*angArray[24]+rad[ip+2*ng]*75*angArray[69]-rad[ip+3*ng]*18*angArray[146]+rad[ip+4*ng]*angArray[263];
         bas[14] = rad[ip]*24*angArray[4]-rad[ip+ng]*168*angArray[25]+rad[ip+2*ng]*123*angArray[70]-rad[ip+3*ng]*22*angArray[147]+rad[ip+4*ng]*angArray[264];
         bas[15] = rad[ip+2*ng]*3*angArray[71]-rad[ip+3*ng]*6*angArray[150]+rad[ip+4*ng]*angArray[269];
         bas[16] = rad[ip+2*ng]*15*angArray[72]-rad[ip+3*ng]*10*angArray[151]+rad[ip+4*ng]*angArray[270];
         bas[17] = -rad[ip+ng]*12*angArray[26]+rad[ip+2*ng]*39*angArray[73]-rad[ip+3*ng]*14*angArray[152]+rad[ip+4*ng]*angArray[271];
         bas[18] = -rad[ip+ng]*60*angArray[27]+rad[ip+2*ng]*75*angArray[74]-rad[ip+3*ng]*18*angArray[153]+rad[ip+4*ng]*angArray[272];
         bas[19] = rad[ip]*24*angArray[5]-rad[ip+ng]*168*angArray[28]+rad[ip+2*ng]*123*angArray[75]-rad[ip+3*ng]*22*angArray[154]+rad[ip+4*ng]*angArray[273];
         bas[20] = rad[ip]*120*angArray[6]-rad[ip+ng]*360*angArray[29]+rad[ip+2*ng]*183*angArray[76]-rad[ip+3*ng]*26*angArray[155]+rad[ip+4*ng]*angArray[274];
         bas[21] = rad[ip+2*ng]*3*angArray[77]-rad[ip+3*ng]*6*angArray[158]+rad[ip+4*ng]*angArray[279];
         bas[22] = rad[ip+2*ng]*15*angArray[78]-rad[ip+3*ng]*10*angArray[159]+rad[ip+4*ng]*angArray[280];
         bas[23] = -rad[ip+ng]*12*angArray[30]+rad[ip+2*ng]*39*angArray[79]-rad[ip+3*ng]*14*angArray[160]+rad[ip+4*ng]*angArray[281];
         bas[24] = -rad[ip+ng]*60*angArray[31]+rad[ip+2*ng]*75*angArray[80]-rad[ip+3*ng]*18*angArray[161]+rad[ip+4*ng]*angArray[282];
         bas[25] = rad[ip]*24*angArray[7]-rad[ip+ng]*168*angArray[32]+rad[ip+2*ng]*123*angArray[81]-rad[ip+3*ng]*22*angArray[162]+rad[ip+4*ng]*angArray[283];
         bas[26] = rad[ip]*120*angArray[8]-rad[ip+ng]*360*angArray[33]+rad[ip+2*ng]*183*angArray[82]-rad[ip+3*ng]*26*angArray[163]+rad[ip+4*ng]*angArray[284];
         bas[27] = rad[ip]*360*angArray[9]-rad[ip+ng]*660*angArray[34]+rad[ip+2*ng]*255*angArray[83]-rad[ip+3*ng]*30*angArray[164]+rad[ip+4*ng]*angArray[285];
      }

   }


}


