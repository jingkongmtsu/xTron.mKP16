//
// Jan. 2014
// this program is used to set the correct 
// shell pair sequence for any two given 
// shell pairs
// author fenglai liu
//

#include<stdlib.h>
#include "libgen.h"
#include "shellprop.h"
using namespace shellprop;

bool switchShellPair(Int ang1, Int ang2);

//
// the code inside is actually taken from the shellprop.h in the 
// shell include folder, from the function of switchShellPair
// if it returns true, it means that ang1 should be ahead of ang2
// in the codeArray, see main function
//
bool switchShellPair(Int ang1, Int ang2)
{
	// decode 
	Int braLmin1 = -1;
	Int braLmax1 = -1;
	Int braLmin2 = -1;
	Int braLmax2 = -1;
	decodeL(ang1,braLmin1,braLmax1,braLmin2,braLmax2);

	// decode 
	Int ketLmin1 = -1;
	Int ketLmax1 = -1;
	Int ketLmin2 = -1;
	Int ketLmax2 = -1;
	decodeL(ang2,ketLmin1,ketLmax1,ketLmin2,ketLmax2);

	// consider rule 3
	if (ketLmax1+ketLmax2>braLmax1+braLmax2) return true;

	// consider rule 4 and 5
	if (ketLmax1+ketLmax2==braLmax1+braLmax2) {

		// rule 4
		// compare bra2 and ket2
		// for bra2/ket2, we will take the bigger one as in the bra side
		if (braLmax2 < ketLmax2) return true;
		if (braLmax2 > ketLmax2) return false;

		// now we have braLmax2 == ketLmax2
		// compare the lmin
		if (braLmin2 < ketLmin2) return true;
		if (braLmin2 > ketLmin2) return false;

		// rule 5
		// now in this case, we must have:
		// braLmax1 == ketLmax1
		// braLmax2 == ketLmax2
		// braLmin2 == ketLmin2
		if (braLmin1 < ketLmin1) return true;
	}

	//
	// all of other cases
	//
	return false;
}

Int main(Int argc, char *argv[]) 
{
	// test the number of args
	if (argc != 2) {
		printf("you need to give the maximum angular momentum value for the shell pair array generation\n");
		printf("for example: \"./test 4\" means you want to max_l up to the G shell\n");
		printf("therefore all of shell pair only have shell data up to G shell\n");
		exit(0);
	}
	Int max_L = atoi(argv[1]);

	// additional check
	if (max_L <= 0 || max_L >= SHELL_ANG_MOM_CODE[MAX_SHELL_TYPES-1]) {
		printf("invalid max_L value you provided through command line: %d\n", max_L);
		printf("it should be >0 and less than the highest boundary value defined in SHELL_ANG_MOM_CODE\n");
		exit(0);
	}

	Int maxL1= max_L;
	Int maxL2= max_L;
	Int maxLPair = maxL1 + maxL2;

	vector<Int> angCodes;
	angCodes.reserve(MAX_SHELL_TYPES*MAX_SHELL_TYPES);
	for(Int totalLPair=0; totalLPair<=maxLPair; totalLPair++) {

		//
		// firstly, let's form the shell pair code
		// for each total L sum
		// on the other hand, for each shell pair
		// the two shell positions are fixed.
		// that is to say, for the shell pair (i,j|:
		// shell i's position must be equal to 
		// or higher than shell j's position
		// the shell position information is 
		// defined in the SHELL_ANG_MOM_CODE
		// see the shellprop.h for more information
		//
		vector<Int> codeArray;
		for(Int i=0; i<MAX_SHELL_TYPES; i++) {
			for(Int j=i; j<MAX_SHELL_TYPES; j++) {

				// get the angular momentum code
				Int L1Code = SHELL_ANG_MOM_CODE[j];
				Int L2Code = SHELL_ANG_MOM_CODE[i];

				// decode 
				Int L1Min = -1;
				Int L1Max = -1;
				Int L2Min = -1;
				Int L2Max = -1;
				decodeL(L1Code,L1Min,L1Max);
				decodeL(L2Code,L2Min,L2Max);

				// L1+L2 = totalLPair
				if (L1Max+L2Max != totalLPair) continue;

				// should not exceed the maxL value
				if (L1Max>maxL1) continue;
				if (L2Max>maxL2) continue;

				// We should have L(bra1) >= L(bra2)
				if (L1Max<L2Max) continue;
				if (L1Max==L2Max && L1Min<L2Min) continue;

				// now let's push in the new code
				// make sure that there's no two LCodes 
				// are same
				Int ang = codeL(L1Min,L1Max,L2Min,L2Max);

				//
				// add in
				//
				if (codeArray.size() == 0) {
					codeArray.push_back(ang);
				} else {

					// comapre and inseret
					for(Int iCode=0; iCode<codeArray.size(); iCode++) {

						// check whether we have same code?
						if (codeArray[iCode] == ang) {
							cout << codeArray[iCode] << endl;
							cout << ang << endl;
							cout << "we should not have same ang appear in the array" << endl;
							exit(1);
						}

						// now find the proper position for insert
						Int ang1 = ang;
						Int ang2 = codeArray[iCode];
						if (switchShellPair(ang1,ang2)) {
							vector<Int>::iterator it = codeArray.begin() + iCode;
							codeArray.insert(it,ang);
							break;
						}

						// finally, we will push it back to the end
						codeArray.push_back(ang);
						break;
					}
				}
			}
		}

		// now let's push the code Array into the angCodes
		for(Int iCode=0; iCode<codeArray.size(); iCode++) {
			Int ang = codeArray[iCode];
			angCodes.push_back(ang);
		}
	}

	// now let's take a look of the ang codes
	// or just print out the results
	cout << "const Int MAX_SHELL_PAIR_NUMBER = " << angCodes.size() << ";" << endl;
	cout << "const Int SHELL_PAIR_ORDER_ARRAY[ ]  =  { " << endl;
	for(Int iCode=0; iCode<angCodes.size(); iCode++) {
		Int ang = angCodes[iCode];
		Int L1Min = -1;
		Int L1Max = -1;
		Int L2Min = -1;
		Int L2Max = -1;
		decodeL(ang,L1Min,L1Max,L2Min,L2Max);
		printf("%6d  %-4s  %-2d  %-2d  %-2d  %-2d\n", ang, ", //", L1Min,L1Max,L2Min,L2Max);
	}
	cout << "}" << endl;

	return 0;

}
