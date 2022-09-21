//
// note:
// set up integral engine information reading test
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<string> 
#include<vector>
#include<cmath>
#include "libgen.h"
#include "gintsinfor.h"
#include "integralinfor.h"
using namespace std;
using namespace gintsinfor;
using namespace integralinfor;

void integralinfor_test()
{
	UInt oper = COULOMB;
	UInt order = 2;
	IntegralInfor infor(oper,order);
	infor.print();
}
