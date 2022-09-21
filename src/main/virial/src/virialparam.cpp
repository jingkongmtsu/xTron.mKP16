#include "textread.h"
#include "parameterparsing.h"
#include "cluster.h"
#include "virialparam.h"
using namespace textread;
using namespace parameterparsing;
using namespace cluster;
using namespace virialparam;

Int VirialParam::findNumMonomers(const string& input) const
{

	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		crash(true, "Can not open input in VirialParam class");
	}

	// find how many sections we have
	Int n = 0;
	string line;
	WordConvert w;
	inf.seekg(0,ios::beg);
	string secKey = "molecule";
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isSec() && w.compare(l.getSecName(), secKey)) {
			n++;
		}
	}

	// close the input file
	inf.close();

	// double check the n we find - should be at least 1
	if (n<1) {
		cout << "The monomer number we found: " << n << endl;
		crash(true, "Something wrong in findNumMonomers of VirialParam class");
	}

	// now return what we find
	return n;
}

VirialParam::VirialParam(const string& input):nMonomers(0),method(PERTURBED_SCF_WAY_VIRIAL)
{
	//
	// firstly, build the culster from given input file
	//
	Cluster cluster(input);

	//
	// figure out how many monomers  
	//
	nMonomers = findNumMonomers(input);

	//
	// now deal with input key words  
	// 
	Int section = cluster.getSec();
	ParameterParsing pp(input,"virial",section);
	if (pp.hasAnyParameters()) {

		// Virial coefficient calculation methods
		string key = "virial_method";
		WordConvert w;
		if (pp.hasKeyDefined(key)) {
			string value = pp.getValue(key);
			w.capitalize(value);
			if (value == "PERTURBED_SCF") {
				method = PERTURBED_SCF_WAY_VIRIAL;
			}else if (value == "FULL_SCF") {
				method = FULL_SCF_WAY_VIRIAL;
			}else{
				crash(true, "In VirialParam virial_method choice is invalid.");
			}
		}
	}
}


