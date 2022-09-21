/**
 * \file    basissetnames.h
 * \brief   parse and read in basis set files name
 * \author  Fenglai Liu and Jing Kong
 *
 * defined by the env variable "EMUL_BASIS_SET_HOME", it's a folder that contains the 
 * standard basis set files data. All of the standard defined basis sets data, like
 * Pople basis sets 6-31G, or Dunning basis sets cc-pvdz etc. are all well defined
 * in that folder. Therefore the user can give a basis set name string, then we 
 * will know what kind of basis set the user want to use.
 *
 * In that folder we also have an important file, basis.conf; which is generated from
 * a perl script in that folder (basisgeneration.pl). It will generate the basis set
 * names for the basis set data file, so each line we will have two fields of information:
 *
 * basis set name  basis set file name   
 *
 * Hence through this way we set up one to one mapping relation between the basis set
 * name, the basis set file name. 
 *
 * In this file we defined several functions. Based on the configuration file basis.conf
 * we are able to process the basis set name(used by the user in the input file), and 
 * the real basis set data file (so that we can locate the data).
 *
 * For the basis set file name, all of files are with suffix of ".bas". However, here
 * the basis set name are not with this suffix. Please make sure of that.
 */

#ifndef BASISSETNAMES_H
#define BASISSETNAMES_H
#include "libgen.h"
#include <string>

namespace basissetnames {

	using namespace std;

	///
	/// this is to return the basis set code and the file name from the input basis set name
	///
	string getBasisSetFileNameFromInput(const string& inputBasisSetName);

	///
	/// this is to return whether the given name is a basis set name?
	///
	bool isValidBasisSetName(const string& name);
}

#endif
