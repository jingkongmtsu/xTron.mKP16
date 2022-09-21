/**
 * \file    parseshell.h
 * \brief   classes for parsing the shell data 
 * \author  Fenglai Liu and Jing Kong
 * 
 * for all of the classes in this module, we just simply use public.
 * this is for easy access to the data.
 *
 * Basically, the shell data read in is just following the natrual way:
 * shell -> atomShell -> molecule shell
 * just as what we do in the shell class. 
 *
 * Currently, the shell data format we support to read in is Gaussian
 * format and NWChem library format (the format that NWChem program uses).
 * The NWChem library format is directly used in NWChem program, therefore
 * you can copy the basis set data file directly from NWChem program
 * to the program here. We can use it directly.
 *
 * For the basis set data, there are two ways you can use it. The first
 * way is that you define (or standard defined) shell data and put them
 * into a single file (with suffix of .bas) to use. In that way you need
 * to access the env variable EMUL_BASIS_SET_HOME, all of the basis set
 * files are stored in that directory.
 *
 * The other way to use a basis set data is you defined it in the user
 * input file, in the %basis  %end  section. Here we support both of the 
 * two ways.
 */

#ifndef PARSESHELL_H
#define PARSESHELL_H
#include "libgen.h"
#include <vector>
#include <string>
#include "scalablevec.h"
#include "molecule.h"

namespace parseshell {

	using namespace std;
	using namespace molecule;

	///
	/// this is the default invalid format of shell data, used to indicate
	/// that something wrong happens
	///
	/// however, if the user uses the standard definded shell data; then 
	/// in the input file the basis set input is also "invalid". We use
	/// this status to indicate that the shell data is not in the user
	/// input file
	///
	const UInt INVALID_BASIS_INPUT  = 0;

	///
	/// this is the constant to define the Gaussian format shell data as input
	///
	const UInt GAUSSIAN_BASIS_INPUT = 1;

	///
	/// this is the constant to define the NWChem format shell data as input
	///
	const UInt NWCHEM_BASIS_INPUT   = 2;

	///
	/// this is a utility function to find out the basis set style in the input
	/// file
	///
	/// if the input loc = 0, the search will begin from the beginning of the 
	/// file. Else it will get to the location, and search the basis set 
	/// style after the location.
	///
	UInt findBasisSetStyle(const string& input, const UInt& loc);

	///
	/// define the data structure of raw atom shell
	///
	class RawShell;
	typedef std::vector<RawShell,tbb::scalable_allocator<RawShell> > RawShellVec;

	/*******************************************************************
	 *                       Raw Shell Data                            *
	 *******************************************************************/

	/**
	 * \class RawShell
	 * \brief This class is used to describe the shell data read in
	 */
	class RawShell {

		public:

			bool    isPure;       ///< whether this shell is pure or cartesian 
			UInt     ngau;        ///< number of gaussian primitives in the shell
			Double  scale;        ///< scaled factor in the data section
			string  shellType;    ///< shell type, S, SP, P, D etc. 
			DoubleVec exp;        ///< read in exponents for the shell
			DoubleVec coe;        ///< read in contraction coefficients 
			
			///
			/// default constructor, we need a couple of information for initialization
			///
			RawShell(bool pureState, const Double& scalef, const UInt& nTypeGau, 
					const UInt& ng, const string& type):isPure(pureState),
			ngau(ng),scale(scalef),shellType(type),exp(ng,ZERO),coe(ng*nTypeGau,ZERO) { 
				// additional check for SP shell etc.
				// for S and P, Cartesian and spherical are same
				if(shellType == "S" || shellType == "P" || shellType == "SP") {
					isPure = false;
				}
		  	};

			///
			/// destructor
			///
			~RawShell() { };

			///
			/// split the input SP shell data into a "S" shell or a "P" shell
			///
			void splitSPShellData(const RawShell& s);
	};


	/**
	 * \class RawAtomShell
	 * \brief This class is used to describe the atom shell data read in
	 */
	class RawAtomShell {

		public:

			UInt atomic;                 ///< atom types
			RawShellVec shells;          ///< one atom's shell data

			/**
			 * this is the driver function to read in the atom shell data 
			 *
			 * \param input   this is the input data file, where we can read the raw atom shell data
			 * \param dataLoc if user defined the shell data in the input file, then we need the data location
			 *                related to the molecule section. So this is what dataLoc means. If this is the 
			 *                shell data file(*.bas file), then this dataLoc is 0.
			 */
			void readRawAtomShellData(const string& input, const UInt& dataLoc);

			/**
			 * read in Gaussian format of data
			 *
			 * \param input   this is the input data file, where we can read the raw atom shell data
			 * \param dataLoc if user defined the shell data in the input file, then we need the data location
			 *                related to the molecule section. So this is what dataLoc means. If this is the 
			 *                shell data file(*.bas file), then this dataLoc is 0.
			 */
			void readGaussianShellData(const string& input, const UInt& dataLoc);

			/**
			 * read in NWChem format of data
			 *
			 * \param input   this is the input data file, where we can read the raw atom shell data
			 * \param dataLoc if user defined the shell data in the input file, then we need the data location
			 *                related to the molecule section. So this is what dataLoc means. If this is the 
			 *                shell data file(*.bas file), then this dataLoc is 0.
			 */
			void readNWChemShellData(const string& input, const UInt& dataLoc);

			/**
			 * get number of shells
			 */
			UInt getNShell() const { return shells.size(); };

			/**
			 * get the shell element
			 */
			const RawShell& getShell(UInt i) const { return shells[i]; };

			/**
			 * constructor 
			 */
			RawAtomShell(const string& input, const UInt& dataLoc, const UInt& atom0):atomic(atom0) { 
				readRawAtomShellData(input,dataLoc);
			};

			/**
			 * destructor
			 */
			~RawAtomShell() { };
	};

	/**
	 * \class RawMolShell
	 * \brief raw molecule shell data
	 */
	class RawMolShell {

		public:

			/**
			 * The whole shell data read from file for the given molecule
			 * the length is number of type atoms
			 */ 
			vector<RawAtomShell> molShells;    

			/**
			 * constructor
			 * \param input  input file
			 * \param mol    molecular geometry
			 * \param part   which part of basis set are we going to form? must >= 0!
			 */
			RawMolShell(const string& input, const Molecule& mol, UInt part, bool more_debug);

			/**
			 */
			~RawMolShell(){ };

			/**
			 * return the shell data beginning position for given input file,
			 * section and part number
			 */
			UInt parsingLoc(const string& input, UInt section, UInt part) const; 

			///
			/// through the input file and given location, we will read the basis set
			/// data information per each element, or maybe basis set data for all of 
			/// elements
			///
			void parsingBasisSetInfor(const string& input, const UInt& loc, 
					vector<string>& basisSetNames, UIntVec& elementNames) const;

			/**
			 * get the raw atom shell element by given atomic number
			 */
			const RawAtomShell& getAtomShell(UInt atomic) const;

			/**
			 * get number of shells
			 */
			UInt getNAtomShell() const { return molShells.size(); };

			///
			/// whether it has SP shell data
			///
			bool hasSPShellData() const;

			///
			/// whether the shell data has angular momentum part > D?
			///
			bool hasHighAngShellData() const;

			/**
			 * print out the shell details after read in-for debug purpose
			 */
			void print() const;

	};


}

#endif
