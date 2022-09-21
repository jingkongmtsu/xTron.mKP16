/**
 * \file    virialparam.h
 * \brief   parsing parameters for virial coefficients calculation
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef VIRIALPARAM_H
#define VIRIALPARAM_H
#include "libgen.h"

namespace virialparam {

	//
	// We have several ways to perform virial coefficients
	// calculation. 
	//
	const Int PERTURBED_SCF_WAY_VIRIAL = 0;
	const Int FULL_SCF_WAY_VIRIAL      = 1;

	///
	/// parsing the parameter for virial coefficients calculation
	/// all of the input data should be attached to cluster section
	/// therefore, you should have cluster section of data defined
	/// in the given input file
	///
	class VirialParam {

		private:

			Int nMonomers;       ///< number of monomers
			Int method;          ///< virial coefficients calculation choice

			///
			/// from the given input file, find how many monombers 
			/// (which is included in the molecule section) for 
			/// Virial calculation
			///
			Int findNumMonomers(const string& input) const;

		public:

			/**
			 * constructor 
			 */
			VirialParam(const string& input);

			/**
			 * destructor
			 */
			~VirialParam() { };

			/**
			 * return the number of monomers
			 */
			Int getNMono() const { return nMonomers; };

			/**
			 * return the virial coefficient calculation method
			 */
			Int getMethod() const { return method; };

	};

}

#endif
