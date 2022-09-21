/**
 * \file   parameterparsing.h
 * \author Fenglai Liu 
 */

#ifndef PARAMETERPARSING_H
#define PARAMETERPARSING_H
#include "libgen.h"
#include<string>
#include<map>

namespace parameterparsing {

	using namespace std;

	/**
	 * In program, one class or several classes would share a parameter section
	 * defined by the user, which is usually from the input file. This class is used
	 * to parse the parameters defined in the input file for a given class.
	 *
	 * We note, that since we support multiple geometry sections, therefore each 
	 * geometry section will has its own shell data and calculation process - so,
	 * we need to know which calculation section it belongs to.
	 *
	 * on the other hand, there's global information that shared by all of sections.
	 * For example, the xcfunc class, and the parallel settings. In this case, the 
	 * section value should be set to negative value so we know that the information
	 * going to be parsed is global.
	 *
	 */
	class ParameterParsing {

		private:

			map<string,string> parameters;   ///< parameters and their values from the input

		public:

			/**
			 * constructor
			 * \param input      input file that we can find the definition of paramter
			 * \param className  which class we want to get it's parameter section?
			 * \param sec        which section the definition in? 
			 *                   0   : cluster section or global data section
			 *                   1-n : molecule section
			 * \param isClusterSec if the sec is given as 0, whether it's corresponding to
			 *                     cluster data section?
			 *
			 * For more information related to section definition, please refer to "geom"
			 * module.
			 */
			ParameterParsing(const string& input, const string& className, 
					const UInt& sec, bool isClusterSec = true);

			/**
			 * destructor
			 */
			~ParameterParsing() { };

			/**
			 * whether we have this key word defined? 
			 */
			bool hasKeyDefined(string oriKey) const;

			/**
			 * if the key is defined, we return it's value
			 * if the key is not defined, we return "NONE"
			 */
			string getValue(string oriKey) const;

			/**
			 * for the input index, let's get the corresponding key defined
			 */
			string getKey(const UInt& index) const;

			/**
			 * sometimes we may not have the entire class defined, just use it's default
			 * definition. So here we want to know whether it's defined or not
			 */
			bool hasAnyParameters() const {
				if (parameters.size() == 0) return false;
				return true;
			};

			/**
			 * return the key number of this map
			 */
			UInt getKeyNumber() const { return parameters.size(); };

			///
			/// debug printing
			///
			void print() const;
	};

}


#endif

