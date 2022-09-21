/**
 * \file    virialcoe.h
 * \brief   monitoring virial coefficient calculation
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef VIRIALCOE_H
#define VIRIALCOE_H
#include "libgen.h"
#include "mo.h"
#include "virialparam.h"

namespace globalinfor {
	class GlobalInfor;
}

namespace virialcoe {

	using namespace mo;
	using namespace globalinfor;
	using namespace virialparam;

	/**
	 * this class mornitors the Virial coefficients calculation
	 * Basically we have two ways to perform the compuation:
	 * -  Purturbed way,
	 * -  Viariational way
	 *
	 *
	 */
	class VirialCoe {

		private:

			Double eMonomerSum;  ///< sum of e for monomers
			Double eCluster;     ///< e for cluster
			VirialParam param;   ///< virial coefficients parameter 
			vector<MOs> moslist; ///< the unique molecular orbital data for monomers

			///////////////////////////////////////
			//          member functions         //
			///////////////////////////////////////

			/**
			 * calculate eMonomerSum
			 */
			void monomerCalculation(const GlobalInfor& infor);

			/**
			 * calculate the cluster energy
			 */
			void clusterCalculation(const GlobalInfor& infor);

		public:

			/**
			 * constructor to build the energy
			 */
			VirialCoe(const string& input);

			/**
			 * destructor
			 */
			~VirialCoe() { };

			/**
			 * return the energy of monomer sum
			 */
			Double getDeltaE() const { return (eCluster-eMonomerSum)*27.212*23.061; };

	};

}

#endif
