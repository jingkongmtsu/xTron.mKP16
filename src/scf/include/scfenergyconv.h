/**
 * \file    scfenergyconv.h
 * \brief   this is a coarse method for converging SCF, only judging from the energy
 * \author  Fenglai Liu 
 * \note
 *
 * This is a simple method to judge the SCF convergence. We only counting on the energy
 * changes from the previous iterations. If the energy changes is less than the threshold
 * value for the given iterations, we will consider that the SCF is converged.
 *
 * This method should be the last method we use for converging the SCF.
 */
#ifndef SCFENERGYCONV_H
#define SCFENERGYCONV_H
#include "libgen.h"

namespace scfconv {
	class SCFConv;
}

namespace molecule {
	class Molecule;
}

namespace globalinfor {
	class GlobalInfor;
}

namespace scfenergyconv {

	using namespace globalinfor;
	using namespace molecule;
	using namespace scfconv;

	/**
	 * this class is used to obtain the SCF DIIS result
	 */
	class SCFEnergyConvController {

		private:

			UInt nIterConvTest;     ///< number of scf iterations that SCF is converged while energy diff is below thresh
			Double convThresh;      ///< the threshold value for comparing energy difference for SCF convergence test

		public:

			/**
			 * constructor - we just simply read in the parameters
			 */
			SCFEnergyConvController(const GlobalInfor& infor, const Molecule& mol);

			/**
			 * destructor
			 */
			~SCFEnergyConvController() { };

			/**
			 * we need the previous energy result, and here we will judge that whether it's converged?
			 */
			bool hasEnergyConv(const SCFConv& conv) const;

	};
}

#endif

