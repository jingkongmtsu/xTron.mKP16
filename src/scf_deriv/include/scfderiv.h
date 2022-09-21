/**
 * \file    scfderiv.h
 * \brief   class mornitoring the 1st and 2ed scf derivatives calculation
 * \author  Fenglai Liu 
 * \note
 */
#ifndef SCFDERIV_H
#define SCFDERIV_H
#include "libgen.h"

namespace shell {
	class MolShell;
}

namespace molecule {
	class Molecule;
}

namespace denmtrx {
	class DenMtrx;
}

namespace scfparam {
	class SCFParam;
}

namespace scfderiv {

	using namespace molecule;
	using namespace shell;
	using namespace denmtrx;
	using namespace scfparam;

	class SCFDeriv {

		private:

			UInt jobOrder;      ///< the deriv order
			Mtrx derivData;     ///< the deriv results

		public:

			/**
			 * constructor for SCF deriv calculation process
			 * basically it uses the molecule data to build the 
			 * result
			 */
			SCFDeriv(const Molecule& mol, UInt order);

			/**
			 * default destructor
			 */
			~SCFDeriv() { };

			/**
			 * perform SCF derivatives calculation
			 * \param infor      global information for the job
			 * \param mol        molecule data
			 * \param ms         shell data
			 * \param den        input density matrix for SCF derivative
			 *
			 * the density matrix here inside actually does not change
			 * however, for analytical integral derivatives we may need
			 * to transform it from Pure basis set order to Cartesian 
			 * (so called "P2C"), therefore it's a non-constant reference pass in
			 */
			void doDeriv(const SCFParam& param, const Molecule& mol, 
					const MolShell& ms, DenMtrx& den);

	};
}

#endif

