/**
 * \file    fock.h
 * \brief   describe the Fock matrix, how to build it and how to handle it in SCF
 * \author  Fenglai Liu 
 */
#ifndef FOCK_H
#define FOCK_H
#include "libgen.h"
#include "spinmatrix.h"
using namespace spinmatrix;

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

namespace oneemtrx {
	class OneEMtrx;
}

namespace scfintscontroller {
	class SCFIntsController;
}

namespace fock {

	using namespace shell;
	using namespace molecule;
	using namespace denmtrx;
	using namespace scfparam;
	using namespace oneemtrx;
	using namespace scfintscontroller; 

	/**
	 * this class is forming and updating Fock matrix in SCF process
	 */
	class Fock : public SpinMatrix {

		private:

			UInt nRow;                  ///< keep a copy of Fock's row dimension
			UInt nCol;                  ///< keep a copy of Fock's col dimension

		public:

			/**
			 * update the Fock matrix with J/K matrix part
			 * return the JK energy
			 */
			Double addJKMtrx(const SCFParam& param, const SCFIntsController& intController,
					const Molecule& mol, const MolShell& ms, DenMtrx& denMtrx);

			/**
			 * update the Fock matrix with XC matrix part
			 * return the XC energy
			 */
			Double addXCMtrx(const SCFParam& param, const SCFIntsController& intController,
					const Molecule& mol, const MolShell& ms, const DenMtrx& denMtrx);

			/**
			 * update the Fock matrix with core Hamiltonian part
			 * return the core energy
			 */
			Double addCoreMtrx(const SCFParam& param, const OneEMtrx& oneEMtrx, 
					const MolShell& ms, const DenMtrx& denMtrx);

			/**
			 * get the nRow
			 */
			UInt getRow() const { return nRow; };

			/**
			 * get the nCol
			 */
			UInt getCol() const { return nCol; };

			/**
			 * constructor for prepare the Fock matrix
			 * \param  param     SCF parameters
			 * \param  ms        shell data
			 */
			Fock(const SCFParam& param, const MolShell& ms);

			/**
			 * destructor
			 */
			~Fock() { };

	};
}

#endif

