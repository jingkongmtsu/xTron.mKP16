/**
 * \file    oneemtrx.h
 * \brief   class for creating/fetching the one electron matrix
 * \author  Fenglai Liu
 */
#ifndef ONEEMTRX_H
#define ONEEMTRX_H
#include "libgen.h"
#include "histdataman.h"

namespace matrix {
	class Mtrx;
}

namespace shell {
	class MolShell;
}

namespace molecule {
	class Molecule;
}

namespace globalinfor {
	class GlobalInfor;
}

namespace gintsinfor {
	class GIntsInfor;
}

namespace oneemtrx {

	using namespace histdataman;
	using namespace matrix;
	using namespace molecule;
	using namespace shell;
	using namespace gintsinfor;
	using namespace globalinfor;

	/**
	 * \brief one electron matrix
	 *
	 * this type of matrix has a feature that the integrals
	 * only has two dimension, therefore it's able to pre-computed
	 * on the other hand, they does not have spin character
	 * and they share the same dimension
	 *
	 * this class will compute these matrices and store them for future use
	 */
	class OneEMtrx : public HistDataMan {

		private:

			UIntVec  jobList;   ///< matrix list for one electron matrix

		public:

			/**
			 * constructor to build the one electron matrix
			 * \param infor       global job information
			 * \param section     get from molecule
			 * \param matrixList  integral matrix name list for computing here
			 * \param useFile     do we use file to store data?
			 */
			OneEMtrx(const GlobalInfor& infor, const UInt& section, 
					const UIntVec& matrixList, bool useFile);

			/**
			 * constructor to build the one electron matrix with default job assignment
			 * This is particularly used for SCF part
			 * \param infor       global job information
			 * \param section     get from molecule
			 * \param useFile     do we use file to store data?
			 */
			OneEMtrx(const GlobalInfor& infor, const UInt& section, bool useFile);

			/**
			 * destructor
			 */
			~OneEMtrx() { };

			/**
			 * form the data matrix 
			 */
			void formMtrx(const GIntsInfor& intInfor, const Molecule& mol, 
					const MolShell& rs, const MolShell& cs, bool printTiming, bool printMatrix);

			/**
			 * get the data matrix
			 *
			 * \param name:  the name of integral matrix
			 * \param intoFullMatrix: do we make the matrix into full matrix(upper and lower)?
			 * \return M:    the result data matrix
			 */
			void getM(const UInt& name, Mtrx& M, bool intoFullMatrix = false) const;

			/**
			 * get the dimension for the given matrix
			 */
			void getDim(const UInt& name, UInt& nRow, UInt& nCol) const;
	};

}

#endif
