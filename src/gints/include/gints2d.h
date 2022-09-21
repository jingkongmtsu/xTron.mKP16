/**
 * \file   gints2d.h
 * \brief  working class to handle the 2D analytical integrals calculation
 * \author Fenglai Liu 
 */
#ifndef GINTS2D_H
#define GINTS2D_H
#include "tbb/tbb.h"
#include "libgen.h"
#include "scalablevec.h"
#include "gintsinfor.h"
#include "sigshellpairinfor.h"
#include "matrix.h"

namespace molecule {
	class Molecule;
};

namespace shell {
	class MolShell;
};

namespace gints2d {

	using namespace tbb;
	using namespace shell;
	using namespace matrix;
	using namespace molecule;
	using namespace gintsinfor;
	using namespace sigshellpairinfor;

	/**
	 * \class TBB_GInts2DIntMatrix
	 *
	 * This is a working class to do the normal 2D analytical integral calculation. 
	 * The output result is 2D AO integral matrix.
	 *
	 * This working class supports the following GInts job:
	 * - two body overlap;
	 * - two body kinetic;
	 * - two body nuclear-electron attraction;
	 * - core Hamiltonian (KI+NAI)
	 */
	class TBB_GInts2DIntMatrix {

		private:

			//
			// reference for the input and output
			//
			const Molecule&     mol;              ///< geometry information
			const MolShell&     rs;               ///< row shell information
			const MolShell&     cs;               ///< column shell information
			const SigMolShellPairInfor& infor;    ///< shell pair information
			const GIntJobInfor& ginfor;           ///< information center for real calculation
			Mtrx& intMtrx;                        ///< output result 

			//
			// geometry data for using in NAI etc.
			//
			DoubleVec coord;                      ///< xyz coordiates derived from molecule 
			UIntVec Z;                            ///< atomic data 

		public:

			/**
			 * constructor 
			 */
			TBB_GInts2DIntMatrix(const GIntJobInfor& ginfor0, const SigMolShellPairInfor& shellPairInfor, 
					const MolShell& rowShell, const MolShell& colShell, const Molecule& mol0, Mtrx& result);

			/**
			 * destructor
			 */
			~TBB_GInts2DIntMatrix() { };

			/**
			 * functional operator to generate the integral results
			 * this is used as parallel_for
			 */
			void operator()(const blocked_range<UInt>& r) const; 
	};

	/**
	 * \class TBB_GInts2DMultiIntMatrix
	 *
	 * This is a working class to do the normal 2D analytical 
	 * integral calculation in terms of multiple integral
	 * matrices. For example, the momentum integrals such as
	 * (1|x,y,z|2) outputs the momentum integral matrix for 
	 * x, y and z.
	 */
	class TBB_GInts2DMultiIntMatrix {

		private:

			//
			// reference for the input and output
			//
			const MolShell&     rs;               ///< row shell information
			const MolShell&     cs;               ///< column shell information
			const SigMolShellPairInfor& infor;    ///< shell pair information
			const GIntJobInfor& ginfor;           ///< information center for real calculation
			Double C[3];                          ///< the possible new center for momentum integrals
			MtrxVec& intMtrxVec;                  ///< output result matrix list

		public:

			/**
			 * constructor - for this constructor the center of C is all set to zero
			 */
			TBB_GInts2DMultiIntMatrix(const GIntJobInfor& ginfor0, const SigMolShellPairInfor& shellPairInfor, 
					const MolShell& rowShell, const MolShell& colShell, MtrxVec& result);

			/**
			 * constructor - for this constructor the center of C is set
			 */
			TBB_GInts2DMultiIntMatrix(const GIntJobInfor& ginfor0, const SigMolShellPairInfor& shellPairInfor, 
					const MolShell& rowShell, const MolShell& colShell, const Double* angCenter, MtrxVec& result);

			/**
			 * destructor
			 */
			~TBB_GInts2DMultiIntMatrix() { };

			/**
			 * functional operator to generate the integral results
			 * this is used as parallel_for
			 */
			void operator()(const blocked_range<UInt>& r) const; 
	};

	/**
	 * \class GInts2D
	 * \brief top level class to mornitor 2D analytical integral calculation
	 */
	class GInts2D {

		private:

			GIntJobInfor ginfor;           ///< information center for real calculation
			SigMolShellPairInfor infor;    ///< shell pair information

		public:

			////////////////////////////////////////////////////////////
			//                  constructor etc.                      //
			////////////////////////////////////////////////////////////

			/**
			 * constructor - build the int job inside
			 */
			GInts2D(const MolShell& rs, const MolShell& cs, const GIntsInfor& infor0, const UInt& job, 
					const UInt& order):ginfor(infor0,job,order),infor(rs,cs,ginfor) { };

			/**
			 * constructor - directly pass in the int job 
			 */
			GInts2D(const MolShell& rs, const MolShell& cs,
					const GIntJobInfor& infor0):ginfor(infor0),infor(rs,cs,ginfor) { };

			/**
			 * destructor
			 */
			~GInts2D() { };

			///
			/// form one electron integral matrix based on two given shells
			///
			/// \param  rS           input shell data corresponding to row of matrix
			/// \param  cS           input shell data corresponding to col of matrix   
			/// \param  mol          geometry data
			/// \return intMtrx      the result integral matrix in matrix form
			///
			void doMtrx(const MolShell& rs, const MolShell& cs, const Molecule& mol, 
					Mtrx& intMtrx, bool printTiming) const;

			///
			/// form the integral matrix on two given shells
			/// the input is the multiple matrix result, matrix vector
			///
			void doMultiMtrx(const MolShell& rs, const MolShell& cs,
					MtrxVec& intMtrxList, bool printTiming) const;
	};
}

#endif

