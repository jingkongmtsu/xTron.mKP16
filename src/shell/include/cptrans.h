/**
 * \file    cptrans.h
 *
 * This file is to do transformation between pure and Cartesian shells etc. 
 * namely, the following work is done here:
 *   - C2P transformation (cartesian to pure);
 *   - P2C transformation (pure to cartesian);
 *   - basis set normalization factor scale (for D, F etc. high L shells)
 *
 * \author  Fenglai Liu 
 */
#ifndef CPTRANS_H
#define CPTRANS_H
#include "libgen.h"
#include "matrix.h"
#include "blockmatrix.h"
#include "globalinfor.h"
#include "tbb/tbb.h"
#include "scalablevec.h"
using namespace matrix;
using namespace blockmatrix;

namespace shell {
	class AtomShell;
	class MolShell;
}

namespace cptrans {

	using namespace shell;
	using namespace globalinfor;
	using namespace tbb;
	using namespace std;

	////////////////////////////////////////////////////////////////////////
	//          @@@@ first section: constant data definition              //
	////////////////////////////////////////////////////////////////////////

	//
	// firstly, let's define the matrix data property
	// whether is C2P or P2C?
	// what kind of cartesian normalization factor it use?
	// for normalization with L00(lx=L,Ly=0,Lz=0) or 
	// for normalization with LxLyLz(lx,Ly,Lz)? 
	//
	// for scaling cartesian normalization factors,
	// do we do the scaling work? or undo the previous
	// scaling work?
	//
	const UInt  NO_TRANS     = 1;  ///< do not do transformation work (neither c2p nor p2c)
	const UInt  C2P_WITH_L00 = 2;  ///< do c2p work with cartesian basis set (lx=L,Ly=0,Lz=0)
	const UInt  C2P_WITH_XYZ = 3;  ///< do c2p work with cartesian basis set (lx,Ly,Lz)
	const UInt  P2C_WITH_L00 = 4;  ///< do p2c work with cartesian basis set (lx=L,Ly=0,Lz=0)
	const UInt  P2C_WITH_XYZ = 5;  ///< do p2c work with cartesian basis set (lx,Ly,Lz)
	const UInt  NO_SCALE     = 6;  ///< do not do basis set scaling work 
	const UInt  DO_SCALE     = 7;  ///< do basis set scaling for D, F etc.
	const UInt  UNDO_SCALE   = 8;  ///< undo the previous scaling work

	//
	// secondly, we note that the matrix could be working 
	// in two ways. in the first situation, the matrix 
	// data is directly appiled to the transformation process;
	// for example, in C2P situation the matrix with C2P_WITH_L00
	// state is directly appiled.
	//
	// In the second case, it's the transpose of matrix data 
	// is appiled in the transformation process. For example,
	// for a matrix you want to do C2P_WITH_L00, but with
	// matrix transpose; then in fact this matrix is from
	// nPurBas dimension to nCarBas dimension. This is useful
	// when in J/K calculation, see gints4d part
	//
	// Therefore, which choice you want to take? User must
	// specify this.
	//
	//
	const UInt  WITH_MATRIX_ITSELF     = 100;  ///< just use matrix itself
	const UInt  WITH_MATRIX_TRANSPOSE  = 200;  ///< use matrix transpose

	//
	// finally, the transformation/scaling work is appiled to
	// row, col or row and col? 
	//
	const UInt  CP_WITH_ROW            = 300;  ///< just work on row
	const UInt  CP_WITH_COL            = 400;  ///< just work on col
	const UInt  CP_WITH_ROW_COL        = 500;  ///< work on row and col

	//
	// this is the largest dimension that currently our cptrans
	// code could handle maximum angular momentum
	//
	const UInt  CPTRANS_MAX_SHELL      = 6;    ///< maximum angular momentum is I shell

	/**
	 * \class cptransdata
	 * \brief this class stores the data used in cp transformation as well
	 * as basis set normalization
	 */
	class CPTransData {

		private:

			UInt maxL;             ///< maximum L for setting up matrix/vector data

			//
			// data section for basis set normalization 
			//
			DoubleVec scaleVec_D;  ///< normalization vector for D shell
			DoubleVec scaleVec_F;  ///< normalization vector for F shell
			DoubleVec scaleVec_G;  ///< normalization vector for G shell
			DoubleVec scaleVec_H;  ///< normalization vector for H shell
			DoubleVec scaleVec_I;  ///< normalization vector for I shell

			//
			// data section for CP transformation
			//
			Mtrx convertMtrx_D;    ///< C2P/P2C transformation matrix for D shell
			Mtrx convertMtrx_F;    ///< C2P/P2C transformation matrix for F shell
			Mtrx convertMtrx_G;    ///< C2P/P2C transformation matrix for G shell
			Mtrx convertMtrx_H;    ///< C2P/P2C transformation matrix for H shell
			Mtrx convertMtrx_I;    ///< C2P/P2C transformation matrix for I shell

			///
			/// form the converting matrix 
			/// \param L           : angular momentum
			/// \param transWork   : transformation work assignment
			/// \param matrixStatus: with matrix itself or transpose?
			/// \return convert    : the final converting matrix
			///
			void formConvertMatrix(const UInt& L, const UInt& transWork, 
					const UInt& matrixStatus, Mtrx& convert);

			///
			/// form the converting vector 
			/// \param L           : angular momentum
			/// \param scaleWork   : basis set scaling work assignment
			/// \return convert    : the final converting vector
			///
			void formConvertVector(const UInt& L, const UInt& scaleWork, 
					DoubleVec& convert);

		public:

			///
			/// constructor 
			/// \param  maxL: maximum L 
			/// \param  scaleWork: DOSCALE, NOSCALE or UNDO_SCALE
			/// \param  transWork: see the macro above
			/// \param  mtrxstatus: with matrix itlsef or it's transpose?
			///
			CPTransData(const UInt& maxL, const UInt& transWork, const UInt& scaleWork, 
					const UInt& mtrxStatus);

			///
			/// simpler contructor, build the content later
			///
			CPTransData(const UInt& maxL0):maxL(maxL0) { };

			///
			/// destructor
			///
			~CPTransData() { };

			///
			/// initialize the scaling data, for DO_SCALE work
			///
			void initScaleData();

			/**
			 * return the  C2P/P2C matrix
			 */
			const Mtrx& getConvertMatrix(const UInt& L) const; 

			/**
			 * get the normalization scale vector the given L 
			 */
			const DoubleVec& getConvertVec(const UInt& L) const;
	};

	/**
	 * \class c2ptransatomshell
	 * \brief this is the class to perform c2p/p2c transformation work
	 * based on atom-shell-pair/atom-shell data
	 *	
	 *	for the source type and target type, we note that 
	 * it's able to derive the source matrix and target matrix status
	 * from the transWork. The value of source/target type
	 * must be in TYPE_CART, TYPE_NORM or TYPE_PURE
	 * see shell.h
	 */
	class CPTransAtomShell {

		private:

			UInt sourceType;  ///< source matrix type, cart or norm?
			UInt targetType;  ///< targrt matrix type, cart or norm?
			CPTransData c2p;  ///< keep a copy of CP transformation work data

		public:

			///
			/// constructor 
			/// since we do not perform basis set scaling work here, therefore
			/// wew do not have scaleWork as parameter here
			/// \param  maxL: maximum L 
			/// \param  transWork: C2P or P2C work assignment?
			/// \param  mtrxstatus: with matrix itlsef or it's transpose?
			///
			CPTransAtomShell(const UInt& maxL, const UInt& transWork, const UInt& mtrxStatus);

			///
			/// destructor
			///
			~CPTransAtomShell() { };

			///
			/// perform C2P/P2C work based on atom shell data
			/// this must be combined with the work either on row, or on col
			/// for the transformation work on both row and col, it will
			/// be performed on shell pairs data
			/// \param work: CP_WITH_ROW or CP_WITH_COL?
			/// \param s   : shell data
			/// \param S   : input data matrix, source matrix
			/// \return T  : output dataa matrix, target matrix
			///
			void cpTransformOnAtomShell(const UInt& work, 
					const AtomShell& s, const Mtrx& S, Mtrx& T) const;

			///
			/// perform C2P/P2C work based on atom shell pair data
			/// both of the matrix row and col will be transformed
			///
			/// alert!!! Here we note that the input data "S" matrix should be 
			/// "full" even if the atom shell rs and cs are same, we do not
			/// check that whether this is some symmetrical matrix and only 
			/// use half of the data. The user should be caution about this!!!
			///
			/// \param rs,cs    : row and col shell data
			/// \param S        : input data matrix, source matrix
			/// \param tmp      : input scratch matrix, used in BLAS3 operation
			/// \return T       : output dataa matrix, target matrix
			///
			///
			void cpTransformOnAtomShellPair(const AtomShell& rs, 
					const AtomShell& cs, const Mtrx& S, Mtrx& tmp, Mtrx& T) const; 

			///
			/// return the source type
			///
			UInt getSourceType() const { return sourceType; };

			///
			/// return the target type
			///
			UInt getTargetType() const { return targetType; };
	};

	/**
	 * \class normatomshelldata
	 * \brief this is the class to perform basis set normalization work
	 * for the row/column data based on atom shell dimension
	 */
	class NormAtomShellData {

		private:

			DoubleVec   rowVec; ///< for transformation on row, we need a scratch vector
			CPTransData scale;  ///< keep a copy of scale transformation work data

		public:

			///
			/// constructor 
			/// \param  maxL     :  set up the scale work
			/// \param  maxNBasis:  this is used to set up the scratch data
			/// \param  scaleWork:  what kind of scale work it will do?
			///
			NormAtomShellData(const UInt& maxL, const UInt maxNBasis, 
					const UInt& scaleWork):rowVec(maxNBasis),
			scale(maxL,NO_TRANS,scaleWork,WITH_MATRIX_ITSELF) { };

			///
			/// destructor
			///
			~NormAtomShellData() { };

			///
			/// normalize the row side data for matrix T with use of atom shell data s
			/// since we need to use the scratch vector, this is not a constant function
			///
			void normRowData(const AtomShell& s, Mtrx& T); 

			///
			/// normalize the col side data for matrix T with use of atom shell data s
			///
			void normColData(const AtomShell& s, Mtrx& T) const; 
	};

	///
	/// this is the functor to perform C2P/P2C transformation
	/// for matrix type of data
	///
	class TBB_CPTransMatrix {

		private:

			UInt rowColStatus;            ///< whetherr the work is on row, col or both of them?
			CPTransAtomShell  cptrans;    ///< perform C2P/P2C transformation work for atom shell pair etc.
			const MolShell&   rs;         ///< row shell for the mtrx 
			const MolShell&   cs;         ///< col shell for the mtrx 
			const Mtrx&       M0;         ///< input matrix
			Mtrx&             M1;         ///< output matrix

			///
			/// this is the work routine to perform transformation on both row and col;
			///
			void doCPTransOnRowAndCol(const blocked_range<UInt>& r) const; 

			///
			/// this is the work routine to perform transformation on either row or col;
			///
			void doCPTransOnRowOrCol(const blocked_range<UInt>& r) const;

		public:

			///
			/// constructor
			/// this is work on both of the shells for row and col
			///
			TBB_CPTransMatrix(const MolShell& rowShell, const MolShell& colShell,
					const UInt& transWork, const UInt& matrixStatus,
					const Mtrx& source, Mtrx& target);

			///
			/// constructor
			/// this is work on either row or col
			///
			TBB_CPTransMatrix(const MolShell& s, const UInt& transWork, 
					const UInt& matrixStatus, const UInt& rowColStatus0,
					const Mtrx& source, Mtrx& target);

			///
			/// default destructor
			///
			~TBB_CPTransMatrix() { };

			///
			/// operator to do the parallel work
			///
			void operator()(const blocked_range<UInt>& r) const {
				if (rowColStatus == CP_WITH_ROW_COL) {
					doCPTransOnRowAndCol(r); 
				}else{
					doCPTransOnRowOrCol(r); 
				}
			}; 
	};

	///
	/// this is the functor to perform basis set scale 
	/// the scale is for the whole matrix
	///
	class TBB_NormBasisMatrix {

		private:

			//
			// input data and output data source
			//
			bool doRowWork;                 ///< do we do row transformation?
			bool doColWork;                 ///< do we do row transformation?
			const DoubleVec&   rowNormVec;  ///< normalization vector for row shell
			const DoubleVec&   colNormVec;  ///< normalization vector for col shell
			Mtrx&              target;      ///< the reference of target matrix   

		public:

			///
			/// constructor
			/// just copy the input reference 
			///
			TBB_NormBasisMatrix(bool doRow, bool doCol, const DoubleVec& rowVec, 
					const DoubleVec& colVec, Mtrx& target0):doRowWork(doRow),
			doColWork(doCol),rowNormVec(rowVec),colNormVec(colVec),target(target0) { };

			///
			/// default destructor
			///
			~TBB_NormBasisMatrix() { };

			///
			/// operator to do the parallel work
			///
			void operator()(const blocked_range<UInt>& r) const; 
	};

	///
	/// this is the top interface class to perform CP transformation
	/// and/or basis set normalization work
	///
	class CPTransBasisNorm {

		private:

			//
			// the global infor
			//
			GlobalInfor infor;  ///< information 

			//
			// work assignment
			//
			UInt  transWork;    ///< trans work status (use the above macro)
			UInt  scaleWork;    ///< scale work status (use the above macro)
			UInt  matrixStatus; ///< matrix status, with itself or with transpose?
			UInt  rowColStatus; ///< whether the work is on row, col or both column and row?

			///
			/// internal work function to perform basis set normation for 
			/// the matrix type data
			///
			void normBasis(const MolShell& rs, const MolShell& cs, Mtrx& M) const;

			///
			/// this is the working function to perform C2P/P2C transformation
			/// for matrix type data
			///
			void CPTransform(const MolShell& rs, const MolShell& cs, Mtrx& M) const;

		public:

			///
			/// constructor 
			/// \param  infor0   : global information to control use of parallel
			/// \param  scaleWork: DOSCALE, NOSCALE or UNDO_SCALE
			/// \param  transWork: see the macro above
			/// \param  mtrxstatus: with matrix itlsef or it's transpose?
			/// we will check the parameters internally
			///
			CPTransBasisNorm(const GlobalInfor& infor0, const UInt& transWork0, 
					const UInt& scaleWork0, const UInt& mtrxStatus0, 
					const UInt& rowColStatus0);

			///
			/// destructor
			///
			~CPTransBasisNorm() { };

			///
			/// top interface for doing C2P/P2C conversion job and/or basis set 
			/// normalization
			/// if the conversion work is carried out for either row/col,
			/// then only the rs/cs shell data is really used.
			/// \param  rs, cs:  input row shell and column shell for M
			/// \param  M     :  original data matrix going to be converted
			///
			void transform(const MolShell& rs, const MolShell& cs, Mtrx& M) const{

				// do C2P/P2C transformation
				CPTransform(rs,cs,M);

				// do normalization work
				normBasis(rs,cs,M);
			};

	};

}
#endif

