/**
 * \file    xdm.h
 * \brief   class used to calculate VDW by XDM model
 * \author  Fenglai Liu and Jing Kong
 * \note
 *
 * XDM references:
 * 
 * E. Johnson and A. Becke, J. Chem. Phys. 123, 024101 (2005).
 *
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 122, 154104 (2005).
 *
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 123, 154101 (2005).
 * 
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 124, 014104 (2006).
 * 
 * E. R. Johnson and A. D. Becke, J. Chem. Phys. 124, 174104 (2006).
 * 
 * E. R. Johnson and A. D. Becke, Chem. Phys. Lett. 432, 600 (2006).
 * 
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 124108 (2007).
 * 
 * A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).
 * 
 * F. O. Kannemann and A. D. Becke, J. Chem. Theory Comput. 5, 719 (2009).
 * 
 * F. O. Kannemann and A. D. Becke, J. Chem. Theory Comput. 6, 1081 (2010).
 * 
 * A. D. Becke, A. Arabi, and F. O. Kannemann, Can. J. Chem. 88, 1057 (2010).
 * 
 * E. Johnson, J. Chem. Phys. 135, 234109 (2011).
 * 
 * A. Arabi and A. Becke, J. Chem. Phys. 137, 014104 (2012).
 * 
 * A. Otero-de-la Roza and E. R. Johnson, J. Chem. Phys. 138, 054103 (2013).
 *
 */
#ifndef XDM_H
#define XDM_H
#include "tbb/tbb.h"
#include "libgen.h"
#include "scalablevec.h"
#include "xcintsinfor.h"
#include "xcfunc.h"
#include "xcvar.h"
#include "molegrids.h"
#include "shell.h"
#include "shellsize.h"
#include "matrix.h"
using namespace xcintsinfor;
using namespace xcfunc;
using namespace xcvar;
using namespace molegrids;
using namespace shell;
using namespace shellsize;
using namespace matrix;

namespace vdwinfor {
	class VDWInfor;
}

namespace denmtrx {
	class DenMtrx;
}

namespace hirshfeld {
	class Hirshfeld;
}

namespace atomdenmtrx {
	class AtomDenMtrx;
}

namespace sigatombasis {
	class SigAtomBasis;
}

namespace batchbasis {
	class BatchBasis;
}

namespace batchgrid {
	class BatchGrid;
}

namespace batchvar {
	class BatchVar;
}

namespace batchfunc {
	class BatchFunc;
}

namespace xdm {

	using namespace vdwinfor;
	using namespace denmtrx;
	using namespace hirshfeld;
	using namespace atomdenmtrx; 
	using namespace sigatombasis;
	using namespace batchbasis;
	using namespace batchgrid;
	using namespace batchvar;
	using namespace batchfunc;


	/**
	 * \class TBB_XDM
	 *
	 * this is the working class to calculate the atom volumn and 
	 * square multiple. The results are calculated in DFT way
	 *
	 * the implementation of the class is based on the following reference:
	 * E. R. Johnson and A. D. Becke, J. Chem. Phys. 124, 174104 (2006).
	 */
	class TBB_XDM {

		private:

			//
			// reference for the input 
			//
			UInt nAtoms;                        ///< number of atoms for the whole molecule
			const MolShell& ms;                 ///< input shell
			const MoleGrids& molGrids;          ///< molecular grids
			const MolShellSize& molShellSize;   ///< molecular shell size data
			const XCFunc& xcfunc;               ///< xc functional information
			const XCVar& xcvar;                 ///< xc variable information
			const XCIntJobInfor& infor;         ///< information center for real calculation
			const DenMtrx& denMtrx;             ///< density matrix
			const AtomDenMtrx& atomDenMatrices; ///< density matrices data for free atom

			//
			// the output data - local copy of results
			// we note that the cols in squareMultiple is always 3 in this scheme
			//
			Mtrx squareMultiple;                ///< this is the term of \f$ M^{2}_{l}\f$ defined in equation 9/10
			DoubleVec effectiveAtomVol;         ///< this term is defined in equation 17
			DoubleVec freeAtomVol;              ///< this term is defined in equation 18

			///
			/// update the free atom vol in batch
			///
			void updateFreeAtomVol(const SigAtomBasis& sigList, const Hirshfeld& hirshfeld);

			///
			/// compute the square multiple and atom vol in batch
			///
			void formSquareMultipleAndAtomVol(const MolShell& ms, const SigAtomBasis& sigList, 
					const BatchGrid& bg, const BatchVar& bvar, const BatchFunc& bfunc, 
					const Hirshfeld& hirfeld);

		public:

			/**
			 * constructor 
			 *
			 * note that the columns in squareMultiple is always 3
			 * this is because we need C6, C8, C10 always
			 */
			TBB_XDM(const UInt& nAtoms0,const MolShell& ms0,const MoleGrids& molGrids0, 
					const MolShellSize& molShellSize0, 
					const XCFunc& xcfunc0, const XCVar& xcvar0, 
					const XCIntJobInfor& infor0, const DenMtrx& den, 
					const AtomDenMtrx& atomDenMatrices0):nAtoms(nAtoms0),
			ms(ms0),molGrids(molGrids0),molShellSize(molShellSize0),xcfunc(xcfunc0),xcvar(xcvar0),
			infor(infor0),denMtrx(den),atomDenMatrices(atomDenMatrices0),
			squareMultiple(nAtoms,3),effectiveAtomVol(nAtoms,ZERO),freeAtomVol(nAtoms,ZERO) { }

			/**
			 * destructor
			 */
			~TBB_XDM() { };

			/**
			 * functional operator to generate the results
			 * this is used as parallel_reduce
			 */
			void operator()(const blocked_range<UInt>& r); 

			/**
			 * splitting constructor used by TBB
			 * we copy all of reference data except result
			 * which should be always empty initially
			 */
			TBB_XDM(const TBB_XDM& tbb_xdm, split):nAtoms(tbb_xdm.nAtoms),
			ms(tbb_xdm.ms),molGrids(tbb_xdm.molGrids),
			molShellSize(tbb_xdm.molShellSize),xcfunc(tbb_xdm.xcfunc),xcvar(tbb_xdm.xcvar),
			infor(tbb_xdm.infor),denMtrx(tbb_xdm.denMtrx),atomDenMatrices(tbb_xdm.atomDenMatrices),
			squareMultiple(nAtoms,tbb_xdm.squareMultiple.getCol()),
			effectiveAtomVol(nAtoms,ZERO),freeAtomVol(nAtoms,ZERO) { }

			/**
			 * for assemble results in threads together
			 */
			void join(TBB_XDM& tbb_xdm) {  
				squareMultiple.add(tbb_xdm.squareMultiple);
				for(UInt i=0; i<nAtoms; i++) {
					effectiveAtomVol[i] += tbb_xdm.effectiveAtomVol[i];
					freeAtomVol[i] += tbb_xdm.freeAtomVol[i];
				}
			};

			/**
			 * return the square multiple
			 */
			const Mtrx& getSquareMultiple() const { return squareMultiple; };

			/**
			 * return the effective atom vol
			 */
			const DoubleVec& getEffectiveAtomVol() const { return effectiveAtomVol; };

			/**
			 * return free atom vol
			 */
			const DoubleVec& getFreeAtomVol() const { return freeAtomVol;};
	};

	/**
	 * top interface class for doing XDM work
	 *
	 * right now we only have energy calculation
	 * so we do not have jobOrder, it's set to 0
	 */
	class XDM {

		protected:

			//
			// class members for perform xcints job
			//
			XCIntJobInfor infor;    ///< information of for doing xcints job
			XCFunc      xcfunc;     ///< xcfunctional information
			XCVar       xcvar;      ///< xc variable information
			MolShellSize size;      ///< molecular shell size data
			MoleGrids   molGrids;   ///< atom and molecule grids setting

			//
			// member data for XDM itself
			//
			Double E6;              ///< energy for R6 term
			Double E8;              ///< energy for R8 term
			Double E10;             ///< energy for R10 term
			Double nonAdditiveE3;   ///< non-additive three body energy term
			Double E;               ///< the total energy for VDW

			///
			/// making C6, C8 C10 parameters for the atom i and atom j pair
			///
			void makeCijForIJAtomPair(const UInt& iAtom, const UInt& jAtom, const Mtrx& squareMultiple,
					const DoubleVec& effectiveAtomPol, Double& C6ij, Double& C8ij, Double& C10ij) const;

			///
			/// making Becke damping factors in terms of C6, C8 etc.
			///
			Double makeBeckeDampingFac(const VDWInfor& vdwInfor, const Double& C6ij, 
					const Double& C8ij, const Double& C10ij) const;

			///
			/// driver function to perform XDM energy calculation
			///
			void xdmEnergy(const MolShell& ms, const VDWInfor& xdmInfor, 
					const Mtrx& squareMultiple, const DoubleVec& effectiveAtomVol, 
					const DoubleVec& freeAtomVol); 

		public:

			/**
			 * construtor for the XDM - only for information processing
			 * \param ms        shell data
			 * \param xcinfor   a copy of xc information center
			 */
			XDM(const VDWInfor& vdwInfor, const MolShell& ms, const XCIntJobInfor& xcinfor);

			/**
			 * destructor
			 */
			~XDM() { };

			/**
			 * top interface to form XDM energy calculation
			 * \param   molShell  shell data
			 * \param   denMtrx   density matrix
			 */
			void energy(const VDWInfor& xdminfor, const MolShell& ms, const Molecule& mol,
					const DenMtrx& denMtrx, bool printTiming);

			/**
			 * return the E
			 */
			Double getEnergy() const { 
				return E;
			};
	};

}

#endif

