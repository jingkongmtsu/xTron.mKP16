/**
 * \file   espints.h
 * \brief  working class to handle the integrals calculation based on ESP integrals
 * \author Fenglai Liu 
 */
#ifndef ESPINTS_H
#define ESPINTS_H
#include "libgen.h"
#include "scalablevec.h"
#include "gintsinfor.h"

namespace shell {
	class MolShell;
};

namespace spinmatrix {
	class SpinMatrix;
};

namespace sigshellpairinfor {
	class SigMolShellPairInfor;
};

namespace xcintsinfor {
	class XCIntJobInfor;
};

namespace espints {

	using namespace shell;
	using namespace spinmatrix;
	using namespace gintsinfor;
	using namespace xcintsinfor; 
	using namespace sigshellpairinfor;

	/**
	 * \class  ESPInts
	 * \brief top level class to mornitor numerical integral calculation based on esp integrals
	 */
	class ESPInts {

		private:

			GIntJobInfor ginfor;           ///< information center for real calculation

			///
			/// this is the working function to do ESP integral work and digest 
			/// the raw integrals into result of halfExRho, halfJRhoVec and halfJRhoVecMtrx
			///
			/// the half exrho is the half of exchange energy density (Eq.6 in the ex_density
			/// paper), and is defined as:
			///
			/// \f$ ex_{\nu}(r)  = \sum_{i}^{occ} \phi_{i}(r)
			/// \int \frac{\phi_{i}(r^{'})\phi_{\nu}(r^{'})}{|r-r^{'}|} dr^{'}\f$
			///
			/// i is the occupied orbital, \f$ nu \f$ is the basis set index. This is done
			/// on each given grid point of r. \f$ \phi(r) \f$ is the basis set value
			/// on the given grid point. It's needed for forming the exchange energy density,
			/// as well as forming the exchange matrix.
			///
			/// the denPhi, which is the combination between the density matrix
			/// and \f$ \phi(r) \f$ value; is defined as:
			/// \f$ p_{nu}(r) = \sum_{\lambda\nu}P_{\nu\lambda}\phi_{\lambda}(r)\f$
			///
			/// P is the density matrix, and \f$\phi_{\lambda}(r)\f$ is the basis set
			/// value on grid point.
			///
			/// for computing Coulomb matrix, we need two parts:
			///
			/// \f$ J_{\mu\nu}_{1}(r) = \phi_{\mu}(r)\phi_{\nu}(r)\sum_{\eta\lambda}
			///     P_{\eta\lambda}\int dr^{'}\frac{\phi_{\eta}(r^{'})\phi_{\lambda}(r^{'})}{|r-r^{'}|} \f$
			/// 
			/// For this part, we compute the 
			///
			/// \f$ ec(r) = \sum_{\eta\lambda}P_{\eta\lambda}
			/// \int dr^{'} \frac{\phi_{\eta}(r^{'})\phi_{\lambda}(r^{'})}{|r-r^{'}|} \f$
			///
			/// this is the halfJRhoVec, since it's a vector on the grid points.
			///
			/// Another part for calculating Coulomb matrix is:
			///
			/// \f$ J_{\mu\nu}_{2}(r) = (\sum_{\eta\lambda}P_{\eta\lambda}\phi_{\eta}(r)\phi_{\lambda}(r))
			///     \int dr^{'}\frac{\phi_{\mu}(r^{'})\phi_{\nu}(r^{'})}{|r-r^{'}|} \f$
			/// 
			/// it can also be written as:
			/// 
			/// \f$ J_{\mu\nu}_{2}(r) = \rho(r)
			///     \int dr^{'}\frac{\phi_{\mu}(r^{'})\phi_{\nu}(r^{'})}{|r-r^{'}|} \f$
			///
			/// this part is just the halfJRhoVecMtrx. By combining the \f$ rho(r) \f$ with the 
			/// weights, this part directly form the contribution to the J matrix.
			///
			/// also for the cartDen, both of alpha and beta has been combined together so that
			/// to form the spin-resolved density
			///
			/// we note that both of the input cartDen, and the output halfJRhoVecMtrx; are in
			/// shell pair vector form. For more details about the form of this vector, please
			/// see the sigshellpairinfor.h function of convertMatrixIntoSigSPVec.
			///
			/// \param grids       : the coordinates of the grid points
			/// \param ms          : the shell information
			/// \param cartDen     : the shell pair vector form density matrix in Cartesian order
			/// \param denPhi      : combination of the density matrix and phi value (nCarBas,nGrids)
			/// \param infor       : the shell pair information
			/// \param wrho        : the density on this batch grid together with weights
			/// \return halfExRho  : half of exchange energy density
			/// \return halfJRhoVec: half of the Coulomb energy density vector in J1 part
			/// \return halfJRhoVecMtrx: half of the Coulomb energy density vector form matrix in J2 part
			///
			void doESPIntWork(const DoubleVec& grids, const MolShell& ms, const SpinMatrix& denPhi,
					  const DoubleVec& cartDen, const SigMolShellPairInfor& infor, const DoubleVec& wRho,
                                          SpinMatrix& halfExRho, DoubleVec& halfJRhoVec, DoubleVec& halfJRhoVecMtrx) const;


			void doESPIntWorkXHole(const DoubleVec& grids, const MolShell& ms, const SpinMatrix& denPhi, 
					const DoubleVec& cartDen, const SigMolShellPairInfor& infor, const DoubleVec& wRho, 
					  SpinMatrix& halfExRho, DoubleVec& halfJRhoVec, DoubleVec& halfJRhoVecMtrx, Double sValue) const;//yw

		public:

			/**
			 * constructor - build the int job inside
			 */
			ESPInts(const GIntsInfor& infor0, UInt job):ginfor(infor0,job,0) { };

			/**
			 * destructor
			 */
			~ESPInts() { };

			///
			///
			/// \param  printTiming     whether to print out the timing
			/// \param  ms              input shell data for both row and col
			/// \param  spData          the significant shell pair information
			/// \param  grid            the grid data for esp integrals
			/// \param  denPhi          combination of the density matrix and phi value (nCarBas,nGrids)
			/// \param  cartDen         the shell pair vector form density matrix in Cartesian order
			/// \param  wrho            the density on this batch grid together with weights
			/// \return halfExRho       the result half exchange rho (nGrids,nBas)
			/// \return halfJRhoVec     half of the Coulomb energy density vector in J1 part
			/// \return halfJRhoVecMtrx half of the Coulomb energy density vector form matrix in J2 part
			///
			void doMtrx(const XCIntJobInfor& xcInfor, const MolShell& ms,
				    const SigMolShellPairInfor& spData, const DoubleVec& grid,
				    const SpinMatrix& denPhi, const DoubleVec& cartDen, const DoubleVec& wRho,
				    SpinMatrix& halfExRhoResult, DoubleVec& halfJRhoVec,
                                    DoubleVec& halfJRhoVecMtrx) const;

			void doMtrxXHole(const XCIntJobInfor& xcInfor, const MolShell& ms, 
					const SigMolShellPairInfor& spData, const DoubleVec& grid, 
					const SpinMatrix& denPhi, const DoubleVec& cartDen, const DoubleVec& wRho,
					SpinMatrix& halfExRhoResult, DoubleVec& halfJRhoVec, 
				    DoubleVec& halfJRhoVecMtrx, Double sValue) const;//yw
	};
}

#endif

