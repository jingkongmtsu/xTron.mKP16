/**
 * \file    scf.h
 * \brief   class mornitoring the whole SCF calculation
 * \author  Fenglai Liu 
 * \note
 */
#ifndef SCF_H
#define SCF_H
#include "libgen.h"
#include "scfparam.h"
#include "scfintscontroller.h"
#include "mo.h"
#include "oneemtrx.h"

namespace shell {
	class MolShell;
}

namespace molecule {
	class Molecule;
}

namespace denmtrx {
	class DenMtrx;
}

namespace fock {
	class Fock;
}

namespace scfconv {
	class SCFConv;
}

namespace scf {

	using namespace molecule;
	using namespace shell;
	using namespace denmtrx;
	using namespace scfparam;
	using namespace scfintscontroller;
	using namespace oneemtrx;
	using namespace mo;
	using namespace fock;
	using namespace scfconv;

	/**
	 * this class is used to monitor the SCF calculation
	 */
	class SCF {

		private:

			//
			//  enegy section
			//
			Double     eNulRep;                ///< nuclear repulsion energy from molecule 
			Double     eCore;                  ///< core Hamiltonian energy
			Double     eJK;                    ///< energy for Coulomb/exchange
			Double     eXC;                    ///< energy for XC
			Double     eTotal;                 ///< energy account for the total energy 

			//
			// scf parameters
			//
			SCFParam   param;                  ///< scf parameters center

			//
			// scf integrals controlling center
			//
			SCFIntsController intsController;  ///< controlling scf integrals calculation

			//
			// one electron matrix
			// keep it for use in recycle way
			//
			OneEMtrx   oneEMtrx;               ///< a copy of one electron matrix

			//
			// scf result
			//
			MO         mo;                     ///< the final result MO

		public:

			/**
			 * constructor for SCF process
			 * preparation for the practical SCF calculation
			 * \param  infor     global information for the job
			 * \param  mol       molecule data
			 * \param  s         shell data
			 * 
			 * inside this constructor we do not initilize the mo first
			 * it will be done until the end cycle of SCF is reached
			 */
			SCF(const GlobalInfor& infor, const Molecule& mol, const MolShell& s);

			/**
			 * constructor for SCF process
			 * preparation for the practical SCF calculation
			 * \param  infor     global information for the job
			 * \param  param     the input SCF parameter center
			 * \param  mol       molecule data
			 * \param  s         shell data
			 * 
			 * compare with the constructor above, this one will use the passed
			 * in scfparam
			 */
			SCF(const GlobalInfor& infor, const SCFParam& param0, const Molecule& mol, 
					const MolShell& s);

			/**
			 * default destructor
			 */
			~SCF() { };

			/**
			 * perform SCF calculation
			 * \param infor      global information for the job
			 * \param mol        molecule data
			 * \param ms         shell data
			 * \param den        input density matrix for SCF, will be renewed inside
			 */
			void doSCF(const GlobalInfor& infor, const Molecule& mol, 
					const MolShell& ms, DenMtrx& den);

			void doSCFConv(const GlobalInfor& infor, const Molecule& mol, 
					const MolShell& ms, DenMtrx& den, SCFConv& conv);


			/**
			 * perform post SCF calculation
			 *
			 * compared with the SCF process, the density matrix now is constant
			 * so we do not change it inside
			 */
			void doPostSCF(const GlobalInfor& infor, const Molecule& mol, 
					const MolShell& ms, const DenMtrx& den) const;

			/**
			 * return the result MO
			 */
			const MO& getMO() const { return mo; };

			/**
			 * update current MO by input data
			 */
			void updateMO(const MO& mos); 

			/**
			 * return the total energy
			 */
			Double energy() const { return eTotal; };

			/**
			 * return the oneEMtrx
			 */
			const OneEMtrx& getOneEMtrx() const { return oneEMtrx; };

			/**
			 * return the scfparam
			 */
			const SCFParam& getSCFParam() const { return param; };

	};
}

#endif

