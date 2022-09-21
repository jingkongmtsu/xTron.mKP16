/**
 * \file    gdm.h
 */
#ifndef GDM_H
#define GDM_H
#include "libgen.h"
#include "matrix.h"
#include "shell.h"
#include "globalinfor.h"
#include "scfparam.h"
#include "fock.h"
#include "mo.h"
using namespace matrix;
using namespace globalinfor;
using namespace shell;
using namespace scfparam;
using namespace fock;
using namespace mo;

//
// initialize the GDM file path
//
extern void initGDMFilePath(const char* path0);

#define  rungdm   rungdm_

extern "C" {

	extern void rungdm(double* CA, double* CB, double* theta, double* dEdT,
			Int* iter, double* FA, double* FB, double* E, Int* NBasis,  
			Int* nOrb, Int* nOA, Int* nOB);
}

namespace gdm {

	class GDM {

		private:

			//
			//  data section
			//
			Int iter;         
			Int nBasis;
			Int nOrb;
			Int nOA;
			Int nOB;
			Int nSpin;
			Int nVA;
			Int nVB;
			DoubleVec theta;
			DoubleVec dEdT;
			Mtrx CA;
			Mtrx FA;
			Mtrx CB;
			Mtrx FB;

			///
			/// reset the number of orbitals
			///
			void resetOrbInfor(const Int& nOrb0, const Int& iSpin) {
				nOrb = nOrb0;
				bool reform = true;
				if (iSpin == 0) {
					CA.reset(nBasis,nOrb,reform);
				}else{
					CB.reset(nBasis,nOrb,reform);
				}
			};

		public:

			/**
			 * constructor
			 */
			GDM(const SCFParam& param, const Molecule& mol, const MolShell& ms):iter(0),nBasis((Int)ms.getNBas()),
			nOrb((Int)ms.getNBas()),nOA((Int)mol.getNAlpha()),nOB((Int)mol.getNBeta()),nSpin((Int)param.getNSpin()),
			nVA(nBasis-nOA),nVB(nBasis-nOB),theta(nOA*nVA+nOB*nVB,ZERO),dEdT(nOA*nVA+nOB*nVB,ZERO),CA(nBasis,nBasis),
			FA(nBasis,nBasis),CB(nBasis,nBasis),FB(nBasis,nBasis) { 

				// we need to initialize the path
				const GlobalInfor& infor = param.getGlobalInfor();
				string path(infor.getScratchDir()); 
				initGDMFilePath(path.c_str());
			};

			/**
			 * destructor
			 */
			~GDM() { };

			/**
			 * do the GDM calculation
			 */
			void doGDM(const UInt& scfIter, MO& mos, Fock& fock, Double scfEnergy) {

				// check whether we need to reset MO dimension?
				for(Int iSpin=0; iSpin<nSpin; iSpin++) {
					Int nMO = (Int)mos.getNOrb(iSpin);
					if (nMO != nOrb) {
						resetOrbInfor(nMO,iSpin);
					}
				}

				// set the iter
				if (scfIter == 0) {
					iter = 0;
				}else{
					iter = 1;
				}

				// copy the MO data
				CA = mos.getMtrx(0);
				FA = fock.getMtrx(0);
				FA.copyLowToUpper();
				//CA.print("MO A");
				//FA.print("FOCK A");
				if (nSpin == 2) {
					CB = mos.getMtrx(1);
					FB = fock.getMtrx(1);
					FB.copyLowToUpper();
					//CB.print("MO B");
					//FB.print("FOCK B");
				}
				Double E = scfEnergy;

				// now do calculation
				rungdm_(CA.getPtr(),CB.getPtr(),&theta.front(),&dEdT.front(),
						&iter,FA.getPtr(),FB.getPtr(),&E,&nBasis,&nOrb,
						&nOA,&nOB);

				// finally copy the result back
				Mtrx& MA = mos.getMtrx(0);
				Mtrx& FockA = fock.getMtrx(0);
				MA = CA;
				FockA = FA;
				if (nSpin == 2) {
					Mtrx& MB = mos.getMtrx(1);
					Mtrx& FockB = fock.getMtrx(1);
					MB = CB;
					FockB = FB;
				}
			};

	};
}

#endif

