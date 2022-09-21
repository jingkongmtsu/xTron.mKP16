//
// set up testing file for exrho
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string> 
#include<vector>
#include<cmath>
#include<boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include "excep.h"
#include "globalinfor.h"
#include "molecule.h"
#include "shell.h"
#include "blas1.h"
#include "gintsinfor.h"
#include "integraljobs.h"
#include "sigshellpairinfor.h"
using namespace std;
using namespace excep;
using namespace globalinfor;
using namespace molecule;
using namespace shell;
using namespace blas;
using namespace gintsinfor;
using namespace integraljobs;
using namespace sigshellpairinfor;

void espints_test()
{

	/*

	// set input file
	string inf = "espints/input.in";

	// set up grid points
	UInt nGrids = MAXNG_PHI_THREAD;
	DoubleVec grid(3*nGrids);
	for(UInt i=0; i<nGrids; i++) {
		Double t = (50.0E0/nGrids)*i;
		grid[3*i+0] = 0.1E0+t;
		grid[3*i+1] = 0.2E0+t;
		grid[3*i+2] = 0.3E0+t;
	}

	// set up molecule and shell data
	GlobalInfor globalInfor(inf);

	// loop over molecules
	for(UInt iMol=1; iMol<=2; iMol++) {
		Molecule molecule(inf,iMol);
		MolShell s(inf,molecule);

		// now set up the compact shell pair infor
		// we just use arbitrary data, large enough to hold the memory
		CompactShellPairData_Phi spData; 
		setCompactShellPairData(1000,4000,100,spData);

		// let's pick up two arbitrary atom shell and form the compact shell pair data
		UInt maxNBas = 0;
		UInt job = EXCHANGE;

		// let's pick up the atom index
		UInt rowAtomIndex = 1;
		UInt colAtomIndex = 0;
		{
			// form the molshell data array
			MolShellDataArray msArray(s);

			// set up shell pair information data
			GIntsInfor ginfor(globalInfor,molecule);
			GIntJobInfor jobInfor(ginfor,job,0);
			SigMolShellPairInfor sp(s,s,jobInfor); 
			CompactSigMolShellPairInforArray comSP(sp);
			const CompactSigShellPairInfor*     spArray  = comSP.getCompactSigShellPairInforArray();
			const CompactSigAtomShellPairInfor* aspArray = comSP.getCompactSigAtomShellPairInforArray(); 

			// set the atom center
			Double A[3];
			Double B[3];
			const AtomShell& rowAS = s.getAtomShell(rowAtomIndex);
			const AtomShell& colAS = s.getAtomShell(colAtomIndex);
			const Double* xyz = rowAS.getXYZ(); 
			A[0] = xyz[0];
			A[1] = xyz[1];
			A[2] = xyz[2];
			xyz  = colAS.getXYZ(); 
			B[0] = xyz[0];
			B[1] = xyz[1];
			B[2] = xyz[2];
			UInt rowAtomic = rowAS.getAtomic();
			UInt colAtomic = colAS.getAtomic();
			if (rowAS.getNCarBas() > maxNBas) maxNBas = rowAS.getNCarBas();
			if (colAS.getNCarBas() > maxNBas) maxNBas = colAS.getNCarBas();

			// get the row and column atom
			UInt status;
			const AtomShellTypeDataInfor& rowAtomShellData = msArray.getAtomShellTypeDataInfor(rowAtomic,status);
			if (status != 0) {
				printf("can not get the rowAtomShellData\n");
				exit(1);
			}
			const AtomShellTypeDataInfor& colAtomShellData = msArray.getAtomShellTypeDataInfor(colAtomic,status);
			if (status != 0) {
				printf("can not get the colAtomShellData\n");
				exit(1);
			}

			// get the caspInfor
			const CompactSigAtomShellPairInfor* caspInfor = NULL;
			for(UInt j=0; j<comSP.getLenCompactSigAtomShellPairInforArray();j++) {
				UInt rowAtomIndex0 = compactsigshellpairinfor::getRowAtomShellIndex(aspArray[j]);
				UInt colAtomIndex0 = compactsigshellpairinfor::getColAtomShellIndex(aspArray[j]);
				//printf("row atom index %d, col %d\n", rowAtomIndex0, colAtomIndex0);
				if (rowAtomIndex0 == rowAtomIndex && colAtomIndex  == colAtomIndex0) {
					caspInfor = &aspArray[j];
					break;
				}
			}
			if (caspInfor == NULL) {
				printf("can not get the corresponding caspInfor\n");
				exit(1);
			}

			// build the compact atom shell pair data
			Double thresh = 1.0E-12;
			compactshellpairdata_phi::init(thresh,*caspInfor,spArray,A,B,rowAtomShellData,colAtomShellData,spData);

			// let's arbitraryly set the C2 to be all one
			setC2ArrayVal(ONE,4000,spData);
		}

		// set up the exchange den phi
		Mtrx aBraAtomBlockDenPhi(maxNBas,nGrids);
		aBraAtomBlockDenPhi.set(ONE);
		Mtrx aKetAtomBlockDenPhi(maxNBas,nGrids);
		aKetAtomBlockDenPhi.set(ONE);
		Mtrx bBraAtomBlockDenPhi(maxNBas,nGrids);
		bBraAtomBlockDenPhi.set(ONE);
		Mtrx bKetAtomBlockDenPhi(maxNBas,nGrids);
		bKetAtomBlockDenPhi.set(ONE);

		// set the result
		Mtrx aBraAtomBlockExRho0(maxNBas,nGrids);
		Mtrx aKetAtomBlockExRho0(maxNBas,nGrids);
		Mtrx bBraAtomBlockExRho0(maxNBas,nGrids);
		Mtrx bKetAtomBlockExRho0(maxNBas,nGrids);
		Mtrx aBraAtomBlockExRho1(maxNBas,nGrids);
		Mtrx aKetAtomBlockExRho1(maxNBas,nGrids);
		Mtrx bBraAtomBlockExRho1(maxNBas,nGrids);
		Mtrx bKetAtomBlockExRho1(maxNBas,nGrids);

		// set the par
		ESPIntsPar1 par1; 
		par1.intJob        = job;
		par1.nTolCarBas    = 0;
		par1.nSpin         = molecule.getNSpin();
		par1.nPts          = nGrids;
		par1.start         = 0;
		par1.sameAtomShell = 0;
		par1.intThresh     = 1.0E-12;

		// set the par2
		ESPIntsPar2 par2; 
		par2.rowNBas = maxNBas;
		par2.colNBas = maxNBas;
		par2.rowAtomShellIndex = rowAtomIndex;
		par2.colAtomShellIndex = colAtomIndex;
		par2.rowCartBasOffset  = 0;
		par2.colCartBasOffset  = 0;

		// set the coefficient array
		// it's always one
		DoubleVec C2(4000,ONE);

		// set the esp raw integrals
		DoubleVec esp(900*nGrids);

		// set up the L Code for testing
		UInt nAngCode = compactshellpairdata_phi::getNAngCode(spData);
		for(UInt iAng=0; iAng<nAngCode; iAng++) {

			// clear the results
			aBraAtomBlockExRho0.set(ZERO);
			aKetAtomBlockExRho0.set(ZERO);
			bBraAtomBlockExRho0.set(ZERO);
			bKetAtomBlockExRho0.set(ZERO);
			aBraAtomBlockExRho1.set(ZERO);
			aKetAtomBlockExRho1.set(ZERO);
			bBraAtomBlockExRho1.set(ZERO);
			bKetAtomBlockExRho1.set(ZERO);

			// get the pointer
			const Double* aBraAtomBlockDenPhi_ptr = aBraAtomBlockDenPhi.getPtr(0,0);
			const Double* bBraAtomBlockDenPhi_ptr = bBraAtomBlockDenPhi.getPtr(0,0);
			const Double* aKetAtomBlockDenPhi_ptr = aKetAtomBlockDenPhi.getPtr(0,0);
			const Double* bKetAtomBlockDenPhi_ptr = bKetAtomBlockDenPhi.getPtr(0,0);
			Double* aBraAtomBlockExRho_ptr = aBraAtomBlockExRho0.getPtr(0,0);
			Double* bBraAtomBlockExRho_ptr = bBraAtomBlockExRho0.getPtr(0,0);
			Double* aKetAtomBlockExRho_ptr = aKetAtomBlockExRho0.getPtr(0,0);
			Double* bKetAtomBlockExRho_ptr = bKetAtomBlockExRho0.getPtr(0,0);

			// now call the mic version
			UInt code = compactshellpairdata_phi::getAngCode(spData,iAng);
			switch(code) {
				case 0:
					hgp_os_esp_s_s_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 1:
					hgp_os_esp_p_s_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 100:
					hgp_os_esp_sp_s_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 1001:
					hgp_os_esp_p_p_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 100100:
					hgp_os_esp_sp_sp_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 100001:
					hgp_os_esp_p_sp_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 2:
					hgp_os_esp_d_s_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 1002:
					hgp_os_esp_d_p_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 100002:
					hgp_os_esp_d_sp_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 2002:
					hgp_os_esp_d_d_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 3:
					hgp_os_esp_f_s_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 1003:
					hgp_os_esp_f_p_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 100003:
					hgp_os_esp_f_sp_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 2003:
					hgp_os_esp_f_d_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				case 3003:
					hgp_os_esp_f_f_mic(iAng,spData,par1,par2,&grid.front(),
							aBraAtomBlockDenPhi_ptr,aBraAtomBlockExRho_ptr,
							bBraAtomBlockDenPhi_ptr,bBraAtomBlockExRho_ptr, 
							aKetAtomBlockDenPhi_ptr,aKetAtomBlockExRho_ptr,
							bKetAtomBlockDenPhi_ptr,bKetAtomBlockExRho_ptr);
					break;
				default:
					{
						printf("%s %d\n","Un-recognized LCode in the esp integrals calculation on phi", (Int)code);
						assert(0);
						break;
					}
			}

			// next call the normal codes
			UInt nSP = compactshellpairdata_phi::getNSP(spData,iAng);
			UInt spOffset = compactshellpairdata_phi::getSPOffset(spData,iAng);
			for(UInt iSP=spOffset; iSP<spOffset+nSP; iSP++) {
				UInt np2          = compactshellpairdata_phi::getNP2(spData,iSP);
				const Double* e2  = compactshellpairdata_phi::getE2(spData,iSP);
				const Double* P   = compactshellpairdata_phi::getP(spData,iSP);
				const Double* fac = compactshellpairdata_phi::getFac(spData,iSP);
				const Double* c2  = &C2.front();

				// center
				const Double* A = compactshellpairdata_phi::getA(spData);
				const Double* B = compactshellpairdata_phi::getB(spData);
				UInt spStatus = compactshellpairdata_phi::getShellStatus(spData,iSP); 
				if (spStatus == WITH_SHELL_SWITCH) {
					A = compactshellpairdata_phi::getB(spData);
					B = compactshellpairdata_phi::getA(spData);
				}

				// do the raw integral
				UInt iLmin,iLmax,jLmin,jLmax;
				decodeL(code,iLmin,iLmax,jLmin,jLmax);
				UInt iNCarBas = getCartBas(iLmin,iLmax);
				UInt jNCarBas = getCartBas(jLmin,jLmax);
				UInt nPts   = nGrids;
				UInt intLen = iNCarBas*jNCarBas*nPts;
				vset(&esp.front(),ZERO,intLen);
				hgp_os_esp(code,np2,nPts,c2,e2,fac,P,A,B,&grid[0],&esp.front()); 

				// set the local shell pair offset
				UInt iLocOffset,jLocOffset;  
				compactshellpairdata_phi::getLocCarBasOffset(spData,iSP,iLocOffset,jLocOffset); 

				// set up the tmp result
				Double shellBlockDenPhi[40];
				Double shellBlockExRho[40];

				// digestion
				for(UInt iPt=0; iPt<nPts; iPt++) {

					// get the esp data
					const Double* espData = &esp[iNCarBas*jNCarBas*iPt];

					// now let's do digestion
					// firstly it's exchange part
					if (doK(job)) {

						//
						// this is alpha spin part
						//

						// firstly let's do the bra side denphi digestion
						// this will produce the ket side result
						// iLocOffset and iNCarBas are all associated with the bra shell block index
						// get the denphi
						if (spStatus == WITH_SHELL_SWITCH) {
							vcopy(aKetAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
						}else{
							vcopy(aBraAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
						}

						// digest data to form the result, result is in ket side shell block
						// jLocOffset and jNCarBas are all associated with the ket shell block index
						for(UInt j=0; j<jNCarBas; j++) {
							shellBlockExRho[j] = MINUS_ONE*vdot(&espData[j*iNCarBas],&shellBlockDenPhi[0],iNCarBas);
						}

						// let's put the result back to atom shell block
						// the result is in ket shell block
						if (spStatus == WITH_SHELL_SWITCH) {
							vaxpy(&shellBlockExRho[0],aBraAtomBlockExRho1.getPtr(jLocOffset,iPt),ONE,jNCarBas);
						}else{
							vaxpy(&shellBlockExRho[0],aKetAtomBlockExRho1.getPtr(jLocOffset,iPt),ONE,jNCarBas);
						}


						// read in denphi
						// this is the ket side denphi digestion
						// jLocOffset and jNCarBas are all associated with the ket shell block index
						if (spStatus == WITH_SHELL_SWITCH) {
							vcopy(aBraAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
						}else{
							vcopy(aKetAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
						}

						// form results
						// result is in bra side shell block
						// iLocOffset and iNCarBas are all associated with the ket shell block index
						for(UInt i=0; i<iNCarBas; i++) shellBlockExRho[i] = ZERO;
						for(UInt i=0; i<iNCarBas; i++) {
							for(UInt j=0; j<jNCarBas; j++) {
								shellBlockExRho[i] += MINUS_ONE*espData[i+j*iNCarBas]*shellBlockDenPhi[j];
							}
						}

						// let's put the result into the atom block result
						// remember that the result is in bra side
						// iLocOffset and iNCarBas are all associated with the bra shell block index
						if (spStatus == WITH_SHELL_SWITCH) {
							vaxpy(&shellBlockExRho[0],aKetAtomBlockExRho1.getPtr(iLocOffset,iPt),ONE,iNCarBas);
						}else{
							vaxpy(&shellBlockExRho[0],aBraAtomBlockExRho1.getPtr(iLocOffset,iPt),ONE,iNCarBas);
						}

						//
						// this is beta spin part
						//
						if (molecule.getNSpin() == 2) {

							// firstly let's do the bra side denphi digestion
							// this will produce the ket side result
							// iLocOffset and iNCarBas are all associated with the bra shell block index
							// get the denphi
							if (spStatus == WITH_SHELL_SWITCH) {
								vcopy(bKetAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
							}else{
								vcopy(bBraAtomBlockDenPhi.getPtr(iLocOffset,iPt),&shellBlockDenPhi[0],iNCarBas);
							}

							// digest data to form the result, result is in ket side shell block
							// jLocOffset and jNCarBas are all associated with the ket shell block index
							for(UInt j=0; j<jNCarBas; j++) {
								shellBlockExRho[j] = MINUS_ONE*vdot(&espData[j*iNCarBas],&shellBlockDenPhi[0],iNCarBas);
							}

							// let's put the result back to atom shell block
							// the result is in ket shell block
							if (spStatus == WITH_SHELL_SWITCH) {
								vaxpy(&shellBlockExRho[0],bBraAtomBlockExRho1.getPtr(jLocOffset,iPt),ONE,jNCarBas);
							}else{
								vaxpy(&shellBlockExRho[0],bKetAtomBlockExRho1.getPtr(jLocOffset,iPt),ONE,jNCarBas);
							}


							// read in denphi
							// this is the ket side denphi digestion
							// jLocOffset and jNCarBas are all associated with the ket shell block index
							if (spStatus == WITH_SHELL_SWITCH) {
								vcopy(bBraAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
							}else{
								vcopy(bKetAtomBlockDenPhi.getPtr(jLocOffset,iPt),&shellBlockDenPhi[0],jNCarBas);
							}

							// form results
							// result is in bra side shell block
							// iLocOffset and iNCarBas are all associated with the ket shell block index
							for(UInt i=0; i<iNCarBas; i++) shellBlockExRho[i] = ZERO;
							for(UInt i=0; i<iNCarBas; i++) {
								for(UInt j=0; j<jNCarBas; j++) {
									shellBlockExRho[i] += MINUS_ONE*espData[i+j*iNCarBas]*shellBlockDenPhi[j];
								}
							}

							// let's put the result into the atom block result
							// remember that the result is in bra side
							// iLocOffset and iNCarBas are all associated with the bra shell block index
							if (spStatus == WITH_SHELL_SWITCH) {
								vaxpy(&shellBlockExRho[0],bKetAtomBlockExRho1.getPtr(iLocOffset,iPt),ONE,iNCarBas);
							}else{
								vaxpy(&shellBlockExRho[0],bBraAtomBlockExRho1.getPtr(iLocOffset,iPt),ONE,iNCarBas);
							}
						}
					}
				}
			}

			// now let's compare the two matrix
			// alpha first
			Double th = 1.0E-10;
			if (aBraAtomBlockExRho0.diffComp(aBraAtomBlockExRho1,th,false,true)) {
				printf("for the %d code the aBraAtomBlockExRho0 testing is failed\n", code);
			}
			if (aKetAtomBlockExRho0.diffComp(aKetAtomBlockExRho1,th,false,true)) {
				printf("for the %d code the aKetAtomBlockExRho0 testing is failed\n", code);
			}
			if (molecule.getNSpin() == 2) {
				if (bBraAtomBlockExRho0.diffComp(bBraAtomBlockExRho1,th,false,true)) {
					printf("for the %d code the bBraAtomBlockExRho0 testing is failed\n", code);
				}
				if (bKetAtomBlockExRho0.diffComp(bKetAtomBlockExRho1,th,false,true)) {
					printf("for the %d code the bKetAtomBlockExRho0 testing is failed\n", code);
				}
			}

			// debug code
			bool doDebug = false;
			if (doDebug) {
				if (code == 3003) {
					aBraAtomBlockExRho0.print("matrix from compact form");
					aBraAtomBlockExRho1.print("matrix from normal  form");
				}
			}
		}

		// finally free the memory
		freeCompactShellPairData(spData);
	}

*/
}
