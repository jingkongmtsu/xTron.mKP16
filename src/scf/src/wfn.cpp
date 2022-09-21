/**
 * \file    wfn.cpp
 * \brief   create the wfn file for the given set of data
 * \author  Fenglai Liu and Jing Kong
 */
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include "parameterparsing.h"
#include "filerw.h"
#include "shell.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "element.h"
#include "molecule.h"
#include "textread.h"
#include "cptrans.h"
#include "matrix.h"
#include "blas.h"
#include "mo.h"
#include "wfn.h"
using namespace boost::filesystem;
using namespace parameterparsing;
using namespace filerw;
using namespace shell;
using namespace shellprop;
using namespace matrix;
using namespace element;
using namespace molecule;
using namespace textread;
using namespace cptrans;
using namespace blas;
using namespace mo;
using namespace wfn;
using namespace std;

WFN::WFN(const GlobalInfor& infor, const Molecule& mol):outputFile("data.wfx"),section(-1)
{ 
	string input = infor.getInputFile();
	ParameterParsing pp(input,"wfn_data",mol.getSec());
	if (pp.hasAnyParameters()) {

		// now assign section number
		section = mol.getSec();

		// wfn file path
		string key = "wfn_path";
		if (pp.hasKeyDefined(key)) {

			// whether we use the default path?
			string p(pp.getValue(key));
			if (p == "default") {
				// if no path defined, we just issue a warning message
				string infor = "in default, the result wfn file will be generated in scratch_dir/section_number/wfn/data.wfx";
				Excep excep("WFN","constructor",EXCEPTION_USE_DEFAULT_FILE_PATH,infor);
				handleExcep(excep);
			}else{
				UInt begin = p.length()-4;
				string sp = p.substr(begin);
				if (!(sp == ".wfx")) {
					string infor = "the path you give should have the file name, with suffix of .wfx";
					Excep excep("WFN","constructor",EXCEPTION_IMPROPER_FILE_NAME,infor);
					handleExcep(excep);
				}
				outputFile = p;
			}
		}else{
			string infor = "You did not define where to store the wfn data file";
			Excep excep("WFN","constructor",EXCEPTION_IMPROPER_FILE_NAME,infor);
			handleExcep(excep);
		}
	}
}

void WFN::dump(const Double& totalEnergy, const GlobalInfor& infor, const Molecule& mol, 
		const MolShell& ms, const MO& mo) const
{
	// check that whether we do wfn file
	if (! doWFNFile()) return; 

	// get the path 
	// if no dir information, and the same default
	// file name used, we just go to default dir
	path wfnName(outputFile.c_str());
	string file = wfnName.string();
	if (outputFile == "data.wfx") {

		// construct the dir for wfn file
		const string& scratchDir = infor.getScratchDir();
		path p(scratchDir.c_str());
		string sec = boost::lexical_cast<string>(section);
		path sect(sec.c_str());
		string wfn = "wfn";
		path wfnPath(wfn.c_str());
		p /= sect;
		p /= wfnPath;

		// remove old data if there's same name folder exists
		// else create new one
		if (exists(p)) {
			remove_all(p);
		}
		create_directories(p);

		// add in file name
		p /= wfnName;
		file = p.string();
	}

	// set up molecule information
	UInt nAtoms = mol.getNAtoms();
	DoubleVec coords(3*nAtoms);
	UIntVec atomicVec(nAtoms);
	for(UInt i=0; i<nAtoms; i++) {
		const Double* xyz = mol.getXYZ(i);
		coords[3*i+0] = xyz[0];
		coords[3*i+1] = xyz[1];
		coords[3*i+2] = xyz[2];
		atomicVec[i]  = mol.getAtomic(i);
	}

	// set up the data
	// firstly calculate how many primitive functions
	UInt nTotPrims = 0;
	for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {
		const AtomShell& atomShell = ms.getAtomShell(iAtom);
		for(UInt iShell=0; iShell<atomShell.getNShell(); iShell++) {
			const Shell& s = atomShell.getShell(iShell);
			UInt lmin  = s.getLmin();
			UInt lmax  = s.getLmax();
			UInt nBas  = getCartBas(lmin,lmax);
			UInt nPrim = nBas*s.getNPrim();
			nTotPrims += nPrim;
		}
	}   

	// count the mo information
	UInt nAlphaOcc   = mo.getNOcc(0);
	UInt nBetaOcc    = mo.getNOcc(1);
	UInt nTolBas     = ms.getNBas();
	UInt nTolCarBas  = ms.getNCarBas();
	UInt nCol = nAlphaOcc >= nBetaOcc ? nAlphaOcc : nBetaOcc;
	DoubleVec moEnergy(mo.getNSpin()*nCol);
	DoubleVec moOCCNumber(mo.getNSpin()*nCol);
	SpinMatrix moData(mo.getNSpin(),nTolCarBas,nCol);
	{
		// decide the occupation number and associated MO energy
		if (mo.getNSpin() == 1) {
			const DoubleVec& E = mo.getEnergyVec(0);
			for(UInt i=0; i<nAlphaOcc; i++) {
				moOCCNumber[i] = TWO;
				moEnergy[i]    = E[i];
			}
		}else{
			const DoubleVec& aE = mo.getEnergyVec(0);
			const DoubleVec& bE = mo.getEnergyVec(1);
			for(UInt i=0; i<nAlphaOcc+nBetaOcc; i++) {
				moOCCNumber[i] = ONE;
				if (i<nAlphaOcc) {
					moEnergy[i] = aE[i];
				}else{
					moEnergy[i] = bE[i-nAlphaOcc];
				}
			}
		}

		// now let's create the occupied mo data
		for(UInt i=0; i<mo.getNSpin(); i++) {
			Mtrx& M = moData.getMtrx(i);
			const Mtrx& M0 = mo.getMtrx(i);
			UInt nC = nAlphaOcc;
			if (i == 1) nC = nBetaOcc;
			M.reset(nTolBas,nC,false);
			M.copyFromMatrix(0,0,0,0,nTolBas,nC,M0);
		}

		// now let's transform the data
		// the mo should be appended with reverse of C2P matrix, where
		// C is lxlylz to spherical
		if (! ms.allPure() && ! ms.allCart()) {
			string infor = "for generating the wfn file the shell data should be either all spherical shells, "
			"or all Cartesian shell data. The mixing between the two is not supported";
			Excep excep("WFN","constructor",EXCEPTION_INVALID_SHELL_DATA,infor);
			handleExcep(excep);
		}
		if (ms.allPure()) {
			CPTransBasisNorm p2c(infor,C2P_WITH_XYZ,NO_SCALE,WITH_MATRIX_TRANSPOSE,CP_WITH_ROW);
			for(UInt i=0; i<mo.getNSpin(); i++) {
				Mtrx& M = moData.getMtrx(i);
				p2c.transform(ms,ms,M);
			}
		}
	}


	// now reserve data
	// one is recording each primitive center
	// and the primitive angular types (see the http://aim.tkgristmill.com/wfxformat.html#primitivetypes)
	// one is for primitive exponents
	// one is for original coefficients
	// one is for the final mo data per each primitive
	UIntVec   primCenters(nTotPrims);
	UIntVec   primTypes(nTotPrims);
	DoubleVec primExp(nTotPrims);
	DoubleVec primCoe(nTotPrims);

	// now let's fetch the exp data and coefficient data
	// we note, that originally the coeff is normalized with X^{L}Y^{0}Z^{0}
	// now we need to load in the data so that each basis set normalized 
	// to it's own x^{l}y^{m}z^{n}
	UInt offset = 0;
	for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {
		const AtomShell& atomShell = ms.getAtomShell(iAtom);
		for(UInt iShell=0; iShell<atomShell.getNShell(); iShell++) {
			const Shell& s = atomShell.getShell(iShell);
			UInt lmin  = s.getLmin();
			UInt lmax  = s.getLmax();
			UInt nP    = s.getNPrim();
			const Double* e = s.getExp();
			for(UInt iL=lmin; iL<=lmax; iL++) {
				UInt nBas = getCartBas(iL,iL);
				const Double* c = s.getCoe(iL);

				// set up the possible normalization
				// for D, F shell etc. we need to re-scale the 
				// normalization from x^{L}y^{0}Z^{0} to x^{l}y^{m}Z^{n}
				DoubleVec normVec(nBas,ONE);
				if (iL>=2) {
					Double baseNorm = normCartBasisSets(iL,0,0,nP,c,e);
					for(UInt iBas = 0; iBas<nBas; iBas++) {
						UInt l,m,n;
						getlmn(iL,iL,iBas,l,m,n);
						Double norm = normCartBasisSets(l,m,n,nP,c,e);
						normVec[iBas] = norm/baseNorm;
					}
				}

				// now let's assign data
				// we note that the type returned from the function starts from 0
				// however for AIMPAC order it starts from 1
				// so we add 1
				for(UInt iBas = 0; iBas<nBas; iBas++) {
					Double n  = normVec[iBas];
					UInt type = shellprop::getBasIndexForTheBasisSetOrder(iL,iBas,AIMPAC_BASIS_SET_ORDER);
					for(UInt i = 0; i<nP; i++) {
						primExp[offset+i] = e[i];
						primCoe[offset+i] = c[i]*n;
						primCenters[offset+i] = iAtom+1;
						primTypes[offset+i]   = type+1;
					}
					offset += nP;
				}
			}
		}
	}   

	// for the mo coefficients, we need to do P2C into lxlylz??
	
	// open the output stream
	const char * output_file = file.c_str();
	FILE* fp = NULL;
	if ((fp = fopen(output_file,"w")) == NULL) {
		cout << "wfn file name and path is " << file << endl;
		string infor = "can not open the given wfn data file for writing";
		Excep excep("WFN","dump",EXCEPTION_FILE_WRITE_FAIL,infor);
		handleExcep(excep);
	}

	// print out the title
	fprintf(fp, "<Title>\n");
	fprintf(fp, "Currently no title information available\n");
	fprintf(fp, "</Title>\n\n");

	// keywords
	fprintf(fp, "<Keywords>\n");
	fprintf(fp, "GTO \n");
	fprintf(fp, "</Keywords>\n\n");

	// number of nuclei
	fprintf(fp, "<Number of Nuclei>\n");
	fprintf(fp, "%d \n", (Int)nAtoms);
	fprintf(fp, "</Number of Nuclei>\n\n");

	// number of primitives
	fprintf(fp, "<Number of Primitives>\n");
	fprintf(fp, "%d \n", (Int)nTotPrims);
	fprintf(fp, "</Number of Primitives>\n\n");

	// number of occupied orbitals
	fprintf(fp, "<Number of Occupied Molecular Orbitals>\n");
	if (mo.getNSpin() == 1) {
		fprintf(fp, "%d \n", (Int)nAlphaOcc);
	}else{
		fprintf(fp, "%d \n", (Int)(nAlphaOcc+nBetaOcc));
	}
	fprintf(fp, "</Number of Occupied Molecular Orbitals>\n\n");

	// pertubation section
	fprintf(fp, "<Number of Perturbations>\n");
	fprintf(fp, "0 \n");
	fprintf(fp, "</Number of Perturbations>\n\n");

	// nuclei name
	fprintf(fp, "<Nuclear Names>\n");
	for(UInt i=0; i<nAtoms; i++) {
		UInt atomic = atomicVec[i];
		string name = getAtomicSymbol(atomic);
		name = name + boost::lexical_cast<string>(i+1);
		fprintf(fp, "%-8s\n", name.c_str());
	}
	fprintf(fp, "</Nuclear Names>\n\n");

	// atomic number
	fprintf(fp, "<Atomic Numbers>\n");
	for(UInt i=0; i<nAtoms; i++) {
		UInt atomic = atomicVec[i];
		fprintf(fp, "%d \n", (Int)atomic);
	}
	fprintf(fp, "</Atomic Numbers>\n\n");

	// nuclei charge
	fprintf(fp, "<Nuclear Charges>\n");
	for(UInt i=0; i<nAtoms; i++) {
		UInt atomic = atomicVec[i];
		fprintf(fp, "%-12.6f \n", (Double)atomic);
	}
	fprintf(fp, "</Nuclear Charges>\n\n");

	// coordinates
	fprintf(fp, "<Nuclear Cartesian Coordinates>\n");
	for(UInt i=0; i<nAtoms; i++) {
		Double x = coords[3*i+0];
		Double y = coords[3*i+1];
		Double z = coords[3*i+2];
		fprintf(fp, "%-14.7f  %-14.7f  %-14.7f\n", x,y,z);
	}
	fprintf(fp, "</Nuclear Cartesian Coordinates>\n\n");

	// net charge
	fprintf(fp, "<Net Charge>\n");
	fprintf(fp, "%d \n", (Int)mol.getCharge());
	fprintf(fp, "</Net Charge>\n\n");

	// total number of electrons
	fprintf(fp, "<Number of Electrons>\n");
	fprintf(fp, "%d \n", (Int)mol.getTotalEle());
	fprintf(fp, "</Number of Electrons>\n\n");

	// alpha electrons
	fprintf(fp, "<Number of Alpha Electrons>\n");
	fprintf(fp, "%d \n", (Int)mol.getNEle(0));
	fprintf(fp, "</Number of Alpha Electrons>\n\n");

	// beta electrons
	fprintf(fp, "<Number of Beta Electrons>\n");
	fprintf(fp, "%d \n", (Int)mol.getNEle(1));
	fprintf(fp, "</Number of Beta Electrons>\n\n");

	// multiplicity
	fprintf(fp, "<Electronic Spin Multiplicity>\n");
	fprintf(fp, "%d \n", (Int)mol.getMult());
	fprintf(fp, "</Electronic Spin Multiplicity>\n\n");

	// model
	fprintf(fp, "<Model>\n");
	if (mol.getMult() == 1) {
		fprintf(fp, "Restricted DFT\n");
	}else{
		fprintf(fp, "Unrestricted DFT\n");
	}
	fprintf(fp, "</Model>\n\n");

	// now print out primitive centers
	fprintf(fp, "<Primitive Centers>\n");
	{
		const UInt nDataPerLine = 10;
		UInt nLeft   = nTotPrims%nDataPerLine;
		UInt nRounds = (nTotPrims-nLeft)/nDataPerLine;
		if (nRounds>0) {
			for(UInt i=0; i<nRounds; i++) {
				const UInt* c = &primCenters[i*nDataPerLine];
				fprintf(fp, "%-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d\n",
						(Int)c[0],(Int)c[1],(Int)c[2],(Int)c[3],(Int)c[4],(Int)c[5],
						(Int)c[6],(Int)c[7],(Int)c[8],(Int)c[9]);
			}
		}
		if (nLeft > 0) {
			const UInt* c = &primCenters[(nRounds)*nDataPerLine];
			for(UInt i=0; i<nLeft; i++) {
				fprintf(fp, "%-6d ",(Int)c[i]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "</Primitive Centers>\n\n");

	// now print out primitive centers
	fprintf(fp, "<Primitive Types>\n");
	{
		const UInt nDataPerLine = 10;
		UInt nLeft   = nTotPrims%nDataPerLine;
		UInt nRounds = (nTotPrims-nLeft)/nDataPerLine;
		if (nRounds>0) {
			for(UInt i=0; i<nRounds; i++) {
				const UInt* c = &primTypes[i*nDataPerLine];
				fprintf(fp, "%-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d %-6d\n",
						(Int)c[0],(Int)c[1],(Int)c[2],(Int)c[3],(Int)c[4],(Int)c[5],
						(Int)c[6],(Int)c[7],(Int)c[8],(Int)c[9]);
			}
		}
		if (nLeft > 0) {
			const UInt* c = &primTypes[(nRounds)*nDataPerLine];
			for(UInt i=0; i<nLeft; i++) {
				fprintf(fp, "%-6d ",(Int)c[i]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "</Primitive Types>\n\n");

	// now print out primitive centers
	fprintf(fp, "<Primitive Exponents>\n");
	{
		const UInt nDataPerLine = 3;
		UInt nLeft   = nTotPrims%nDataPerLine;
		UInt nRounds = (nTotPrims-nLeft)/nDataPerLine;
		if (nRounds>0) {
			for(UInt i=0; i<nRounds; i++) {
				const Double* e = &primExp[i*nDataPerLine];
				fprintf(fp, "%-16.10f %-16.10f %-16.10f\n", e[0],e[1],e[2]);
			}
		}
		if (nLeft > 0) {
			const Double* e = &primExp[(nRounds)*nDataPerLine];
			for(UInt i=0; i<nLeft; i++) {
				fprintf(fp, "%-16.10f ",e[i]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "</Primitive Exponents>\n\n");

	// mo occupation number
	UInt nSpin = mo.getNSpin();
	fprintf(fp, "<Molecular Orbital Occupation Numbers>\n");
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		UInt nMO = nAlphaOcc;
		if (iSpin == 1) nMO = nBetaOcc;
		for(UInt i=0; i<nMO; i++) {
			fprintf(fp, "%-16.10f\n", moOCCNumber[iSpin*nAlphaOcc+i]);
		}
	}
	fprintf(fp, "</Molecular Orbital Occupation Numbers>\n\n");
	
	// mo energy
	fprintf(fp, "<Molecular Orbital Energies>\n");
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		UInt nMO = nAlphaOcc;
		if (iSpin == 1) nMO = nBetaOcc;
		for(UInt i=0; i<nMO; i++) {
			fprintf(fp, "%-16.10f\n", moEnergy[iSpin*nAlphaOcc+i]);
		}
	}
	fprintf(fp, "</Molecular Orbital Energies>\n\n");
	
	// mo status
	fprintf(fp, "<Molecular Orbital Spin Types>\n");
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {
		string status = "Alpha and Beta";
		if (nSpin == 2) {
			if (iSpin == 0) {
				status = "Alpha";
			}else{
				status = "Beta";
			}
		}
		UInt nMO = nAlphaOcc;
		if (iSpin == 1) nMO = nBetaOcc;
		for(UInt i=0; i<nMO; i++) {
			fprintf(fp, "%-s\n", status.c_str());
		}
	}
	fprintf(fp, "</Molecular Orbital Spin Types>\n\n");
	
	// molecular orbital dependent information
	// here we need to loop over all possible mo
	DoubleVec primMO(nTotPrims);
	UInt MOIndex = 1;
	fprintf(fp, "<Molecular Orbital Primitive Coefficients>\n");
	for(UInt iSpin=0; iSpin<nSpin; iSpin++) {

		// mo coefficient based on the primitive functions
		{
			const Mtrx& M = moData.getMtrx(iSpin);
			for(UInt iMO=0; iMO<M.getCol(); iMO++) {

				// print out MO label
				fprintf(fp, "<MO Number>\n");
				fprintf(fp, "%-6d\n", (Int)MOIndex);
				fprintf(fp, "</MO Number>\n");
				MOIndex++;

				// form the mo data
				const Double* v = M.getPtr(0,iMO);
				UInt offset     = 0;
				UInt moOffset   = 0;
				for(UInt iAtom=0; iAtom<ms.getNAtomShells(); iAtom++) {
					const AtomShell& atomShell = ms.getAtomShell(iAtom);
					for(UInt iShell=0; iShell<atomShell.getNShell(); iShell++) {
						const Shell& s = atomShell.getShell(iShell);
						UInt lmin  = s.getLmin();
						UInt lmax  = s.getLmax();
						UInt nP    = s.getNPrim();
						for(UInt iL=lmin; iL<=lmax; iL++) {
							UInt nBas  = getCartBas(iL,iL);
							for(UInt iBas = 0; iBas<nBas; iBas++) {
								Double moCoeff = v[moOffset];
								for(UInt i = 0; i<nP; i++) {
									primMO[offset+i] = primCoe[offset+i]*moCoeff;
								}
								offset += nP;
								moOffset++;
							}
						}
					}
				}

				// now print out the mo coefficient
				const UInt nDataPerLine = 3;
				UInt nLeft   = nTotPrims%nDataPerLine;
				UInt nRounds = (nTotPrims-nLeft)/nDataPerLine;
				if (nRounds>0) {
					for(UInt i=0; i<nRounds; i++) {
						const Double* e = &primMO[i*nDataPerLine];
						fprintf(fp, "%-16.10f %-16.10f %-16.10f\n", e[0],e[1],e[2]);
					}
				}
				if (nLeft > 0) {
					const Double* e = &primMO[(nRounds)*nDataPerLine];
					for(UInt i=0; i<nLeft; i++) {
						fprintf(fp, "%-16.10f ",e[i]);
					}
					fprintf(fp, "\n");
				}
			}   
		}
	}
	fprintf(fp, "</Molecular Orbital Primitive Coefficients>\n\n");

	// system energy
	fprintf(fp, "<Energy = T + Vne + Vee + Vnn>\n");
	fprintf(fp, "%-16.10f\n", totalEnergy);
	fprintf(fp, "</Energy = T + Vne + Vee + Vnn>\n\n");
	
	// system energy
	fprintf(fp, "<Virial Ratio (-V/T)>\n");
	fprintf(fp, "%-16.10f\n", TWO);
	fprintf(fp, "</Virial Ratio (-V/T)>\n\n");

	// finally close the whole file 
	fclose(fp);
}

