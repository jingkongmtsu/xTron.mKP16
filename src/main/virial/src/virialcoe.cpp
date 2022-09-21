/**
 * cpp file corresponding to virialcoe.h
 */
#include "globalinfor.h"
#include "molecule.h"
#include "cluster.h"
#include "shell.h"
#include "denmtrx.h"
#include "scf.h"
#include "virialcoe.h"
using namespace globalinfor;
using namespace molecule;
using namespace cluster;
using namespace shell;
using namespace denmtrx;
using namespace scf;
using namespace virialcoe;

void VirialCoe::monomerCalculation(const GlobalInfor& infor)
{
	//
	// do we initilize the data
	//
	if (param.getMethod() == PERTURBED_SCF_WAY_VIRIAL) {
		moslist.reserve(param.getNMono());
	}

	//
	// now go into the fragment calculation
	//
	string inf = infor.getInputFile();
	for(Int sec=1; sec<=param.getNMono(); sec++) {

		// form geom and shell
		Molecule mol(inf,sec);
		MolShell ms(inf,mol);

		// perform SCF calculation
		SCF scf(infor,mol,ms);
		DenMtrx den(ms,ms,scf.getNSpin());
		den.scfGuess(scf.getSCFParam(),scf.getOneEMtrx(),ms,mol);
		scf.doSCF(infor,mol,ms,den);

		// push the result MO back to the moslist
		// prepare this data for purturbed calculation
		if (param.getMethod() == PERTURBED_SCF_WAY_VIRIAL) {
			// finally, let's note that density matrix is only 
			// related to the occupied mos. In the following 
			// density matrix we will do a transformation to make 
			// the MO idempotent, therefore only occupied MO
			// is involved
			MOs mos(scf.getMOs());
			for(Int iSpin=0; iSpin<mos.getNSpin(); iSpin++) {
				MO& mo = mos.getMO(iSpin);
				Int nocc = mo.getNOcc();
				mo.resetMO(nocc);
			}
			moslist.push_back(mos);
		}

		// add e into result
		Double e = scf.totalEnergy();
		eMonomerSum += e;
		cout << "monomer energy is " << e << endl;
	}
}

void VirialCoe::clusterCalculation(const GlobalInfor& infor)
{

	// obtain the cluster and its shell data
	string input = infor.getInputFile();
	Cluster cluster(input);
	MolShell ms(input,cluster);

	// get the density matrix first
	SCF scf(infor,cluster,ms);
	DenMtrx den(ms,ms,scf.getNSpin());
	Int method = param.getMethod();
	if (method == PERTURBED_SCF_WAY_VIRIAL) {
		den.monomerFragmentGuess(scf.getOneEMtrx(),ms,cluster,moslist);
	}else if (method == FULL_SCF_WAY_VIRIAL){
		den.scfGuess(scf.getSCFParam(),scf.getOneEMtrx(),ms,cluster);
	}

	// perform SCF calculation
	scf.doSCF(infor,cluster,ms,den);
	eCluster = scf.totalEnergy();
	printf("total energy plus vdw: %-16.12f\n", eCluster);
}

VirialCoe::VirialCoe(const string& input):eMonomerSum(ZERO),eCluster(ZERO),param(input)
{
	GlobalInfor infor(input);
	monomerCalculation(infor);
	clusterCalculation(infor);
}
