#ifndef FUNCTIONALLIST_H
#define FUNCTIONALLIST_H
#include "libgen.h"

#define functional_s           functional_s_
#define functional_vwn5        functional_vwn5_
#define functional_vwn1rpa     functional_vwn1rpa_
#define functional_pw92c       functional_pw92c_

#define functional_pw86x       functional_pw86x_
#define functional_pw86c       functional_pw86c_
#define functional_becke88     functional_becke88_
#define functional_pbex        functional_pbex_
#define functional_pbec        functional_pbec_
#define functional_pw91x       functional_pw91x_
#define functional_pw91c       functional_pw91c_
#define functional_lyp         functional_lyp_
#define functional_tpssx       functional_tpssx_
#define functional_tpssc       functional_tpssc_
#define functional_scanx       functional_scanx_
#define functional_scanc       functional_scanc_

#define br89b2_vdwx            br89b2_vdwx_ 
#define br94coor_op            br94coor_op_ 
#define br94coor_par           br94coor_par_ 
#define br89hole               br89hole_ 
#define br89xx                 br89xx_ 
#define becke05_ndop           becke05_ndop_
#define becke05_ndpar          becke05_ndpar_
#define becke05_odd_electron   becke05_odd_electron_
#define beck_rsc               beck_rsc_
#define b13coor_opp            b13coor_opp_
#define b13coor_par            b13coor_par_
#define b13strong_ac2          b13strong_ac2_
#define b13strong_ac3          b13strong_ac3_
#define kp14ec                 kp14ec_
#define kp14exc                kp14exc_
#define functional_psts        functional_psts_

// For hole functions.
#define br89xhole_s            br89xhole_s_
#define hirao00xhole_s         hirao00xhole_s_


#define find_diff_check_deriv1_var    find_diff_check_deriv1_var_
#define find_diff_check_deriv1_array  find_diff_check_deriv1_array_

/**
 * C wrapper function for all of fortran files doing derivatives calculation
 */
extern "C" {

	/**
	 * the files which calculates the functional value and 1st derivatives
	 */

	/**
	 * LDA
	 */
	void functional_s(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, 
			double* F, double* D1F);

	void functional_vwn5(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, 
			double* F, double* D1F);

	void functional_vwn1rpa(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, 
			double* F, double* D1F);

	void functional_pw92c(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, 
			double* F, double* D1F);

	/**
	 * GGA
	 */
	void functional_pw86x(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, const double* DRA, 
			const double* DRB, double* F, double* D1F);

	void functional_pw86c(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* rhoA, const double* rhoB, const double* GAA, 
			const double* GAB, const double* GBB, double* F, double* D1F);

	void functional_becke88(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, const double* DRA, 
			const double* DRB, double* F, double* D1F);

	void functional_pbex(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, const double* DRA, 
			const double* DRB, double* F, double* D1F);

	void functional_pbec(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, const double* DRA, 
			const double* DRB, double* F, double* D1F);

	void functional_pw91x(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, const double* DRA, 
			const double* DRB, double* F, double* D1F);

	void functional_lyp(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, const double* DRA, 
			const double* DRB, double* F, double* D1F);

	void functional_pw91c(const Int* infor, const Int* NG, const Int* NDEN, 
			const double* TOL, const double* rhoA, const double* rhoB, const double* DRA, 
			const double* DRB, double* F, double* D1F);


	/**
	 * META-GGA
	 */
	void br89b2_vdwx(double* TF, double* TD1F,const Int* infor, const double* rhoA,
			const double* rhoB, const double* DRA, const double* DRB,
			const double* LA, const double* LB, const double* TA,
			const double* TB, const Int* ng, const Int* nDen);

	void br94coor_op(const Int* INFOR, const Int* NDEN,const Int* NA,const Int* NB,const Int* NG,
			const double* THRESH, 
                        const double* VAL_P, const double* HIRSHWTS, 
                        const double* RhoA,const double* RhoB, const double* DRA,
			const double* DRB, const double* TauA,const double* TauB,const double* LapA,
			const double* LapB, const double* EXA, const double* EXB, double* F,double* D1F);

	void br94coor_par(const Int* INFOR, const Int* NDEN,const Int* NA,const Int* NB,const Int* NG,
			const double* THRESH, 
                        const double* VAL_P, const double* VAL_Q, const double* HIRSHWTS, 
                        const double* RhoA,const double* RhoB, const double* DRA,
			const double* DRB, const double* TauA,const double* TauB,const double* LapA,
			const double* LapB, const double* EXA, const double* EXB, double* F,double* D1F);

	void br89hole(const double* THRESH, const double* Rho, const double* Gam, const double* Tau, 
			const double* Lap, double* U, double* D1R, double* D1G, double* D1T, double* D1L);

	void br89xx(double* F, double* D1F, const Int* INFOR, const double* RhoA, const double* RhoB,
			const double* DRA, const double* DRB, const double* LapA, const double* LapB, 
			const double* TauA, const double* TauB, const Int* ND, const Int* NDEN);

	void functional_tpssx(const Int* INFOR, const Int* NG, const Int* NDEN, const double* THRESH, 
			const double* RhoA,const double* RhoB, const double* DRA,
			const double* DRB, const double* TauA,const double* TauB, 
			double* F,double* D1F);

	void functional_tpssc(const Int* INFOR, const Int* NG, const Int* NDEN, const double* THRESH, 
			const Int* USE_LAMBDA, const double* RhoA,const double* RhoB, const double* DRA,
			const double* DRB, const double* TauA,const double* TauB, 
			double* F,double* D1F);

	void functional_scanx(const Int* INFOR, const Int* NG, const Int* NDEN, const double* THRESH, 
			const double* RhoA,const double* RhoB, const double* DRA,
			const double* DRB, const double* TauA,const double* TauB, 
			double* F,double* D1F);

	void functional_scanc(const Int* INFOR, const Int* NG, const Int* NDEN, const double* THRESH, 
			const Int* USE_LAMBDA, const double* RhoA,const double* RhoB, const double* DRA,
			const double* DRB, const double* TauA,const double* TauB, 
			double* F,double* D1F);

	/**
	 * HYPER-FUNCTIONAL
	 * THESE WHO CONTAINS THE EXCHANGE ENERGY DENSITY
	 */
	void becke05_odd_electron(const Int* NG,const Int* NA, const Int* NB, const double* THRESH, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* LapA, const double* LapB, 
			const double* EXA, const double* EXB, const double* a1, const double* a2, double* F);

	void becke05_ndop(const Int* INFOR,const Int* NDEN,const Int* NG,const double* THRESH, 
			const double* VAL_P, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* LapA, const double* LapB, 
			const double* EXA, const double* EXB, double* F, double* D1F);

	void becke05_ndpar(const Int* INFOR,const Int* NDEN,const Int* NG, const Int* NA,
			const Int* NB, const double* THRESH, 
			const double* VAL_P, const double* VAL_Q, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* LapA, const double* LapB, 
			const double* EXA, const double* EXB, double* F, double* D1F);

   void beck_rsc(double* F, double* D1F, const double* D1FEXCH, const double* W, 
			const double* RA, const double* RB, const double* D1RA, const double* D1RB, 
			const double* RLapA, const double* RLapB, const double* TA, const double* TB, 
			const double* EX_HF_DENSA, const double* EX_HF_DENSB, const double* ACOP, const double* ACPAR, 
			const double* ODEL, const double* ODELW, const double* ODDSUM, const Int* NGrid, const Int* iterSCF, 
			const Int* ICEN, const Int* NA, const Int* NB);

   void kp14exc(const Int* iFunc, double* F, double* D1F, const double* D1FEXCH, const double* W, 
			const double* RA, const double* RB, const double* D1RA, const double* D1RB, 
			const double* RLapA, const double* RLapB, const double* TA, const double* TB, 
			const double* EX_HF_DENSA, const double* EX_HF_DENSB, 
			const double* ODEL, const double* ODELW, const double* ODDSUM, const Int* NGrid, const Int* iterSCF, 
			const Int* ICEN, const Int* NA, const Int* NB);

	void kp14ec(const Int* NDPARMETHOD, const Int* INFOR, const Int* NDEN, const Int* NG, const Int* NA, const Int* NB, const double* b,
			const double* alpha,const double* c_ndpar, const double* c_ndpar_cap, const double* THRESH, 
			const double* VAL_P, const double* VAL_Q, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB,
			const double* DRA, const double* DRB, const double* TauA, const double* TauB,
			const double* LapA, const double* LapB, const double* EXA, const double* EXB,
			double* F, double* T1F, double* D1F);

	void b13coor_par(const Int* method, const Int* INFOR,const Int* NDEN,
			const Int* NA, const Int* NB, const Int* NG, const double* THRESH, 
			const double* VAL_P, const double* VAL_Q, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* LapA, const double* LapB, 
			const double* EXA, const double* EXB, double* F, double* D1F);

	void b13coor_opp(const Int* method, const Int* INFOR,const Int* NDEN,
			const Int* NA, const Int* NB, const Int* NG, const double* THRESH, 
			const double* VAL_P, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* LapA, const double* LapB, 
			const double* EXA, const double* EXB, double* F, double* D1F);

	void b13strong_ac2(const Int* INFOR,const Int* NDEN,
			const Int* NG, const Int* NA, const Int* NB, const double* THRESH, 
			const double* VAL_P, const double* VAL_Q, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* LapA, const double* LapB, 
			const double* EXA, const double* EXB, double* F, double* D1F);

	void b13strong_ac3(const Int* INFOR,const Int* NDEN,
			const Int* NG, const Int* NA, const Int* NB, const double* THRESH, 
			const double* VAL_P, const double* VAL_Q, const double* HIRSHWTS,
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* LapA, const double* LapB, 
			const double* EXA, const double* EXB, double* F, double* D1F);

	void functional_psts(const Int* INFOR,const Int* NG,const Int* NDEN,const double* THRESH, 
			const double* RhoA, const double* RhoB, const double* DRA, const double* DRB, 
			const double* TauA, const double* TauB, const double* EXA, const double* EXB, 
			double* F, double* D1F);


	/*
	 * Hole functions.
	 */
	void br89xhole_s(double* hx, const double* sVal, const double* RhoA, 
	    const double* RhoB, const double* DRA, const double* DRB, const double* LapA, 
	    const double* LapB, const double* TauA, const double* TauB, const Int* Ng, 
	    const Int* NDEN);
	void hirao00xhole_s(double* hx, const double* sVal, const double* RhoA, 
	    const double* RhoB, const double* DRA, const double* DRB, const double* LapA, 
	    const double* LapB, const double* TauA, const double* TauB, const Int* Ng, 
	    const Int* NDEN);


	/*
	 * this is used only for debugging purpose
	 */
	void find_diff_check_deriv1_var(const Int* INFOR, const Int* NG, 
			const double* THRESH, const double* DELTA, const double* CRITERIA, 
			const double* RhoA, const double* RhoB,
			const double* DRA, const double* DRB, const double* TauA, const double* TauB,
			const double* LapA, const double* LapB, const double* EXA, const double* EXB);

	void find_diff_check_deriv1_array(const Int* INFOR, const Int* NG, 
			const double* THRESH, const double* DELTA, const double* CRITERIA, 
			const double* RhoA, const double* RhoB,
			const double* DRA, const double* DRB, const double* TauA, const double* TauB,
			const double* LapA, const double* LapB, const double* EXA, const double* EXB);

}

#endif
