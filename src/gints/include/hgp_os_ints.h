/**
 * \file   hgp_os_ints.h
 * \brief  integral function interface for HGP_OS engine (see engine/hgp_os)
 * \author Fenglai Liu 
 */
#ifndef HGP_OS_INTS_H
#define HGP_OS_INTS_H
#include "libgen.h"

namespace localmemscr {
	class LocalMemScr;
}

using namespace localmemscr;

////////////////////////////////////
//       energy section           //
////////////////////////////////////

///
/// interface function for ESP calculation
///
extern void hgp_os_esp(const LInt& LCode, const UInt& inp2, const UInt& nGrids, 
		const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* R, Double* abcd);

///
/// interface function for kinetic calculation
///
extern void hgp_os_kinetic(const LInt& LCode, const UInt& inp2, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, Double* abcd);

///
/// interface function for NAI calculation
///
extern void hgp_os_nai(const LInt& LCode, const UInt& inp2, const UInt& nAtoms, 
		const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

///
/// interface function for two body overlap calculation
///
extern void hgp_os_twobodyoverlap(const LInt& LCode, const UInt& inp2, const Double* icoe, 
		const Double* iexp, const Double* ifac, const Double* P, const Double* A, 
		const Double* B, Double* abcd);

///
/// interface function for ERI calculation
///
extern void hgp_os_eri(const LInt& LCode, const UInt& maxL, const UInt& inp2, const UInt& jnp2, 
		const Double& pMax, const Double& omega, const Double* icoe, 
		const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, 
		const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, 
		const Double* D, Double* abcd, LocalMemScr& scr);

///
/// interface function for the MOM integral in terms of P shell
/// this is actually for diple moment calculation
///
extern void hgp_os_mom_p(const LInt& LCode, const UInt& inp2, const Double* icoe, const Double* iexp, 
		const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

////////////////////////////////////
//    first derivative section    //
////////////////////////////////////

///
/// 1st order derivatives for ESP
///
extern void hgp_os_esp_d1(const LInt& LCode, const UInt& inp2, const UInt& nGrids, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, 
		const Double* B, const Double* R, Double* abcd);

///
/// 1st order derivatives for kinetic
///
extern void hgp_os_kinetic_d1(const LInt& LCode, const UInt& inp2, const Double* icoe, const Double* iexp, 
		const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, 
		Double* abcd);

///
/// 1st order derivatives for NAI
///
extern void hgp_os_nai_d1(const LInt& LCode, const UInt& inp2, const UInt& nAtoms, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, 
		const Double* B, const Double* N, const UInt* Z, Double* abcd);

///
/// 1st order derivatives for two body overlap
///
extern void hgp_os_twobodyoverlap_d1(const LInt& LCode, const UInt& inp2, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, Double* abcd);

///
/// 1st order derivatives for ERI
///
extern void hgp_os_eri_d1(const LInt& LCode, const UInt& maxL, const UInt& inp2, const UInt& jnp2, 
		const Double& pMax, const Double& omega, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jexpdiff, 
		const Double* jfac, const Double* Q, const Double* C, const Double* D, 
		Double* abcd, LocalMemScr& scr);

////////////////////////////////////
//   second derivative section    //
////////////////////////////////////

///
/// 2ed order derivatives for ESP
///
extern void hgp_os_esp_d2(const LInt& LCode, const UInt& inp2, const UInt& nGrids, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, 
		const Double* B, const Double* R, Double* abcd);

///
/// 2ed order derivatives for kinetic
///
extern void hgp_os_kinetic_d2(const LInt& LCode, const UInt& inp2, const Double* icoe, const Double* iexp, 
		const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, 
		const Double* B, Double* abcd);

///
/// 2ed order derivatives for NAI
///
extern void hgp_os_nai_d2(const LInt& LCode, const UInt& inp2, const UInt& nAtoms, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, 
		const Double* B, const Double* N, const UInt* Z, Double* abcd);

///
/// 2ed order derivatives for two body overlap
///
extern void hgp_os_twobodyoverlap_d2(const LInt& LCode, const UInt& inp2, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, Double* abcd);

///
/// 2ed order derivatives for ERI
///
extern void hgp_os_eri_d2(const LInt& LCode, const UInt& maxL, const UInt& inp2, const UInt& jnp2, 
		const Double& pMax, const Double& omega, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* jcoe, const Double* jexp, 
		const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, 
		const Double* D, Double* abcd, LocalMemScr& scr);

#endif

