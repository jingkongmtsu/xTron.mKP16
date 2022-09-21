#ifndef  HGP_OS_MOM_P_H
#define HGP_OS_MOM_P_H

#include <cstddef>
#include <cassert>
#include <cstdio>

typedef int             Int; 
typedef long long       LInt; 
typedef size_t          UInt; 
#ifdef WITH_SINGLE_PRECISION 
typedef float           Double;
#else
typedef double          Double;
#endif

void hgp_os_mom_g_g_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_s_s_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_g_s_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_f_sp_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_f_d_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_g_p_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_d_d_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_f_p_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_p_p_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_sp_s_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_d_sp_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_f_f_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_g_f_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_p_sp_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_d_s_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_d_p_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_g_d_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_f_s_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_g_sp_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_sp_sp_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

void hgp_os_mom_p_s_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd);

#endif
