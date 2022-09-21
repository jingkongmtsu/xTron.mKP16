#ifndef  HGP_OS_TWOBODYOVERLAP_H
#define HGP_OS_TWOBODYOVERLAP_H

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

void hgp_os_twobodyoverlap_g_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_d_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_f_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_s_s(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_g_d(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_sp_s(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_h_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_h_d(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_h_h(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_f_d(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_f_f(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_h_f(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_h_s(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_g_f(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_d_s(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_f_sp(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_p_s(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_h_g(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_f_s(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_d_sp(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_h_sp(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_d_d(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_g_s(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_p_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_g_sp(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_p_sp(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_g_g(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

void hgp_os_twobodyoverlap_sp_sp(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);

#endif
