#ifndef  HGP_OS_NAI_H
#define HGP_OS_NAI_H

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

void hgp_os_nai_p_sp(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_h_d(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_h_s(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_d_s(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_g_sp(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_d_sp(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_h_h(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_p_p(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_f_s(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_h_p(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_g_d(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_sp_sp(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_h_f(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_f_sp(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_f_p(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_f_f(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_h_g(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_f_d(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_g_s(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_g_f(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_h_sp(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_sp_s(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_p_s(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_d_p(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_g_g(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_s_s(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_d_d(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

void hgp_os_nai_g_p(const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

#endif
