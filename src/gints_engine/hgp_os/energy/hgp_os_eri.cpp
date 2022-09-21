#include "hgp_os_eri.h"
#include <cstdio>
#include <cassert>

void hgp_os_eri(const LInt& LCode, const UInt& maxL, const UInt& inp2, const UInt& jnp2, const Double& pMax, const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd, LocalMemScr& scr)
{
  if (maxL<=1) {
    switch(LCode) {
    case 1000001:
      hgp_os_eri_p_s_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100:
      hgp_os_eri_sp_s_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001001001:
      hgp_os_eri_p_p_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001001:
      hgp_os_eri_p_p_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1:
      hgp_os_eri_p_s_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100001:
      hgp_os_eri_p_sp_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 0:
      hgp_os_eri_s_s_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001001:
      hgp_os_eri_p_p_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001:
      hgp_os_eri_p_sp_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100100001:
      hgp_os_eri_p_sp_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001001001:
      hgp_os_eri_p_p_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100:
      hgp_os_eri_sp_sp_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100100:
      hgp_os_eri_sp_sp_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001100001:
      hgp_os_eri_p_sp_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100000001:
      hgp_os_eri_p_s_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100100100:
      hgp_os_eri_sp_sp_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100001001:
      hgp_os_eri_p_p_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100000100:
      hgp_os_eri_sp_s_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1100001:
      hgp_os_eri_p_sp_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001:
      hgp_os_eri_p_p_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1100100:
      hgp_os_eri_sp_sp_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
    }

  }else if (maxL==2) {
    switch(LCode) {
    case 100001002002:
      hgp_os_eri_d_d_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001100002:
      hgp_os_eri_d_sp_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2000002:
      hgp_os_eri_d_s_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100001002:
      hgp_os_eri_d_p_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100002002:
      hgp_os_eri_d_d_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001002002:
      hgp_os_eri_d_d_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001001002:
      hgp_os_eri_d_p_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002100002:
      hgp_os_eri_d_sp_d_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1000002:
      hgp_os_eri_d_s_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001002:
      hgp_os_eri_d_p_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2100100:
      hgp_os_eri_sp_sp_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001100002:
      hgp_os_eri_d_sp_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002002:
      hgp_os_eri_d_d_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2001002:
      hgp_os_eri_d_p_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1100002:
      hgp_os_eri_d_sp_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002:
      hgp_os_eri_d_sp_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2:
      hgp_os_eri_d_s_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2100002:
      hgp_os_eri_d_sp_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002:
      hgp_os_eri_d_p_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002002:
      hgp_os_eri_d_d_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002001002:
      hgp_os_eri_d_p_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002001002:
      hgp_os_eri_d_p_d_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002002002:
      hgp_os_eri_d_d_d_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2001001:
      hgp_os_eri_p_p_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100000002:
      hgp_os_eri_d_s_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002002002:
      hgp_os_eri_d_d_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002002:
      hgp_os_eri_d_d_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2100001:
      hgp_os_eri_p_sp_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002002002:
      hgp_os_eri_d_d_d_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100002:
      hgp_os_eri_d_sp_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002:
      hgp_os_eri_d_d_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001002:
      hgp_os_eri_d_p_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100100002:
      hgp_os_eri_d_sp_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001001002:
      hgp_os_eri_d_p_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
    }

  }else if (maxL==3) {
    switch(LCode) {
    case 100001100003:
      hgp_os_eri_f_sp_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100003:
      hgp_os_eri_f_sp_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002001003:
      hgp_os_eri_f_p_d_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001003003:
      hgp_os_eri_f_f_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3100003:
      hgp_os_eri_f_sp_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003:
      hgp_os_eri_f_d_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003003:
      hgp_os_eri_f_f_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003003003:
      hgp_os_eri_f_f_f_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3:
      hgp_os_eri_f_s_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100100003:
      hgp_os_eri_f_sp_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001002003:
      hgp_os_eri_f_d_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100000003:
      hgp_os_eri_f_s_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003003:
      hgp_os_eri_f_f_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002001003:
      hgp_os_eri_f_p_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003002002:
      hgp_os_eri_d_d_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002003:
      hgp_os_eri_f_d_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003003003:
      hgp_os_eri_f_f_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3000003:
      hgp_os_eri_f_s_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2001003:
      hgp_os_eri_f_p_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002003003:
      hgp_os_eri_f_f_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3001002:
      hgp_os_eri_d_p_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2100003:
      hgp_os_eri_f_sp_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002002003:
      hgp_os_eri_f_d_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001000003:
      hgp_os_eri_f_s_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001003:
      hgp_os_eri_f_p_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100003002002:
      hgp_os_eri_d_d_f_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3003:
      hgp_os_eri_f_f_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3003003003:
      hgp_os_eri_f_f_f_f(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002100003:
      hgp_os_eri_f_sp_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3002003:
      hgp_os_eri_f_d_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003001003:
      hgp_os_eri_f_p_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001003003:
      hgp_os_eri_f_f_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002100003:
      hgp_os_eri_f_sp_d_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100003001003:
      hgp_os_eri_f_p_f_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3002002:
      hgp_os_eri_d_d_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002003:
      hgp_os_eri_f_d_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100003002003:
      hgp_os_eri_f_d_f_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100003100003:
      hgp_os_eri_f_sp_f_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002003:
      hgp_os_eri_f_d_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100003003:
      hgp_os_eri_f_f_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001001003:
      hgp_os_eri_f_p_p_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100001003:
      hgp_os_eri_f_p_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001002003:
      hgp_os_eri_f_d_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001001003:
      hgp_os_eri_f_p_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1100003:
      hgp_os_eri_f_sp_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003002003:
      hgp_os_eri_f_d_f_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100003003003:
      hgp_os_eri_f_f_f_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002002003:
      hgp_os_eri_f_d_d_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100000003:
      hgp_os_eri_f_s_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2000003:
      hgp_os_eri_f_s_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100002003:
      hgp_os_eri_f_d_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3001003:
      hgp_os_eri_f_p_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001000003:
      hgp_os_eri_f_s_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100003003:
      hgp_os_eri_f_f_sp_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002003003:
      hgp_os_eri_f_f_d_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100002003003:
      hgp_os_eri_f_f_d_sp(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001100003:
      hgp_os_eri_f_sp_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002002003:
      hgp_os_eri_f_d_d_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003:
      hgp_os_eri_f_p_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100001003:
      hgp_os_eri_f_p_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 100100003:
      hgp_os_eri_f_sp_sp_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3003003:
      hgp_os_eri_f_f_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3100002:
      hgp_os_eri_d_sp_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003002003:
      hgp_os_eri_f_d_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1000003:
      hgp_os_eri_f_s_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
    }

  }else if (maxL==4) {
    switch(LCode) {
    case 1002002004:
      hgp_os_eri_g_d_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3004004004:
      hgp_os_eri_g_g_g_f(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 1001004:
      hgp_os_eri_g_p_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003003004:
      hgp_os_eri_g_f_f_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003002004:
      hgp_os_eri_g_d_f_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002001004:
      hgp_os_eri_g_p_d_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003001004:
      hgp_os_eri_g_p_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004:
      hgp_os_eri_g_p_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001000004:
      hgp_os_eri_g_s_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001003004:
      hgp_os_eri_g_f_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4002003:
      hgp_os_eri_f_d_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2001004:
      hgp_os_eri_g_p_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4004:
      hgp_os_eri_g_g_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001004004:
      hgp_os_eri_g_g_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4004004004:
      hgp_os_eri_g_g_g_g(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 2004004004:
      hgp_os_eri_g_g_g_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 4003004:
      hgp_os_eri_g_f_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004001004:
      hgp_os_eri_g_p_g_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3002004:
      hgp_os_eri_g_d_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4:
      hgp_os_eri_g_s_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003002004:
      hgp_os_eri_g_d_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3003004004:
      hgp_os_eri_g_g_f_f(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 1002000004:
      hgp_os_eri_g_s_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002003004:
      hgp_os_eri_g_f_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2000004:
      hgp_os_eri_g_s_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2004003004:
      hgp_os_eri_g_f_g_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 2004004:
      hgp_os_eri_g_g_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002004004:
      hgp_os_eri_g_g_d_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001002004:
      hgp_os_eri_g_d_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3004:
      hgp_os_eri_g_f_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2004:
      hgp_os_eri_g_d_s_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4002002:
      hgp_os_eri_d_d_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002004:
      hgp_os_eri_g_d_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002002004:
      hgp_os_eri_g_d_d_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4000004:
      hgp_os_eri_g_s_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4002004:
      hgp_os_eri_g_d_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004002003:
      hgp_os_eri_f_d_g_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003004004:
      hgp_os_eri_g_g_f_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 3004004:
      hgp_os_eri_g_g_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2003004:
      hgp_os_eri_g_f_d_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004004004:
      hgp_os_eri_g_g_g_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004002004:
      hgp_os_eri_g_d_g_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002004004:
      hgp_os_eri_g_g_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004003004:
      hgp_os_eri_g_f_g_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004003003:
      hgp_os_eri_f_f_g_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2004002004:
      hgp_os_eri_g_d_g_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4001004:
      hgp_os_eri_g_p_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1000004:
      hgp_os_eri_g_s_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3003003004:
      hgp_os_eri_g_f_f_f(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 1003003004:
      hgp_os_eri_g_f_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4001003:
      hgp_os_eri_f_p_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003004:
      hgp_os_eri_g_f_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4003003:
      hgp_os_eri_f_f_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 4004004:
      hgp_os_eri_g_g_g_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3004003004:
      hgp_os_eri_g_f_g_f(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd,scr);
      break;
    case 3001004:
      hgp_os_eri_g_p_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3003004:
      hgp_os_eri_g_f_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2004003003:
      hgp_os_eri_f_f_g_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002001004:
      hgp_os_eri_g_p_d_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 3000004:
      hgp_os_eri_g_s_f_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1003004004:
      hgp_os_eri_g_g_f_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1002004:
      hgp_os_eri_g_d_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 2002003004:
      hgp_os_eri_g_f_d_d(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1004004:
      hgp_os_eri_g_g_p_s(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    case 1001001004:
      hgp_os_eri_g_p_p_p(inp2,jnp2,pMax,omega,icoe,iexp,ifac,P,A,B,jcoe,jexp,jfac,Q,C,D,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
    }

  }
}
