#include "hgp_os_nai.h"
#include <cstdio>
#include <cassert>

void hgp_os_nai(const LInt& LCode, const UInt& inp2, const UInt& nAtoms, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd)
{

  switch (LCode) {
    case 100001:
      hgp_os_nai_p_sp(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 2005:
      hgp_os_nai_h_d(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 5:
      hgp_os_nai_h_s(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 2:
      hgp_os_nai_d_s(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 100004:
      hgp_os_nai_g_sp(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 100002:
      hgp_os_nai_d_sp(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 5005:
      hgp_os_nai_h_h(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 1001:
      hgp_os_nai_p_p(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 3:
      hgp_os_nai_f_s(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 1005:
      hgp_os_nai_h_p(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 2004:
      hgp_os_nai_g_d(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 100100:
      hgp_os_nai_sp_sp(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 3005:
      hgp_os_nai_h_f(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 100003:
      hgp_os_nai_f_sp(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 1003:
      hgp_os_nai_f_p(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 3003:
      hgp_os_nai_f_f(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 4005:
      hgp_os_nai_h_g(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 2003:
      hgp_os_nai_f_d(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 4:
      hgp_os_nai_g_s(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 3004:
      hgp_os_nai_g_f(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 100005:
      hgp_os_nai_h_sp(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 100:
      hgp_os_nai_sp_s(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 1:
      hgp_os_nai_p_s(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 1002:
      hgp_os_nai_d_p(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 4004:
      hgp_os_nai_g_g(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 0:
      hgp_os_nai_s_s(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 2002:
      hgp_os_nai_d_d(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    case 1004:
      hgp_os_nai_g_p(inp2,nAtoms,icoe,iexp,ifac,P,A,B,N,Z,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
  }

}
