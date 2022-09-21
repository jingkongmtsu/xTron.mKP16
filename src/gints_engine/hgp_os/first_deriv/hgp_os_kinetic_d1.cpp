#include "hgp_os_kinetic_d1.h"
#include <cstdio>
#include <cassert>

void hgp_os_kinetic_d1(const LInt& LCode, const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{

  switch (LCode) {
    case 2003:
      hgp_os_kinetic_f_d_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 3003:
      hgp_os_kinetic_f_f_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 2004:
      hgp_os_kinetic_g_d_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 100005:
      hgp_os_kinetic_h_sp_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 4005:
      hgp_os_kinetic_h_g_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 0:
      hgp_os_kinetic_s_s_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 100002:
      hgp_os_kinetic_d_sp_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 100004:
      hgp_os_kinetic_g_sp_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 100003:
      hgp_os_kinetic_f_sp_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 1:
      hgp_os_kinetic_p_s_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 1003:
      hgp_os_kinetic_f_p_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 4004:
      hgp_os_kinetic_g_g_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 5:
      hgp_os_kinetic_h_s_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 1005:
      hgp_os_kinetic_h_p_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 5005:
      hgp_os_kinetic_h_h_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 3004:
      hgp_os_kinetic_g_f_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 3:
      hgp_os_kinetic_f_s_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 2005:
      hgp_os_kinetic_h_d_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 4:
      hgp_os_kinetic_g_s_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 100001:
      hgp_os_kinetic_p_sp_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 3005:
      hgp_os_kinetic_h_f_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 100:
      hgp_os_kinetic_sp_s_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 1002:
      hgp_os_kinetic_d_p_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 2002:
      hgp_os_kinetic_d_d_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 100100:
      hgp_os_kinetic_sp_sp_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 2:
      hgp_os_kinetic_d_s_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 1001:
      hgp_os_kinetic_p_p_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    case 1004:
      hgp_os_kinetic_g_p_d1(inp2,icoe,iexp,iexpdiff,ifac,P,A,B,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
  }

}
