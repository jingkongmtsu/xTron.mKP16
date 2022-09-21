#include "hgp_os_mom_p.h"
#include <cstdio>
#include <cassert>

void hgp_os_mom_p(const LInt& LCode, const UInt& inp2, const Double* icoe, const Double* iexp, 
		const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* C, Double* abcd)
{
  switch (LCode) {
    case 1004004:
      hgp_os_mom_g_g_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1000000:
      hgp_os_mom_s_s_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1000004:
      hgp_os_mom_g_s_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1100003:
      hgp_os_mom_f_sp_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1002003:
      hgp_os_mom_f_d_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1001004:
      hgp_os_mom_g_p_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1002002:
      hgp_os_mom_d_d_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1001003:
      hgp_os_mom_f_p_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1001001:
      hgp_os_mom_p_p_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1000100:
      hgp_os_mom_sp_s_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1100002:
      hgp_os_mom_d_sp_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1003003:
      hgp_os_mom_f_f_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1003004:
      hgp_os_mom_g_f_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1100001:
      hgp_os_mom_p_sp_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1000002:
      hgp_os_mom_d_s_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1001002:
      hgp_os_mom_d_p_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1002004:
      hgp_os_mom_g_d_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1000003:
      hgp_os_mom_f_s_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1100004:
      hgp_os_mom_g_sp_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1100100:
      hgp_os_mom_sp_sp_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    case 1000001:
      hgp_os_mom_p_s_p(inp2,icoe,iexp,ifac,P,A,B,C,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
  }

}
