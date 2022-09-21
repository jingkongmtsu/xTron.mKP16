#include "hgp_os_esp_d2.h"
#include <cstdio>
#include <cassert>

void hgp_os_esp_d2(const LInt& LCode, const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{

  switch (LCode) {
    case 2003:
      hgp_os_esp_f_d_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 0:
      hgp_os_esp_s_s_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 100004:
      hgp_os_esp_g_sp_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 3:
      hgp_os_esp_f_s_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 3004:
      hgp_os_esp_g_f_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 4:
      hgp_os_esp_g_s_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 1002:
      hgp_os_esp_d_p_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 1004:
      hgp_os_esp_g_p_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 2004:
      hgp_os_esp_g_d_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 3003:
      hgp_os_esp_f_f_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 5005:
      hgp_os_esp_h_h_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 100001:
      hgp_os_esp_p_sp_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 2005:
      hgp_os_esp_h_d_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 100100:
      hgp_os_esp_sp_sp_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 100002:
      hgp_os_esp_d_sp_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 4004:
      hgp_os_esp_g_g_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 1003:
      hgp_os_esp_f_p_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 1005:
      hgp_os_esp_h_p_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 2002:
      hgp_os_esp_d_d_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 100005:
      hgp_os_esp_h_sp_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 100:
      hgp_os_esp_sp_s_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 100003:
      hgp_os_esp_f_sp_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 3005:
      hgp_os_esp_h_f_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 2:
      hgp_os_esp_d_s_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 1:
      hgp_os_esp_p_s_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 4005:
      hgp_os_esp_h_g_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 1001:
      hgp_os_esp_p_p_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    case 5:
      hgp_os_esp_h_s_d2(inp2,nGrids,icoe,iexp,iexpdiff,ifac,P,A,B,R,abcd);
      break;
    default:
      printf("%s %lld\n","Un-recognized LCode in the integrals calculation ", LCode);
      assert(0);
      break;
  }

}
