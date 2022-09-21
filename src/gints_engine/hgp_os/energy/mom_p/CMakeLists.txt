#
# build the analytical integral engine for energy of mom
#
set(SRC_GINTS_ENGINE_MOM_P_D0
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/hgp_os_mom_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_d_d_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_d_p_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_d_s_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_d_sp_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_f_d_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_f_f_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_f_p_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_f_s_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_f_sp_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_g_d_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_g_f_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_g_g_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_g_p_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_g_s_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_g_sp_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_p_p_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_p_s_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_p_sp_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_s_s_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_sp_s_p.cpp
   ${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy/mom_p/hgp_os_mom_sp_sp_p.cpp
   )
add_library(gints_engine_mom_p_d0 ${SRC_GINTS_ENGINE_MOM_P_D0})
