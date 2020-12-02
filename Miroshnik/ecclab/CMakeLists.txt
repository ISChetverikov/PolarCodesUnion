cmake_minimum_required(VERSION 3.13)
project(ecclib)

set(SIM_BG simulators/main_txt.c simulators/sim_bg.c common/ui_txt.c common/spf_par.c)
set(SIM_SRM_BG common/srm_utils.c ${SIM_BG})

add_compile_options(-DUSE_STDLIB_RND)

add_executable(ca_polar_scl_bg polar_scl/ca_polar_scl_main.c polar_scl/polar_scl_inner.c ${SIM_SRM_BG})
target_compile_options(ca_polar_scl_bg PUBLIC -DLISTFLIPPINGFAST -DCUSTOM_SELECTION_METRIC -DDEC_NEEDS_SIGMA)
target_link_libraries(ca_polar_scl_bg m)

INSTALL(TARGETS ca_polar_scl_bg DESTINATION ${CMAKE_SOURCE_DIR}/work)
