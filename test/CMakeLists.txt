find_package(Pagmo REQUIRED)

add_library(kep3_test STATIC catch_main.cpp)
target_compile_options(kep3_test PRIVATE
  "$<$<CONFIG:Debug>:${kep3_CXX_FLAGS_DEBUG}>"
  "$<$<CONFIG:Release>:${kep3_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:RelWithDebInfo>:${kep3_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:MinSizeRel>:${kep3_CXX_FLAGS_RELEASE}>"
)

# Setup the C++ standard.
target_compile_features(kep3_test PRIVATE cxx_std_23)
set_property(TARGET kep3_test PROPERTY CXX_EXTENSIONS NO)
target_link_libraries(kep3_test PRIVATE kep3)

function(ADD_kep3_TESTCASE arg1)
  add_executable(${arg1} ${arg1}.cpp)
  target_link_libraries(${arg1} PRIVATE kep3_test kep3 xtensor xtensor-blas Boost::boost Pagmo::pagmo)
  target_compile_options(${arg1} PRIVATE
    "$<$<CONFIG:Debug>:${kep3_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${kep3_CXX_FLAGS_RELEASE}>"
  )

  if (YACMA_COMPILER_IS_GNUCXX)
      target_compile_options(${arg1} PRIVATE -Wno-array-bounds)
      target_compile_options(${arg1} PRIVATE -Wno-stringop-overflow)
  endif()

  # Setup the C++ standard.
  target_compile_features(${arg1} PRIVATE cxx_std_23)
  set_property(TARGET ${arg1} PROPERTY CXX_EXTENSIONS NO)
  add_test(${arg1} ${arg1})
endfunction()

ADD_kep3_TESTCASE(convert_anomalies_test)
ADD_kep3_TESTCASE(epoch_test)
ADD_kep3_TESTCASE(planet_test)
ADD_kep3_TESTCASE(udpla_keplerian_test)
ADD_kep3_TESTCASE(udpla_jpl_lp_test)
ADD_kep3_TESTCASE(udpla_vsop2013_test)
ADD_kep3_TESTCASE(stm_test)
ADD_kep3_TESTCASE(ic2par2ic_test)
ADD_kep3_TESTCASE(ic2mee2ic_test)
ADD_kep3_TESTCASE(mee2par2mee_test)
ADD_kep3_TESTCASE(propagate_lagrangian_test)
ADD_kep3_TESTCASE(propagate_keplerian_test)
ADD_kep3_TESTCASE(lambert_problem_test)
ADD_kep3_TESTCASE(leg_sims_flanagan_test)
ADD_kep3_TESTCASE(leg_sims_flanagan_alpha_test)
ADD_kep3_TESTCASE(leg_zoh_test)
ADD_kep3_TESTCASE(ta_kep_test)
ADD_kep3_TESTCASE(ta_zoh_kep_test)
ADD_kep3_TESTCASE(ta_zoh_eq_test)
ADD_kep3_TESTCASE(ta_zoh_cr3bp_test)
ADD_kep3_TESTCASE(ta_zoh_ss_test)
ADD_kep3_TESTCASE(ta_cr3bp_test)
ADD_kep3_TESTCASE(ta_bcp_test)
ADD_kep3_TESTCASE(ta_pontryagin_cartesian_test)
ADD_kep3_TESTCASE(ta_pontryagin_equinoctial_test)
ADD_kep3_TESTCASE(flyby_test)
ADD_kep3_TESTCASE(encodings_test)
ADD_kep3_TESTCASE(mima_test)
ADD_kep3_TESTCASE(basic_transfers_test)
