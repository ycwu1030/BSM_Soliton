
add_executable(LocalMinimaTest cxSM_Z2_LocalMinimaTest.cpp)
add_executable(THDMLocalMinimaTest THDM_SCPV_LocalMinimaTest.cpp)
add_executable(BiasedLocalMinimaTest cxSM_Z2_biased_LocalMinimaTest.cpp)
add_executable(CPLocalMinimaTest cxSM_CP_LocalMinimaTest.cpp)
add_executable(CPReducedTest cxSM_CP_reduced_test.cpp)
add_executable(CPReducedA1Test cxSM_CP_reduced_a1_test.cpp)
add_executable(cxSMZ2DWTest RelaxationDW_cxSMZ2.cpp)
add_executable(cxSMZ2BiasedDWTest RelaxationDW_cxSMZ2_biased.cpp)
add_executable(cxSMCPreducedDWTest RelaxationDW_cxSMCP_reduced.cpp)
add_executable(cxSMCPreducedPDDWTest DeformationDW_cxSMCP_reduced.cpp)
add_executable(cxSMCPreducedA1PDDWTest DeformationDW_cxSMCP_reduced_a1.cpp)
add_executable(THDMSCPVPDDWTest DeformationDW_THDM_SCPV.cpp)
add_executable(PDDWTemplate DeformationDW_template.cpp)
add_executable(QuarticSolver QuarticSolver_test.cpp)
add_executable(Kink1D Kink1D_test.cpp)
add_executable(Kink1DMem Kink1D_memtest.cpp)
add_executable(KinknD KinknD_test.cpp)
add_executable(KinknDMem KinknD_memtest.cpp)
add_executable(SplineTest Spline_Inter_Potential_test.cpp)
add_executable(SplineMemTest Spline_Inter_Potential_memtest.cpp)
add_executable(SplineInterMemTest Spline_Inter_memtest.cpp)
add_executable(CMMTest CMMonopole_test.cpp)
add_executable(CMMUVTest CMMonopoleUV_test.cpp)
add_executable(THDMCMMTest THDMCMM_test.cpp)
add_executable(THDMCMMUVTest THDMCMMUV_test.cpp)
# add_executable(RKTest RungeKutta_Test.cpp)
# add_executable(ShootingTest Shooting_Test.cpp)
# add_executable(VortexTest NOVortex_Test.cpp)
# add_executable(VortexRKTest NOVortex_RK_Test.cpp)
# add_executable(RelaxationTest RelaxationTest.cpp)
add_executable(RelaxationVortexTest RelaxationVortex.cpp)
# add_executable(RelaxationVortexTheta RelaxationVortex_theta.cpp)

target_link_libraries(CMMTest ${cxSMSoliton_LIBRARY})
target_link_libraries(CMMUVTest ${cxSMSoliton_LIBRARY})
target_link_libraries(THDMCMMTest ${cxSMSoliton_LIBRARY})
target_link_libraries(THDMCMMUVTest ${cxSMSoliton_LIBRARY})
target_link_libraries(LocalMinimaTest ${cxSMSoliton_LIBRARY})
target_link_libraries(THDMLocalMinimaTest ${cxSMSoliton_LIBRARY})
target_link_libraries(BiasedLocalMinimaTest ${cxSMSoliton_LIBRARY})
target_link_libraries(CPLocalMinimaTest ${cxSMSoliton_LIBRARY})
target_link_libraries(CPReducedTest ${cxSMSoliton_LIBRARY})
target_link_libraries(CPReducedA1Test ${cxSMSoliton_LIBRARY})
target_link_libraries(cxSMZ2DWTest ${cxSMSoliton_LIBRARY})
target_link_libraries(cxSMZ2BiasedDWTest ${cxSMSoliton_LIBRARY})
target_link_libraries(cxSMCPreducedDWTest ${cxSMSoliton_LIBRARY})
target_link_libraries(cxSMCPreducedPDDWTest ${cxSMSoliton_LIBRARY})
target_link_libraries(cxSMCPreducedA1PDDWTest ${cxSMSoliton_LIBRARY})
target_link_libraries(THDMSCPVPDDWTest ${cxSMSoliton_LIBRARY})
target_link_libraries(PDDWTemplate ${cxSMSoliton_LIBRARY})
target_link_libraries(QuarticSolver ${cxSMSoliton_LIBRARY})
target_link_libraries(Kink1D ${cxSMSoliton_LIBRARY})
target_link_libraries(Kink1DMem ${cxSMSoliton_LIBRARY})
target_link_libraries(KinknD ${cxSMSoliton_LIBRARY})
target_link_libraries(KinknDMem ${cxSMSoliton_LIBRARY})
target_link_libraries(SplineTest ${cxSMSoliton_LIBRARY})
target_link_libraries(SplineMemTest ${cxSMSoliton_LIBRARY})
target_link_libraries(SplineInterMemTest ${cxSMSoliton_LIBRARY})
# target_link_libraries(RKTest ${cxSMSoliton_LIBRARY})
# target_link_libraries(ShootingTest ${cxSMSoliton_LIBRARY})
# target_link_libraries(VortexTest ${cxSMSoliton_LIBRARY})
# target_link_libraries(VortexRKTest ${cxSMSoliton_LIBRARY})
# target_link_libraries(RelaxationTest ${cxSMSoliton_LIBRARY})
target_link_libraries(RelaxationVortexTest ${cxSMSoliton_LIBRARY})
# target_link_libraries(RelaxationVortexTheta ${cxSMSoliton_LIBRARY})