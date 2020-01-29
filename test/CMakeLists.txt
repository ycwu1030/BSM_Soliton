
add_executable(LocalMinimaTest cxSM_LocalMinimaTest.cpp)
add_executable(SolitonTest cxSM_SolitonTest.cpp)

target_link_libraries(LocalMinimaTest ${cxSMSoliton_LIBRARY})
target_link_libraries(SolitonTest ${cxSMSoliton_LIBRARY})