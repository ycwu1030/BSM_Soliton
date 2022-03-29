
add_executable(ReduceTest Reduce_to_Zero_test.cpp)
add_executable(RelaxationTest RelaxationTest.cpp)
add_executable(PDDWToyZ4Scan DeformationDW_ToyZ4_Scan.cpp)
add_executable(PDDWToyZ3Scan DeformationDW_ToyZ3_Scan.cpp)
add_executable(PDDWToyZ4Test DeformationDW_ToyZ4_Test.cpp)

target_link_libraries(ReduceTest ${cxSMSoliton_LIBRARY})
target_link_libraries(RelaxationTest ${cxSMSoliton_LIBRARY})
target_link_libraries(PDDWToyZ4Scan ${cxSMSoliton_LIBRARY})
target_link_libraries(PDDWToyZ3Scan ${cxSMSoliton_LIBRARY})
target_link_libraries(PDDWToyZ4Test ${cxSMSoliton_LIBRARY})
