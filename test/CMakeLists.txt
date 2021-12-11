
add_executable(ReduceTest Reduce_to_Zero_test.cpp)
add_executable(RelaxationTest RelaxationTest.cpp)

target_link_libraries(ReduceTest ${cxSMSoliton_LIBRARY})
target_link_libraries(RelaxationTest ${cxSMSoliton_LIBRARY})
