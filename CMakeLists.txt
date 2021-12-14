CMAKE_MINIMUM_REQUIRED(VERSION 3.1 FATAL_ERROR)
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
PROJECT(cxSMSoliton CXX)


SET(cxSMSoliton_VERSION "0.2.2")

FIND_PACKAGE(GSL 2.1 REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR} ${GSL_INCLUDE_DIRS} include)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
# SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wno-unknown-warning-option")

## Collecting Source Files
FILE(GLOB lib_SOURCES src/*.cpp)

## Build the Libraries
# For .a library
# ADD_LIBRARY(cxSMSoliton_static STATIC ${lib_SOURCES})
# SET_TARGET_PROPERTIES(cxSMSoliton_static PROPERTIES
#    LINKER_LANGUAGE CXX
#    VERSION ${cxSMSoliton_VERSION}
#    OUTPUT_NAME cxSMSoliton
#)
# TARGET_LINK_LIBRARIES(cxSMSoliton_static ${GSL_LIBRARIES})

# For .so/.dylib library
ADD_LIBRARY(cxSMSoliton_shared SHARED ${lib_SOURCES})
SET_TARGET_PROPERTIES(cxSMSoliton_shared PROPERTIES
    LINKER_LANGUAGE CXX
    VERSION ${cxSMSoliton_VERSION}
    OUTPUT_NAME cxSMSoliton
)
TARGET_LINK_LIBRARIES(cxSMSoliton_shared ${GSL_LIBRARIES})

SET(cxSMSoliton_LIBRARY cxSMSoliton_shared)

add_executable(ToyZ4 models/ToyZ4/ToyZ4.cpp)
add_executable(ToyKink models/ToyKink/ToyKink.cpp)
# add_executable(Z2Scan scanner/cxSM_Z2_Scan.cpp)
# add_executable(Z2BiasedScan scanner/cxSM_Z2_Biased_Scan.cpp)
# add_executable(CPParamScan scanner/cxSM_CP_reduced_Scan_param_only.cpp)
# add_executable(CPParamScanCons scanner/cxSM_CP_reduced_Scan_param_only_cons.cpp)
# add_executable(CPScan scanner/cxSM_CP_reduced_Scan.cpp)
# add_executable(THDMScan scanner/THDM_SCPV_Scan.cpp)
# add_executable(cxSMSCPVScan scanner/cxSM_CP_reduced_a1_Scan_cons.cpp)
# add_executable(cxSMSCPVScanFM scanner/cxSM_CP_reduced_a1_Scan_cons_fix_m2.cpp)
target_link_libraries(ToyZ4 ${cxSMSoliton_LIBRARY})
target_link_libraries(ToyKink ${cxSMSoliton_LIBRARY})
# target_link_libraries(Z2Scan ${cxSMSoliton_LIBRARY})
# target_link_libraries(Z2BiasedScan ${cxSMSoliton_LIBRARY})
# target_link_libraries(CPScan ${cxSMSoliton_LIBRARY})
# target_link_libraries(CPParamScan ${cxSMSoliton_LIBRARY})
# target_link_libraries(CPParamScanCons ${cxSMSoliton_LIBRARY})
# target_link_libraries(THDMScan ${cxSMSoliton_LIBRARY})
# target_link_libraries(cxSMSCPVScan ${cxSMSoliton_LIBRARY})
# target_link_libraries(cxSMSCPVScanFM ${cxSMSoliton_LIBRARY})

add_subdirectory(test)
