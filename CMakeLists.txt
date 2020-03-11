CMAKE_MINIMUM_REQUIRED(VERSION 3.1 FATAL_ERROR)
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
PROJECT(cxSMSoliton CXX)


SET(cxSMSoliton_VERSION "0.1.2")

FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR} include)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")

## Collecting Source Files
FILE(GLOB lib_SOURCES src/*.cpp)

## Build the Libraries
# For .a library
ADD_LIBRARY(cxSMSoliton_static STATIC ${lib_SOURCES})
SET_TARGET_PROPERTIES(cxSMSoliton_static PROPERTIES
    LINKER_LANGUAGE CXX
    VERSION ${cxSMSoliton_VERSION}
    OUTPUT_NAME cxSMSoliton
)

# For .so/.dylib library
ADD_LIBRARY(cxSMSoliton_shared SHARED ${lib_SOURCES})
SET_TARGET_PROPERTIES(cxSMSoliton_shared PROPERTIES
    LINKER_LANGUAGE CXX
    VERSION ${cxSMSoliton_VERSION}
    OUTPUT_NAME cxSMSoliton
)

SET(cxSMSoliton_LIBRARY cxSMSoliton_shared)

add_executable(Scan cxSM_Scan.cpp)
target_link_libraries(Scan ${cxSMSoliton_LIBRARY})

add_subdirectory(test)