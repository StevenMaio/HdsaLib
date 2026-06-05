cmake_minimum_required(VERSION 3.6.1)
project(OedLib) # some languages stuff?

set(CMAKE_CXX_STANDARD 20)

add_library(OedLib INTERFACE)

target_include_directories(OedLib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/bayesian_inversion/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/optimal_experimental_design/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/oed/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/algorithms/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/linear_algebra/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/vectors/
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
