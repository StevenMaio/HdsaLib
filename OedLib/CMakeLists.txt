cmake_minimum_required(VERSION 3.6.1)
project(OedLib
    LANGUAGES NONE
) # some languages stuff?

set(CMAKE_CXX_STANDARD 20)

add_library(OedLib INTERFACE)

target_include_directories(OedLib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/bayesian_inversion/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/optimal_experimental_design/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/oed/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/algorithms/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/linear_algebra/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/vectors/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/impl/bayesian_inversion/model_error/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

IF(MrHyDE_ENABLE_OED)
MESSAGE("OedLib::CMakeLists.txt MrHyDE enabling OED")
target_include_directories(OedLib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interfaces/mrhyde/>
)
ENDIF()
