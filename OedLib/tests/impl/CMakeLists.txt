cmake_minimum_required(VERSION 3.6.1)

project(OedLibImpl)

add_library(OedLibImpl INTERFACE)

target_link_libraries(OedLibImpl INTERFACE OedLib)
add_subdirectory(../../ oed-lib/)  # Include the src directory -- maybe add conditional to avoid this

target_include_directories(OedLibImpl INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/base/bayesian_inversion
        ${CMAKE_CURRENT_SOURCE_DIR}/base/vectors
        ${CMAKE_CURRENT_SOURCE_DIR}/base/linear_algebra
)
