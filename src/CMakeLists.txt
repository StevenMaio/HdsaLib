# Create an interface library for HDSALib
add_library(HDSALib INTERFACE)

# Specify include directories for the interface library
target_include_directories(HDSALib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/interfaces/hyperparameter_interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/interfaces/problem_interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/interfaces/ouu_interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/interfaces/transient_prior_interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/interfaces/u_prior_interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/interfaces/z_prior_interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/analysis/prior>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/analysis/posterior>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/model_discrepancy/analysis/auxillary>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/parametric_continuation/interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/parametric_continuation/analysis>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/base>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/trilinos_adapter/adapters>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/trilinos_adapter/core>
    ${Trilinos_INCLUDE_DIRS}
)

# Optionally, link against Trilinos libraries if needed
# target_link_libraries(HDSALib INTERFACE ${Trilinos_LIBRARIES})