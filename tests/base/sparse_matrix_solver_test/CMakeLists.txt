# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(sparse_matrix_solver_test_example.exe example.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(sparse_matrix_solver_test_example.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(sparse_matrix_solver_test_example.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/base/sparse_matrix_solver_test
)

