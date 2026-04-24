# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(sparse_matrix_sqrt_test_example.exe example.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(sparse_matrix_sqrt_test_example.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(sparse_matrix_sqrt_test_example.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/base/sparse_matrix_sqrt_test
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/tests/base/sparse_matrix_sqrt_test)

# Copy the required files to the output directory
add_custom_command(TARGET sparse_matrix_sqrt_test_example.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/run_test.py 
        ${TARGET_DIR}/run_test.py 
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/error.txt
        ${TARGET_DIR}/reference_output/error.txt
)