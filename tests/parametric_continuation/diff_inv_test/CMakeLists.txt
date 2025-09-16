# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(diff_inv_test_example_sensitivity.exe example_sensitivity.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(diff_inv_test_example_sensitivity.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(diff_inv_test_example_sensitivity.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/parametric_continuation/diff_inv_test
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/tests/parametric_continuation/diff_inv_test)

# Copy the required files to the output directory
add_custom_command(TARGET diff_inv_test_example_sensitivity.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Plotting_Script.m
        ${TARGET_DIR}/Plotting_Script.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/obs_data.txt
        ${TARGET_DIR}/obs_data.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/z_bar.txt
        ${TARGET_DIR}/z_bar.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/theta_bar.txt
        ${TARGET_DIR}/theta_bar.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/theta_star.txt
        ${TARGET_DIR}/theta_star.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/random_numbers.txt 
        ${TARGET_DIR}/random_numbers.txt 
    COMMENT "Copying files to the output directory"
)
