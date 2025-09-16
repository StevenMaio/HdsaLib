# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(synthetic_test_transient_multi_state_example_sensitivity.exe example_sensitivity.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(synthetic_test_transient_multi_state_example_sensitivity.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(synthetic_test_transient_multi_state_example_sensitivity.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/model_discrepancy/synthetic_test_transient_multi_state
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/tests/model_discrepancy/synthetic_test_transient_multi_state)

# Copy the required files to the output directory
add_custom_command(TARGET synthetic_test_transient_multi_state_example_sensitivity.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Plotting_Script.m
        ${TARGET_DIR}/Plotting_Script.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/run_test.py 
        ${TARGET_DIR}/run_test.py 
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/random_numbers.txt
        ${TARGET_DIR}/random_numbers.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Sabl_Output.mat
        ${TARGET_DIR}/Sabl_Output.mat
    # Copy the specified files from the reference_output directory
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/prior_discrepancy_sample_1/Vector_3_time_9.txt
        ${TARGET_DIR}/reference_output/prior_discrepancy_sample_1/Vector_3_time_9.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/prior_discrepancy_evaluated_at_z_opt/Vector_10_time_4.txt
        ${TARGET_DIR}/reference_output/prior_discrepancy_evaluated_at_z_opt/Vector_10_time_4.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/posterior_discrepancy_mean_1_time_5.txt
        ${TARGET_DIR}/reference_output/posterior_discrepancy_mean_1_time_5.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/posterior_discrepancy_mean_3_time_2.txt
        ${TARGET_DIR}/reference_output/posterior_discrepancy_mean_3_time_2.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/posterior_discrepancy_samples_3/Vector_20_time_7.txt
        ${TARGET_DIR}/reference_output/posterior_discrepancy_samples_3/Vector_20_time_7.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/posterior_update_mean.txt 
        ${TARGET_DIR}/reference_output/posterior_update_mean.txt 
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/posterior_update_samples/Vector_40.txt
        ${TARGET_DIR}/reference_output/posterior_update_samples/Vector_40.txt
    COMMENT "Copying files to the output directory"
)
