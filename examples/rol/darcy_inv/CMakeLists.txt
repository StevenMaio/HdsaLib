add_subdirectory(Data_Generation)

# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(darcy_inv_example_sensitivity.exe example_sensitivity.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(darcy_inv_example_sensitivity.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(darcy_inv_example_sensitivity.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/darcy_inv
)

add_executable(darcy_inv_example_inversion.exe example_inversion.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(darcy_inv_example_inversion.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(darcy_inv_example_inversion.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/darcy_inv
)

add_executable(darcy_inv_example_perturbed_inversion.exe example_perturbed_inversion.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(darcy_inv_example_perturbed_inversion.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(darcy_inv_example_perturbed_inversion.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/darcy_inv
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/examples/rol/darcy_inv)

# Copy the required files to the output directory
add_custom_command(TARGET darcy_inv_example_sensitivity.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/plotresults.m
        ${TARGET_DIR}/plotresults.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/input.xml
        ${TARGET_DIR}/input.xml
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Sensitivity_input.xml
        ${TARGET_DIR}/Sensitivity_input.xml
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Convert_results_to_mat_file.m
        ${TARGET_DIR}/Convert_results_to_mat_file.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Read_Inputs_from_Xml.m
        ${TARGET_DIR}/Read_Inputs_from_Xml.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/plot_sensitivity_results.m
        ${TARGET_DIR}/plot_sensitivity_results.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Read_Results.m 
        ${TARGET_DIR}/Read_Results.m 
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Plot_cost_per_timestep_comparison.m
        ${TARGET_DIR}/Plot_cost_per_timestep_comparison.m
    COMMENT "Copying files to the output directory"
)




    
