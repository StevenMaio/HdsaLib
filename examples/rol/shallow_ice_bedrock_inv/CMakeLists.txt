add_subdirectory(Data_Generation)

# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(shallow_ice_bedrock_inv_example_sensitivity.exe example_sensitivity.cpp)

# Link the test executable against the HDSALib interface library and any necessary libraries
target_link_libraries(shallow_ice_bedrock_inv_example_sensitivity.exe HDSALib HDSALibROL ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(shallow_ice_bedrock_inv_example_sensitivity.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/shallow_ice_bedrock_inv
)

add_executable(shallow_ice_bedrock_inv_example_inversion.exe example_inversion.cpp)

# Link the test executable against the HDSALib interface library and any necessary libraries
target_link_libraries(shallow_ice_bedrock_inv_example_inversion.exe HDSALib HDSALibROL ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(shallow_ice_bedrock_inv_example_inversion.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/shallow_ice_bedrock_inv
)

add_executable(shallow_ice_bedrock_inv_example_reoptimization.exe example_reoptimization.cpp)

# Link the test executable against the HDSALib interface library and any necessary libraries
target_link_libraries(shallow_ice_bedrock_inv_example_reoptimization.exe HDSALib HDSALibROL ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(shallow_ice_bedrock_inv_example_reoptimization.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/shallow_ice_bedrock_inv
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/examples/rol/shallow_ice_bedrock_inv)

# Copy the required files to the output directory
add_custom_command(TARGET shallow_ice_bedrock_inv_example_sensitivity.exe POST_BUILD
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
        ${SOURCE_DIR}/copy_results.py
        ${TARGET_DIR}/copy_results.py
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/initial_bedrock.m
        ${TARGET_DIR}/initial_bedrock.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/add_noise.m
        ${TARGET_DIR}/add_noise.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/prior_mean.m
        ${TARGET_DIR}/prior_mean.m
    COMMENT "Copying files to the output directory"
)



    
