# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(shallow_ice_bedrock_inv_example_data_generation.exe shallow_ice_inversion_data_generation.cpp)

# Link the test executable against the HDSALib interface library and any necessary libraries
target_link_libraries(shallow_ice_bedrock_inv_example_data_generation.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(shallow_ice_bedrock_inv_example_data_generation.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/shallow_ice_bedrock_inv/Data_Generation
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/examples/rol/shallow_ice_bedrock_inv/Data_Generation)

# Copy the required files to the output directory
add_custom_command(TARGET shallow_ice_bedrock_inv_example_data_generation.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/plotresults.m
        ${TARGET_DIR}/plotresults.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/input_data_generation.xml
        ${TARGET_DIR}/input_data_generation.xml
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/plot_states.m
        ${TARGET_DIR}/plot_states.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/sparsify_data.m
        ${TARGET_DIR}/sparsify_data.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/copy_results.py
        ${TARGET_DIR}/copy_results.py
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/push_data_to_ip.py
        ${TARGET_DIR}/push_data_to_ip.py
    COMMENT "Copying files to the output directory"
)


