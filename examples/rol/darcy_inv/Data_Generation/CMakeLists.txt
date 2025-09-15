# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(darcy_inv_example_data_generation.exe example_data_generation.cpp)

# Link the test executable against the HDSALib interface library and any necessary libraries
target_link_libraries(darcy_inv_example_data_generation.exe HDSALib HDSALibROL ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(darcy_inv_example_data_generation.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/darcy_inv/Data_Generation
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/examples/rol/darcy_inv/Data_Generation)

# Copy the required files to the output directory
add_custom_command(TARGET darcy_inv_example_data_generation.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/plotresults.m
        ${TARGET_DIR}/plotresults.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Generate_Obs_Data.m
        ${TARGET_DIR}/Generate_Obs_Data.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Design_Permeability_Field.m
        ${TARGET_DIR}/Design_Permeability_Field.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Design_Forcing.m
        ${TARGET_DIR}/Design_Forcing.m
    COMMENT "Copying files to the output directory"
)


