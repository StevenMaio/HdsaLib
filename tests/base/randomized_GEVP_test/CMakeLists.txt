# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(randomized_GEVP_test_example.exe example_gevp.cpp)

# Link the test executable against the HdsaLib interface library and any necessary libraries
target_link_libraries(randomized_GEVP_test_example.exe HdsaLib ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(randomized_GEVP_test_example.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests/base/randomized_GEVP_test
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/tests/base/randomized_GEVP_test)

# Copy the required files to the output directory
add_custom_command(TARGET randomized_GEVP_test_example.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Plotting_Script.m
        ${TARGET_DIR}/Plotting_Script.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/random_numbers.txt
        ${TARGET_DIR}/random_numbers.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Sabl_Output.mat
        ${TARGET_DIR}/Sabl_Output.mat
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/evals.txt
        ${TARGET_DIR}/reference_output/evals.txt
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/reference_output/Evec/Vector_10.txt
        ${TARGET_DIR}/reference_output/Evec/Vector_10.txt
    COMMENT "Copying files to the output directory"
)
