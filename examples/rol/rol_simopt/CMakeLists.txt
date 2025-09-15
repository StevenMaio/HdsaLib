# Include directories for the test
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Include the current directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../..)  # Include the src directory

# Add the test executable
add_executable(rol_simopt_example_sensitivity.exe example_sensitivity.cpp)

# Link the test executable against the HDSALib interface library and any necessary libraries
target_link_libraries(rol_simopt_example_sensitivity.exe HDSALib HDSALibROL ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(rol_simopt_example_sensitivity.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/rol_simopt
)

add_executable(rol_simopt_example_optimization.exe example_optimization.cpp)

# Link the test executable against the HDSALib interface library and any necessary libraries
target_link_libraries(rol_simopt_example_optimization.exe HDSALib HDSALibROL ${Trilinos_LIBRARIES} ${MPI_LIBRARIES})

# Set the output directory for the executable
set_target_properties(rol_simopt_example_optimization.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/rol/rol_simopt
)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/examples/rol/rol_simopt)

# Copy the required files to the output directory
add_custom_command(TARGET rol_simopt_example_sensitivity.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/Plotting_Script.m
        ${TARGET_DIR}/Plotting_Script.m
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/input.xml
        ${TARGET_DIR}/input.xml
    COMMENT "Copying files to the output directory"
)
