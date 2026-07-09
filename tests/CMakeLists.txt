add_subdirectory(base)
add_subdirectory(model_discrepancy)
add_subdirectory(pseudo_time_continuation)
add_subdirectory(bifidelity_correction)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/tests/)

# Copy the run_tests.py file to the build directory above the test directories
add_custom_command(OUTPUT ${TARGET_DIR}/execute_all_tests
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/execute_all_tests
        ${TARGET_DIR}/execute_all_tests
    COMMENT "Copying execute_all_tests to the build directory"
)

# Ensure the custom command runs after the subdirectories are processed
add_custom_target(copy_execute_all_tests ALL DEPENDS ${TARGET_DIR}/execute_all_tests)