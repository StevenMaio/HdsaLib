add_subdirectory(rosenbrock_test)
add_subdirectory(diff_inv_test)

# Define the source directory and the target directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_DIR ${CMAKE_BINARY_DIR}/tests/pseudo_time_continuation/)

# Copy the run_tests.py file to the build directory above the test directories
add_custom_command(OUTPUT ${TARGET_DIR}/run_tests.py
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/run_tests.py
        ${TARGET_DIR}/run_tests.py
    COMMENT "Copying run_tests.py to the build directory"
)

# Copy the tester.py file to the build directory above the test directories
add_custom_command(OUTPUT ${TARGET_DIR}/tester.py
    COMMAND ${CMAKE_COMMAND} -E copy
        ${SOURCE_DIR}/tester.py
        ${TARGET_DIR}/tester.py
    COMMENT "Copying tester.py to the build directory"
)

# Ensure the custom command runs after the subdirectories are processed
add_custom_target(copy_run_tests_pseudo_time_continuation ALL DEPENDS ${TARGET_DIR}/run_tests.py)

# Ensure the custom command runs after the subdirectories are processed
add_custom_target(copy_tester_pseudo_time_continuation ALL DEPENDS ${TARGET_DIR}/tester.py)