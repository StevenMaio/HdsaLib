INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)	

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_sensitivity
    SOURCES example_sensitivity.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

 TRIBITS_COPY_FILES_TO_BINARY_DIR(
    SyntheticTestTransientDataCopy
    SOURCE_FILES
	Plotting_Script.m Sabl_Output.mat random_numbers.txt run_test.py 
	reference_output/prior_discrepancy_sample_1/Vector_3.txt
	reference_output/prior_discrepancy_evaluated_at_z_opt/Vector_10.txt
	reference_output/posterior_discrepancy_mean_1.txt
	reference_output/posterior_discrepancy_mean_3.txt
	reference_output/posterior_discrepancy_samples_3/Vector_20.txt
	reference_output/posterior_update_mean.txt 
	reference_output/posterior_update_samples/Vector_40.txt
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

# Specify the source directory containing the files you want to copy
#set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Specify the binary directory where you want to copy the files
#set(BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

# Use file(GLOB ...) to get all files in the source directory
#file(GLOB ALL_FILES "${SOURCE_DIR}/reference_output/*.txt")

# Copy all files to the binary directory
#foreach(FILE ${ALL_FILES})
#    get_filename_component(FILENAME ${FILE} NAME)
#    TRIBITS_COPY_FILES_TO_BINARY_DIR("${FILENAME}" "${FILE}" "${BINARY_DIR}")
#endforeach()
