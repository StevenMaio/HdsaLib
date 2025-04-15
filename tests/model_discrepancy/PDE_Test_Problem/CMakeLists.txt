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
    PDETestProblemDataCopy
    SOURCE_FILES
	u_opt.txt z_opt.txt Z.txt D.txt random_numbers.txt Sabl_Output.mat Plotting_Script.m run_test.py 
	reference_output/random_numbers.txt
	reference_output/u_opt.txt
	reference_output/z_opt.txt
	reference_output/Z.txt
	reference_output/D.txt
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
