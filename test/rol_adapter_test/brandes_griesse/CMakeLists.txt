INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)	

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_control
    SOURCES example_control.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

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
    BrandesGriesseHDSARegTestDataCopy
    SOURCE_FILES
      input.xml Sensitivity_input.xml plotresults.m control_read.txt ref_control.txt ref_state.txt ref_Sensitivity_Indices_1.txt  ref_Singular_Values_1.txt  ref_theta_Sample_1.txt  ref_theta_Singular_Vector_1.txt  ref_z_Singular_Vector_1.txt regression_test
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
