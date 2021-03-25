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
    CDRModelErrorTestHDSADataCopy
    SOURCE_FILES
      Sensitivity_input.xml Driver.m  Model_Error_HDSA_Kronecker.m  Model_Error_HDSA.m  Model_Error_HDSA_React_Diff.m  Write_Data_for_Test_Problem.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
