ADD_SUBDIRECTORY(brandes_griesse)
ADD_SUBDIRECTORY(darcy_flow)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegressionTestRolAdaptersDataCopy
    SOURCE_FILES
     regression_test	
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
