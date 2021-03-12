ADD_SUBDIRECTORY(reduced_space)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegressionTestSenSolversRegUpdateTrueDirectDataCopy
    SOURCE_FILES
     regression_test	
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
