ADD_SUBDIRECTORY(regularization_update_true)
ADD_SUBDIRECTORY(regularization_update_false)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegressionTestHessianEvalSolversDataCopy
    SOURCE_FILES
     regression_test	
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
