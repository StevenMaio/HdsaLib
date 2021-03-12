ADD_SUBDIRECTORY(double_pass)
ADD_SUBDIRECTORY(nystrom)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegressionTestHessianEvalSolversRegUpdateTrueRandDataCopy
    SOURCE_FILES
     regression_test	
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
