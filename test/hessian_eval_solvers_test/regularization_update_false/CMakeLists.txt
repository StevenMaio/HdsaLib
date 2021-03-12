ADD_SUBDIRECTORY(synthetic_test_source_code)
ADD_SUBDIRECTORY(randomized_hessian_EVP)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegressionTestHessianEvalSolversRegUpdateFalseDataCopy
    SOURCE_FILES
     regression_test	
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
