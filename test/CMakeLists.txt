ADD_SUBDIRECTORY(sensitivity_solvers_test)
ADD_SUBDIRECTORY(hessian_eval_solvers_test)
ADD_SUBDIRECTORY(rol_adapter_test)
ADD_SUBDIRECTORY(LIS_solvers_test)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegressionTestDataCopy
    SOURCE_FILES
     regression_test	
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
