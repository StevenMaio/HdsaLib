ADD_SUBDIRECTORY(synthetic_test)
ADD_SUBDIRECTORY(synthetic_test_elliptic_prior)
ADD_SUBDIRECTORY(synthetic_test_hessian_gevp)
ADD_SUBDIRECTORY(synthetic_test_hyperparam_1D)
ADD_SUBDIRECTORY(PDE_Test_Problem)

 TRIBITS_COPY_FILES_TO_BINARY_DIR(
    TestDataCopy
    SOURCE_FILES
	run_tests.py
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )