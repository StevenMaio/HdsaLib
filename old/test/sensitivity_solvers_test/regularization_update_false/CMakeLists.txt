ADD_SUBDIRECTORY(synthetic_test_source_code)
ADD_SUBDIRECTORY(randomized_GSVD)
ADD_SUBDIRECTORY(randomized_GEVP_JW)
ADD_SUBDIRECTORY(randomized_GEVP_gram)
ADD_SUBDIRECTORY(anasazi_GEVP_JW)
ADD_SUBDIRECTORY(anasazi_GEVP_gram)
ADD_SUBDIRECTORY(direct_sensitivities)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegressionTestSenSolversRegUpdateFalseDataCopy
    SOURCE_FILES
     regression_test	
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
