INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)	

 TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegTestRegUpdateTrueRandGevpGramRsConserveMemoryTrueEnforceZerosFalseDataCopy
    SOURCE_FILES
      Sensitivity_input.xml ref_Sensitivity_Indices_1.txt  ref_Singular_Values_1.txt  ref_theta_Sample_1.txt  ref_theta_Singular_Vector_1.txt Second_Order_Update_evals.txt Second_Order_Update_evecs.txt
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
