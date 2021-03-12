INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)	

 TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegTestRegUpdateTrueLISGEVPNystromConserveMemoryTrueEnforceZerosTrueDataCopy
    SOURCE_FILES
        ref_eigenvector_weightmat_products_1.txt  ref_LIS_evals_1.txt  ref_LIS_evecs_1.txt  Sensitivity_input.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
