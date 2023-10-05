INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)	

 TRIBITS_COPY_FILES_TO_BINARY_DIR(
    RegTestRegUpdateFalseDirectRsEnforceZerosTrueDataCopy
    SOURCE_FILES
      ref_Sensitivity_Operator_Column_13.txt  ref_Sensitivity_Operator_Column_18.txt  ref_Sensitivity_Operator_Column_3.txt  ref_Sensitivity_Operator_Column_8.txt ref_Sensitivity_Indices.txt ref_Sensitivity_Operator_Column_14.txt  ref_Sensitivity_Operator_Column_19.txt  ref_Sensitivity_Operator_Column_4.txt  ref_Sensitivity_Operator_Column_9.txt ref_Sensitivity_Operator_Column_10.txt  ref_Sensitivity_Operator_Column_15.txt  ref_Sensitivity_Operator_Column_1.txt   ref_Sensitivity_Operator_Column_5.txt ref_Sensitivity_Operator_Column_11.txt  ref_Sensitivity_Operator_Column_16.txt  ref_Sensitivity_Operator_Column_20.txt  ref_Sensitivity_Operator_Column_6.txt  Sensitivity_input.xml ref_Sensitivity_Operator_Column_12.txt  ref_Sensitivity_Operator_Column_17.txt  ref_Sensitivity_Operator_Column_2.txt   ref_Sensitivity_Operator_Column_7.txt
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
