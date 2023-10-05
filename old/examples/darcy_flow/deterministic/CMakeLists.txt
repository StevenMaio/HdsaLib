INCLUDE(TribitsAddExecutableAndTest)

IF(${PROJECT_NAME}_ENABLE_Intrepid AND 
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  # Need ROL_TpetraTeuchosBatchManager.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)


  TRIBITS_ADD_EXECUTABLE(
    example_control
    SOURCES example_control.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_sensitivity
    SOURCES example_sensitivity.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    DarcyFlowHDSADataCopy
    SOURCE_FILES
      input.xml Sensitivity_input.xml plotresults.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
