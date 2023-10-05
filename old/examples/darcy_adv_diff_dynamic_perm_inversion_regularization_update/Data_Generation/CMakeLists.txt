INCLUDE(TribitsAddExecutableAndTest)

IF(${PROJECT_NAME}_ENABLE_Intrepid AND 
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  TRIBITS_ADD_EXECUTABLE(
    perm_inversion_data_generation
    SOURCES perm_inversion_data_generation.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    DarcyAdvDiffPermInversionDataGenerationHDSADataCopy
    SOURCE_FILES
      input_data_generation.xml plotresults.m Sparcify_Data.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
