ADD_SUBDIRECTORY(Data_Generation)

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
    example_inversion
    SOURCES example_inversion.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_perturbed_inversion
    SOURCES example_perturbed_inversion.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_sensitivity
    SOURCES example_sensitivity.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    DarcyInvDataCopy
    SOURCE_FILES
      input.xml plotresults.m Sensitivity_input.xml Convert_results_to_mat_file.m Read_Inputs_from_Xml.m plot_sensitivity_results.m Read_Results.m Plot_cost_per_timestep_comparison.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
