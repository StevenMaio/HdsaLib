cmake_minimum_required(VERSION 3.6.1)
project(OedLib
    LANGUAGES NONE
) # some languages stuff?

set(CMAKE_CXX_STANDARD 20)

add_library(OedLib INTERFACE)

target_include_directories(OedLib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/bayesian_inversion/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/optimal_experimental_design/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/oed/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/algorithms/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/linear_algebra/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/base/vectors/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/core/utils/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/impl/bayesian_inversion/model_error/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/impl/bayesian_inversion/observation_operator/>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

IF(MrHyDE_ENABLE_OED)
    # add MrHyDE interfaces -- there's got to be a better way to do this
    MESSAGE("OedLib::CMakeLists.txt MrHyDE enabling OED")
    target_include_directories(OedLib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interfaces/mrhyde/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/adapter/trilinos/adapter/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/adapter/trilinos/core/>
    )

    # load TRILINOS
    SET(TRILINOS_DIR ${Trilinos_INSTALL_DIR})
    SET(Trilinos_SRC_DIR ${Trilinos_SRC_DIR})

    set(CMAKE_PREFIX_PATH ${TRILINOS_DIR} ${CMAKE_PREFIX_PATH})

    find_package(Trilinos REQUIRED)

    IF(NOT Trilinos_FOUND)
    MESSAGE(FATAL_ERROR "Could not find Trilinos!")
    ENDIF()

    # Echo trilinos build info just for fun
    MESSAGE("\nFound Trilinos!  Here are the details: ")
    MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
    MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
    MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
    MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
    MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
    MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
    MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
    MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
    MESSAGE("   Trilinos_CXX_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
    IF(Trilinos_VERSION LESS 14.0) # GH: if a version is not detected, CMake defaults to 0, so I'm making <14 the default behavior
    MESSAGE("   Trilinos_VERSION < 14 detected, using old names for Kokkos::Compat objects")
    ELSE()
    MESSAGE("   Trilinos_VERSION >= 14 detected, using new names for Tpetra::KokkosCompat objects")
    ADD_DEFINITIONS(-DHdsaLib_HAVE_TRILINOS14)
    ENDIF()
    MESSAGE("End of Trilinos details\n")

    MESSAGE("End of OedLib configuration details")
    MESSAGE("=======================================================")
    MESSAGE("\n")

    # Make sure to use same compilers and flags as Trilinos
    SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
    SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
    SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

    SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
    SET(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
    SET(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

    enable_language(C)
    enable_language(CXX)
    if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    endif()

    include(CTest)

    # MPI check 
    MESSAGE("-- Checking if MPI is enabled in Trilinos:")
    LIST(FIND Trilinos_TPL_LIST MPI MPI_List_ID)
    IF (MPI_List_ID GREATER -1)
    MESSAGE("-- Checking if MPI is enabled in Trilinos: MPI ENABLED")
    SET(main_MPI TRUE)
    ADD_DEFINITIONS(-Dmain_MPI)
    ELSE()
    MESSAGE("-- Checking if MPI is enabled in Trilinos: MPI NOT ENABLED")
    SET(main_MPI FALSE)
    ENDIF()

    MESSAGE("   CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

ENDIF()
