cmake_minimum_required(VERSION 3.11)

project(raven VERSION 1.5.3
              LANGUAGES CXX
              DESCRIPTION "Raven is a de novo genome assembler for long uncorrected reads.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(raven_main_project ON)
endif ()
option(raven_build_tests "Build unit tests" ${raven_main_project})

find_package(bioparser 3.0.13 QUIET)
if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.13)

  FetchContent_GetProperties(bioparser)
  if (NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(
      ${bioparser_SOURCE_DIR}
      ${bioparser_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(cereal 1.3.0 QUIET)
if (NOT cereal_FOUND)
  FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal
    GIT_TAG v1.3.0)

  FetchContent_GetProperties(cereal)
  if (NOT cereal_POPULATED)
    FetchContent_Populate(cereal)
    add_subdirectory(
      ${cereal_SOURCE_DIR}
      ${cereal_BINARY_DIR}
      EXCLUDE_FROM_ALL)
    add_library(cereal::cereal ALIAS cereal)
  endif ()
endif ()

find_package(racon 3.0.2 QUIET)
if (NOT racon_FOUND)
  FetchContent_Declare(
    racon
    GIT_REPOSITORY https://github.com/lbcb-sci/racon
    GIT_TAG library)

  FetchContent_GetProperties(racon)
  if (NOT racon_POPULATED)
    FetchContent_Populate(racon)
    add_subdirectory(
      ${racon_SOURCE_DIR}
      ${racon_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

if (raven_build_tests)
  find_package(GTest 1.10.0 QUIET)
  if (NOT GTest_FOUND)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG release-1.10.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()
endif ()

add_library(raven
  src/graph.cpp
  src/pile.cpp)

target_link_libraries(raven
  cereal::cereal
  racon::racon)

target_include_directories(raven PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

add_executable(raven_exe
  src/main.cpp)

target_link_libraries(raven_exe
  raven
  bioparser::bioparser)

if (racon_enable_cuda)
  target_compile_definitions(raven_exe PRIVATE CUDA_ENABLED)
endif ()

target_compile_definitions(raven_exe PRIVATE VERSION="${PROJECT_VERSION}")
set_property(TARGET raven_exe PROPERTY OUTPUT_NAME raven)

install(TARGETS raven_exe DESTINATION ${CMAKE_INSTALL_BINDIR})

if (raven_build_tests)
  add_executable(raven_test
    test/raven_test.cpp)

  target_link_libraries(raven_test
    raven
    bioparser::bioparser
    GTest::Main)

  target_compile_definitions(raven_test
    PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/")
endif ()
