cmake_minimum_required(VERSION 3.9)

project(raven VERSION 1.2.2
              LANGUAGES CXX
              DESCRIPTION "Raven is a de novo genome assembler for long uncorrected reads.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (NOT TARGET bioparser)
  add_subdirectory(vendor/racon/vendor/ram/vendor/bioparser EXCLUDE_FROM_ALL)
endif ()
if (NOT TARGET racon)
  add_subdirectory(vendor/racon EXCLUDE_FROM_ALL)
endif ()
if (NOT TARGET cereal)
  add_subdirectory(vendor/cereal EXCLUDE_FROM_ALL)
endif ()

add_executable(${PROJECT_NAME}
  src/graph.cpp
  src/main.cpp
  src/pile.cpp)
target_link_libraries(${PROJECT_NAME}
  bioparser
  cereal
  racon)
target_compile_definitions(${PROJECT_NAME} PRIVATE
  RAVEN_VERSION="v${PROJECT_VERSION}")
if (racon_enable_cuda)
  target_compile_definitions(${PROJECT_NAME} PRIVATE CUDA_ENABLED)
endif ()

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_BINDIR})

option(raven_build_tests "Build raven unit tests" OFF)
if (raven_build_tests)
  find_package(GTest REQUIRED)
  include_directories(${PROJECT_SOURCE_DIR}/src)
  add_executable(${PROJECT_NAME}_test
    src/graph.cpp
    src/pile.cpp
    test/raven_test.cpp)
  target_link_libraries(${PROJECT_NAME}_test
    bioparser
    cereal
    GTest::Main
    racon)
  target_compile_definitions(${PROJECT_NAME}_test PRIVATE
    RAVEN_DATA_PATH="${PROJECT_SOURCE_DIR}/test/data/")
endif ()
