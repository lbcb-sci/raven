cmake_minimum_required(VERSION 3.11)

project(raven VERSION 1.7.0
        LANGUAGES CXX
        DESCRIPTION "Raven is a de novo genome assembler for long uncorrected reads.")


set(CXX_VISIBILITY_PRESET)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -fvisibility=hidden")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(raven_main_project ON)
endif ()

option(raven_build_tests "Build unit tests" ${raven_main_project})

include(${CMAKE_CURRENT_LIST_DIR}/Raven.deps.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/RavenLib/RavenLib.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/Demo/Demo.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/RavenTest/RavenTest.cmake)


if (RAVEN_BUILD_PYTHON)
  include(${CMAKE_CURRENT_LIST_DIR}/PythonLib/RavenPy.cmake)
endif()

