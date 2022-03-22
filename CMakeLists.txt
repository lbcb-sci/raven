cmake_minimum_required(VERSION 3.11)
project(raven 
          VERSION 1.8.1
          LANGUAGES CXX
          DESCRIPTION "Raven is a de novo genome assembler for long uncorrected reads.")

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(RAVEN_MAIN_PROJECT TRUE)
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT DEFINED CMAKE_CXX_VISIBLITY_PRESET AND
    NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
  set(CMAKE_CXX_VISIBLITY_PRESET hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif ()

option(RAVEN_BUILD_TESTS "Build unit tests" OFF)
option(RAVEN_BUILD_PYTHON "Build python bindings to raven library using pybind11" OFF)
option(RAVEN_BUILD_SHARED_LIBS "Build raven and it's dependencies as a shared libraries" OFF)
option(RAVEN_BUILD_EXE "Build raven executable" OFF)

if (RAVEN_BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
  message(STATUS "[raven] building shared libs")
endif ()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(${CMAKE_CURRENT_LIST_DIR}/Raven.deps.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/RavenLib/RavenLib.cmake)

if (RAVEN_BUILD_TESTS)
  include(${CMAKE_CURRENT_LIST_DIR}/RavenTest/RavenTest.cmake)
endif ()

if (RAVEN_BUILD_EXE)
  include(${CMAKE_CURRENT_LIST_DIR}/RavenExe/Exe.cmake)
endif ()

if (RAVEN_BUILD_PYTHON)
  include(${CMAKE_CURRENT_LIST_DIR}/PythonLib/RavenPy.cmake)
endif ()

if (RAVEN_MAIN_PROJECT)
  include(${CMAKE_CURRENT_LIST_DIR}/Packing/Packing.cmake) 
endif ()
