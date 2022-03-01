cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT DEFINED CMAKE_CXX_VISIBLITY_PRESET AND
    NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
  set(CMAKE_CXX_VISIBLITY_PRESET hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif()

project(raven VERSION 1.8
        LANGUAGES CXX
        DESCRIPTION "Raven is a de novo genome assembler for long uncorrected reads.")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

option(RAVEN_BUILD_TESTS "Build unit tests" ${RAVEN_MAIN_PROJECT})
option(RAVEN_BUILD_PYTHON "Build python bindings to raven library using pybind11" OFF)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(RAVEN_MAIN_PROJECT ON)
  message(STATUS "[raven] main project")

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "[raven] compling with g++ and using gold linker") 
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold")
  endif()
endif ()

if (RAVEN_BUILD_PYTHON)
  if (UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()

include(${CMAKE_CURRENT_LIST_DIR}/Raven.deps.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/RavenLib/RavenLib.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/RavenTest/RavenTest.cmake)

if (RAVEN_MAIN_PROJECT)
  include(${CMAKE_CURRENT_LIST_DIR}/RavenExe/Exe.cmake)
endif()

if (RAVEN_BUILD_PYTHON)
  include(${CMAKE_CURRENT_LIST_DIR}/PythonLib/RavenPy.cmake)
endif()
