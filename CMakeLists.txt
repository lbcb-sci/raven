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

option(RAVEN_BUILD_TESTS "Build unit tests" OFF)
option(RAVEN_BUILD_PYTHON "Build python bindings to raven library using pybind11" OFF)
option(RAVEN_BUILD_EXE "Build raven executable" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (RAVEN_BUILD_PYTHON)
  if (UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()

include(${CMAKE_CURRENT_LIST_DIR}/Raven.deps.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/RavenLib/RavenLib.cmake)

if (RAVEN_BUILD_TESTS)
  include(${CMAKE_CURRENT_LIST_DIR}/RavenTest/RavenTest.cmake)
endif()

if (RAVEN_BUILD_EXE)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "[raven] compling with g++ and using gold linker") 
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold")
  endif()
  include(${CMAKE_CURRENT_LIST_DIR}/RavenExe/Exe.cmake)
endif()

if (RAVEN_BUILD_PYTHON)
  include(${CMAKE_CURRENT_LIST_DIR}/PythonLib/RavenPy.cmake)
endif()
