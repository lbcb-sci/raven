cmake_minimum_required(VERSION 3.2)
project(raven)

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(raven_build_tests "Build raven unit tests" OFF)

add_executable(raven
    src/graph.cpp
    src/main.cpp
    src/pile.cpp)

if (NOT TARGET bioparser)
    add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
endif()
if (NOT TARGET logger)
    add_subdirectory(vendor/logger EXCLUDE_FROM_ALL)
endif()
if (NOT TARGET ram)
    add_subdirectory(vendor/ram EXCLUDE_FROM_ALL)
endif()

target_link_libraries(raven bioparser logger ram)

install(TARGETS raven DESTINATION ${CMAKE_INSTALL_BINDIR})

if (raven_build_tests)
endif()
