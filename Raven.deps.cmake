include(FetchContent)
include(GNUInstallDirs)


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

find_package(racon 3.0.4 QUIET)
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

if (RAVEN_BUILD_PYTHON)
    find_package(pybind11 QUIET)
    if (NOT pybind11_FOUND)
        FetchContent_Declare(
                pybind11
                GIT_REPOSITORY https://github.com/pybind/pybind11
                GIT_TAG v2.9.1)


        FetchContent_GetProperties(pybind11)
        if (NOT pybind11_POPULATED)
            FetchContent_Populate(pybind11)
            add_subdirectory(
                    ${pybind11_SOURCE_DIR}
                    ${pybind11_BINARY_DIR}
                    EXCLUDE_FROM_ALL)
        endif()
    endif()
endif()

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
