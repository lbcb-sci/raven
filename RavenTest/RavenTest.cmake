if (NOT TARGET RavenTest)

    include(${CMAKE_CURRENT_LIST_DIR}/RavenTest.srcs.cmake)

    add_executable(RavenTest ${SOURCES})
    target_link_libraries(RavenTest RavenLib bioparser::bioparser GTest::Main)

    target_compile_definitions(RavenTest  PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/RavenTest/data/")

endif()
