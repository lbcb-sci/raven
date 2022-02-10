if (NOT TARGET RavenPy)

    include(${CMAKE_CURRENT_LIST_DIR}/RavenPy.srcs.cmake)

    pybind11_add_module(RavenPy ${SOURCES})
    target_link_libraries(RavenPy PRIVATE RavenLib bioparser::bioparser)

endif()
