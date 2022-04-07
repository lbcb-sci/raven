if (NOT TARGET ravenpy)

    include(${CMAKE_CURRENT_LIST_DIR}/RavenPy.srcs.cmake)

    pybind11_add_module(ravenpy ${SOURCES})
    target_link_libraries(ravenpy PRIVATE raven bioparser::bioparser)
endif()
