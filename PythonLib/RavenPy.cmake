if (NOT TARGET ravenpy)

    include(${CMAKE_CURRENT_LIST_DIR}/RavenPy.srcs.cmake)

    pybind11_add_module(ravenpy ${SOURCES})
    target_link_libraries(ravenpy PRIVATE RavenLib bioparser::bioparser)

    if (RAVENPY_EXT_DIR)
      install(TARGETS ravenpy
              LIBRARY DESTINATION ${RAVENPY_EXT_DIR})
    endif()
endif()
