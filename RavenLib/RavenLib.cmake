if (NOT TARGET RavenLib)

    include(${CMAKE_CURRENT_LIST_DIR}/RavenLib.srcs.cmake)

    add_library(RavenLib STATIC ${SOURCES})

    target_include_directories(RavenLib INTERFACE ${CMAKE_CURRENT_LIST_DIR}/Source)
    target_link_libraries(RavenLib PUBLIC cereal::cereal racon::racon)

endif()
