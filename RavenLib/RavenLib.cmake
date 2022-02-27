if (NOT TARGET RavenLib)

    include(${CMAKE_CURRENT_LIST_DIR}/RavenLib.srcs.cmake)

    add_library(RavenLib STATIC ${SOURCES})

    if (raven_main_project)

      if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(ASAN_FLAGS -fno-omit-frame-pointer -fsanitize=address)
        set(WARNINGS -Wall -Wextra -pedantic -fvisibility=hidden)

        add_compile_options("${ASAN_FLAGS} ${WARNINGS}")
        add_link_options(${ASAN_FLAGS})
      endif()
    endif()

    target_include_directories(RavenLib INTERFACE ${CMAKE_CURRENT_LIST_DIR}/Source)
    target_link_libraries(RavenLib PUBLIC cereal::cereal racon::racon)

endif()
