if (NOT TARGET RavenLib)

  include(${CMAKE_CURRENT_LIST_DIR}/RavenLib.srcs.cmake)
  add_library(RavenLib STATIC ${SOURCES})

  if (RAVEN_MAIN_PROJECT OR RAVEN_BUILD_PYTHON)
    set(ASAN_FLAGS -fno-omit-frame-pointer -fsanitize=address)
    set(DEBUG_WANINGS -Wall -Wextra -pedantic)

    target_compile_options(RavenLib PUBLIC ${DEBUG_WANINGS})
    target_compile_options(RavenLib PUBLIC 
      "$<$<CONFIG:Debug>:${ASAN_FLAGS}>"
      "$<$<CONFIG:RelWithDebInfo>:${ASAN_FLAGS}>")

    target_link_options(RavenLib PUBLIC
      "$<$<CONFIG:Debug>:${ASAN_FLAGS}>"
      "$<$<CONFIG:RelWithDebInfo>:${ASAN_FLAGS}>")
  endif()

  target_include_directories(RavenLib PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>)

  target_link_libraries(RavenLib 
    PUBLIC 
      bioparser::bioparser
      cereal::cereal 
      racon::racon)

endif()
