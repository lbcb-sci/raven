if (NOT TARGET RavenLib)

  include(GenerateExportHeader)
  include(${CMAKE_CURRENT_LIST_DIR}/RavenLib.srcs.cmake)

  add_library(raven ${SOURCES})
  add_library(${PROJECT_NAME}::raven ALIAS raven)
  set_target_properties(raven PROPERTIES
    VERSION ${raven_VERSION} 
    SOVERSION ${raven_VERSION_MAJOR})

  if (RAVEN_MAIN_PROJECT OR RAVEN_BUILD_PYTHON)
    set(ASAN_FLAGS -fno-omit-frame-pointer -fsanitize=address)
    set(DEBUG_WANINGS -Wall -Wextra -pedantic)

    target_compile_options(raven PUBLIC ${DEBUG_WANINGS})
    target_compile_options(raven PUBLIC 
      "$<$<CONFIG:Debug>:${ASAN_FLAGS}>"
      "$<$<CONFIG:RelWithDebInfo>:${ASAN_FLAGS}>")

    target_link_options(raven PUBLIC
      "$<$<CONFIG:Debug>:${ASAN_FLAGS}>"
      "$<$<CONFIG:RelWithDebInfo>:${ASAN_FLAGS}>")
  endif()

  target_include_directories(raven PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>)

  generate_export_header(raven EXPORT_FILE_NAME include/raven/export.h)
  target_compile_definitions(raven PUBLIC
    "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:RAVEN_STATIC_DEFINE>")
  target_include_directories(raven PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")

  target_link_libraries(raven 
    PUBLIC 
      bioparser::bioparser
      cereal::cereal 
      ram::ram
      racon::racon)
endif ()
