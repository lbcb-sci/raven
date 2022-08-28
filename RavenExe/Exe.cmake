if (NOT TARGET RavenExe)

    include(${CMAKE_CURRENT_LIST_DIR}/Exe.srcs.cmake)

    add_executable(raven_exe ${SOURCES})
    target_link_libraries(raven_exe PRIVATE cxxopts::cxxopts raven)

    if (racon_enable_cuda)
        target_compile_definitions(raven_exe PRIVATE CUDA_ENABLED)
    endif ()

    target_compile_definitions(raven_exe PRIVATE VERSION="${PROJECT_VERSION}")
    set_property(TARGET raven_exe PROPERTY OUTPUT_NAME raven)

    configure_file(
      "${CMAKE_CURRENT_LIST_DIR}/src/raven_cfg.h.in" raven_cfg.h)
    target_include_directories(raven_exe PRIVATE ${PROJECT_BINARY_DIR})
endif()
