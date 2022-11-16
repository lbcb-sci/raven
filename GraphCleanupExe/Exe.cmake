if (NOT TARGET GraphCleanupExe)

    include(${CMAKE_CURRENT_LIST_DIR}/Exe.srcs.cmake)

    add_executable(graph_cleanup_exe ${SOURCES})
    target_link_libraries(graph_cleanup_exe PRIVATE raven)

    if (racon_enable_cuda)
        target_compile_definitions(graph_cleanup_exe PRIVATE CUDA_ENABLED)
    endif ()

    target_compile_definitions(graph_cleanup_exe PRIVATE VERSION="${PROJECT_VERSION}")
    set_property(TARGET graph_cleanup_exe PROPERTY OUTPUT_NAME graph_cleanup)

endif()
