if (NOT TARGET RavenExe)

    include(${CMAKE_CURRENT_LIST_DIR}/Exe.srcs.cmake)

    add_executable(raven ${SOURCES})
    target_link_libraries(raven PRIVATE RavenLib)

    if (racon_enable_cuda)
        target_compile_definitions(raven PRIVATE CUDA_ENABLED)
    endif ()

    target_compile_definitions(raven PRIVATE VERSION="${PROJECT_VERSION}")
    install(TARGETS raven DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
