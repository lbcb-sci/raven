if (NOT TARGET Demo)

    include(${CMAKE_CURRENT_LIST_DIR}/Demo.srcs.cmake)

    add_executable(Raven ${SOURCES})
    target_link_libraries(Raven PRIVATE RavenLib bioparser::bioparser)

    if (racon_enable_cuda)
        target_compile_definitions(Raven PRIVATE CUDA_ENABLED)
    endif ()

    target_compile_definitions(Raven PRIVATE VERSION="${PROJECT_VERSION}")
    install(TARGETS Raven DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
