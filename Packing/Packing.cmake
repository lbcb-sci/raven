include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (NOT DEFINED raven_INSTALL_CMAKEDIR) 
  set(raven_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/raven"
    CACHE STRING "Path to raven cmake files")
endif ()

set(RAVEN_TARGET_LIST 
  raven 
    bioparser
    biosoup 
    racon 
    ram
    spoa
    thread_pool)

if (RAVEN_BUILD_EXE)
  list(APPEND RAVEN_TARGET_LIST raven_exe)
endif ()

install(TARGETS ${RAVEN_TARGET_LIST}  EXPORT raven_TARGETS
        RUNTIME COMPONENT raven_Runtime
        LIBRARY COMPONENT raven_Runtime
        NAMELINK_COMPONENT raven_Development
        ARCHIVE COMPONENT raven_Development
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(DIRECTORY "${raven_SOURCE_DIR}/RavenLib/include/" "${raven_BINARY_DIR}/include/"
        TYPE INCLUDE
        COMPONENT raven_Development)

if (BUILD_SHARED_LIBS)
  set(type shared)
else ()
  set(tpe static)
endif ()

install(EXPORT raven_TARGETS 
        DESTINATION "${raven_INSTALL_CMAKEDIR}"
        NAMESPACE raven
        FILE raven-${type}-targets.cmake
        COMPONENT raven_Development)

write_basic_package_version_file(
  ravenConfigVersion.cmake 
  COMPATIBILITY SameMajorVersion)

install(FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/ravenConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/ravenConfigVersion.cmake"
        DESTINATION "${raven_INSTALL_CMAKEDIR}"
        COMPONENT raven_Development)

include(CPack)
