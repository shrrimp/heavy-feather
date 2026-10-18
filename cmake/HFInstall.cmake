# Install rules and the CMake package, so a host can do:
#
#   find_package(HeavyFeather 0.1 REQUIRED)
#   target_link_libraries(game PRIVATE HeavyFeather::HeavyFeather)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(HF_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/HeavyFeather")

install(TARGETS heavy_feather
  EXPORT HeavyFeatherTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# The hand-written headers, then the generated ones (Version.hpp, Export.hpp).
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/hf"
                  "${PROJECT_BINARY_DIR}/include/hf"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.hpp")

install(EXPORT HeavyFeatherTargets
  NAMESPACE HeavyFeather::
  DESTINATION ${HF_CMAKE_DIR})

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/HeavyFeatherConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/HeavyFeatherConfig.cmake"
  INSTALL_DESTINATION ${HF_CMAKE_DIR})

# Before 1.0 a minor version may break the API, so 0.1 only accepts 0.1.x.
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/HeavyFeatherConfigVersion.cmake"
  COMPATIBILITY SameMinorVersion)

install(FILES
  "${PROJECT_BINARY_DIR}/HeavyFeatherConfig.cmake"
  "${PROJECT_BINARY_DIR}/HeavyFeatherConfigVersion.cmake"
  DESTINATION ${HF_CMAKE_DIR})
