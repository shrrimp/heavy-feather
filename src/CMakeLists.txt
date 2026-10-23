# The library. Public headers live in include/hf/; everything else lives here, with a private
# header beside the .cpp that owns it.
#
# Every source compiles ONCE, into hf_objects. The shipped library is built from those objects,
# and the unit tests link them directly: a shared build hides everything not marked HF_API, and
# the tests have to reach the internals (the LCP, the articulation) all the same.

include(GNUInstallDirs)

set(HF_SOURCES
  Version.cpp
  contact/Contacts.cpp
  contact/DirectContacts.cpp
  contact/MachineSides.cpp
  dynamics/Articulation.cpp
  math/Dense.cpp
  solver/BoxLcp.cpp
  solver/Machine.cpp
  world/Faults.cpp
  world/Groups.cpp
  world/MachineRows.cpp
  world/Machines.cpp
  world/Sleep.cpp
  world/Step.cpp
  world/World.cpp)

add_library(hf_objects OBJECT ${HF_SOURCES})
add_library(heavy_feather $<TARGET_OBJECTS:hf_objects>)
add_library(HeavyFeather::HeavyFeather ALIAS heavy_feather)

# ---- generated headers: hf/Version.hpp, and hf/Export.hpp for HF_API ----
set(HF_GENERATED_INCLUDE "${PROJECT_BINARY_DIR}/include")
configure_file("${PROJECT_SOURCE_DIR}/cmake/Version.hpp.in"
               "${HF_GENERATED_INCLUDE}/hf/Version.hpp" @ONLY)

include(GenerateExportHeader)
set_target_properties(heavy_feather PROPERTIES DEFINE_SYMBOL HF_BUILDING)
generate_export_header(heavy_feather
  BASE_NAME HF
  EXPORT_MACRO_NAME HF_API
  EXPORT_FILE_NAME "${HF_GENERATED_INCLUDE}/hf/Export.hpp")

# ---- the objects: how every source compiles ----
target_include_directories(hf_objects PUBLIC
  "${PROJECT_SOURCE_DIR}/include" "${HF_GENERATED_INCLUDE}" "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(hf_objects PUBLIC cxx_std_20)
target_link_libraries(hf_objects PRIVATE hf_warnings)

# Symbols are hidden unless marked HF_API, so a shared build exports the API and nothing else,
# and Linux behaves the way Windows has to.
set_target_properties(hf_objects PROPERTIES
  CXX_EXTENSIONS OFF
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  POSITION_INDEPENDENT_CODE "${BUILD_SHARED_LIBS}")   # quoted: unset unless a host sets it

# HF_API is dllexport while the objects compile for a shared library, and nothing when they are
# linked straight into something (a static library, the tests).
if(BUILD_SHARED_LIBS)
  target_compile_definitions(hf_objects PRIVATE HF_BUILDING INTERFACE HF_STATIC_DEFINE)
else()
  target_compile_definitions(hf_objects PUBLIC HF_STATIC_DEFINE)
  target_compile_definitions(heavy_feather PUBLIC HF_STATIC_DEFINE)
endif()

# ---- the shipped library: what a host sees ----
target_include_directories(heavy_feather PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${HF_GENERATED_INCLUDE}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(heavy_feather PUBLIC cxx_std_20)

# Pre-1.0, every minor version may break the ABI.
set_target_properties(heavy_feather PROPERTIES
  OUTPUT_NAME heavy-feather
  EXPORT_NAME HeavyFeather   # installed as HeavyFeather::HeavyFeather, same as the alias
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
