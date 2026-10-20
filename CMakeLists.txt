cmake_minimum_required(VERSION 3.25)

# The version this tree is heading to.
project(HeavyFeather
  VERSION 0.3.0
  DESCRIPTION "Rigid-body dynamics for heavy machinery and complex mechanisms"
  HOMEPAGE_URL "https://github.com/shrrimp/heavy-feather"
  LANGUAGES CXX)

# Tests and the warning wall are for working ON Heavy Feather. A host that pulls it in with
# add_subdirectory or FetchContent gets the library alone, built with its own flags.
if(PROJECT_IS_TOP_LEVEL)
  set(HF_DEV_DEFAULT ON)
else()
  set(HF_DEV_DEFAULT OFF)
endif()

option(HF_BUILD_TESTS "Build the unit tests" ${HF_DEV_DEFAULT})
option(HF_WARNINGS "Build with the full warning set" ${HF_DEV_DEFAULT})
option(HF_WARNINGS_AS_ERRORS "Turn those warnings into errors (CI)" OFF)
option(HF_INSTALL "Generate the install and package rules" ON)
# BUILD_SHARED_LIBS (CMake's own) picks a static or a shared library; static by default.

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(HFWarnings)

add_subdirectory(src)

if(HF_INSTALL)
  include(HFInstall)
endif()

if(HF_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
