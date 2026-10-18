include(FetchContent)

set(DOCTEST_NO_INSTALL ON CACHE BOOL "" FORCE)
FetchContent_Declare(doctest
  URL      https://github.com/doctest/doctest/archive/refs/tags/v2.5.3.tar.gz
  URL_HASH SHA256=174ebc4e769928959614789c5b4e9c3d0a0f81a62bb608756b127bfebfb21331
  SYSTEM)
FetchContent_MakeAvailable(doctest)

# A small game around the library: collision, scenes, and an allocation counter. An object library, so the
# counting operator new is always linked in.
add_library(hf-test-host OBJECT
  host/Allocations.cpp
  host/Collision.cpp
  host/Host.cpp
  host/Scenes.cpp)
target_include_directories(hf-test-host PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hf-test-host PUBLIC hf_objects PRIVATE hf_warnings)

set(HF_TEST_SOURCES
  ArticulationTests.cpp
  BoxLcpTests.cpp
  DenseTests.cpp
  MachineTests.cpp
  SpatialTests.cpp
  VersionTests.cpp
  WorldTests.cpp)

# hf_objects rather than the shipped library: a shared build hides the internals these test.
add_executable(hf-tests ${HF_TEST_SOURCES})
target_link_libraries(hf-tests PRIVATE hf-test-host hf_objects doctest::doctest_with_main hf_warnings)
add_test(NAME hf-tests COMMAND hf-tests)
