include(FetchContent)

set(DOCTEST_NO_INSTALL ON CACHE BOOL "" FORCE)
FetchContent_Declare(doctest
  URL      https://github.com/doctest/doctest/archive/refs/tags/v2.5.3.tar.gz
  URL_HASH SHA256=174ebc4e769928959614789c5b4e9c3d0a0f81a62bb608756b127bfebfb21331
  SYSTEM)
FetchContent_MakeAvailable(doctest)

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
target_link_libraries(hf-tests PRIVATE hf_objects doctest::doctest_with_main hf_warnings)

add_test(NAME hf-tests COMMAND hf-tests)
