# Unit tests: one doctest executable, registered with CTest. doctest is fetched here and only
# here, so it is a dependency of the tests and never of the library.

include(FetchContent)

set(DOCTEST_NO_INSTALL ON CACHE BOOL "" FORCE)
FetchContent_Declare(doctest
  URL      https://github.com/doctest/doctest/archive/refs/tags/v2.5.3.tar.gz
  URL_HASH SHA256=174ebc4e769928959614789c5b4e9c3d0a0f81a62bb608756b127bfebfb21331
  SYSTEM)
FetchContent_MakeAvailable(doctest)

# hf_objects, not the shipped library: the tests reach internals that a shared build hides.
# The package check (tests/package/) is what exercises the library the way a host links it.
add_executable(hf-tests
  version_test.cpp)
target_link_libraries(hf-tests PRIVATE hf_objects doctest::doctest_with_main hf_warnings)

add_test(NAME hf-tests COMMAND hf-tests)
