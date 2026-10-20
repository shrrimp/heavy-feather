# Building

You need C++20 and CMake 3.25+. There are no dependencies. CI builds with GCC, Clang and MSVC.

## Working on it

```sh
cmake --workflow --preset debug   # configure, build, test
```

Each preset builds into `build/<preset>/`:

- `debug`: for day-to-day work, and the build clangd reads
- `release`: optimised
- `asan`: AddressSanitizer + UBSan
- `shared`: a shared library, to check that only `HF_API` symbols get exported

The options are `HF_BUILD_TESTS`, `HF_WARNINGS`, `HF_WARNINGS_AS_ERRORS`, `HF_INSTALL`, and
CMake's own `BUILD_SHARED_LIBS` (static by default). Tests and warnings are only on when Heavy
Feather is the top-level project, so a game that pulls it in just gets the library.

## Using it

```cmake
find_package(HeavyFeather 0.3 REQUIRED)   # or FetchContent / add_subdirectory
target_link_libraries(game PRIVATE HeavyFeather::HeavyFeather)
```

```cpp
#include <hf/hf.hpp>
```

Until 1.0, asking for 0.3 gets you any 0.3.x and nothing else.

To install, run `cmake --install build/release --prefix <dir>`. `tests/package/` is a tiny host
that builds against an install, and CI runs it.

## Under the hood

The sources compile once, into an object library. The shipped library is built from it, and
the tests link it directly, so they can reach internals that a shared build hides.
`hf/Version.hpp` and `hf/Export.hpp` are generated into the build directory and installed with
the rest.
