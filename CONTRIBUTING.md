# Contributing

For anything bigger than a small fix, open an issue first so we can talk it through.

## Before a PR

```sh
cmake --workflow --preset debug
cmake --workflow --preset asan
```

Run `clang-format` on what you touched. CI builds with GCC, Clang and MSVC, and treats
warnings as errors.

## Rules

- SI units, right-handed, no built-in up axis. Quaternions are `(w, x, y, z)`, body to world.
- Double precision everywhere.
- Real masses. Never clamp a mass ratio or inflate an inertia to make a solve behave.
- The library doesn't throw, and a step doesn't allocate once it's warmed up.
- No dependencies in the library. Tests can have them.

## Style

`.clang-format` handles layout. For naming:

- `PascalCase` for types and files, `camelCase` for functions and variables
- `member_` for private members, `kName` for constants, `HF_` for macros
- everything public goes in namespace `hf`; tests go in `tests/<module>_test.cpp`

Comments say why, not what.
