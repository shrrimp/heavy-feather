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
- everything goes in namespace `hf`

Comments are rare: a short one on public functions and types when the name isn't enough. The code
should read on its own.

## Tests

One doctest file per topic, `tests/<Topic>Tests.cpp`, with each check saying what it shows:
`CHECK_MESSAGE(condition, "what it shows")`. Expected values come from a hand derivation or an independent
referee (another algorithm, a finite difference), never from running the code under test.

`tests/host/` is a small game around the library: collision for boxes, spheres and planes, mixed stress scenes,
and a count of every allocation. The stress, memory and determinism tests run on it, and so does `hf-bench`, which
times the large scenes; build it in release and run it with `--check` to hold each scene to its budget.
