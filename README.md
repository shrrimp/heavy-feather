<p align="center">
  <img src="assets/heavy-feather.png" alt="Heavy Feather" width="180">
</p>

<h1 align="center">Heavy Feather</h1>

<p align="center">
  A physics library for dynamic heavy machinery and complex mechanisms.
</p>

---

Heavy Feather is a C++ rigid-body dynamics library for things built from parts that push hard on each other: linkages, hydraulic arms, cranes, presses, drivetrains. Each machine is solved as one [Featherstone](https://en.wikipedia.org/wiki/Featherstone%27s_algorithm) articulation, so its joints hold exactly under load and every part keeps its real mass, however lopsided the masses get.

It doesn't depend on any engine and for now the plan is for it to do no collision detection. You give it bodies, joints and contacts, and it does the dynamics.

## Why

**Heavy Metal**, my voxel physics game full of heavy machinery, is what drove me to make it, because I found no open-source physics engine that could do that and fast.

<p align="center">
  <a href="https://www.youtube.com/watch?v=VFehGWf1e1U">
    <img src="https://img.youtube.com/vi/VFehGWf1e1U/maxresdefault.jpg" alt="Heavy Metal - Voxel Physics Demo" width="560">
  </a>
</p>

## Where it's at

There's no physics in the repo yet, just the project around it: the build, the tests and the docs. A prototype of the design is loosely implemented in Heavy Metal, and it worked well enough that I decided to turn it into a standalone library, written from scratch.

## Roadmap
(vague prediction)

- **v0.1:** Machines, motors, contacts, friction, and a simple API
- **v0.2:** Enhanced contacts for mechanisms: cams, gears and other contact-driven parts
- **v0.3:** CCD, sub-steps, simulation islands and sleeping
- **v0.4:** Break stuff: breakable joints, and load analysis for rigid-body rupture
- **v0.5:** Performance and multi-threading
- **v1.0:** Clean API, docs, and a basic host for testing and demos

## Building

C++20 and CMake 3.25+, no dependencies.

```sh
cmake --workflow --preset debug   # configure, build, test
```

To use it from your own project, `find_package(HeavyFeather)` or FetchContent, then link
`HeavyFeather::HeavyFeather`. [docs/building.md](docs/building.md) has the details.

## Docs

- [Architecture](docs/architecture.md): how it works, and why
- [Building](docs/building.md): presets, options, and using it from your project
- [Contributing](CONTRIBUTING.md)

## License

[MIT](LICENSE)
