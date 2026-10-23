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

v0.3: machines built from bodies and joints and solved as articulations, with motors, pistons and servos that work like an operator's levers, all behind the `hf::World` API.

Contacts come with friction and rolling resistance, and can be reported before they touch, so fast bodies stop on thin walls instead of passing through. A contact marked as part of a mechanism, a cam on its follower or a gear on its gear, is solved exactly with its machine. Bodies that stay still fall asleep and cost nothing until something wakes them.

Steps are deterministic and allocate nothing once warm.

It started as a prototype loosely implemented in Heavy Metal, which worked well enough that I decided to turn it into a standalone library. v0.1 was that prototype's machine solver, ported and cleaned up, and the versions since have built on it.

## Roadmap
(vague prediction)

- **v0.1:** Machines, motors, contacts, friction, and a simple API (done)
- **v0.2:** Enhanced contacts for mechanisms: cams, gears and other contact-driven parts (done)
- **v0.3:** CCD, sub-steps, simulation islands and sleeping (done)
- **v0.4:** Robust, and moved by the host: a fault tripwire that rolls an island back before a NaN leaves it, inputs checked at the API, and one call for each way a host moves a body: kinematic bodies, whole machines moved rigidly, anchors that carry what hangs on them
- **v0.5:** Performance and multi-threading: each island solved on its own budget, lean contact rows, resting contacts kept in the library, no whole-world rebuilds; joints posed directly and a mouse grab; then islands solved in parallel
- **v0.6:** Break stuff: body frames apart from the centre of mass, breakable joints, and load analysis for rigid-body rupture
- **v1.0:** Clean API, docs, and a basic host for testing and demos
- **Later:** A follower-on-profile joint, for cams that are exact without any collision geometry

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
