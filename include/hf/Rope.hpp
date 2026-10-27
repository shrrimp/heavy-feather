#pragma once

#include <cstdint>
#include <limits>

#include "hf/Types.hpp"

namespace hf {

    // A point a rope is routed through: one of its ends, a sheave it runs round, an eye it slides through.
    struct RopePoint {
        BodyId body;     // invalid: the world, and position is then in world coordinates
        Vec3 position{}; // from the body's origin, in its axes: an end's point, or a sheave's or drum's centre
        // A sheave the rope runs round, or a drum it leaves at an end: the axle, in the body's axes, and its radius.
        // Radius 0 is a point, an eye the rope slides through.
        Vec3 axis{1.0, 0.0, 0.0};
        double radius = 0.0;
        // Which way round the axle the rope goes: +1 anticlockwise about it, -1 clockwise. 0: the way the route bends
        // there when the rope is added.
        int side = 0;
    };

    enum class RopeSpool : std::uint8_t { None, Start, End, Both };

    // A rope, a cable, a chain: it only pulls, and only once it is taut, along its whole route at once. Pulleys
    // between its ends turn its pull back on itself, so a block on four strands of it is pulled with four times its
    // tension. Its tension is the same all along it: its sheaves turn freely.
    struct RopeDesc {
        double length = 1.0; // m: paid out between its ends
        // N: its axial stiffness, EA; infinite, it doesn't stretch at all.
        double stiffness = std::numeric_limits<double>::infinity();
        double damping = 1.0; // a stretchy rope's damping ratio
        // N: it breaks once its tension goes above this, and pulls nothing from then on. See World::brokenRopes.
        double strength = std::numeric_limits<double>::infinity();

        // A spool at an end winds the rope in and pays it out, worked by World::setSpoolLever: its motor pulls no
        // harder than maxForce, and idle its brake holds up to brakeForce. Either one overloaded, the rope pays out
        // under the load instead of stretching.
        RopeSpool spool = RopeSpool::None;
        double speed = 0.5;                                          // m/s each spool, at a full lever
        double maxForce = std::numeric_limits<double>::infinity();   // N
        double brakeForce = std::numeric_limits<double>::infinity(); // N
        // What the spools hold of it: never wound in shorter than minLength, and paid out to maxLength it is tied to
        // the drum and holds whatever pulls.
        double minLength = 0.0;
        double maxLength = std::numeric_limits<double>::infinity();
    };

    // Where a rope ended up as the last step solved it.
    struct RopeInfo {
        int machine = -1;      // held exactly with that machine's joints; -1: swept between bodies or machines
        bool taut = false;     // pulling
        bool slipping = false; // a spool overloaded: paying out under the load
        bool broken = false;
    };

    // The rope at one point of its route, as the last step left it: where it arrives and where it leaves, which are
    // the point itself unless it runs round a sheave. A sheave the rope has come off doesn't touch it.
    struct RopeBend {
        Vec3 in{}, out{};
        bool touching = true;
    };

} // namespace hf
