#pragma once

#include <cstdint>
#include <limits>

#include "hf/Types.hpp"

namespace hf {

    // A guide a belt runs round: the body a hinge turns, the belt riding on it about the hinge's axle.
    struct BeltGuide {
        JointId hinge;         // a Hinge, Motor, AngularSpring or Servo; the guide is the body it turns, its b
        double radius = 0.1;   // m: where the belt rides on it; for one that meshes, its pitch radius
        int side = 1;          // +1 the belt goes round its outside, -1 its inside
        bool meshes = false;   // toothed: it turns with the belt, never slipping
        double friction = 0.6; // against the belt, for one that doesn't mesh
    };

    // A closed loop of fixed length round its guides: a track round its rollers, a chain round its sprockets, a
    // conveyor. Its one freedom is its speed along the loop, positive the way its guides are listed. It turns the
    // guides that mesh with it and the others as far as their friction on it allows, so a motor on one guide drives
    // the rest through it.
    struct BeltDesc {
        BodyId carrier;      // what its own motion is measured against; invalid: the world, as a fixed conveyor's
        double length = 1.0; // m round the loop
        double mass = 1.0;   // kg
        // m between its teeth. A meshing guide with a whole number of teeth is kept in step with them, so a sprocket
        // never drifts against the links it carries. 0: no teeth to keep.
        double pitch = 0.0;
        double stiffness = std::numeric_limits<double>::infinity(); // N: its EA; infinite, it never stretches
        double damping = 1.0;                                        // a stretchy belt's damping ratio
        double strength = std::numeric_limits<double>::infinity();  // N: past it the belt breaks
    };

    // Where a belt ended up as the last step solved it.
    struct BeltInfo {
        int machine = -1;  // the machine whose solve holds it; -1: nothing it touches moves
        bool taut = false; // pulling
        int slipping = 0;  // guides slipping on it
        bool broken = false;
    };

} // namespace hf
