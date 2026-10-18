#pragma once

#include <cstdint>

#include "hf/Types.hpp"

namespace hf {

    // One contact point, found by the host's collision detection and handed to the world for the next step.
    struct Contact {
        BodyId a, b;
        Vec3 point;                // world
        Vec3 normal;               // world, pointing from a to b
        double depth = 0.0;        // penetration along the normal; negative is a gap, closed this step and no faster
        std::uint64_t feature = 0; // identifies this point across steps, for warm starting
        double friction = 0.6;
        double restitution = 0.0;
        // Resistance to rolling across the surface and to spinning about the normal: the lever arm, in m, of the
        // normal force that resists them. A tyre's is a few millimetres; 0 leaves them free.
        double rolling = 0.0;
        double spinning = 0.0;
    };

} // namespace hf
