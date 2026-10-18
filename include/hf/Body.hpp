#pragma once

#include "hf/math/Mat3.hpp"
#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

namespace hf {

    struct BodyDesc {
        Vec3 position{}; // of the centre of mass, world
        Quat orientation{};
        Vec3 velocity{}; // of the centre of mass
        Vec3 angularVelocity{};

        double mass = 1.0;
        Mat3 inertia{1.0}; // about the centre of mass, body axes

        // Bounding radius about the centre of mass, which caps how far a position correction may move the body in
        // one sub-step. 0: estimated from the inertia.
        double radius = 0.0;

        bool isStatic = false;
    };

} // namespace hf
