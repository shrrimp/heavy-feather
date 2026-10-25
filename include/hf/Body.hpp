#pragma once

#include "hf/math/Mat3.hpp"
#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

namespace hf {

    // A body's pose is its origin's: the point in it the host measures from, whose joints' frames and grabs are
    // placed from it, and which an edit of its mass leaves where it is. Its centre of mass sits at centerOfMass from
    // there; its velocity is that centre's.
    struct BodyDesc {
        Vec3 position{}; // of the origin, world
        Quat orientation{};
        Vec3 velocity{}; // of the centre of mass
        Vec3 angularVelocity{};

        double mass = 1.0;
        Mat3 inertia{1.0};   // about the centre of mass, body axes
        Vec3 centerOfMass{}; // from the origin, body axes

        // Bounding radius about the centre of mass, which caps how far a position correction may move the body in
        // one sub-step. 0: estimated from the inertia.
        double radius = 0.0;

        bool isStatic = false;
        // Never turns: contacts, joints and torques push it about, but its orientation stays as set, and a pose write
        // is the only thing that changes it. A character standing upright. A weld group holding one turns with none.
        // Not held as a machine's link: there its joints turn it like any other.
        bool fixedRotation = false;
    };

    // A body's pose and motion, read at once.
    struct BodyState {
        Vec3 position{}; // of the origin
        Quat orientation{};
        Vec3 velocity{}; // of the centre of mass
        Vec3 angularVelocity{};
    };

} // namespace hf
