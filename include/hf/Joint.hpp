#pragma once

#include <cstdint>
#include <limits>

#include "hf/Types.hpp"

namespace hf {

    enum class JointType : std::uint8_t {
        Fixed,         // a weld: the two bodies move as one
        Ball,          // 3 rotations free
        Hinge,         // 1 rotation free, about the axis
        Motor,         // a hinge driven to an angular velocity
        Slider,        // 1 translation free along the axis, within [min, max]
        Piston,        // a slider driven to a position
        Spring,        // a soft distance between the two anchors
        AngularSpring, // a hinge pulled softly toward angle 0
    };

    struct JointDesc {
        JointType type = JointType::Hinge;
        BodyId a, b;          // one may be invalid: the world, whose frame is then in world coordinates
        Frame frameA, frameB; // where the joint sits on each body; they coincide at rest, and X is the axis

        // Drives. Motor: angular velocity in rad/s. Piston: position along the axis in m.
        double target = 0.0;
        double maxForce = std::numeric_limits<double>::infinity();  // N, or N·m for a motor
        double forceSlew = std::numeric_limits<double>::infinity(); // how fast the drive force may grow, per second
        double speed = 1.0;                                         // Piston: how fast the target may move, m/s

        double min = -1.0, max = 1.0; // Slider and Piston: the course, in m

        double hertz = 0.0, damping = 1.0; // Spring and AngularSpring
        double restLength = 0.0;           // Spring
    };

    // Where a joint ended up once the world built its machines.
    struct JointInfo {
        int machine = -1;       // -1: a weld, or a joint between two static bodies
        bool tree = false;      // held exactly by the articulation's tree; false: a loop's cut, held by rows
        bool duplicate = false; // a second hinge on an axle another joint already holds
    };

} // namespace hf
