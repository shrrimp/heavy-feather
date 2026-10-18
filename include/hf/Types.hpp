#pragma once

#include <cstdint>

#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

namespace hf {

    // Handles: a slot and the generation it was issued in, so a handle to something removed stays invalid.
    struct BodyId {
        std::uint32_t index = ~0u;
        std::uint32_t generation = 0;

        bool valid() const { return index != ~0u; }
        bool operator==(const BodyId &) const = default;
    };

    struct JointId {
        std::uint32_t index = ~0u;
        std::uint32_t generation = 0;

        bool valid() const { return index != ~0u; }
        bool operator==(const JointId &) const = default;
    };

    // A frame on a body, relative to its centre of mass in body axes. A joint's axis is its frame's X.
    struct Frame {
        Vec3 position{};
        Quat rotation{};
    };

} // namespace hf
