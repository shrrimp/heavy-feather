#pragma once

#include <cmath>

#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

// What a value has to be for the world to take it.
namespace hf {

    inline bool finite(const Vec3 &v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

    // A rotation the world can normalise: finite, and not of almost no length.
    inline bool rotation(const Quat &q) {
        return std::isfinite(q.w) && std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z) &&
               length(q) >= 1e-6;
    }

} // namespace hf
