#pragma once

#include <cmath>

#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

// What a value has to be for the world to take it.
namespace hf {

    inline bool finite(const Vec3 &v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

    // At or above zero, and a number; infinite only where infinite means no limit.
    inline bool nonNegative(double x, bool infinite = false) { return x >= 0.0 && (infinite || std::isfinite(x)); }
    inline bool positive(double x) { return x > 0.0 && std::isfinite(x); }

    // A rotation the world can normalise: finite, and not of almost no length.
    inline bool rotation(const Quat &q) {
        return std::isfinite(q.w) && std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z) &&
               length(q) >= 1e-6;
    }

} // namespace hf
