#pragma once

#include <cmath>

#include "hf/math/Vec3.hpp"

// What a contact is, wherever it's solved: swept beside the machines, or inside one machine's direct solve.
namespace hf {

    constexpr double kContactHertz = 120.0;
    constexpr double kContactDamping = 4.0; // a larger ratio makes the row stiffer: 4 is the quiet point
    constexpr double kSlop = 0.002;         // m of penetration left alone
    constexpr double kGapReach = 2.0;       // a gap is solved once its approach could close it this many times over
    constexpr double kMaxBias = 3.0;        // m/s: the fastest penetration is pushed out
    constexpr double kRestitutionMinSpeed = 0.5;

    // A stable tangent pair for any unit normal; (n, t1, t2) is right-handed.
    inline void contactBasis(const Vec3 &n, Vec3 &t1, Vec3 &t2) {
        const double s = std::copysign(1.0, n.z);
        const double a = -1.0 / (s + n.z);
        const double b = n.x * n.y * a;
        t1 = {1.0 + s * n.x * n.x * a, s * b, -s * n.x};
        t2 = {b, s + n.y * n.y * a, -n.y};
    }

} // namespace hf
