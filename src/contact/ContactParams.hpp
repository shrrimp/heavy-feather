#pragma once

#include <algorithm>
#include <cmath>

#include "hf/Contact.hpp"
#include "hf/math/Vec3.hpp"

// What a contact is, wherever it's solved: swept beside the machines, or inside one machine's direct solve.
namespace hf {

    constexpr double kContactHertz = 120.0;
    constexpr double kContactDamping = 4.0; // a larger ratio makes the row stiffer: 4 is the quiet point
    constexpr double kSlop = 0.002;         // m of penetration left alone
    constexpr double kGapReach = 2.0;       // a gap is solved once its approach could close it this many times over
    constexpr double kMaxBias = 3.0;        // m/s: the fastest penetration is pushed out
    constexpr double kRestitutionMinSpeed = 0.5;

    inline bool compliant(const Contact &k) { return k.stiffness > 0.0 && !k.mechanism; }
    inline double sinkReach(const Contact &k) { return compliant(k) ? k.sinkLimit : 0.0; }

    // A row's spring, from a frequency and a damping ratio.
    struct SoftParams {
        double biasRate = 0.0;
        double massCoeff = 1.0;
        double impulseCoeff = 0.0;

        static SoftParams rigid() { return {}; }
        static SoftParams fromHertz(double hertz, double zeta, double h);
        // A spring of stiffness k, N/m, on a row of effective inverse mass kn, anchored at no push at all.
        static SoftParams fromSpring(double k, double zeta, double h, double kn);
    };

    // A compliant contact's rigid surface is its sink limit: within it the row is a spring, short of it a speculative
    // gap, as any other contact is short of its touch.
    struct RowKind {
        double excess;
        bool spring;
        bool gap;
    };

    inline RowKind rowKind(const Contact &k, double depth) {
        const double excess = depth - sinkReach(k);
        const bool spring = compliant(k) && depth >= 0.0 && excess <= 0.0;
        return {excess, spring, excess < 0.0 && !spring};
    }

    // Penetration is recovered in the position pass, where it moves the pose without leaving momentum. A spring's push
    // is momentum, so it is solved with the velocities and never re-anchored.
    struct RowTargets {
        SoftParams soft;
        double vtarget;
        double posBias;
    };

    inline RowTargets rowTargets(const RowKind &kind, const Contact &k, double depth, double h, double kn,
                                 const SoftParams &resting) {
        if (kind.spring) {
            const SoftParams s = SoftParams::fromSpring(k.stiffness, k.damping, h, kn);
            return {s, s.biasRate * depth, 0.0};
        }
        if (kind.gap) return {SoftParams::rigid(), kind.excess / h, 0.0};
        return {resting, 0.0, std::min(kMaxBias, resting.biasRate * std::max(0.0, kind.excess - kSlop))};
    }

    // Relaxed to no approach, a spring would stop the very sink it carries the load through: it keeps its push.
    inline double relaxedPush(double ln, double vn, double vtarget, double kn, bool spring) {
        return spring ? ln : std::max(0.0, ln + (vtarget - vn) / kn);
    }

    // A stable tangent pair for any unit normal; (n, t1, t2) is right-handed.
    inline void contactBasis(const Vec3 &n, Vec3 &t1, Vec3 &t2) {
        const double s = std::copysign(1.0, n.z);
        const double a = -1.0 / (s + n.z);
        const double b = n.x * n.y * a;
        t1 = {1.0 + s * n.x * n.x * a, s * b, -s * n.x};
        t2 = {b, s + n.y * n.y * a, -n.y};
    }

} // namespace hf
