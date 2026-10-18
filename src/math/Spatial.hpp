#pragma once

#include <cmath>

#include "hf/math/Mat3.hpp"
#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

// Spatial algebra in world axes, every vector taken about one reference point O. With one frame for
// everything, Featherstone's transforms are the identity and composite inertias simply add.

namespace hf {

    // A twist: angular velocity, and the velocity of the body point at O.
    struct Motion {
        Vec3 w{}, v{};
    };

    // A wrench: moment about O, and the force.
    struct Force {
        Vec3 n{}, f{};
    };

    inline Motion operator+(const Motion &a, const Motion &b) { return {a.w + b.w, a.v + b.v}; }
    inline Motion operator-(const Motion &a, const Motion &b) { return {a.w - b.w, a.v - b.v}; }
    inline Motion operator*(double s, const Motion &a) { return {s * a.w, s * a.v}; }
    inline Motion &operator+=(Motion &a, const Motion &b) {
        a.w += b.w;
        a.v += b.v;
        return a;
    }

    inline Force operator+(const Force &a, const Force &b) { return {a.n + b.n, a.f + b.f}; }
    inline Force operator-(const Force &a, const Force &b) { return {a.n - b.n, a.f - b.f}; }
    inline Force operator*(double s, const Force &a) { return {s * a.n, s * a.f}; }
    inline Force &operator+=(Force &a, const Force &b) {
        a.n += b.n;
        a.f += b.f;
        return a;
    }

    // Power: a twist against a wrench, both about the same point.
    inline double dot(const Motion &m, const Force &f) { return dot(m.w, f.n) + dot(m.v, f.f); }

    // a × b: how fast twist b changes when carried by a frame moving with a.
    inline Motion crossM(const Motion &a, const Motion &b) {
        return {cross(a.w, b.w), cross(a.w, b.v) + cross(a.v, b.w)};
    }

    // a ×* f: the same for a wrench.
    inline Force crossF(const Motion &a, const Force &f) {
        return {cross(a.w, f.n) + cross(a.v, f.f), cross(a.w, f.f)};
    }

    // Velocity of the material point at x (relative to O) on a body moving with twist V.
    inline Vec3 pointVelocity(const Motion &V, const Vec3 &x) { return V.v + cross(V.w, x); }

    // A force f applied at x, as a wrench about O.
    inline Force pointForce(const Vec3 &x, const Vec3 &f) { return {cross(x, f), f}; }

    // Spatial inertia about O: mass, first moment h = m c, rotational inertia about O.
    struct Inertia {
        double m = 0.0;
        Vec3 h{};
        Mat3 I{};

        // A body of mass m, centre of mass c (relative to O), inertia Ic about that centre in world axes.
        static Inertia body(double m, const Vec3 &c, const Mat3 &Ic) {
            Inertia r;
            r.m = m;
            r.h = m * c;
            r.I = Ic + m * (dot(c, c) * Mat3(1.0) - outerProduct(c, c));
            return r;
        }

        // Momentum of a body moving with twist x.
        Force operator*(const Motion &x) const { return {I * x.w + cross(h, x.v), m * x.v - cross(h, x.w)}; }

        Inertia &operator+=(const Inertia &o) {
            m += o.m;
            h += o.h;
            I = I + o.I;
            return *this;
        }
    };

    // exp([phi]×) as a quaternion, exact for any angle.
    inline Quat rotationVector(const Vec3 &phi) {
        const double a = length(phi);
        if (a < 1e-12) return normalize(Quat(1.0, 0.5 * phi.x, 0.5 * phi.y, 0.5 * phi.z));
        const double s = std::sin(0.5 * a) / a;
        return {std::cos(0.5 * a), s * phi.x, s * phi.y, s * phi.z};
    }

    // Its inverse: the rotation vector of q, the shorter way round.
    inline Vec3 rotationVectorOf(const Quat &q) {
        const double sign = q.w < 0.0 ? -1.0 : 1.0;
        const Vec3 v{sign * q.x, sign * q.y, sign * q.z};
        const double s = length(v);
        if (s < 1e-12) return 2.0 * v;
        return (2.0 * std::atan2(s, sign * q.w) / s) * v;
    }

} // namespace hf
