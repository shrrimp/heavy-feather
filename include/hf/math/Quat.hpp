#pragma once

#include <cmath>

#include "hf/math/Mat3.hpp"
#include "hf/math/Vec3.hpp"

namespace hf {

    // Unit quaternions are rotations. q1 * q2 applies q2 first.
    struct Quat {
        double w = 1.0, x = 0.0, y = 0.0, z = 0.0;

        constexpr Quat() = default;
        constexpr Quat(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}
        constexpr Quat(double w_, const Vec3 &v) : w(w_), x(v.x), y(v.y), z(v.z) {}

        constexpr bool operator==(const Quat &) const = default;
    };

    constexpr Quat operator*(const Quat &p, const Quat &q) {
        return {
            p.w * q.w - p.x * q.x - p.y * q.y - p.z * q.z,
            p.w * q.x + p.x * q.w + p.y * q.z - p.z * q.y,
            p.w * q.y + p.y * q.w + p.z * q.x - p.x * q.z,
            p.w * q.z + p.z * q.w + p.x * q.y - p.y * q.x,
        };
    }

    constexpr Quat operator+(const Quat &p, const Quat &q) { return {p.w + q.w, p.x + q.x, p.y + q.y, p.z + q.z}; }

    // Rotates v by q.
    constexpr Vec3 operator*(const Quat &q, const Vec3 &v) {
        const Vec3 u{q.x, q.y, q.z};
        const Vec3 uv = cross(u, v);
        const Vec3 uuv = cross(u, uv);
        return v + ((uv * q.w) + uuv) * 2.0;
    }

    // The inverse rotation, for a unit quaternion.
    constexpr Quat conjugate(const Quat &q) { return {q.w, -q.x, -q.y, -q.z}; }

    constexpr double dot(const Quat &a, const Quat &b) { return (a.w * b.w + a.x * b.x) + (a.y * b.y + a.z * b.z); }

    inline double length(const Quat &q) { return std::sqrt(dot(q, q)); }

    // A zero quaternion normalises to the identity.
    inline Quat normalize(const Quat &q) {
        const double len = length(q);
        if (len <= 0.0) return {};
        const double k = 1.0 / len;
        return {q.w * k, q.x * k, q.y * k, q.z * k};
    }

    // The rotation of `angle` about `axis`. The axis is used as given, not normalised.
    inline Quat angleAxis(double angle, const Vec3 &axis) {
        return {std::cos(angle * 0.5), axis * std::sin(angle * 0.5)};
    }

    constexpr Mat3 toMat3(const Quat &q) {
        const double xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
        const double xz = q.x * q.z, xy = q.x * q.y, yz = q.y * q.z;
        const double wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;

        return {
            Vec3{1.0 - 2.0 * (yy + zz), 2.0 * (xy + wz), 2.0 * (xz - wy)},
            Vec3{2.0 * (xy - wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz + wx)},
            Vec3{2.0 * (xz + wy), 2.0 * (yz - wx), 1.0 - 2.0 * (xx + yy)},
        };
    }

} // namespace hf
