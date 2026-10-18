#pragma once

#include <cmath>

#include "hf/math/Vec3.hpp"

namespace hf {

    // Stored by columns: m[c][r] is row r of column c.
    struct Mat3 {
        Vec3 col[3]{};

        constexpr Mat3() = default;
        constexpr explicit Mat3(double d) : col{{d, 0.0, 0.0}, {0.0, d, 0.0}, {0.0, 0.0, d}} {}
        constexpr Mat3(const Vec3 &c0, const Vec3 &c1, const Vec3 &c2) : col{c0, c1, c2} {}
        constexpr Mat3(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2,
                       double z2)
            : col{{x0, y0, z0}, {x1, y1, z1}, {x2, y2, z2}} {}

        constexpr Vec3 &operator[](int c) { return col[c]; }
        constexpr const Vec3 &operator[](int c) const { return col[c]; }

        constexpr bool operator==(const Mat3 &) const = default;
    };

    constexpr Vec3 operator*(const Mat3 &m, const Vec3 &v) { return m[0] * v.x + m[1] * v.y + m[2] * v.z; }
    constexpr Mat3 operator*(const Mat3 &a, const Mat3 &b) { return {a * b[0], a * b[1], a * b[2]}; }
    constexpr Mat3 operator*(const Mat3 &m, double s) { return {m[0] * s, m[1] * s, m[2] * s}; }
    constexpr Mat3 operator*(double s, const Mat3 &m) { return m * s; }
    constexpr Mat3 operator+(const Mat3 &a, const Mat3 &b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
    constexpr Mat3 operator-(const Mat3 &a, const Mat3 &b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
    constexpr Mat3 operator-(const Mat3 &m) { return {-m[0], -m[1], -m[2]}; }

    constexpr Mat3 transpose(const Mat3 &m) {
        return {m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]};
    }

    // a bᵀ
    constexpr Mat3 outerProduct(const Vec3 &a, const Vec3 &b) { return {a * b[0], a * b[1], a * b[2]}; }

    constexpr double determinant(const Mat3 &m) {
        return +m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - m[1][0] * (m[0][1] * m[2][2] - m[2][1] * m[0][2]) +
               m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]);
    }

    // The rotation of `angle` about `axis`. The axis is normalised.
    inline Mat3 rotationMatrix(double angle, const Vec3 &axis) {
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        const Vec3 a = normalize(axis);
        const Vec3 t = (1.0 - c) * a;

        return {
            c + t[0] * a[0],        t[0] * a[1] + s * a[2], t[0] * a[2] - s * a[1],
            t[1] * a[0] - s * a[2], c + t[1] * a[1],        t[1] * a[2] + s * a[0],
            t[2] * a[0] + s * a[1], t[2] * a[1] - s * a[0], c + t[2] * a[2],
        };
    }

    // No singularity check: callers test the determinant first.
    constexpr Mat3 inverse(const Mat3 &m) {
        const double k = 1.0 / determinant(m);
        Mat3 r;
        r[0][0] = +(m[1][1] * m[2][2] - m[2][1] * m[1][2]) * k;
        r[1][0] = -(m[1][0] * m[2][2] - m[2][0] * m[1][2]) * k;
        r[2][0] = +(m[1][0] * m[2][1] - m[2][0] * m[1][1]) * k;
        r[0][1] = -(m[0][1] * m[2][2] - m[2][1] * m[0][2]) * k;
        r[1][1] = +(m[0][0] * m[2][2] - m[2][0] * m[0][2]) * k;
        r[2][1] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]) * k;
        r[0][2] = +(m[0][1] * m[1][2] - m[1][1] * m[0][2]) * k;
        r[1][2] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]) * k;
        r[2][2] = +(m[0][0] * m[1][1] - m[1][0] * m[0][1]) * k;
        return r;
    }

} // namespace hf
