#pragma once

#include <cmath>

namespace hf {

    struct Vec3 {
        double x = 0.0, y = 0.0, z = 0.0;

        constexpr Vec3() = default;
        constexpr explicit Vec3(double s) : x(s), y(s), z(s) {}
        constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

        constexpr double &operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
        constexpr double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }

        constexpr bool operator==(const Vec3 &) const = default;

        constexpr Vec3 operator-() const { return {-x, -y, -z}; }

        constexpr Vec3 &operator+=(const Vec3 &v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        constexpr Vec3 &operator-=(const Vec3 &v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }

        constexpr Vec3 &operator*=(double s) {
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }

        constexpr Vec3 &operator/=(double s) {
            x /= s;
            y /= s;
            z /= s;
            return *this;
        }
    };

    constexpr Vec3 operator+(Vec3 a, const Vec3 &b) { return a += b; }
    constexpr Vec3 operator-(Vec3 a, const Vec3 &b) { return a -= b; }
    constexpr Vec3 operator*(Vec3 v, double s) { return v *= s; }
    constexpr Vec3 operator*(double s, Vec3 v) { return v *= s; }
    constexpr Vec3 operator/(Vec3 v, double s) { return v /= s; }

    constexpr double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

    constexpr Vec3 cross(const Vec3 &a, const Vec3 &b) {
        return {a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y};
    }

    inline double length(const Vec3 &v) { return std::sqrt(dot(v, v)); }

    // No zero check: a zero vector gives NaNs.
    inline Vec3 normalize(const Vec3 &v) { return v * (1.0 / std::sqrt(dot(v, v))); }

} // namespace hf
