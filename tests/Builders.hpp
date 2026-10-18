#pragma once

#include <cmath>

#include "hf/World.hpp"

// Small builders the world and drive tests share: bodies by their ends, joints by world points.
namespace build {

    using namespace hf;

    constexpr double kDt = 1.0 / 60.0;

    inline Mat3 boxInertia(double m, const Vec3 &half) {
        const double x = 4 * half.x * half.x, y = 4 * half.y * half.y, z = 4 * half.z * half.z;
        return Mat3(m / 12.0 * (y + z), 0, 0, 0, m / 12.0 * (x + z), 0, 0, 0, m / 12.0 * (x + y));
    }

    inline void run(World &w, double seconds) {
        for (int s = 0; s < int(std::lround(seconds / kDt)); ++s) w.step(kDt);
    }

    // The rotation taking +X onto axis.
    inline Quat xTo(const Vec3 &axis) {
        const Vec3 a = normalize(axis);
        return normalize(Quat{1.0 + a.x, cross(Vec3{1, 0, 0}, a)});
    }

    // A joint at a world point about a world axis, with the bodies where they are: their frames coincide.
    inline JointDesc jointAt(World &w, JointType type, BodyId a, BodyId b, const Vec3 &point, const Vec3 &axis) {
        const auto local = [&](BodyId body) -> Frame {
            if (!body.valid()) return {point, xTo(axis)};
            const Quat q = w.orientation(body);
            return {conjugate(q) * (point - w.position(body)), normalize(conjugate(q) * xTo(axis))};
        };
        JointDesc j;
        j.type = type;
        j.a = a;
        j.b = b;
        j.frameA = local(a);
        j.frameB = local(b);
        return j;
    }

    inline BodyId rodBetween(World &w, const Vec3 &p, const Vec3 &q, double mass) {
        const Vec3 d = q - p;
        const double half = 0.5 * length(d);
        BodyDesc b;
        b.position = 0.5 * (p + q);
        b.orientation = xTo(d);
        b.mass = mass;
        b.inertia = boxInertia(mass, Vec3{half, 0.1, 0.1});
        return w.addBody(b);
    }

} // namespace build
