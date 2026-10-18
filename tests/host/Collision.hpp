#pragma once

#include <cstdint>
#include <vector>

#include "hf/Contact.hpp"
#include "hf/math/Mat3.hpp"
#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

// The test host's collision detection: what a game would do, kept to boxes, spheres and planes.
namespace host {

    enum class ShapeKind { Box, Sphere, Plane };

    struct Shape {
        hf::BodyId body;
        ShapeKind kind = ShapeKind::Box;
        hf::Vec3 half{0.5, 0.5, 0.5}; // box
        double radius = 0.5;          // sphere
        hf::Vec3 offset{};            // from the body's centre of mass, in its frame
        double friction = 0.6;
        double restitution = 0.0;
        double rolling = 0.0, spinning = 0.0; // lever arms, m
    };

    struct Pose {
        hf::Vec3 c;
        hf::Mat3 R;
    };

    // Appends the contacts between two shapes, normal from a to b, within margin of touching. featureBase keeps
    // the features of different pairs apart.
    void collide(const Shape &a, const Pose &pa, const Shape &b, const Pose &pb, double margin,
                 std::uint64_t featureBase, std::vector<hf::Contact> &out);

    double boundingRadius(const Shape &s);

} // namespace host
