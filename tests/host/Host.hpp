#pragma once

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "hf/World.hpp"
#include "host/Collision.hpp"

namespace host {

    // What a game wraps around the library: bodies with shapes, joints authored at world points, and collision
    // detection feeding the world its contacts every step.
    class Host {
    public:
        hf::World world;
        std::vector<Shape> shapes;
        std::vector<hf::Contact> contacts; // the last detection's
        double margin = 0.05;              // how far apart a pair is still reported

        explicit Host(const hf::WorldSettings &settings = {});

        hf::BodyId ground() const { return ground_; }
        // A mass of 0 makes a static body.
        hf::BodyId addBox(const hf::Vec3 &position, const hf::Quat &orientation, const hf::Vec3 &half, double mass,
                          double friction = 0.6);
        hf::BodyId addSphere(const hf::Vec3 &position, double radius, double mass, double restitution = 0.0);
        void removeBody(hf::BodyId body);

        // A joint at a world point, about or along a world axis, with its bodies where they are now.
        hf::JointDesc jointAt(hf::JointType type, hf::BodyId a, hf::BodyId b, const hf::Vec3 &point,
                              const hf::Vec3 &axis) const;
        // Jointed bodies don't collide unless asked.
        hf::JointId addJoint(const hf::JointDesc &desc, bool collide = false);
        void ignore(hf::BodyId a, hf::BodyId b);

        void detect();
        void step(double dt);

    private:
        hf::BodyId ground_;
        std::unordered_set<std::uint64_t> ignored_;
        struct Span {
            double lo, hi;
            int shape;
        };
        std::vector<Span> spans_;
        std::vector<Pose> poses_;
        std::vector<int> planes_;
    };

    // The rotation taking +X onto axis.
    hf::Quat xTo(const hf::Vec3 &axis);
    hf::Mat3 boxInertia(double mass, const hf::Vec3 &half);

} // namespace host
