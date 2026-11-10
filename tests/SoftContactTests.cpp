#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kHalf{0.5, 0.5, 0.5};
    constexpr double kMass = 100.0;
    constexpr double kG = 9.81;

    BodyId addGround(World &w) {
        BodyDesc d;
        d.isStatic = true;
        return w.addBody(d);
    }

    BodyId addBox(World &w, const Vec3 &position, double mass = kMass) {
        BodyDesc d;
        d.position = position;
        d.mass = mass;
        d.inertia = boxInertia(mass, kHalf);
        return w.addBody(d);
    }

    struct Bed {
        double stiffness = 50000.0;
        double damping = 0.3;
        double sinkLimit = 0.05;
        double friction = 0.6;
        double restitution = 0.0;
    };

    // The box's bottom corners on a soft plane at y = 0, handed from a reach above it so a fall lands on springs.
    void onBed(World &w, BodyId ground, BodyId box, const Bed &bed, bool fresh = true) {
        if (fresh) w.clearContacts();
        const Vec3 c = w.position(box);
        const Quat q = w.orientation(box);
        std::uint64_t feature = 0;
        for (double sx : {-1.0, 1.0})
            for (double sz : {-1.0, 1.0}) {
                const Vec3 p = c + q * Vec3{sx * kHalf.x, -kHalf.y, sz * kHalf.z};
                ++feature;
                if (p.y > 0.1) continue;
                Contact k;
                k.a = ground;
                k.b = box;
                k.point = p;
                k.normal = Vec3{0, 1, 0};
                k.depth = -p.y;
                k.feature = feature;
                k.friction = bed.friction;
                k.restitution = bed.restitution;
                k.stiffness = bed.stiffness;
                k.damping = bed.damping;
                k.sinkLimit = bed.sinkLimit;
                REQUIRE(w.addContact(k));
            }
    }

    double sinkOf(const World &w, BodyId box) { return kHalf.y - w.position(box).y; }

    bool within(double value, double expected, double fraction) {
        return std::fabs(value - expected) <= fraction * std::fabs(expected);
    }

    void settle(World &w, BodyId ground, BodyId box, const Bed &bed, double seconds) {
        for (int s = 0; s < int(std::lround(seconds / kDt)); ++s) {
            onBed(w, ground, box, bed);
            w.step(kDt);
        }
    }

    Vec3 pushOn(const World &w, BodyId body) {
        Vec3 sum{};
        for (const ContactForce &f : w.contactForces()) {
            if (f.b == body) sum += f.force;
            if (f.a == body) sum -= f.force;
        }
        return sum;
    }

} // namespace

TEST_CASE("a block on a soft bed sinks to its weight over the bed's stiffness") {
    // Awake: a sleeper is frozen where its swing slowed, a millimetre or so short of rest.
    WorldSettings settings;
    settings.sleeping = false;
    World w(settings);
    const BodyId ground = addGround(w);
    const BodyId box = addBox(w, Vec3{0, kHalf.y, 0});
    Bed bed;
    bed.stiffness = 5000.0;
    bed.sinkLimit = 0.3;
    settle(w, ground, box, bed, 3.0);
    const double expected = kMass * kG / (4.0 * bed.stiffness);
    CHECK(within(sinkOf(w, box), expected, 0.01));
    CHECK(within(pushOn(w, box).y, kMass * kG, 0.01));
}

TEST_CASE("a hinged pair on a soft bed sinks the same, solved with its machine") {
    World w;
    const BodyId ground = addGround(w);
    const BodyId a = addBox(w, Vec3{0, kHalf.y, 0});
    const BodyId b = addBox(w, Vec3{1.0, kHalf.y, 0});
    w.addJoint(jointAt(w, JointType::Hinge, a, b, Vec3{0.5, kHalf.y, 0}, Vec3{0, 0, 1}));
    const Bed bed;
    for (int s = 0; s < 180; ++s) {
        onBed(w, ground, a, bed);
        onBed(w, ground, b, bed, false);
        w.step(kDt);
    }
    const double expected = kMass * kG / (4.0 * bed.stiffness);
    CHECK(within(sinkOf(w, a), expected, 0.02));
    CHECK(within(sinkOf(w, b), expected, 0.02));
}

TEST_CASE("a block let go on a soft bed rings down at the bed's damping ratio") {
    World w;
    const BodyId ground = addGround(w);
    const BodyId box = addBox(w, Vec3{0, kHalf.y, 0});
    Bed bed;
    bed.stiffness = 5000.0;
    bed.damping = 0.1;
    bed.sinkLimit = 0.3;
    const double rest = kMass * kG / (4.0 * bed.stiffness);

    std::vector<double> peaks;
    double last = 0.0, before = 0.0;
    for (int s = 0; s < 240; ++s) {
        onBed(w, ground, box, bed);
        w.step(kDt);
        const double x = sinkOf(w, box) - rest;
        if (s >= 2 && last > before && last > x && last > 0.0) peaks.push_back(last);
        before = last;
        last = x;
    }
    REQUIRE(peaks.size() >= 2);
    const double decrement = std::log(peaks[0] / peaks[1]);
    const double zeta = decrement / std::sqrt(4.0 * std::numbers::pi * std::numbers::pi + decrement * decrement);
    CHECK(within(zeta, bed.damping, 0.1));
}

TEST_CASE("a block dropped on a soft bed bounces off its springs, with no restitution") {
    const auto rebound = [](double stiffness) {
        World w;
        const BodyId ground = addGround(w);
        const BodyId box = addBox(w, Vec3{0, kHalf.y + 0.3, 0});
        Bed bed;
        bed.stiffness = stiffness;
        bed.damping = 0.0;
        bed.sinkLimit = 0.2;
        double top = 0.0;
        bool landed = false;
        for (int s = 0; s < 90; ++s) {
            onBed(w, ground, box, bed);
            w.step(kDt);
            const double y = w.position(box).y - kHalf.y;
            if (y < 0.0) landed = true;
            if (landed) top = std::max(top, y);
        }
        return top;
    };
    CHECK_MESSAGE(rebound(50000.0) > 0.2, "a spring gives back what it took");
    CHECK_MESSAGE(rebound(0.0) < 0.01, "a rigid contact without restitution keeps it");
}

TEST_CASE("a soft bed ten times overloaded bottoms out at its sink limit") {
    World w;
    const BodyId ground = addGround(w);
    const BodyId box = addBox(w, Vec3{0, kHalf.y, 0});
    Bed bed;
    bed.stiffness = kMass * kG / (4.0 * 10.0 * bed.sinkLimit);
    settle(w, ground, box, bed, 3.0);
    CHECK(sinkOf(w, box) < bed.sinkLimit + 0.004);
    CHECK(sinkOf(w, box) > bed.sinkLimit - 0.004);
    CHECK(std::fabs(w.velocity(box).y) < 1e-3);
}

TEST_CASE("a block on a soft bed holds within its friction and slides past it") {
    const auto slide = [](double pushOverGrip) {
        World w;
        const BodyId ground = addGround(w);
        const BodyId box = addBox(w, Vec3{0, kHalf.y, 0});
        const Bed bed;
        settle(w, ground, box, bed, 1.0);
        const double x0 = w.position(box).x;
        for (int s = 0; s < 60; ++s) {
            onBed(w, ground, box, bed);
            w.applyForce(box, Vec3{pushOverGrip * bed.friction * kMass * kG, 0, 0}, w.position(box));
            w.step(kDt);
        }
        return w.position(box).x - x0;
    };
    CHECK(std::fabs(slide(0.5)) < 0.002);
    CHECK(slide(1.5) > 0.5);
}
