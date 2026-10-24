#include <cmath>
#include <numbers>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};
    constexpr double kPi = std::numbers::pi;

    // How far a joint's two frames, each placed from its body's origin, have come apart.
    double separation(World &w, JointId j) {
        const JointDesc d = w.jointDesc(j);
        const auto at = [&](BodyId b, const Frame &f) {
            return b.valid() ? w.position(b) + w.orientation(b) * f.position : f.position;
        };
        return length(at(d.a, d.frameA) - at(d.b, d.frameB));
    }

    WorldSettings weightless() {
        WorldSettings s;
        s.gravity = Vec3{};
        return s;
    }

} // namespace

TEST_CASE("a body's pose is its origin's, its centre of mass where it says from there") {
    World w;
    BodyDesc d;
    d.position = Vec3{1, 2, 3};
    d.orientation = angleAxis(kPi / 2, kZ);
    d.centerOfMass = Vec3{0.5, 0, 0};
    const BodyId b = w.addBody(d);
    CHECK_MESSAGE(length(w.position(b) - Vec3{1, 2, 3}) < 1e-12, "the origin, as given");
    CHECK_MESSAGE(length(w.centerOfMass(b) - Vec3{1, 2.5, 3}) < 1e-12,
                  "the centre half a metre along its x, turned to y");
    CHECK_MESSAGE(length(w.state(b).position - w.position(b)) < 1e-12, "state reads the origin too");
}

TEST_CASE("a body spinning about its centre of mass swings its origin round it") {
    World w(weightless());
    BodyDesc d;
    d.centerOfMass = Vec3{1, 0, 0};
    d.angularVelocity = Vec3{0, 0, kPi / 2}; // a quarter turn a second
    const BodyId b = w.addBody(d);
    run(w, 1.0);
    CHECK_MESSAGE(length(w.centerOfMass(b) - Vec3{1, 0, 0}) < 1e-9, "the centre stays put");
    // From one metre short of the centre along -x, a quarter turn about z puts it one metre along -y.
    CHECK_MESSAGE(length(w.position(b) - Vec3{1, -1, 0}) < 1e-3, "the origin goes round it");
}

TEST_CASE("a pendulum hinged at its origin swings with the period its centre of mass gives it") {
    World w;
    const double m = 2.0, d = 1.0, I = 0.02;
    BodyDesc b;
    b.position = Vec3{0, 3, 0};
    b.mass = m;
    b.inertia = Mat3(I);
    b.centerOfMass = Vec3{d * std::sin(0.1), -d * std::cos(0.1), 0}; // let go a tenth of a radian off
    const BodyId bob = w.addBody(b);
    JointDesc j;
    j.type = JointType::Hinge;
    j.b = bob;
    j.frameA = Frame{Vec3{0, 3, 0}, xTo(kZ)};
    j.frameB = Frame{Vec3{}, xTo(kZ)}; // at its origin
    w.addJoint(j);
    // T = 2π √((I + m d²) / (m g d)) for small swings.
    const double period = 2.0 * kPi * std::sqrt((I + m * d * d) / (m * 9.81 * d));
    double last = -1.0, previous = w.centerOfMass(bob).x;
    int crossings = 0;
    double first = 0.0;
    for (int i = 1; i <= 600; ++i) {
        w.step(kDt);
        const double x = w.centerOfMass(bob).x;
        if (previous > 0.0 && x <= 0.0) {
            const double t = (i - 1 + previous / (previous - x)) * kDt;
            if (crossings == 0) first = t;
            last = t;
            ++crossings;
        }
        previous = x;
    }
    REQUIRE(crossings >= 2);
    CHECK_MESSAGE(std::fabs((last - first) / (crossings - 1) - period) < 0.01 * period,
                  "the period of its centre, not of its origin");
    CHECK_MESSAGE(length(w.position(bob) - Vec3{0, 3, 0}) < 1e-6, "its origin on the pin all the while");
}

TEST_CASE("a machine link's centre of mass moved mid-drive: its joint, its setpoint and its machine carry on") {
    World w;
    const BodyId arm = rodBetween(w, Vec3{0, 2, 0}, Vec3{2, 2, 0}, 50.0);
    JointDesc s = jointAt(w, JointType::Servo, BodyId{}, arm, Vec3{0, 2, 0}, kZ);
    s.maxForce = 1e5;
    s.speed = 0.5;
    s.min = -1.5;
    s.max = 1.5;
    s.target = 1.0;
    const JointId servo = w.addJoint(s);
    run(w, 0.5);
    const double at = w.jointPosition(servo), target = w.target(servo);
    REQUIRE(at > 0.1);
    REQUIRE(at < 0.9);
    const Vec3 origin = w.position(arm);
    // Five centimetres sideways, along the arm's own y: a voxel edit.
    w.setMassProperties(arm, 50.0, boxInertia(50.0, Vec3{1.0, 0.1, 0.1}), Vec3{0, 0.05, 0});
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo) - at) < 1e-12, "the joint is where it was");
    CHECK_MESSAGE(w.target(servo) == target, "its drive's target too");
    CHECK_MESSAGE(length(w.position(arm) - origin) < 1e-12, "and the arm's origin");
    w.step(kDt);
    CHECK_MESSAGE((w.stats().groupRebuilds == 0 && w.stats().machineRebuilds == 0), "nothing rebuilt");
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo) - at) < 0.5 * kDt + 1e-3, "the drive carried on from where it was");
    CHECK_MESSAGE(separation(w, servo) < 1e-6, "the arm still on its pin, not snapped to where its old centre put it");
    run(w, 3.0);
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo) - 1.0) < 1e-3, "and got there");
    CHECK_MESSAGE(separation(w, servo) < 1e-6, "on its pin all the way");
}

TEST_CASE("a grab holds a point measured from the body's origin") {
    World w(weightless());
    BodyDesc d;
    d.centerOfMass = Vec3{0, -1, 0};
    const BodyId b = w.addBody(d);
    GrabDesc g;
    g.maxForce = 100.0;
    const GrabId hand = w.grab(b, Vec3{}, g); // at the origin, a metre from the centre
    w.setGrabTarget(hand, Vec3{2, 0, 0});
    run(w, 3.0);
    CHECK_MESSAGE(length(w.position(b) - Vec3{2, 0, 0}) < 1e-2, "the origin, not the centre, went to the hand");
}
