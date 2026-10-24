#include <cmath>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"
#include "host/Host.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};

    WorldSettings weightless() {
        WorldSettings s;
        s.gravity = Vec3{};
        return s;
    }

    BodyId addCube(World &w, const Vec3 &at, double mass) {
        BodyDesc d;
        d.position = at;
        d.mass = mass;
        d.inertia = boxInertia(mass, Vec3{0.5, 0.5, 0.5});
        return w.addBody(d);
    }

    // How far a joint's two frames have come apart.
    double separation(World &w, JointId j) {
        const JointDesc d = w.jointDesc(j);
        const auto at = [&](BodyId b, const Frame &f) {
            return b.valid() ? w.position(b) + w.orientation(b) * f.position : f.position;
        };
        return length(at(d.a, d.frameA) - at(d.b, d.frameB));
    }

} // namespace

TEST_CASE("a free box grabbed and lifted rises to the hand") {
    World w;
    const BodyId box = addCube(w, Vec3{0, 1, 0}, 1.0);
    GrabDesc d;
    d.maxForce = 100.0;
    const GrabId g = w.grab(box, Vec3{}, d);
    REQUIRE(g.valid());
    REQUIRE(w.setGrabTarget(g, Vec3{0, 3, 0}));
    run(w, 3.0);
    // A 5 Hz spring on a kilogram sags g / ω² under its own weight: a centimetre.
    const double sag = 9.81 / std::pow(2.0 * std::numbers::pi * 5.0, 2);
    CHECK_MESSAGE(std::fabs(w.position(box).y - (3.0 - sag)) < 2e-3,
                  "it hangs at the hand, sagging as a 5 Hz spring does");
    CHECK_MESSAGE(!w.asleep(box), "held, it never sleeps");
    w.release(g);
    run(w, 0.5);
    CHECK_MESSAGE(w.position(box).y < 2.0, "let go, it falls");
}

TEST_CASE("a grab of ten newtons can't lift ten kilograms") {
    host::Host h;
    const BodyId box = h.addBox(Vec3{0, 0.5, 0}, Quat{}, Vec3{0.5, 0.5, 0.5}, 10.0);
    for (int i = 0; i < 30; ++i) h.step(kDt);
    const double y = h.world.position(box).y;
    GrabDesc d;
    d.maxForce = 10.0;
    const GrabId g = h.world.grab(box, Vec3{0, 0.5, 0}, d);
    h.world.setGrabTarget(g, Vec3{0, 3, 0});
    for (int i = 0; i < 120; ++i) h.step(kDt);
    CHECK_MESSAGE(std::fabs(h.world.position(box).y - y) < 2e-3, "a tenth of its weight leaves it on the ground");
}

TEST_CASE("a grab pulls no harder diagonally than straight") {
    // Three rows capped one by one would pull √3 times the cap along a diagonal.
    World w(weightless());
    const BodyId box = addCube(w, Vec3{}, 1.0);
    GrabDesc d;
    d.maxForce = 10.0;
    const GrabId g = w.grab(box, Vec3{}, d);
    w.setGrabTarget(g, Vec3{100, 100, 100});
    run(w, 1.0);
    CHECK_MESSAGE(length(w.velocity(box)) < 10.0 * 1.0 + 1e-6,
                  "ten newtons on a kilogram for a second: 10 m/s at most");
    CHECK_MESSAGE(length(w.velocity(box)) > 9.9, "and all of it, pulling the whole second");
}

TEST_CASE("an arm's tip grabbed and pulled out of reach straightens the arm toward the hand") {
    World w(weightless());
    // Its joints with a little friction, motors held still and capped low, or it would swing about straight forever.
    std::vector<BodyId> rods;
    std::vector<JointId> hinges;
    for (int i = 0; i < 3; ++i) {
        const Vec3 a{0, double(i), 0};
        rods.push_back(rodBetween(w, a, a + Vec3{0, 1, 0}, 10.0));
        JointDesc j = jointAt(w, JointType::Motor, i ? rods[std::size_t(i) - 1] : BodyId{}, rods.back(), a, kZ);
        j.maxForce = 5.0;
        hinges.push_back(w.addJoint(j));
    }
    GrabDesc d;
    d.maxForce = 200.0;
    // The tip of the last rod: half a metre along its own axis from its centre.
    const GrabId g = w.grab(rods.back(), Vec3{0.5, 0, 0}, d);
    w.setGrabTarget(g, Vec3{10, 0, 0});
    double worst = 0.0;
    for (int i = 0; i < 300; ++i) {
        w.step(kDt);
        for (JointId j : hinges) worst = std::max(worst, separation(w, j));
    }
    const Vec3 tip = w.position(rods.back()) + w.orientation(rods.back()) * Vec3{0.5, 0, 0};
    CHECK_MESSAGE(dot(normalize(tip), Vec3{1, 0, 0}) > std::cos(5.0 * std::numbers::pi / 180.0),
                  "the arm points at the hand");
    CHECK_MESSAGE(length(tip) > 2.95, "stretched out straight");
    CHECK_MESSAGE(worst < 1e-3, "every joint held within a millimetre all the way");
}

TEST_CASE("a four-bar's coupler grabbed and pulled: the loop stays closed") {
    World w(weightless());
    const Vec3 p0{0, 0, 0}, p1{0, 1, 0}, p2{2, 1, 0}, p3{2, 0, 0};
    const BodyId crank = rodBetween(w, p0, p1, 5.0), coupler = rodBetween(w, p1, p2, 5.0),
                 rocker = rodBetween(w, p3, p2, 5.0);
    std::vector<JointId> hinges{w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, crank, p0, kZ)),
                                w.addJoint(jointAt(w, JointType::Hinge, crank, coupler, p1, kZ)),
                                w.addJoint(jointAt(w, JointType::Hinge, coupler, rocker, p2, kZ)),
                                w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, rocker, p3, kZ))};
    GrabDesc d;
    d.maxForce = 100.0;
    const GrabId g = w.grab(coupler, Vec3{}, d);
    w.setGrabTarget(g, Vec3{3, 3, 0});
    double worst = 0.0;
    for (int i = 0; i < 180; ++i) {
        w.step(kDt);
        for (JointId j : hinges) worst = std::max(worst, separation(w, j));
    }
    CHECK_MESSAGE(w.position(coupler).x > 1.1, "the coupler moved toward the hand");
    CHECK_MESSAGE(worst < 1e-3, "its loop closed within a millimetre all along");
}

TEST_CASE("a grab turns a box to an orientation, within its torque cap") {
    World w(weightless());
    const BodyId box = addCube(w, Vec3{}, 1.0);
    GrabDesc d;
    d.maxForce = 50.0;
    d.maxTorque = 5.0;
    const GrabId g = w.grab(box, Vec3{}, d);
    const Quat turned = angleAxis(std::numbers::pi / 2, kZ);
    w.setGrabTarget(g, Vec3{}, &turned);
    run(w, 3.0);
    CHECK_MESSAGE(std::fabs(dot(w.orientation(box), turned)) > std::cos(0.5 * 2.0 * std::numbers::pi / 180.0),
                  "turned a quarter turn, to within two degrees");
    CHECK_MESSAGE(length(w.position(box)) < 1e-3, "and held where it was");
}

TEST_CASE("a grab refuses what it can't use, and goes with its body") {
    World w;
    const BodyId box = addCube(w, Vec3{}, 1.0);
    GrabDesc bad;
    bad.maxForce = -1.0;
    CHECK(!w.grab(box, Vec3{}, bad).valid());
    GrabDesc good;
    good.maxForce = 10.0;
    const GrabId g = w.grab(box, Vec3{}, good);
    CHECK(!w.setGrabTarget(g, Vec3{std::nan(""), 0, 0}));
    w.removeBody(box);
    CHECK_MESSAGE(!w.setGrabTarget(g, Vec3{}), "released with its body");
}
