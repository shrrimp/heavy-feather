#include <cmath>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};

    // A 2 m, 50 kg arm out along x from a joint at the origin: gravity turns it down with 490 N·m.
    JointId armOn(World &w, JointType type, double friction, const Vec3 &axis = kZ) {
        const BodyId arm = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 50.0);
        JointDesc d = jointAt(w, type, BodyId{}, arm, Vec3{0, 0, 0}, axis);
        d.friction = friction;
        return w.addJoint(d);
    }

    // A 2 m, 20 kg bar turning about its middle on a joint at the origin: no weight turns it.
    BodyId barOn(World &w, JointDesc &d, JointType type, double friction) {
        const BodyId bar = rodBetween(w, Vec3{-1, 0, 0}, Vec3{1, 0, 0}, 20.0);
        d = jointAt(w, type, BodyId{}, bar, Vec3{0, 0, 0}, kZ);
        d.friction = friction;
        return bar;
    }

    // A 10 kg block on a slider standing up: 98 N pull it down.
    JointId blockOnSlider(World &w, double friction) {
        const BodyId block = rodBetween(w, Vec3{0, -0.2, 0}, Vec3{0, 0.2, 0}, 10.0);
        JointDesc d = jointAt(w, JointType::Slider, BodyId{}, block, Vec3{0, 0, 0}, Vec3{0, 1, 0});
        d.min = -10.0;
        d.max = 10.0;
        d.friction = friction;
        return w.addJoint(d);
    }

} // namespace

TEST_CASE("a hinge holds what its friction can bear and gives to what it can't") {
    World holds, gives;
    const JointId a = armOn(holds, JointType::Hinge, 600.0);
    const JointId b = armOn(gives, JointType::Hinge, 300.0);
    run(holds, 2.0);
    run(gives, 1.0);
    CHECK_MESSAGE(std::fabs(holds.jointPosition(a)) < 1e-3, "600 N·m holds 490 N·m still");
    CHECK_MESSAGE(gives.jointPosition(b) < -0.3, "300 N·m lets it fall");

    World free;
    const JointId c = armOn(free, JointType::Hinge, 0.0);
    run(free, 0.3);
    World slowed;
    const JointId e = armOn(slowed, JointType::Hinge, 300.0);
    run(slowed, 0.3);
    CHECK_MESSAGE(slowed.jointPosition(e) > free.jointPosition(c) + 0.05, "and slower than it falls free");
}

TEST_CASE("a ball holds a swing its friction can bear on each axis") {
    // Its frame's X square to the swing, then along the arm: the swing comes about one of its axes, then another.
    for (const Vec3 axis : {kZ, Vec3{1, 0, 0}}) {
        World holds, gives;
        const JointId a = armOn(holds, JointType::Ball, 600.0, axis);
        const JointId b = armOn(gives, JointType::Ball, 300.0, axis);
        const BodyId armA = holds.jointDesc(a).b, armB = gives.jointDesc(b).b;
        run(holds, 2.0);
        run(gives, 1.0);
        CHECK_MESSAGE(holds.position(armA).y > -1e-3, "600 N·m holds the arm out");
        CHECK_MESSAGE(gives.position(armB).y < -0.3, "300 N·m lets it swing down");
    }
}

TEST_CASE("a slider holds a load its friction can bear and slides under one it can't") {
    World holds, gives;
    const JointId a = blockOnSlider(holds, 150.0);
    const JointId b = blockOnSlider(gives, 50.0);
    run(holds, 1.0);
    run(gives, 0.5);
    CHECK_MESSAGE(std::fabs(holds.jointPosition(a)) < 1e-3, "150 N holds 98 N");
    // 98 N down, 50 N against: 4.8 m/s^2, 0.6 m in half a second.
    CHECK(gives.jointPosition(b) == doctest::Approx(-0.6).epsilon(0.05));
}

TEST_CASE("a bar spun on a bearing slows to a stop by its friction") {
    WorldSettings settings;
    settings.gravity = Vec3{0, 0, 0};
    World w(settings);
    JointDesc d;
    const BodyId bar = barOn(w, d, JointType::Hinge, 10.0);
    const JointId j = w.addJoint(d);
    w.setVelocity(bar, Vec3{0, 0, 0}, Vec3{0, 0, 3});
    // I = 20 * 2^2 / 12 = 6.67 kg m^2: 10 N·m stops 3 rad/s in 2 s.
    run(w, 1.0);
    CHECK(w.jointVelocity(j) == doctest::Approx(1.5).epsilon(0.05));
    run(w, 1.5);
    CHECK(std::fabs(w.jointVelocity(j)) < 1e-3);
}

TEST_CASE("a motor drives through friction it can beat and stalls against friction it can't") {
    WorldSettings settings;
    settings.gravity = Vec3{0, 0, 0};
    World runs(settings), stalls(settings);
    JointDesc d;
    barOn(runs, d, JointType::Motor, 50.0);
    d.target = 2.0;
    d.maxForce = 100.0;
    const JointId a = runs.addJoint(d);
    barOn(stalls, d, JointType::Motor, 200.0);
    d.target = 2.0;
    d.maxForce = 100.0;
    const JointId b = stalls.addJoint(d);
    run(runs, 2.0);
    run(stalls, 2.0);
    CHECK_MESSAGE(runs.jointVelocity(a) == doctest::Approx(2.0).epsilon(0.01), "100 N·m beats 50 to its speed");
    CHECK_MESSAGE(std::fabs(stalls.jointVelocity(b)) < 1e-3, "100 N·m against 200 doesn't turn it");
    CHECK(runs.stats().lcpFallbacks == 0);
    CHECK(stalls.stats().lcpFallbacks == 0);
}

TEST_CASE("a joint with negative friction is refused") {
    World w;
    const BodyId arm = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 50.0);
    JointDesc d = jointAt(w, JointType::Hinge, BodyId{}, arm, Vec3{0, 0, 0}, kZ);
    d.friction = -1.0;
    CHECK_FALSE(w.addJoint(d).valid());
}
