#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"
#include "world/TestHooks.hpp"
#include "world/WorldImpl.hpp"

using namespace build;

namespace {

    constexpr double kPi = std::numbers::pi;
    constexpr Vec3 kZ{0, 0, 1};
    constexpr double kPitch = 0.05;

    double pitchRadius(int teeth) { return teeth * kPitch / (2.0 * kPi); }

    BodyId wheel(World &w, const Vec3 &at, double radius, double mass = 5.0) {
        BodyDesc b;
        b.position = at;
        b.mass = mass;
        const double flat = 0.5 * mass * radius * radius, across = 0.25 * mass * radius * radius;
        b.inertia = Mat3(across, 0, 0, 0, across, 0, 0, 0, flat);
        return w.addBody(b);
    }

    // A wheel turning about Z on the world, on a joint of this type at its centre.
    struct Guide {
        BodyId body;
        JointId joint;
        double radius = 0.0;
    };

    Guide guideAt(World &w, JointType type, const Vec3 &at, double radius) {
        Guide g;
        g.radius = radius;
        g.body = wheel(w, at, radius);
        g.joint = w.addJoint(jointAt(w, type, BodyId{}, g.body, at, kZ));
        return g;
    }

    void drive(World &w, JointId joint, double speed, double maxTorque) {
        JointDesc d = w.jointDesc(joint);
        d.target = speed;
        d.maxForce = maxTorque;
        w.setJoint(joint, d);
    }

    BeltGuide meshing(const Guide &g) {
        BeltGuide b;
        b.hinge = g.joint;
        b.radius = g.radius;
        b.meshes = true;
        return b;
    }

    // A chain round a 12-tooth sprocket and a 30-tooth one a metre apart, both on the world.
    struct ChainDrive {
        World w;
        Guide small, large;
        BeltId chain;

        ChainDrive(JointType smallType, JointType largeType) {
            small = guideAt(w, smallType, Vec3{0, 0, 0}, pitchRadius(12));
            large = guideAt(w, largeType, Vec3{1, 0, 0}, pitchRadius(30));
            BeltDesc d;
            d.length = 3.0;
            d.mass = 2.0;
            d.pitch = kPitch;
            const BeltGuide loop[] = {meshing(small), meshing(large)};
            chain = w.addBelt(d, loop);
        }
    };

    // A joint's angle, followed through its turns.
    struct Turned {
        double total = 0.0, last = 0.0;
        void follow(double angle) {
            total += std::remainder(angle - last, 2.0 * kPi);
            last = angle;
        }
    };

} // namespace

TEST_CASE("a chain carries a motor's turn to a larger sprocket at the ratio of their teeth") {
    ChainDrive c(JointType::Motor, JointType::Hinge);
    REQUIRE(c.w.contains(c.chain));
    drive(c.w, c.small.joint, 2.0, 1000.0);
    run(c.w, 1.0);
    const double small = c.w.jointVelocity(c.small.joint), large = c.w.jointVelocity(c.large.joint);
    CHECK(std::fabs(small - 2.0) < 1e-6);
    CHECK_MESSAGE(std::fabs(large / small - 12.0 / 30.0) < 1e-6, "30 teeth turn 12 for every 30 of the small");
    CHECK_MESSAGE(std::fabs(c.w.beltSpeed(c.chain) - c.small.radius * small) < 1e-6, "the chain runs at its rim");
    CHECK(c.w.beltInfo(c.chain).machine >= 0);
}

TEST_CASE("a chain carries the torque back: a braked large sprocket loads the small one by the same ratio") {
    ChainDrive c(JointType::Motor, JointType::Motor);
    drive(c.w, c.small.joint, 2.0, 100.0);
    drive(c.w, c.large.joint, 0.0, 10.0); // a brake that slips at 10 N·m
    run(c.w, 2.0);
    CHECK_MESSAGE(std::fabs(c.w.jointVelocity(c.large.joint) - 0.8) < 1e-4, "dragged round by the chain");
    CHECK_MESSAGE(std::fabs(std::fabs(c.w.driveForce(c.large.joint)) - 10.0) < 1e-6, "the brake slipping at its cap");
    CHECK_MESSAGE(std::fabs(std::fabs(c.w.driveForce(c.small.joint)) - 4.0) < 1e-3,
                  "12/30 of it at the small sprocket's motor");
}

TEST_CASE("a sprocket keeps its step on the chain through thousands of steps of drive and reversal") {
    ChainDrive c(JointType::Motor, JointType::Hinge);
    Turned small, large;
    small.last = c.w.jointPosition(c.small.joint);
    large.last = c.w.jointPosition(c.large.joint);
    const auto slip = [&](const Turned &t, double radius) { return c.w.beltPosition(c.chain) - radius * t.total; };
    c.w.step(kDt);
    small.follow(c.w.jointPosition(c.small.joint));
    large.follow(c.w.jointPosition(c.large.joint));
    const double small0 = slip(small, c.small.radius), large0 = slip(large, c.large.radius);
    double worst = 0.0;
    for (int s = 0; s < 10000; ++s) {
        if (s % 120 == 0) drive(c.w, c.small.joint, (s / 120) % 2 ? -6.0 : 6.0, 1000.0);
        c.w.step(kDt);
        small.follow(c.w.jointPosition(c.small.joint));
        large.follow(c.w.jointPosition(c.large.joint));
        worst = std::max({worst, std::fabs(slip(small, c.small.radius) - small0),
                          std::fabs(slip(large, c.large.radius) - large0)});
    }
    CHECK_MESSAGE(worst < 0.01 * kPitch, "never a hundredth of a link off its teeth");
}

TEST_CASE("a belt makes one machine of the guides it joins, and they part when it goes") {
    ChainDrive c(JointType::Hinge, JointType::Hinge);
    c.w.step(kDt);
    CHECK(c.w.jointInfo(c.small.joint).machine == c.w.jointInfo(c.large.joint).machine);
    CHECK(c.w.stats().belts == 1);
    c.w.removeBelt(c.chain);
    c.w.step(kDt);
    CHECK(!c.w.contains(c.chain));
    CHECK(c.w.jointInfo(c.small.joint).machine != c.w.jointInfo(c.large.joint).machine);
}

TEST_CASE("a belt goes with a guide's hinge, and with a guide removed") {
    ChainDrive a(JointType::Hinge, JointType::Hinge);
    a.w.removeJoint(a.large.joint);
    CHECK(!a.w.contains(a.chain));
    a.w.step(kDt);

    ChainDrive b(JointType::Hinge, JointType::Hinge);
    b.w.removeBody(b.small.body);
    CHECK(!b.w.contains(b.chain));
    b.w.step(kDt);
}

TEST_CASE("belts the world can't use are refused") {
    World w;
    const Guide a = guideAt(w, JointType::Hinge, Vec3{0, 0, 0}, 0.2), b = guideAt(w, JointType::Hinge, Vec3{1, 0, 0}, 0.2);
    const Guide tilted = [&] {
        Guide g;
        g.radius = 0.2;
        g.body = wheel(w, Vec3{2, 0, 0}, 0.2);
        g.joint = w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, g.body, Vec3{2, 0, 0}, Vec3{1, 0, 0}));
        return g;
    }();
    BeltDesc d;
    const BeltGuide pair[] = {meshing(a), meshing(b)};
    const BeltGuide one[] = {meshing(a)};
    const BeltGuide skew[] = {meshing(a), meshing(tilted)};
    CHECK(!w.addBelt(d, one).valid());
    CHECK(!w.addBelt(d, skew).valid());
    BeltDesc light = d;
    light.mass = 0.0;
    CHECK(!w.addBelt(light, pair).valid());
    BeltGuide flat = meshing(a);
    flat.radius = -0.1;
    const BeltGuide bad[] = {flat, meshing(b)};
    CHECK(!w.addBelt(d, bad).valid());
    CHECK(w.addBelt(d, pair).valid());
}

TEST_CASE("an idle chain sleeps with its sprockets, and a running one keeps them awake") {
    ChainDrive c(JointType::Motor, JointType::Hinge);
    drive(c.w, c.small.joint, 0.0, 1000.0);
    run(c.w, 2.0);
    CHECK_MESSAGE(c.w.asleep(c.small.body), "nothing moving: asleep");
    drive(c.w, c.small.joint, 1.0, 1000.0);
    run(c.w, 2.0);
    CHECK(!c.w.asleep(c.small.body));
    CHECK(c.w.beltSpeed(c.chain) > 0.0);
}

TEST_CASE("an island put back puts its belt back where the step found it") {
    ChainDrive c(JointType::Motor, JointType::Hinge);
    drive(c.w, c.small.joint, 3.0, 1000.0);
    run(c.w, 0.5);
    const double position = c.w.beltPosition(c.chain), speed = c.w.beltSpeed(c.chain);
    testing::phaseHook = [](World::Impl &I, int substep, FaultPhase phase) {
        if (phase != FaultPhase::Internal || substep != 2) return;
        for (MachineData &md : I.machines)
            for (const MachineData::BeltLink &bl : md.belts)
                md.m.art.qdot[std::size_t(md.m.art.links[std::size_t(bl.link)].dof0)] =
                    std::numeric_limits<double>::quiet_NaN();
    };
    c.w.step(kDt);
    testing::phaseHook = nullptr;
    CHECK(c.w.faults().size() == 1);
    CHECK(c.w.beltPosition(c.chain) == position);
    CHECK(c.w.beltSpeed(c.chain) == speed);
    run(c.w, 0.5);
    CHECK_MESSAGE(c.w.beltSpeed(c.chain) > 0.0, "and it runs on");
}

TEST_CASE("a chain on a frame that tumbles keeps its ratio on the frame, and the frame feels its drive") {
    WorldSettings s;
    s.gravity = Vec3{};
    World w(s);
    BodyDesc f;
    f.position = Vec3{0.5, 0, 0};
    f.mass = 20.0;
    f.inertia = boxInertia(20.0, Vec3{0.8, 0.2, 0.1});
    f.angularVelocity = Vec3{0.3, 0.2, 1.0};
    const BodyId frame = w.addBody(f);
    const auto onFrame = [&](JointType type, const Vec3 &at, int teeth) {
        Guide g;
        g.radius = pitchRadius(teeth);
        g.body = wheel(w, at, g.radius);
        g.joint = w.addJoint(jointAt(w, type, frame, g.body, at, kZ));
        return g;
    };
    const Guide small = onFrame(JointType::Motor, Vec3{0, 0, 0}, 12), large = onFrame(JointType::Hinge, Vec3{1, 0, 0}, 30);
    BeltDesc d;
    d.carrier = frame;
    d.length = 3.0;
    d.mass = 2.0;
    d.pitch = kPitch;
    const BeltGuide loop[] = {meshing(small), meshing(large)};
    const BeltId chain = w.addBelt(d, loop);
    REQUIRE(chain.valid());
    drive(w, small.joint, 4.0, 1000.0);
    run(w, 1.0);
    CHECK(std::fabs(w.jointVelocity(small.joint) - 4.0) < 1e-6);
    CHECK_MESSAGE(std::fabs(w.jointVelocity(large.joint) / w.jointVelocity(small.joint) - 0.4) < 1e-6,
                  "on the frame, whatever the frame does");
    CHECK(std::fabs(w.beltSpeed(chain) - small.radius * 4.0) < 1e-6);
}
