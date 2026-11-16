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

namespace {

    // A sprocket on the world, and under it an idler hanging on the belt: on a carriage that slides up and down, the
    // idler turning on it. The belt holds their weight on two strands.
    struct Hanging {
        World w;
        Guide sprocket;
        BodyId carriage, idler;
        JointId turn;
        BeltId belt;
        static constexpr double kR = 0.2, kMass = 5.0, kDrop = 1.0;

        Hanging(JointType idlerJoint, double friction, const BeltDesc &d) {
            sprocket = guideAt(w, JointType::Motor, Vec3{0, 0, 0}, kR);
            BodyDesc c;
            c.position = Vec3{0, -kDrop, 0};
            c.mass = kMass;
            c.inertia = boxInertia(kMass, Vec3{0.1, 0.1, 0.1});
            carriage = w.addBody(c);
            JointDesc slide = jointAt(w, JointType::Slider, BodyId{}, carriage, Vec3{0, -kDrop, 0}, Vec3{0, 1, 0});
            slide.min = -5.0;
            slide.max = 5.0;
            w.addJoint(slide);
            idler = wheel(w, Vec3{0, -kDrop, 0}, kR, kMass);
            turn = w.addJoint(jointAt(w, idlerJoint, carriage, idler, Vec3{0, -kDrop, 0}, kZ));
            BeltGuide free;
            free.hinge = turn;
            free.radius = kR;
            free.friction = friction;
            const BeltGuide loop[] = {meshing(sprocket), free};
            belt = w.addBelt(d, loop);
        }

        double weight() const { return 2.0 * kMass * 9.81; }
    };

    BeltDesc hangingBelt() {
        BeltDesc d;
        d.length = 2.0 * Hanging::kDrop + 2.0 * kPi * Hanging::kR;
        d.mass = 1.0;
        return d;
    }

} // namespace

TEST_CASE("a belt holds what hangs on it, half its weight in each strand") {
    Hanging k(JointType::Hinge, 0.5, hangingBelt());
    REQUIRE(k.w.contains(k.belt));
    run(k.w, 2.0);
    CHECK(k.w.beltInfo(k.belt).taut);
    CHECK(std::fabs(k.w.beltTension(k.belt) / (0.5 * k.weight()) - 1.0) < 1e-3);
    CHECK_MESSAGE(std::fabs(k.w.position(k.carriage).y + Hanging::kDrop) < 1e-4, "where the loop meets its length");
}

TEST_CASE("a guide that doesn't mesh turns with the belt, and slips once it is held harder than the belt grips") {
    // The grip: friction 0.5 × tension 49 N round half a turn, 77 N at the rim, 15.4 N·m.
    Hanging free(JointType::Hinge, 0.5, hangingBelt());
    drive(free.w, free.sprocket.joint, 2.0, 1000.0);
    run(free.w, 2.0);
    CHECK_MESSAGE(std::fabs(free.w.jointVelocity(free.turn) - 2.0) < 1e-4, "free, it turns as the sprocket does");
    CHECK(free.w.beltInfo(free.belt).slipping == 0);

    Hanging light(JointType::Motor, 0.5, hangingBelt());
    drive(light.w, light.turn, 0.0, 5.0);
    drive(light.w, light.sprocket.joint, 2.0, 1000.0);
    run(light.w, 2.0);
    CHECK_MESSAGE(std::fabs(light.w.jointVelocity(light.turn) - 2.0) < 1e-4, "a light brake is dragged round");
    CHECK(std::fabs(std::fabs(light.w.driveForce(light.turn)) - 5.0) < 1e-6);

    Hanging held(JointType::Motor, 0.5, hangingBelt());
    drive(held.w, held.turn, 0.0, 30.0);
    drive(held.w, held.sprocket.joint, 2.0, 1000.0);
    run(held.w, 2.0);
    CHECK_MESSAGE(std::fabs(held.w.jointVelocity(held.turn)) < 1e-3, "a hard brake holds it: the belt slides over it");
    CHECK(held.w.beltInfo(held.belt).slipping == 1);
    CHECK_MESSAGE(std::fabs(held.w.beltSpeed(held.belt) - 2.0 * Hanging::kR) < 1e-4, "while the sprocket runs it on");
}

TEST_CASE("an elastic belt stretches T L / EA under its load") {
    Hanging rigid(JointType::Hinge, 0.5, hangingBelt());
    BeltDesc d = hangingBelt();
    d.stiffness = 1e4;
    Hanging stretchy(JointType::Hinge, 0.5, d);
    run(rigid.w, 3.0);
    run(stretchy.w, 3.0);
    const double stretch = 0.5 * rigid.weight() * d.length / d.stiffness;
    const double drop = rigid.w.position(rigid.carriage).y - stretchy.w.position(stretchy.carriage).y;
    CHECK_MESSAGE(std::fabs(drop / (0.5 * stretch) - 1.0) < 0.02, "both strands give half of it each");
}

TEST_CASE("a belt pulled harder than it is strong breaks, lets go, and is mended by setBelt") {
    BeltDesc d = hangingBelt();
    d.strength = 30.0;
    Hanging k(JointType::Hinge, 0.5, d);
    int steps = 0;
    while (!k.w.beltInfo(k.belt).broken && steps++ < 120) k.w.step(kDt);
    REQUIRE(k.w.beltInfo(k.belt).broken);
    REQUIRE(k.w.brokenBelts().size() == 1);
    CHECK(k.w.brokenBelts()[0] == k.belt);
    k.w.step(kDt);
    CHECK_MESSAGE(k.w.brokenBelts().empty(), "listed the step it broke");
    run(k.w, 0.2);
    CHECK_MESSAGE(k.w.position(k.carriage).y < -Hanging::kDrop - 0.05, "what hung on it falls");
    d.strength = 1000.0;
    CHECK(k.w.setBelt(k.belt, d));
    CHECK(!k.w.beltInfo(k.belt).broken);
}

TEST_CASE("a chain spinning free on free sprockets keeps its energy") {
    WorldSettings s;
    s.gravity = Vec3{};
    World w(s);
    const Guide small = guideAt(w, JointType::Hinge, Vec3{0, 0, 0}, pitchRadius(12));
    const Guide large = guideAt(w, JointType::Hinge, Vec3{1, 0, 0}, pitchRadius(30));
    BeltDesc d;
    d.length = 3.1;
    d.mass = 2.0;
    d.pitch = kPitch;
    const BeltGuide loop[] = {meshing(small), meshing(large)};
    const BeltId chain = w.addBelt(d, loop);
    w.setVelocity(small.body, Vec3{}, Vec3{0, 0, 10});
    w.step(kDt);
    const auto energy = [&] {
        const auto spin = [&](const Guide &g) {
            const double I = 0.5 * 5.0 * g.radius * g.radius, v = w.jointVelocity(g.joint);
            return 0.5 * I * v * v;
        };
        const double v = w.beltSpeed(chain);
        return spin(small) + spin(large) + 0.5 * d.mass * v * v;
    };
    const double start = energy();
    REQUIRE(start > 0.0);
    run(w, 10.0);
    CHECK(std::fabs(energy() / start - 1.0) < 1e-6);
}

TEST_CASE("a frame yawing with a running belt on it is twisted by the belt's turning momentum, m v 2A/L") {
    WorldSettings s;
    s.gravity = Vec3{};
    World w(s);
    const double guideMass = 0.5, beltMass = 4.0, R = 0.1, d = 0.8, yaw = 2.0, spin = 40.0;
    const auto box = [&](double mass, const Vec3 &size) {
        BodyDesc b;
        b.mass = mass;
        b.inertia = boxInertia(mass, size);
        return w.addBody(b);
    };
    // A gimbal: yawed at a steady rate, and held level by a motor whose torque is what the test reads.
    const BodyId yoke = box(5.0, Vec3{0.3, 0.3, 0.3}), frame = box(10.0, Vec3{1.0, 0.3, 0.2});
    const JointId turn = w.addJoint(jointAt(w, JointType::Motor, BodyId{}, yoke, Vec3{}, Vec3{0, 1, 0}));
    const JointId level = w.addJoint(jointAt(w, JointType::Motor, yoke, frame, Vec3{}, Vec3{1, 0, 0}));
    const auto onFrame = [&](JointType type, const Vec3 &at) {
        Guide g;
        g.radius = R;
        g.body = wheel(w, at, R, guideMass);
        g.joint = w.addJoint(jointAt(w, type, frame, g.body, at, kZ));
        return g;
    };
    const Guide driven = onFrame(JointType::Motor, Vec3{-0.5 * d, 0, 0}), idle = onFrame(JointType::Hinge, Vec3{0.5 * d, 0, 0});
    BeltDesc b;
    b.carrier = frame;
    b.length = 2.0 * d + 2.0 * kPi * R;
    b.mass = beltMass;
    const BeltGuide loop[] = {meshing(driven), meshing(idle)};
    const BeltId belt = w.addBelt(b, loop);
    REQUIRE(belt.valid());
    drive(w, driven.joint, spin, 1e4);
    drive(w, turn, yaw, 1e4);
    drive(w, level, 0.0, 1e4);
    run(w, 1.0);
    REQUIRE(std::fabs(w.beltSpeed(belt) - spin * R) < 1e-6);

    const double reach = 2.0 * (2.0 * d * R + kPi * R * R) / b.length;
    const double wheels = 2.0 * 0.5 * guideMass * R * R * spin;
    const double twist = yaw * (beltMass * w.beltSpeed(belt) * reach + wheels);
    CHECK_MESSAGE(std::fabs(std::fabs(w.driveForce(level)) / twist - 1.0) < 0.01,
                  "held level against " << twist << " N·m, read " << w.driveForce(level));
}

TEST_CASE("a belt given another loop runs round its new guides, and lets the old ones go") {
    ChainDrive c(JointType::Motor, JointType::Hinge);
    drive(c.w, c.small.joint, 2.0, 1000.0);
    run(c.w, 0.5);
    const Guide other = guideAt(c.w, JointType::Hinge, Vec3{0.8, 0, 0}, pitchRadius(24));
    const BeltGuide loop[] = {meshing(c.small), meshing(other)};
    const BeltGuide skew[] = {meshing(c.small)};
    CHECK(!c.w.setBeltLoop(c.chain, skew));
    REQUIRE(c.w.setBeltLoop(c.chain, loop));
    run(c.w, 1.0);
    CHECK(std::fabs(c.w.jointVelocity(other.joint) - 1.0) < 1e-6);
    CHECK(c.w.jointInfo(c.small.joint).machine == c.w.jointInfo(other.joint).machine);
    CHECK(c.w.jointInfo(c.small.joint).machine != c.w.jointInfo(c.large.joint).machine);
}

namespace {

    BodyId slab(World &w, const Vec3 &at, const Vec3 &half, double mass, bool fixed = false) {
        BodyDesc b;
        b.position = at;
        b.mass = mass;
        b.inertia = boxInertia(mass, half);
        b.isStatic = fixed;
        return w.addBody(b);
    }

    // The straight from one guide to the next, the way the belt runs along it.
    Vec3 straight(const World &w, BeltId belt, int from) {
        RopeBend path[2];
        w.beltPath(belt, path);
        return normalize(path[(from + 1) % 2].in - path[from].out);
    }

    // A block on a track round a driven sprocket and an idler that doesn't mesh, on the ground: the host's contacts
    // under each guide, two across the track.
    struct Tracked {
        static constexpr double kR = 0.15, kFriction = 0.8;
        World w;
        BodyId ground, block;
        Guide sprocket, idler;
        BeltId track;

        Tracked() {
            ground = slab(w, Vec3{0, -0.5, 0}, Vec3{50, 0.5, 50}, 1.0, true);
            block = slab(w, Vec3{0, kR, 0}, Vec3{0.6, 0.1, 0.3}, 40.0);
            const auto onBlock = [&](JointType type, double x) {
                Guide g;
                g.radius = kR;
                g.body = wheel(w, Vec3{x, kR, 0}, kR, 2.0);
                g.joint = w.addJoint(jointAt(w, type, block, g.body, Vec3{x, kR, 0}, kZ));
                return g;
            };
            sprocket = onBlock(JointType::Motor, -0.5);
            idler = onBlock(JointType::Hinge, 0.5);
            BeltDesc d;
            d.carrier = block;
            d.length = 2.0 + 2.0 * kPi * kR;
            d.mass = 4.0;
            BeltGuide free;
            free.hinge = idler.joint;
            free.radius = kR;
            free.friction = kFriction;
            const BeltGuide loop[] = {meshing(sprocket), free};
            track = w.addBelt(d, loop);
        }

        double weight() const { return (40.0 + 2.0 * 2.0 + 4.0) * 9.81; }
        double mass() const { return 40.0 + 2.0 * 2.0 + 4.0; }

        void touch() {
            w.clearContacts();
            const Vec3 bottom = straight(w, track, 0);
            std::uint64_t feature = 0;
            for (const Guide &g : {sprocket, idler})
                for (double z : {-0.2, 0.2}) {
                    Contact k;
                    k.a = ground;
                    k.b = g.body;
                    k.point = w.position(g.body) + Vec3{0, -kR, z};
                    k.normal = Vec3{0, 1, 0};
                    k.depth = -k.point.y;
                    k.feature = feature++;
                    k.friction = kFriction;
                    k.belt = track;
                    k.beltDirection = bottom;
                    w.addContact(k);
                }
        }

        void drive(double seconds, double resist = 0.0) {
            for (int s = 0; s < int(std::lround(seconds / kDt)); ++s) {
                touch();
                if (resist != 0.0) w.applyForce(block, Vec3{resist, 0, 0}, w.position(block) + Vec3{0, -kR, 0});
                w.step(kDt);
            }
        }
    };

} // namespace

TEST_CASE("a track on the ground carries its block at the track's speed, and turns the idler under it") {
    Tracked t;
    REQUIRE(t.w.contains(t.track));
    REQUIRE(straight(t.w, t.track, 0).x > 0.99);
    t.drive(0.5);
    ::drive(t.w, t.sprocket.joint, 4.0, 1e4);
    t.drive(2.0);
    const double speed = t.w.beltSpeed(t.track);
    CHECK(std::fabs(speed - 4.0 * Tracked::kR) < 1e-6);
    CHECK_MESSAGE(std::fabs(t.w.velocity(t.block).x + speed) < 1e-3, "the bottom stands on the ground while the block goes");
    CHECK_MESSAGE(std::fabs(t.w.jointVelocity(t.idler.joint) - 4.0) < 1e-3, "the ground's push presses the idler on the belt");
    CHECK(std::fabs(t.w.position(t.block).y - Tracked::kR) < 0.01);
}

TEST_CASE("a track pulls with at most its friction times its weight, and spins past it") {
    Tracked held;
    held.drive(0.5);
    ::drive(held.w, held.sprocket.joint, 4.0, 1e4);
    held.drive(2.0, 0.5 * Tracked::kFriction * held.weight());
    // The sweep's crawl under load: driven wheels in this place creep 8 mm/s, the track 2.5.
    CHECK_MESSAGE(std::fabs(held.w.velocity(held.block).x + held.w.beltSpeed(held.track)) < 0.01 * 4.0 * Tracked::kR,
                  "half its grip: it holds");

    Tracked spun;
    spun.drive(0.5);
    ::drive(spun.w, spun.sprocket.joint, 4.0, 1e4);
    const double push = 1.5 * Tracked::kFriction * spun.weight();
    spun.drive(0.5, push);
    const double v0 = spun.w.velocity(spun.block).x;
    spun.drive(1.0, push);
    const double pulled = push - spun.mass() * (spun.w.velocity(spun.block).x - v0);
    CHECK_MESSAGE(std::fabs(pulled / (Tracked::kFriction * spun.weight()) - 1.0) < 0.02, "pulled with " << pulled);
    CHECK_MESSAGE(spun.w.velocity(spun.block).x > 0.0, "pushed back while its track spins forward");
    CHECK(std::fabs(spun.w.beltSpeed(spun.track) - 4.0 * Tracked::kR) < 1e-5);
}

TEST_CASE("a conveyor carries a box along its top at its speed") {
    World w;
    const double R = 0.1;
    const Guide head = guideAt(w, JointType::Motor, Vec3{0, 0, 0}, R), tail = guideAt(w, JointType::Hinge, Vec3{1, 0, 0}, R);
    BeltDesc d;
    d.length = 2.0 + 2.0 * kPi * R;
    d.mass = 3.0;
    const BeltGuide loop[] = {meshing(head), meshing(tail)};
    const BeltId belt = w.addBelt(d, loop);
    const BodyId box = slab(w, Vec3{0.5, R + 0.1, 0}, Vec3{0.1, 0.1, 0.1}, 2.0);
    drive(w, head.joint, 3.0, 1e4);
    const Vec3 top = straight(w, belt, 1);
    REQUIRE(top.x < -0.99);
    for (int s = 0; s < 60; ++s) {
        w.clearContacts();
        const Vec3 c = w.position(box);
        for (int k = 0; k < 4; ++k) {
            Contact t;
            t.a = box;
            t.point = c + Vec3{k & 1 ? 0.1 : -0.1, -0.1, k & 2 ? 0.1 : -0.1};
            t.b = t.point.x < 0.5 ? head.body : tail.body;
            t.normal = Vec3{0, -1, 0};
            t.depth = R - t.point.y;
            t.feature = std::uint64_t(k);
            t.belt = belt;
            t.beltDirection = top;
            w.addContact(t);
        }
        w.step(kDt);
    }
    CHECK(std::fabs(w.beltSpeed(belt) - 3.0 * R) < 1e-6);
    CHECK_MESSAGE(std::fabs(dot(w.velocity(box), top) - w.beltSpeed(belt)) < 1e-3, "riding the belt");
    CHECK(std::fabs(w.position(box).y - (R + 0.1)) < 0.005);
}
