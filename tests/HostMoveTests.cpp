#include <cmath>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"
#include "host/Host.hpp"
#include "world/TestHooks.hpp"
#include "world/WorldImpl.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};

    WorldSettings weightless() {
        WorldSettings s;
        s.gravity = Vec3{};
        return s;
    }

    BodyId addCube(World &w, const Vec3 &at, bool isStatic = false) {
        BodyDesc d;
        d.position = at;
        d.inertia = boxInertia(1.0, Vec3{0.5, 0.5, 0.5});
        d.isStatic = isStatic;
        return w.addBody(d);
    }

    // The rigid move taking pose (fromC, fromQ) to (toC, toQ), applied to a point.
    Vec3 carried(const Vec3 &x, const Vec3 &fromC, const Quat &fromQ, const Vec3 &toC, const Quat &toQ) {
        return toC + normalize(toQ * conjugate(fromQ)) * (x - fromC);
    }

    // A chain of rods from start along step, each hinged to the last about z; the first to anchor, if anchored.
    std::vector<BodyId> chain(World &w, const Vec3 &start, int links, BodyId anchor, bool anchored,
                              std::vector<JointId> &joints, const Vec3 &step = Vec3{1, 0, 0}) {
        std::vector<BodyId> rods;
        for (int i = 0; i < links; ++i) {
            const Vec3 a = start + double(i) * step, b = a + step;
            rods.push_back(rodBetween(w, a, b, 5.0));
            if (i > 0)
                joints.push_back(
                    w.addJoint(jointAt(w, JointType::Hinge, rods[std::size_t(i) - 1], rods.back(), a, kZ)));
            else if (anchored) joints.push_back(w.addJoint(jointAt(w, JointType::Hinge, anchor, rods.back(), a, kZ)));
        }
        return rods;
    }

    std::vector<double> coordinates(World &w, const std::vector<JointId> &joints) {
        std::vector<double> q;
        for (JointId j : joints) q.push_back(w.jointPosition(j));
        return q;
    }

    // The warm start each contact begins the step with, read before any sweep has touched it.
    std::vector<double> warmStarts;
    void readWarmStarts(World::Impl &I, int substep, FaultPhase phase) {
        if (substep != 0 || phase != FaultPhase::Internal) return;
        warmStarts.clear();
        for (const ContactData &c : I.contacts) warmStarts.push_back(c.ln);
    }

} // namespace

TEST_CASE("a static body let go and thrown in the same step leaves at the throw") {
    // The rebuild that letting go asks for used to read every body back from where the world had it, the throw
    // lost: a velocity gizmo in Play did nothing.
    for (int order = 0; order < 2; ++order) {
        World w(weightless());
        const BodyId cube = addCube(w, Vec3{0, 1, 0}, true);
        w.step(kDt);
        if (order == 0) {
            w.setStatic(cube, false);
            w.setVelocity(cube, Vec3{3, 0, 0}, Vec3{0, 2, 0});
        } else {
            w.setVelocity(cube, Vec3{3, 0, 0}, Vec3{0, 2, 0});
            w.setStatic(cube, false);
        }
        w.step(kDt);
        CHECK_MESSAGE((w.velocity(cube) == Vec3{3, 0, 0}), "it leaves at the velocity thrown");
        CHECK_MESSAGE((w.angularVelocity(cube) == Vec3{0, 2, 0}), "and the spin");
        CHECK_MESSAGE(std::fabs(w.position(cube).x - 3 * kDt) < 1e-12, "and has gone that far");
    }
}

TEST_CASE("a static body let go and placed in the same step starts from there") {
    World w(weightless());
    const BodyId cube = addCube(w, Vec3{0, 1, 0}, true);
    w.step(kDt);
    w.setStatic(cube, false);
    w.setPose(cube, Vec3{5, 2, 1}, angleAxis(0.5, kZ));
    w.step(kDt);
    CHECK_MESSAGE((w.position(cube) == Vec3{5, 2, 1}), "where it was put");
    CHECK_MESSAGE(std::fabs(dot(w.orientation(cube), angleAxis(0.5, kZ))) > 1.0 - 1e-15, "turned as it was put");
}

TEST_CASE("a body the host moves starts its contacts cold") {
    World w;
    const BodyId ground = addCube(w, Vec3{0, -0.5, 0}, true);
    const BodyId cube = addCube(w, Vec3{0, 0.5, 0});
    const auto contacts = [&] {
        w.clearContacts();
        const Vec3 c = w.position(cube);
        for (int k = 0; k < 4; ++k) {
            Contact t;
            t.a = ground;
            t.b = cube;
            t.point = c + Vec3{k & 1 ? 0.5 : -0.5, -0.5, k & 2 ? 0.5 : -0.5};
            t.normal = Vec3{0, 1, 0};
            t.depth = -t.point.y;
            t.feature = std::uint64_t(k);
            w.addContact(t);
        }
    };
    for (int i = 0; i < 5; ++i) {
        contacts();
        w.step(kDt);
    }
    testing::phaseHook = readWarmStarts;
    contacts();
    w.step(kDt);
    CHECK_MESSAGE(warmStarts.size() == 4, "four contacts");
    CHECK_MESSAGE(warmStarts[0] > 0.0, "resting, a contact starts with what held it last step");
    w.setPose(cube, w.position(cube) + Vec3{2, 0, 0}, w.orientation(cube));
    contacts();
    w.step(kDt);
    testing::phaseHook = nullptr;
    for (double l : warmStarts) CHECK_MESSAGE(l == 0.0, "moved by the host, it starts with nothing");
}

TEST_CASE("a body's state reads in one call, and the step says what it moved") {
    World w;
    const BodyId ground = addCube(w, Vec3{0, -5, 0}, true);
    const BodyId falling = addCube(w, Vec3{0, 3, 0});
    w.step(kDt);
    const BodyState s = w.state(falling);
    CHECK((s.position == w.position(falling) && s.velocity == w.velocity(falling)));
    CHECK((s.orientation == w.orientation(falling) && s.angularVelocity == w.angularVelocity(falling)));
    REQUIRE(w.moved().size() == 1);
    CHECK_MESSAGE(w.moved()[0] == falling, "the falling body moved; the static one didn't");
    (void)ground;
}

TEST_CASE("a static post moved carries the chain hanging from it, its joints as they were") {
    for (bool asleep : {false, true}) {
        World w;
        const BodyId post = addCube(w, Vec3{0, 3, 0}, true);
        std::vector<JointId> joints;
        const std::vector<BodyId> rods = chain(w, Vec3{0, 2.5, 0}, 3, post, true, joints, Vec3{0, -1, 0});
        run(w, asleep ? 1.0 : 0.05);
        CHECK_MESSAGE(w.asleep(rods[0]) == asleep, "hanging still, the chain fell asleep");
        const std::vector<double> before = coordinates(w, joints);
        std::vector<Vec3> at;
        for (BodyId r : rods) at.push_back(w.position(r));

        const Vec3 from = w.position(post), to = Vec3{4, 5, -1};
        const Quat fromQ = w.orientation(post), toQ = angleAxis(0.4, Vec3{1, 0, 0}) * fromQ;
        w.setPose(post, to, toQ);
        CHECK_MESSAGE(!w.asleep(rods[0]), "the chain wakes");
        for (std::size_t i = 0; i < rods.size(); ++i)
            CHECK_MESSAGE(length(w.position(rods[i]) - carried(at[i], from, fromQ, to, toQ)) < 1e-12,
                          "each link went where the post took it");
        const std::vector<double> after = coordinates(w, joints);
        for (std::size_t i = 0; i < joints.size(); ++i)
            CHECK_MESSAGE(std::fabs(after[i] - before[i]) < 1e-12, "every joint as it was");
        w.step(kDt);
        CHECK_MESSAGE(length(w.position(rods[0]) - carried(at[0], from, fromQ, to, toQ)) < 0.01,
                      "and it hangs on from there");
    }
}

TEST_CASE("a floating machine moved and turned whole keeps its joints, its rates and its next step") {
    World moved(weightless()), control(weightless());
    std::vector<JointId> jm, jc;
    const std::vector<BodyId> rm = chain(moved, Vec3{}, 5, BodyId{}, false, jm);
    const std::vector<BodyId> rc = chain(control, Vec3{}, 5, BodyId{}, false, jc);
    for (World *w : {&moved, &control}) {
        const std::vector<BodyId> &r = w == &moved ? rm : rc;
        for (std::size_t i = 0; i < r.size(); ++i)
            w->setVelocity(r[i], Vec3{0, 0.3 * double(i), 0.1}, Vec3{0, 0, 0.2 * double(i % 2)});
        w->step(kDt);
    }
    const std::vector<double> before = coordinates(moved, jm);

    const Vec3 from = moved.position(rm[2]), to = from + Vec3{3, -2, 7};
    const Quat fromQ = moved.orientation(rm[2]), toQ = normalize(angleAxis(1.1, normalize(Vec3{1, 2, 3})) * fromQ);
    const Vec3 v0 = moved.velocity(rm[4]);
    REQUIRE(moved.moveMachine(rm[2], to, toQ));
    const std::vector<double> after = coordinates(moved, jm);
    for (std::size_t i = 0; i < jm.size(); ++i)
        CHECK_MESSAGE(std::fabs(after[i] - before[i]) < 1e-12, "every joint coordinate as it was");
    CHECK_MESSAGE(length(moved.velocity(rm[4]) - normalize(toQ * conjugate(fromQ)) * v0) < 1e-12,
                  "its velocities turned with it");

    moved.step(kDt);
    control.step(kDt);
    const std::vector<double> m = coordinates(moved, jm), c = coordinates(control, jc);
    for (std::size_t i = 0; i < jm.size(); ++i)
        CHECK_MESSAGE(std::fabs(m[i] - c[i]) < 1e-9, "and steps on exactly as it would have where it was");
}

TEST_CASE("an anchored machine can't be moved away from its anchor") {
    World w;
    const BodyId post = addCube(w, Vec3{0, 3, 0}, true);
    std::vector<JointId> joints;
    const std::vector<BodyId> rods = chain(w, Vec3{0.5, 3, 0}, 3, post, true, joints);
    w.step(kDt);
    const Vec3 at = w.position(rods[1]);
    CHECK_MESSAGE(!w.moveMachine(rods[1], at + Vec3{1, 0, 0}, w.orientation(rods[1])), "refused");
    CHECK_MESSAGE((w.position(rods[1]) == at), "and nothing moved");
    CHECK_MESSAGE(w.moveMachine(rods[1], at, w.orientation(rods[1])), "where it already is, it is no move at all");
    CHECK_MESSAGE(!w.setMachineVelocity(rods[1], Vec3{1, 0, 0}, Vec3{}), "nor thrown");
}

TEST_CASE("a floating machine thrown whole flies as one body") {
    World w(weightless());
    std::vector<JointId> joints;
    const std::vector<BodyId> rods = chain(w, Vec3{}, 4, BodyId{}, false, joints);
    w.step(kDt);
    const std::vector<double> before = coordinates(w, joints);
    const Vec3 v{1, 2, 0}, om{0, 0, 0.5};
    const Vec3 about = w.position(rods[1]);
    REQUIRE(w.setMachineVelocity(rods[1], v, om));
    for (BodyId r : rods)
        CHECK_MESSAGE(length(w.velocity(r) - (v + cross(om, w.position(r) - about))) < 1e-12,
                      "every link moving with the one rigid motion");
    run(w, 1.0);
    const std::vector<double> after = coordinates(w, joints);
    for (std::size_t i = 0; i < joints.size(); ++i)
        CHECK_MESSAGE(std::fabs(after[i] - before[i]) < 1e-6, "and no joint moves");
}

TEST_CASE("a kinematic paddle sweeps a row of boxes; let go, it flies on at the sweep speed") {
    host::Host h;
    const BodyId paddle = h.addBox(Vec3{-1.5, 0.6, 0}, Quat{}, Vec3{0.1, 0.5, 1.0}, 50.0);
    std::vector<BodyId> boxes;
    for (int i = 0; i < 3; ++i)
        boxes.push_back(h.addBox(Vec3{0.5 + 1.2 * i, 0.5, 0}, Quat{}, Vec3{0.5, 0.5, 0.5}, 10.0));
    for (int i = 0; i < 30; ++i) h.step(kDt);
    REQUIRE(h.world.setKinematic(paddle, true));
    const double speed = 2.0, x0 = h.world.position(paddle).x;
    for (int i = 0; i < 90; ++i) {
        const Vec3 p = h.world.position(paddle);
        h.world.setPose(paddle, p + Vec3{speed * kDt, 0, 0}, h.world.orientation(paddle));
        h.step(kDt);
        CHECK(h.world.stats().groupRebuilds == 0);
    }
    CHECK_MESSAGE(std::fabs(h.world.position(paddle).x - (x0 + 90 * speed * kDt)) < 1e-9,
                  "the paddle went where the host put it, the boxes pushing back on nothing");
    CHECK_MESSAGE(h.world.position(boxes[0]).x > 1.0, "and pushed the boxes ahead of it");
    REQUIRE(h.world.setKinematic(paddle, false));
    CHECK_MESSAGE(std::fabs(h.world.velocity(paddle).x - speed) < 1e-9, "let go, it keeps the sweep's speed");
    h.step(kDt);
    CHECK_MESSAGE(h.world.velocity(paddle).x > 0.8 * speed, "and flies on");
    CHECK_MESSAGE(h.world.stats().groupRebuilds == 0, "let go without rebuilding the world");
}

TEST_CASE("a kinematic platform carries a box resting on it, and wakes it to") {
    host::Host h;
    const BodyId platform = h.addBox(Vec3{0, 0.25, 0}, Quat{}, Vec3{3.0, 0.25, 3.0}, 100.0);
    REQUIRE(h.world.setKinematic(platform, true));
    const BodyId box = h.addBox(Vec3{0, 1.0, 0}, Quat{}, Vec3{0.5, 0.5, 0.5}, 10.0);
    run(h.world, 0.0);
    for (int i = 0; i < 600 && !h.world.asleep(box); ++i) h.step(kDt);
    REQUIRE(h.world.asleep(box));
    const double x0 = h.world.position(box).x;
    for (int i = 0; i < 60; ++i) {
        const Vec3 p = h.world.position(platform);
        h.world.setPose(platform, p + Vec3{1.0 * kDt, 0, 0}, h.world.orientation(platform));
        h.step(kDt);
    }
    CHECK_MESSAGE(!h.world.asleep(box), "the platform moving under it woke it");
    CHECK_MESSAGE(std::fabs(h.world.position(box).x - x0 - 1.0) < 0.1, "and carried it the metre it went");
    CHECK_MESSAGE(h.world.position(platform).y == 0.25, "nothing pushed the platform down");
}

TEST_CASE("a machine link can't be made kinematic") {
    World w;
    std::vector<JointId> joints;
    const std::vector<BodyId> rods = chain(w, Vec3{0, 2, 0}, 2, BodyId{}, true, joints);
    CHECK(!w.setKinematic(rods[1], true));
    CHECK(!w.kinematic(rods[1]));
}

TEST_CASE("a servo set to half a radian by hand holds there, every other joint as it was") {
    World w;
    std::vector<BodyId> rods;
    std::vector<JointId> servos;
    for (int i = 0; i < 3; ++i) {
        const Vec3 a{double(i), 2, 0};
        rods.push_back(rodBetween(w, a, a + Vec3{1, 0, 0}, 5.0));
        JointDesc s = jointAt(w, JointType::Servo, i ? rods[std::size_t(i) - 1] : BodyId{}, rods.back(), a, kZ);
        s.maxForce = 1e4;
        s.min = -2.0;
        s.max = 2.0;
        servos.push_back(w.addJoint(s));
    }
    run(w, 0.5);
    const std::vector<double> before = coordinates(w, servos);
    REQUIRE(w.setJointPosition(servos[1], 0.5));
    CHECK_MESSAGE(std::fabs(w.jointPosition(servos[1]) - 0.5) < 1e-12, "it is there at once");
    CHECK_MESSAGE(std::fabs(w.jointPosition(servos[0]) - before[0]) < 1e-12, "the joint before it as it was");
    CHECK_MESSAGE(std::fabs(w.jointPosition(servos[2]) - before[2]) < 1e-12, "and the one after it");
    CHECK_MESSAGE(w.target(servos[1]) == 0.5, "its drive holds there");
    for (int i = 0; i < 300; ++i) {
        w.step(kDt);
        if (std::fabs(w.jointPosition(servos[1]) - 0.5) > 1e-3) break;
    }
    CHECK_MESSAGE(std::fabs(w.jointPosition(servos[1]) - 0.5) < 1e-3, "and holds half a radian five seconds on");
    CHECK_MESSAGE(std::fabs(w.jointPosition(servos[0]) - before[0]) < 1e-3, "every other joint holding where it was");
    CHECK_MESSAGE(std::fabs(w.jointPosition(servos[2]) - before[2]) < 1e-3, "this one too");
}

TEST_CASE("a piston set mid-course by hand holds there within its deadband") {
    World w;
    const BodyId ram = rodBetween(w, Vec3{0, 1, 0}, Vec3{0, 2, 0}, 20.0);
    JointDesc p = jointAt(w, JointType::Piston, BodyId{}, ram, Vec3{0, 1, 0}, Vec3{0, 1, 0});
    p.maxForce = 1e5;
    p.min = -0.5;
    p.max = 0.5;
    const JointId piston = w.addJoint(p);
    run(w, 0.5);
    REQUIRE(w.setJointPosition(piston, 0.3));
    CHECK_MESSAGE(std::fabs(w.position(ram).y - 1.8) < 1e-9, "the ram went 30 cm up its axis");
    run(w, 5.0);
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - 0.3) < 1e-3, "and holds there");
    CHECK_MESSAGE(w.setJointPosition(piston, 9.0), "past the stop, it is clamped");
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - 0.5) < 1e-12, "to the course's end");
}

TEST_CASE("only a tree joint of one coordinate can be set by hand") {
    World w;
    const BodyId a = rodBetween(w, Vec3{0, 2, 0}, Vec3{1, 2, 0}, 5.0),
                 b = rodBetween(w, Vec3{1, 2, 0}, Vec3{2, 2, 0}, 5.0);
    const JointId ball = w.addJoint(jointAt(w, JointType::Ball, BodyId{}, a, Vec3{0, 2, 0}, kZ));
    w.addJoint(jointAt(w, JointType::Hinge, a, b, Vec3{1, 2, 0}, kZ));
    // A third link closing a loop back to the world: one hinge of it is a cut.
    const BodyId c = rodBetween(w, Vec3{2, 2, 0}, Vec3{2, 1, 0}, 5.0);
    const JointId h1 = w.addJoint(jointAt(w, JointType::Hinge, b, c, Vec3{2, 2, 0}, kZ));
    const JointId h2 = w.addJoint(jointAt(w, JointType::Hinge, c, BodyId{}, Vec3{2, 1, 0}, kZ));
    w.step(kDt);
    CHECK_MESSAGE(!w.setJointPosition(ball, 0.1), "a ball joint has three coordinates");
    const JointId cut = w.jointInfo(h1).tree ? h2 : h1;
    REQUIRE(!w.jointInfo(cut).tree);
    CHECK_MESSAGE(!w.setJointPosition(cut, 0.1), "a loop's cut isn't a coordinate at all");
}
