#include <cmath>
#include <vector>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};

    BodyId addCube(World &w, const Vec3 &at, bool isStatic = false) {
        BodyDesc d;
        d.position = at;
        d.inertia = boxInertia(1.0, Vec3{0.5, 0.5, 0.5});
        d.isStatic = isStatic;
        return w.addBody(d);
    }

    // An arm on a servo from the world, held level: x along it, its weight on the servo.
    struct Arm {
        BodyId body;
        JointId servo;
    };

    Arm addArm(World &w, const Vec3 &at) {
        const BodyId arm = rodBetween(w, at, at + Vec3{2, 0, 0}, 50.0);
        JointDesc s = jointAt(w, JointType::Servo, BodyId{}, arm, at, kZ);
        s.maxForce = 1e5;
        s.min = -1.5;
        s.max = 1.5;
        return {arm, w.addJoint(s)};
    }

    WorldStats stepped(World &w) {
        w.step(kDt);
        return w.stats();
    }

} // namespace

TEST_CASE("bodies added one a step beside five machines rebuild nothing") {
    World w;
    std::vector<Arm> arms;
    for (int i = 0; i < 5; ++i) arms.push_back(addArm(w, Vec3{4.0 * i, 3, 0}));
    stepped(w);
    run(w, 0.5);
    int groupRebuilds = 0, machineRebuilds = 0;
    for (int i = 0; i < 1000; ++i) {
        addCube(w, Vec3{double(i % 40) * 2.0, 20.0 + double(i / 40) * 2.0, 10.0}); // apart: none touches another
        const WorldStats s = stepped(w);
        groupRebuilds += s.groupRebuilds;
        machineRebuilds += s.machineRebuilds;
    }
    CHECK_MESSAGE(groupRebuilds == 0, "each body took a group of its own, in place");
    CHECK_MESSAGE(machineRebuilds == 0, "and no machine was built again: their rows stayed warm");
    CHECK_MESSAGE(w.stats().bodies == 1005, "every body is there");
    for (const Arm &a : arms)
        CHECK_MESSAGE(std::fabs(w.jointPosition(a.servo)) < 1e-3, "and every arm still holds level");
}

TEST_CASE("a joint changed on one machine rebuilds that machine alone") {
    World w;
    std::vector<Arm> arms;
    for (int i = 0; i < 5; ++i) arms.push_back(addArm(w, Vec3{4.0 * i, 3, 0}));
    run(w, 0.5);
    // A second link hinged to the end of the third arm.
    const BodyId tip = rodBetween(w, Vec3{10, 3, 0}, Vec3{11, 3, 0}, 5.0);
    stepped(w);
    w.addJoint(jointAt(w, JointType::Hinge, arms[2].body, tip, Vec3{10, 3, 0}, kZ));
    const WorldStats s = stepped(w);
    CHECK_MESSAGE(s.groupRebuilds == 0, "no group was rebuilt");
    CHECK_MESSAGE(s.machineRebuilds == 1, "and one machine: the arm and its new link");
    CHECK_MESSAGE(w.jointInfo(arms[2].servo).machine >= 0, "the arm is still a machine");
    run(w, 0.4);
    CHECK_MESSAGE(w.position(tip).y < 2.8, "the new link swings down from the arm's end");
    CHECK_MESSAGE(std::fabs(length(w.position(tip) - Vec3{10, 3, 0}) - 0.5) < 1e-3, "on its pin");
    for (int i : {0, 1, 3, 4})
        CHECK_MESSAGE(std::fabs(w.jointPosition(arms[std::size_t(i)].servo)) < 1e-3, "the other arms hold on");
}

TEST_CASE("a body removed takes its group with it, and the next one added reuses it") {
    World w;
    const BodyId ground = addCube(w, Vec3{0, -0.5, 0}, true);
    std::vector<BodyId> cubes;
    for (int i = 0; i < 4; ++i) cubes.push_back(addCube(w, Vec3{3.0 * i, 5, 0}));
    stepped(w);
    w.removeBody(cubes[1]);
    const BodyId again = addCube(w, Vec3{3, 8, 0});
    const WorldStats s = stepped(w);
    CHECK_MESSAGE(s.groupRebuilds == 0, "removed and added in place");
    CHECK_MESSAGE(!w.contains(cubes[1]), "the removed body is gone");
    CHECK_MESSAGE(w.velocity(again).y < 0.0, "the new one falls");
    CHECK_MESSAGE(w.velocity(cubes[2]).y < 0.0, "and so do the others");
    (void)ground;
}

TEST_CASE("a machine link removed rebuilds its machine without it") {
    World w;
    const BodyId a = rodBetween(w, Vec3{0, 3, 0}, Vec3{1, 3, 0}, 5.0),
                 b = rodBetween(w, Vec3{1, 3, 0}, Vec3{2, 3, 0}, 5.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, a, Vec3{0, 3, 0}, kZ));
    w.addJoint(jointAt(w, JointType::Hinge, a, b, Vec3{1, 3, 0}, kZ));
    const BodyId lone = rodBetween(w, Vec3{5, 3, 0}, Vec3{6, 3, 0}, 5.0);
    const JointId pin = w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, lone, Vec3{5, 3, 0}, kZ));
    stepped(w);
    w.removeBody(b);
    w.removeBody(lone);
    const WorldStats s = stepped(w);
    CHECK_MESSAGE(s.groupRebuilds == 0, "no group rebuilt");
    CHECK_MESSAGE(s.machineRebuilds == 1, "the pendulum rebuilt without its second link, the lone one gone");
    CHECK_MESSAGE(s.machines == 1, "one machine left");
    CHECK_MESSAGE(!w.contains(pin), "the lone link's joint went with it");
    run(w, 0.5);
    CHECK_MESSAGE(length(w.position(a) - Vec3{0, 3, 0}) < 0.5 + 1e-3, "the first link swings on its pin");
}

TEST_CASE("a body made static and let go again stays in place in the world") {
    World w(WorldSettings{});
    const BodyId cube = addCube(w, Vec3{0, 5, 0});
    stepped(w);
    w.setStatic(cube, true);
    WorldStats s = stepped(w);
    CHECK_MESSAGE(s.groupRebuilds == 0, "made static in place");
    const Vec3 held = w.position(cube);
    s = stepped(w);
    CHECK_MESSAGE((w.position(cube) == held), "and holds");
    w.setStatic(cube, false);
    s = stepped(w);
    CHECK_MESSAGE(s.groupRebuilds == 0, "let go in place");
    CHECK_MESSAGE(w.velocity(cube).y < 0.0, "and falls");
}

TEST_CASE("a machine link's mass changed in place: its drive feels it, nothing rebuilt") {
    World w;
    const Arm arm = addArm(w, Vec3{0, 3, 0});
    run(w, 1.0);
    const double before = w.driveForce(arm.servo);
    w.setMassProperties(arm.body, 100.0, boxInertia(100.0, Vec3{1.0, 0.1, 0.1}));
    const WorldStats s = stepped(w);
    CHECK_MESSAGE((s.groupRebuilds == 0 && s.machineRebuilds == 0), "nothing rebuilt");
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.driveForce(arm.servo) / before - 2.0) < 0.02, "twice the weight, twice the torque");
    CHECK_MESSAGE(std::fabs(w.jointPosition(arm.servo)) < 1e-3, "held level all the same");
}
