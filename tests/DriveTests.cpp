#include <cmath>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"
#include "host/Scenes.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};

    // A 600 kg boom on its pivot, lifted by a cylinder hinged under it: the piston carries the boom's weight.
    JointId boomOnCylinder(World &w, double maxForce = 2e5) {
        const Vec3 pivot{0.5, 1.6, 0}, tip{3.5, 3.6, 0}, base{1.0, 1.3, 0};
        const Vec3 onBoom = pivot + 1.5 * normalize(tip - pivot), mid = 0.5 * (base + onBoom);
        const Vec3 u = normalize(onBoom - base);
        const BodyId boom = rodBetween(w, pivot, tip, 600.0);
        const BodyId barrel = rodBetween(w, base, mid + 0.2 * u, 40.0);
        const BodyId rod = rodBetween(w, mid - 0.2 * u, onBoom, 20.0);
        w.addJoint(jointAt(w, JointType::Hinge, boom, BodyId{}, pivot, kZ));
        w.addJoint(jointAt(w, JointType::Hinge, barrel, BodyId{}, base, kZ));
        w.addJoint(jointAt(w, JointType::Hinge, rod, boom, onBoom, kZ));
        JointDesc p = jointAt(w, JointType::Piston, barrel, rod, mid, u);
        p.maxForce = maxForce;
        p.speed = 0.5;
        p.min = -0.3;
        p.max = 0.2;
        return w.addJoint(p);
    }

    // A 2 m, 50 kg arm out along x from a servo at the origin: gravity pulls it down with 490 N·m.
    JointId armOnServo(World &w, double maxForce) {
        const BodyId arm = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 50.0);
        JointDesc s = jointAt(w, JointType::Servo, BodyId{}, arm, Vec3{0, 0, 0}, kZ);
        s.maxForce = maxForce;
        s.speed = 0.5;
        s.min = -1.5;
        s.max = 1.5;
        return w.addJoint(s);
    }

} // namespace

TEST_CASE("a released lever parks a loaded piston where it got to") {
    World w;
    const JointId piston = boomOnCylinder(w);
    run(w, 0.5);
    w.setLever(piston, 1.0);
    run(w, 0.3);
    w.setLever(piston, 0.0);
    w.step(kDt);
    const double parked = w.jointPosition(piston);
    run(w, 2.0);
    CHECK_MESSAGE(parked > 0.08, "held, the lever moved it");
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - parked) < 1e-3,
                  "let go, it stops where it is instead of creeping on");
}

TEST_CASE("a lever held past the end of the course stops at the end") {
    World w;
    const JointId piston = boomOnCylinder(w);
    w.setLever(piston, 1.0);
    run(w, 1.5);
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - 0.2) < 1e-3, "at its stop, and no further");
    w.setLever(piston, 0.0);
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - 0.2) < 1e-3, "and held there");
}

TEST_CASE("a servo turns to its angle at its speed and holds it against gravity") {
    WorldSettings settings;
    World w(settings);
    const JointId servo = armOnServo(w, 5e3);
    run(w, 0.5);
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo)) < 1e-3, "it holds the arm out where it was built");
    w.setTarget(servo, 1.0);
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo) - 0.5) < 0.05, "half a radian a second");
    run(w, 2.0);
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo) - 1.0) < 1e-3, "it arrives");
    run(w, 2.0);
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo) - 1.0) < 1e-3, "and stays");
}

TEST_CASE("a weak servo can't lift what it can't carry") {
    World w;
    const JointId servo = armOnServo(w, 100.0);
    w.setTarget(servo, 1.0);
    run(w, 3.0);
    CHECK_MESSAGE(w.jointPosition(servo) < 0.05, "100 N·m against 490: the cap is a cap while it drives");
}

TEST_CASE("a servo's lever turns it, and parks it when let go") {
    World w;
    const JointId servo = armOnServo(w, 5e3);
    w.setLever(servo, -1.0);
    run(w, 0.5);
    w.setLever(servo, 0.0);
    w.step(kDt);
    const double parked = w.jointPosition(servo);
    run(w, 2.0);
    CHECK_MESSAGE(parked < -0.15, "it turned while held");
    CHECK_MESSAGE(std::fabs(w.jointPosition(servo) - parked) < 1e-3, "and stays where it was let go");
}

TEST_CASE("a motor's lever sets its speed") {
    WorldSettings settings;
    settings.gravity = Vec3{};
    World w(settings);
    const BodyId wheel = rodBetween(w, Vec3{-0.5, 0, 0}, Vec3{0.5, 0, 0}, 10.0);
    JointDesc m = jointAt(w, JointType::Motor, BodyId{}, wheel, Vec3{0, 0, 0}, kZ);
    m.speed = 4.0;
    const JointId motor = w.addJoint(m);
    w.setLever(motor, 0.5);
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.angularVelocity(wheel).z - 2.0) < 1e-6, "half a lever, half its speed");
}

TEST_CASE("twin cylinders authored a hair apart don't fight when they're born") {
    // Two rigid cylinders side by side under one boom, their targets 0.3 mm apart: two rows on one freedom.
    World w;
    const Vec3 pivot{0.5, 1.6, 0}, tip{3.5, 3.6, 0};
    const BodyId boom = rodBetween(w, pivot, tip, 600.0);
    w.addJoint(jointAt(w, JointType::Hinge, boom, BodyId{}, pivot, kZ));
    JointId pistons[2];
    for (int i = 0; i < 2; ++i) {
        const Vec3 side{0, 0, i ? 0.4 : -0.4};
        const Vec3 base = Vec3{1.0, 1.3, 0} + side, onBoom = pivot + 1.5 * normalize(tip - pivot) + side;
        const Vec3 mid = 0.5 * (base + onBoom), u = normalize(onBoom - base);
        const BodyId barrel = rodBetween(w, base, mid + 0.2 * u, 40.0);
        const BodyId rod = rodBetween(w, mid - 0.2 * u, onBoom, 20.0);
        w.addJoint(jointAt(w, JointType::Hinge, barrel, BodyId{}, base, kZ));
        w.addJoint(jointAt(w, JointType::Hinge, rod, boom, onBoom, kZ));
        JointDesc p = jointAt(w, JointType::Piston, barrel, rod, mid, u);
        p.target = i ? 3e-4 : 0.0;
        p.min = -0.3;
        p.max = 0.2;
        pistons[i] = w.addJoint(p);
    }
    w.step(kDt);
    const double first = std::max(std::fabs(w.driveForce(pistons[0])), std::fabs(w.driveForce(pistons[1])));
    run(w, 1.0);
    CHECK_MESSAGE(first < 1e5, "the pair carries the boom's 20 kN or so, not a fight between them");
    CHECK_MESSAGE(length(w.angularVelocity(boom)) < 1e-3, "and holds it still");
}

TEST_CASE("a held cylinder that gives while another one works comes back") {
    // Working the stick stresses the boom's loop. The boom's held cylinder yields to the loop's closure, as a servo
    // must; scored as if both mattered equally, the closing step was refused forever and the boom stayed down.
    host::Host h;
    const host::Excavator e = host::addExcavator(h, Vec3{});
    const auto run = [&](double seconds) {
        for (int i = 0; i < int(std::lround(seconds / kDt)); ++i) h.step(kDt);
    };
    h.world.setLever(e.boom, 1.0);
    run(0.3);
    h.world.setLever(e.boom, 0.0);
    run(0.5);
    const double held = h.world.jointPosition(e.boom);
    h.world.setLever(e.stick, -1.0);
    run(0.4);
    h.world.setLever(e.stick, 0.0);
    run(4.0);
    CHECK_MESSAGE(std::fabs(h.world.jointPosition(e.boom) - held) < 1e-4, "the boom is back where it was held");
}
