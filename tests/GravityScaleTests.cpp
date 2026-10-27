#include <cmath>
#include <limits>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"

using namespace hf;
using namespace build;

namespace {

    constexpr double kG = 9.81;
    const Vec3 kZ{0, 0, 1};

    WorldSettings weightless() {
        WorldSettings s;
        s.gravity = Vec3{};
        return s;
    }

    BodyId ball(World &w, const Vec3 &at, double mass = 1.0, double gravityScale = 1.0, double damping = 0.0) {
        BodyDesc d;
        d.position = at;
        d.mass = mass;
        d.inertia = Mat3(0.4 * mass * 0.25);
        d.gravityScale = gravityScale;
        d.damping = damping;
        return w.addBody(d);
    }

    double substep(const World &w) { return kDt / w.settings().substeps; }

    double fallAfter(double gravityScale, double seconds) {
        World w;
        const BodyId b = ball(w, Vec3{}, 1.0, gravityScale);
        run(w, seconds);
        return -w.position(b).y;
    }

    // Two rods hinged end to end, free in the air: a machine of two links.
    struct Arm {
        BodyId a, b;
        JointId hinge;
    };

    Arm arm(World &w, JointType type = JointType::Hinge) {
        Arm m;
        m.a = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 2.0);
        m.b = rodBetween(w, Vec3{2, 0, 0}, Vec3{4, 0, 0}, 6.0);
        m.hinge = w.addJoint(jointAt(w, type, m.a, m.b, Vec3{2, 0, 0}, kZ));
        return m;
    }

    Vec3 momentum(const World &w, const Arm &m) { return 2.0 * w.velocity(m.a) + 6.0 * w.velocity(m.b); }

} // namespace

TEST_CASE("a body without gravity floats where it was left, sleeps, and falls once given its weight back") {
    World w;
    const BodyId b = ball(w, Vec3{0, 5, 0}, 1.0, 0.0);
    run(w, 3.0);
    CHECK_MESSAGE((w.position(b) == Vec3{0, 5, 0}), "it doesn't move by a bit");
    CHECK_MESSAGE(w.asleep(b), "still, it sleeps like a body at rest");

    CHECK(w.setGravityScale(b, 1.0));
    CHECK_MESSAGE(!w.asleep(b), "its weight back wakes it");
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.velocity(b).y + kG) < 1e-9, "and it falls as any body does");
}

TEST_CASE("a body at half gravity falls half as far") {
    CHECK(std::fabs(fallAfter(0.5, 1.0) / fallAfter(1.0, 1.0) - 0.5) < 1e-12);
    World w;
    const BodyId up = ball(w, Vec3{}, 1.0, -1.0);
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.velocity(up).y - kG) < 1e-9, "a negative scale lifts it");
}

TEST_CASE("a gravity scale or a damping the world can't use is refused") {
    World w;
    BodyDesc d;
    d.gravityScale = std::nan("");
    CHECK(!w.addBody(d).valid());
    d.gravityScale = 1.0;
    d.damping = -1.0;
    CHECK(!w.addBody(d).valid());
    const BodyId b = ball(w, Vec3{});
    CHECK(!w.setGravityScale(b, std::numeric_limits<double>::infinity()));
    CHECK(!w.setDamping(b, -0.1));
    CHECK(!w.setDamping(b, std::nan("")));
    w.step(kDt);
    CHECK(w.stats().rejected == 5);
    CHECK(w.setDamping(b, 0.0));
}

TEST_CASE("bodies welded together fall at their share of gravity by mass") {
    World w;
    const BodyId light = ball(w, Vec3{0, 0, 0}, 1.0, 1.0);
    const BodyId heavy = ball(w, Vec3{1, 0, 0}, 3.0, 0.0);
    w.addJoint(jointAt(w, JointType::Fixed, light, heavy, Vec3{0.5, 0, 0}, kZ));
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.velocity(heavy).y + kG / 4.0) < 1e-9, "a quarter of the mass feels gravity");
    CHECK(std::fabs(w.velocity(light).y + kG / 4.0) < 1e-9);

    w.setGravityScale(light, 0.0);
    run(w, 1.0);
    CHECK_MESSAGE(std::fabs(w.velocity(heavy).y + kG / 4.0) < 1e-9, "none of it now: it carries on at that speed");
}

TEST_CASE("damping slows a body by 1 + c h every sub-step, and holds a fall at g / c") {
    World spun(weightless());
    BodyDesc d;
    d.velocity = Vec3{2, 0, 0};
    d.angularVelocity = Vec3{0, 3, 0};
    d.damping = 4.0;
    const BodyId c = spun.addBody(d);
    run(spun, 1.0);
    const double k = std::pow(1.0 + 4.0 * substep(spun), -60 * spun.settings().substeps);
    CHECK(std::fabs(spun.velocity(c).x - 2.0 * k) < 1e-12);
    CHECK(std::fabs(spun.angularVelocity(c).y - 3.0 * k) < 1e-12);

    World falling;
    const BodyId f = ball(falling, Vec3{});
    CHECK(falling.setDamping(f, 2.0));
    run(falling, 10.0);
    CHECK_MESSAGE(std::fabs(falling.velocity(f).y + kG / 2.0) < 1e-6, "its terminal speed");
}

TEST_CASE("a machine's links each feel their own share of gravity") {
    World w;
    const Arm m = arm(w);
    w.step(kDt);
    REQUIRE(w.jointInfo(m.hinge).machine >= 0);

    w.setGravityScale(m.a, 0.0);
    const Vec3 p0 = momentum(w, m);
    run(w, 1.0);
    // A swinging arm's momentum drifts by a tenth of a percent without any force on it: 6 g, not 8 g.
    CHECK_MESSAGE(std::fabs(momentum(w, m).y - p0.y + 6.0 * kG) < 1e-2 * 6.0 * kG, "the heavy link's weight alone");

    World still;
    const Arm n = arm(still);
    still.setGravityScale(n.a, 0.0);
    still.setGravityScale(n.b, 0.0);
    const Vec3 at = still.position(n.b);
    run(still, 2.0);
    CHECK_MESSAGE(length(still.position(n.b) - at) < 1e-12, "with none, it floats where it was built");
}

TEST_CASE("a damped machine thrown through the air slows like a damped body") {
    World w(weightless());
    const Arm m = arm(w);
    w.setDamping(m.a, 2.0);
    w.setDamping(m.b, 2.0);
    w.step(kDt);
    REQUIRE(w.setMachineVelocity(m.a, Vec3{3, 0, 0}, Vec3{0, 0, 1}));
    const double p0 = length(momentum(w, m));
    run(w, 1.0);
    CHECK(std::fabs(length(momentum(w, m)) / p0 - std::exp(-2.0)) < 1e-2);
}

TEST_CASE("a motor turns a damped wheel at its speed, or as fast as its cap holds against the damping") {
    // A 10 kg, 2 m rod turning about its end: 13.3 kg m² there, so damping of 5/s takes 200 N·m at 3 rad/s.
    const auto spin = [](double maxForce) {
        World w(weightless());
        const BodyId rod = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 10.0);
        w.setDamping(rod, 5.0);
        JointDesc d = jointAt(w, JointType::Motor, BodyId{}, rod, Vec3{}, kZ);
        d.target = 3.0;
        d.maxForce = maxForce;
        const JointId motor = w.addJoint(d);
        run(w, 2.0);
        return w.jointVelocity(motor);
    };
    CHECK(std::fabs(spin(1e3) - 3.0) < 1e-3);
    const double held = spin(100.0);
    CHECK_MESSAGE(std::fabs(held - 100.0 / (5.0 * 40.0 / 3.0)) < 0.05 * held, "at its cap: 1.5 rad/s");
}
