#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};
    constexpr double kPi = std::numbers::pi;
    constexpr double kInfinity = std::numeric_limits<double>::infinity();

    // A 2 m, 20 kg rod out along x from a joint of `type` at the origin, about z.
    struct Arm {
        BodyId rod;
        JointId joint;
    };

    Arm armOn(World &w, JointType type, double lower, double upper) {
        Arm a;
        a.rod = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 20.0);
        JointDesc j = jointAt(w, type, BodyId{}, a.rod, Vec3{0, 0, 0}, kZ);
        j.lowerAngle = lower;
        j.upperAngle = upper;
        a.joint = w.addJoint(j);
        return a;
    }

    // The angle between two directions.
    double between(const Vec3 &a, const Vec3 &b) { return std::atan2(length(cross(a, b)), dot(a, b)); }

    // A rod hanging 2 m down from a ball at the origin, its joint's axis down the rod, kicked sideways.
    struct Pendulum {
        BodyId rod;
        JointId ball;
        Quat frameB;
    };

    Pendulum pendulum(World &w, double cone, double kick) {
        Pendulum p;
        p.rod = rodBetween(w, Vec3{0, 0, 0}, Vec3{0, -2, 0}, 20.0);
        JointDesc j = jointAt(w, JointType::Ball, BodyId{}, p.rod, Vec3{0, 0, 0}, Vec3{0, -1, 0});
        j.coneAngle = cone;
        p.ball = w.addJoint(j);
        p.frameB = j.frameB.rotation;
        w.setVelocity(p.rod, Vec3{kick, 0, 0}, Vec3{0, 0, kick});
        return p;
    }

    // How far the rod swings off straight down: the ball's axis on the rod against the world's.
    double swing(World &w, const Pendulum &p) {
        return between(w.orientation(p.rod) * p.frameB * Vec3{1, 0, 0}, Vec3{0, -1, 0});
    }

} // namespace

TEST_CASE("a hinge falls onto the end of its window and stays there") {
    World w;
    const Arm free = armOn(w, JointType::Hinge, -kInfinity, kInfinity);
    World stopped;
    const Arm a = armOn(stopped, JointType::Hinge, -0.5, 0.8);
    double lowest = 0.0, lowestFree = 0.0;
    for (int s = 0; s < 180; ++s) {
        w.step(kDt);
        stopped.step(kDt);
        lowestFree = std::min(lowestFree, w.jointPosition(free.joint));
        lowest = std::min(lowest, stopped.jointPosition(a.joint));
    }
    CHECK_MESSAGE(lowestFree < -1.0, "free, the arm swings down past where the window would end");
    CHECK_MESSAGE(lowest > -0.5 - 0.01, "the arm never falls more than 0.01 rad past its stop");
    CHECK_MESSAGE(std::fabs(stopped.jointPosition(a.joint) + 0.5) < 2e-3, "and rests on it");
}

TEST_CASE("a motor driven into the end of its window is held there") {
    World w;
    const Arm a = armOn(w, JointType::Motor, -0.4, 0.6);
    JointDesc d = w.jointDesc(a.joint);
    d.target = 1.0;
    d.maxForce = 1e4;
    w.setJoint(a.joint, d);
    double highest = 0.0;
    for (int s = 0; s < 120; ++s) {
        w.step(kDt);
        highest = std::max(highest, w.jointPosition(a.joint));
    }
    CHECK_MESSAGE(highest < 0.6 + 0.01, "driven up, it stops at its upper stop");
    CHECK_MESSAGE(std::fabs(w.jointPosition(a.joint) - 0.6) < 2e-3, "and is held there while it still drives");
}

TEST_CASE("a window across the half turn keeps out of the gap") {
    World w;
    const double edge = kPi - 0.17;
    const Arm a = armOn(w, JointType::Motor, -edge, edge);
    JointDesc d = w.jointDesc(a.joint);
    d.target = 2.0;
    d.maxForce = 1e4;
    w.setJoint(a.joint, d);
    double highest = 0.0;
    bool wrapped = false;
    for (int s = 0; s < 180; ++s) {
        w.step(kDt);
        const double angle = w.jointPosition(a.joint);
        highest = std::max(highest, angle);
        wrapped = wrapped || (s > 30 && angle < 0.0);
    }
    CHECK_MESSAGE(!wrapped, "turned toward the half turn, it never comes round the other side");
    CHECK_MESSAGE(std::fabs(w.jointPosition(a.joint) - edge) < 2e-3, "it stops 10 degrees short of it");
    CHECK_MESSAGE(highest < edge + 0.01, "and never reaches into the gap");
}

TEST_CASE("a ball's cone keeps a kicked pendulum inside it") {
    World w;
    const Pendulum free = pendulum(w, kInfinity, 2.0);
    World coned;
    const Pendulum p = pendulum(coned, 0.35, 2.0);
    double widest = 0.0, widestFree = 0.0;
    for (int s = 0; s < 240; ++s) {
        w.step(kDt);
        coned.step(kDt);
        widestFree = std::max(widestFree, swing(w, free));
        widest = std::max(widest, swing(coned, p));
    }
    CHECK_MESSAGE(widestFree > 0.6, "free, the kick swings it well past the cone");
    CHECK_MESSAGE(widest < 0.35 + 0.01, "in its cone, it never swings more than 0.01 rad past it");
    CHECK_MESSAGE(widest > 0.34, "and reaches it");
}

TEST_CASE("a cone tilted off the joint's axis lets the pendulum swing further one way than the other") {
    // A 25 degree cone tilted 20 degrees toward +x off straight down: the rod hanging straight down is 20 degrees
    // off the cone's middle. Kicked toward +x it reaches the cone at 45 degrees; toward -x, at 5.
    const auto widest = [](double kick) {
        World w;
        const BodyId rod = rodBetween(w, Vec3{0, 0, 0}, Vec3{0, -2, 0}, 20.0);
        JointDesc j = jointAt(w, JointType::Ball, BodyId{}, rod, Vec3{0, 0, 0}, Vec3{0, -1, 0});
        j.coneAngle = 25.0 * kPi / 180.0;
        const double tilt = 20.0 * kPi / 180.0;
        j.coneAxis = conjugate(j.frameA.rotation) * Vec3{std::sin(tilt), -std::cos(tilt), 0.0};
        w.addJoint(j);
        w.setVelocity(rod, Vec3{kick, 0, 0}, Vec3{0, 0, kick});
        double most = 0.0;
        for (int s = 0; s < 120; ++s) {
            w.step(kDt);
            const Vec3 down = w.orientation(rod) * j.frameB.rotation * Vec3{1, 0, 0};
            const double swing = std::atan2(down.x, -down.y) * 180.0 / kPi;
            most = kick > 0.0 ? std::max(most, swing) : std::min(most, swing);
        }
        return most;
    };
    const double open = widest(2.5), narrow = widest(-2.5);
    CHECK_MESSAGE(open > 44.4, "toward the tilt it swings out to the cone's far side, 45 degrees off straight down");
    CHECK_MESSAGE(open < 45.6, "and no further");
    CHECK_MESSAGE(narrow < -4.4, "the other way it swings to the cone's near side, 5 degrees off straight down");
    CHECK_MESSAGE(narrow > -5.6, "and no further");
}

TEST_CASE("a cone holds between two moving links") {
    World w;
    // An upper arm hinged to the world, swinging; a lower arm on a ball at its end, its axis along the upper's.
    const BodyId upper = rodBetween(w, Vec3{0, 0, 0}, Vec3{1.5, 0, 0}, 15.0);
    const BodyId lower = rodBetween(w, Vec3{1.5, 0, 0}, Vec3{3, 0, 0}, 10.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, upper, Vec3{0, 0, 0}, kZ));
    JointDesc j = jointAt(w, JointType::Ball, upper, lower, Vec3{1.5, 0, 0}, Vec3{1, 0, 0});
    j.coneAngle = 0.26;
    w.addJoint(j);
    w.setVelocity(lower, Vec3{0, 0, 3.0}, Vec3{0, 2.0, 0});
    double widest = 0.0;
    for (int s = 0; s < 240; ++s) {
        w.step(kDt);
        const Vec3 xA = w.orientation(upper) * j.frameA.rotation * Vec3{1, 0, 0};
        const Vec3 xB = w.orientation(lower) * j.frameB.rotation * Vec3{1, 0, 0};
        widest = std::max(widest, between(xA, xB));
    }
    CHECK_MESSAGE(widest < 0.26 + 0.01, "the lower arm never swings more than 0.01 rad past its cone off the upper");
    CHECK_MESSAGE(widest > 0.25, "and reaches it");
}

TEST_CASE("a cone turns with the link it is on") {
    // A double pendulum in one plane: the upper arm falls about z, and the lower arm on its ball swings off it about z
    // too, the way the cone stops it.
    World w;
    const BodyId upper = rodBetween(w, Vec3{0, 0, 0}, Vec3{1.5, 0, 0}, 15.0);
    const BodyId lower = rodBetween(w, Vec3{1.5, 0, 0}, Vec3{3, 0, 0}, 10.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, upper, Vec3{0, 0, 0}, kZ));
    JointDesc j = jointAt(w, JointType::Ball, upper, lower, Vec3{1.5, 0, 0}, Vec3{1, 0, 0});
    j.coneAngle = 0.26;
    w.addJoint(j);
    w.setVelocity(lower, Vec3{0, 2.0, 0}, Vec3{0, 0, 2.5});
    double widest = 0.0;
    for (int s = 0; s < 240; ++s) {
        w.step(kDt);
        const Vec3 xA = w.orientation(upper) * j.frameA.rotation * Vec3{1, 0, 0};
        const Vec3 xB = w.orientation(lower) * j.frameB.rotation * Vec3{1, 0, 0};
        widest = std::max(widest, between(xA, xB));
    }
    CHECK_MESSAGE(widest < 0.26 + 0.01, "the lower arm never swings more than 0.01 rad past its cone, its link turning");
    CHECK_MESSAGE(widest > 0.25, "and reaches it");
}

TEST_CASE("a cone stops a swing inside the machine: the pair keeps its spin") {
    // Two rods afloat on a ball, one spun off the other into the cone: the stop is between them, so what the pair
    // turns with as a whole is what it turned with before.
    WorldSettings settings;
    settings.gravity = Vec3{};
    World w(settings);
    const BodyId a = rodBetween(w, Vec3{0, 0, 0}, Vec3{1.5, 0, 0}, 15.0);
    const BodyId b = rodBetween(w, Vec3{1.5, 0, 0}, Vec3{3, 0, 0}, 10.0);
    JointDesc j = jointAt(w, JointType::Ball, a, b, Vec3{1.5, 0, 0}, Vec3{1, 0, 0});
    j.coneAngle = 0.26;
    w.addJoint(j);
    w.setVelocity(b, Vec3{0, 1.5, 0}, Vec3{0, 0, 2.0});
    const auto spin = [&] {
        const double mass[2] = {15.0, 10.0};
        const BodyId bodies[2] = {a, b};
        Vec3 c{}, p{};
        for (int i = 0; i < 2; ++i) {
            c += mass[i] * w.position(bodies[i]);
            p += mass[i] * w.velocity(bodies[i]);
        }
        c = c / 25.0;
        const Vec3 v = p / 25.0;
        Vec3 L{};
        for (int i = 0; i < 2; ++i) {
            const Mat3 R = toMat3(w.orientation(bodies[i]));
            const Mat3 I = R * boxInertia(mass[i], Vec3{0.75, 0.1, 0.1}) * transpose(R);
            L += I * w.angularVelocity(bodies[i]) + mass[i] * cross(w.position(bodies[i]) - c, w.velocity(bodies[i]) - v);
        }
        return L;
    };
    const Vec3 L0 = spin();
    double widest = 0.0;
    for (int s = 0; s < 120; ++s) {
        w.step(kDt);
        const Vec3 xA = w.orientation(a) * j.frameA.rotation * Vec3{1, 0, 0};
        const Vec3 xB = w.orientation(b) * j.frameB.rotation * Vec3{1, 0, 0};
        widest = std::max(widest, between(xA, xB));
    }
    CHECK_MESSAGE(widest > 0.25, "the swing reaches the cone");
    CHECK_MESSAGE(length(spin() - L0) < 1e-3 * length(L0), "and the pair's spin is what it was, to 0.1%");
}

TEST_CASE("a window on a loop's cut hinge holds the loop") {
    World w;
    // A parallelogram on two ground pivots, leaning over and falling further under gravity: every hinge turns alike.
    const Vec3 top{std::sin(0.35), std::cos(0.35), 0.0}, side{2, 0, 0};
    const BodyId crank = rodBetween(w, Vec3{0, 0, 0}, top, 5.0);
    const BodyId coupler = rodBetween(w, top, top + side, 5.0);
    const BodyId rocker = rodBetween(w, side, top + side, 5.0);
    const JointId joints[4] = {
        w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, crank, Vec3{0, 0, 0}, kZ)),
        w.addJoint(jointAt(w, JointType::Hinge, crank, coupler, top, kZ)),
        w.addJoint(jointAt(w, JointType::Hinge, coupler, rocker, top + side, kZ)),
        w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, rocker, side, kZ)),
    };
    w.step(kDt);
    JointId cut;
    for (const JointId id : joints)
        if (!w.jointInfo(id).tree) cut = id;
    REQUIRE_MESSAGE(cut.valid(), "the loop is closed at one of its hinges");
    JointDesc d = w.jointDesc(cut);
    d.lowerAngle = -0.3;
    d.upperAngle = 0.3;
    w.setJoint(cut, d);
    double widest = 0.0;
    for (int s = 0; s < 180; ++s) {
        w.step(kDt);
        widest = std::max(widest, std::fabs(w.jointPosition(cut)));
    }
    CHECK_MESSAGE(widest < 0.3 + 0.01, "the cut hinge never turns more than 0.01 rad past its window");
    CHECK_MESSAGE(widest > 0.29, "and the loop falls onto it");
}

TEST_CASE("stops that are not a window or a cone are refused") {
    World w;
    const BodyId rod = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 20.0);
    JointDesc j = jointAt(w, JointType::Hinge, BodyId{}, rod, Vec3{0, 0, 0}, kZ);
    j.lowerAngle = 0.5;
    j.upperAngle = -0.5;
    CHECK_MESSAGE(!w.addJoint(j).valid(), "a window that ends before it starts");
    j.lowerAngle = std::nan("");
    j.upperAngle = 0.5;
    CHECK_MESSAGE(!w.addJoint(j).valid(), "a window that isn't a number");
    j.lowerAngle = -0.5;
    j.type = JointType::Ball;
    j.coneAngle = -0.1;
    CHECK_MESSAGE(!w.addJoint(j).valid(), "a cone narrower than nothing");
    j.coneAngle = 0.4;
    const JointId ok = w.addJoint(j);
    CHECK_MESSAGE(ok.valid(), "a cone and a window together on a ball: the window does nothing, and is allowed");
    j.coneAngle = std::nan("");
    CHECK_MESSAGE(!w.setJoint(ok, j), "a cone that isn't a number, set later");
    j.coneAngle = 0.4;
    j.coneAxis = Vec3{};
    CHECK_MESSAGE(!w.setJoint(ok, j), "a cone with no middle to open round");
}

namespace {

    // The arm afloat, no gravity, spun about its hinge toward the upper end of its window at `spin` rad/s: how fast it
    // turns once it has struck, the first step its turn is no longer the way it came, or where it is after 4 s.
    double struckArm(double bounce, double spin) {
        WorldSettings settings;
        settings.gravity = Vec3{};
        World w(settings);
        const Arm a = armOn(w, JointType::Hinge, -0.5, 0.8);
        JointDesc d = w.jointDesc(a.joint);
        d.bounce = bounce;
        w.setJoint(a.joint, d);
        w.setVelocity(a.rod, Vec3{0, spin, 0}, Vec3{0, 0, spin});
        for (int s = 0; s < 240; ++s) {
            w.step(kDt);
            const double rate = w.angularVelocity(a.rod).z;
            if (rate < spin * 0.5) return rate;
        }
        return w.angularVelocity(a.rod).z;
    }

} // namespace

TEST_CASE("a hinge struck hard against the end of its window comes back off it by its bounce") {
    CHECK_MESSAGE(std::fabs(struckArm(0.0, 3.0)) < 1e-3, "no bounce: it stops dead on its stop, as ever");
    CHECK_MESSAGE(std::fabs(struckArm(0.6, 3.0) + 0.6 * 3.0) < 0.03,
                  "bounce 0.6: it leaves the stop at 0.6 of the 3 rad/s it struck at");
    CHECK_MESSAGE(std::fabs(struckArm(1.0, 3.0) + 3.0) < 0.05, "bounce 1: it leaves at the speed it struck at");
}

TEST_CASE("a stop struck slower than kStopBounceSpeed holds without a bounce") {
    WorldSettings settings;
    settings.gravity = Vec3{};
    World w(settings);
    const Arm a = armOn(w, JointType::Hinge, -0.5, 0.8);
    JointDesc d = w.jointDesc(a.joint);
    d.bounce = 1.0;
    w.setJoint(a.joint, d);
    const double spin = 0.6 * kStopBounceSpeed;
    w.setVelocity(a.rod, Vec3{0, spin, 0}, Vec3{0, 0, spin});
    run(w, 3.5);
    CHECK_MESSAGE(std::fabs(w.jointPosition(a.joint) - 0.8) < 2e-3, "it reaches its stop");
    CHECK_MESSAGE(std::fabs(w.angularVelocity(a.rod).z) < 1e-3, "and stays on it, even at a full bounce");
}

TEST_CASE("a stop's bounce is between the two: a free pair keeps its spin, and turns apart by it") {
    WorldSettings settings;
    settings.gravity = Vec3{};
    World w(settings);
    const BodyId a = rodBetween(w, Vec3{0, 0, 0}, Vec3{1.5, 0, 0}, 15.0);
    const BodyId b = rodBetween(w, Vec3{1.5, 0, 0}, Vec3{3, 0, 0}, 10.0);
    JointDesc j = jointAt(w, JointType::Hinge, a, b, Vec3{1.5, 0, 0}, kZ);
    j.lowerAngle = -0.4;
    j.upperAngle = 0.4;
    j.bounce = 0.6;
    const JointId joint = w.addJoint(j);
    w.setVelocity(b, Vec3{0, 0.75 * 2.5, 0}, Vec3{0, 0, 2.5});
    const auto spinOf = [&] {
        const double mass[2] = {15.0, 10.0};
        const BodyId bodies[2] = {a, b};
        Vec3 c{}, p{};
        for (int i = 0; i < 2; ++i) {
            c += mass[i] * w.position(bodies[i]);
            p += mass[i] * w.velocity(bodies[i]);
        }
        c = c / 25.0;
        const Vec3 v = p / 25.0;
        Vec3 L{};
        for (int i = 0; i < 2; ++i) {
            const Mat3 R = toMat3(w.orientation(bodies[i]));
            const Mat3 I = R * boxInertia(mass[i], Vec3{0.75, 0.1, 0.1}) * transpose(R);
            L += I * w.angularVelocity(bodies[i]) + mass[i] * cross(w.position(bodies[i]) - c, w.velocity(bodies[i]) - v);
        }
        return L;
    };
    const Vec3 L0 = spinOf();
    double before = 0.0, after = 0.0, last = w.jointPosition(joint);
    bool struck = false;
    for (int s = 0; s < 120; ++s) {
        w.step(kDt);
        const double now = w.jointPosition(joint), rate = (now - last) / kDt;
        if (!struck && rate > 0.0) before = rate;
        if (!struck && rate < 0.0) struck = true;
        if (struck && after == 0.0 && rate < 0.0 && s > 0) after = w.angularVelocity(b).z - w.angularVelocity(a).z;
        last = now;
    }
    CHECK_MESSAGE(struck, "the pair turns into the stop and back off it");
    CHECK_MESSAGE(length(spinOf() - L0) < 1e-3 * length(L0), "and what it turns with as a whole is what it was, to 0.1%");
    CHECK_MESSAGE(after < 0.0, "turning apart after it");
    CHECK_MESSAGE(std::fabs(-after / before - 0.6) < 0.06, "at 0.6 of what it turned together at");
}

TEST_CASE("a ball's cone and a slider's course bounce too") {
    WorldSettings settings;
    settings.gravity = Vec3{};
    // The cone: a rod on a ball to the world, swung into a 0.35 rad cone at 2 rad/s.
    {
        World w(settings);
        const Pendulum p = pendulum(w, 0.35, 2.0);
        JointDesc d = w.jointDesc(p.ball);
        d.bounce = 0.6;
        w.setJoint(p.ball, d);
        // Stepped until it swings back from the cone, and a few steps on.
        double widest = 0.0, back = 0.0;
        int steps = 0;
        for (; steps < 120 && !(widest > 0.3 && swing(w, p) < widest - 0.02); ++steps) {
            w.step(kDt);
            widest = std::max(widest, swing(w, p));
        }
        back = length(w.angularVelocity(p.rod));
        CHECK_MESSAGE(widest < 0.35 + 0.01, "the cone still holds it");
        CHECK_MESSAGE(widest > 0.34, "it swings into it");
        CHECK_MESSAGE(std::fabs(back - 0.6 * 2.0) < 0.06, "and is sent back at 0.6 of the 2 rad/s it swung in at");
    }
    // The course: a block on a slider to the world, slid at 2 m/s into the end of its 0.5 m course.
    {
        World w(settings);
        BodyDesc bd;
        bd.mass = 10.0;
        bd.inertia = boxInertia(10.0, Vec3{0.2, 0.2, 0.2});
        const BodyId block = w.addBody(bd);
        JointDesc j = jointAt(w, JointType::Slider, BodyId{}, block, Vec3{0, 0, 0}, Vec3{1, 0, 0});
        j.min = -0.5;
        j.max = 0.5;
        j.bounce = 0.6;
        w.addJoint(j);
        w.setVelocity(block, Vec3{2.0, 0, 0}, Vec3{});
        run(w, 0.5);
        CHECK_MESSAGE(std::fabs(w.velocity(block).x + 0.6 * 2.0) < 0.03, "it comes back off its end at 0.6 of 2 m/s");
        CHECK_MESSAGE(w.position(block).x < 0.5 + 0.01, "having never gone past it");
    }
}

TEST_CASE("a bounce outside 0..1 is refused") {
    World w;
    const BodyId rod = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 20.0);
    JointDesc j = jointAt(w, JointType::Hinge, BodyId{}, rod, Vec3{0, 0, 0}, kZ);
    j.lowerAngle = -0.5;
    j.upperAngle = 0.5;
    j.bounce = -0.1;
    CHECK_MESSAGE(!w.addJoint(j).valid(), "less than none");
    j.bounce = 1.5;
    CHECK_MESSAGE(!w.addJoint(j).valid(), "more than all");
    j.bounce = std::nan("");
    CHECK_MESSAGE(!w.addJoint(j).valid(), "not a number");
    j.bounce = 0.6;
    CHECK_MESSAGE(w.addJoint(j).valid(), "0.6 is a bounce");
}
