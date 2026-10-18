#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "hf/World.hpp"

using namespace hf;

namespace {

    constexpr double kPi = std::numbers::pi;
    constexpr double kDt = 1.0 / 60.0;

    Mat3 boxInertia(double m, const Vec3 &half) {
        const double x = 4 * half.x * half.x, y = 4 * half.y * half.y, z = 4 * half.z * half.z;
        return Mat3(m / 12.0 * (y + z), 0, 0, 0, m / 12.0 * (x + z), 0, 0, 0, m / 12.0 * (x + y));
    }

    BodyId addBox(World &w, const Vec3 &position, double mass, const Vec3 &half) {
        BodyDesc d;
        d.position = position;
        d.mass = mass;
        d.inertia = boxInertia(mass, half);
        return w.addBody(d);
    }

    BodyId addGround(World &w) {
        BodyDesc d;
        d.isStatic = true;
        return w.addBody(d);
    }

    // The host's collision detection for these tests: a box's corners against the plane y = 0.
    void boxOnGround(World &w, BodyId ground, BodyId box, const Vec3 &half, double friction = 0.6,
                     double restitution = 0.0) {
        const Vec3 c = w.position(box);
        const Quat q = w.orientation(box);
        std::uint64_t feature = 0;
        for (double sx : {-1.0, 1.0})
            for (double sy : {-1.0, 1.0})
                for (double sz : {-1.0, 1.0}) {
                    const Vec3 p = c + q * Vec3{sx * half.x, sy * half.y, sz * half.z};
                    ++feature;
                    if (p.y > 0.01) continue;
                    Contact k;
                    k.a = ground;
                    k.b = box;
                    k.point = p;
                    k.normal = Vec3{0, 1, 0};
                    k.depth = -p.y;
                    k.feature = feature;
                    k.friction = friction;
                    k.restitution = restitution;
                    w.addContact(k);
                }
    }

    void run(World &w, double seconds) {
        for (int s = 0; s < int(std::lround(seconds / kDt)); ++s) w.step(kDt);
    }

} // namespace

TEST_CASE("free fall is semi-implicit Euler, exactly") {
    World w;
    const BodyId b = addBox(w, Vec3{0, 10, 0}, 1.0, Vec3{0.5, 0.5, 0.5});
    run(w, 1.0);
    // 960 sub-steps of h = 1/960: v = -g, y = y0 - g h² N(N + 1) / 2.
    const double N = 960.0, h = 1.0 / N;
    CHECK_MESSAGE(std::fabs(w.velocity(b).y + 9.81) < 1e-9, "v = -g t");
    CHECK_MESSAGE(std::fabs(w.position(b).y - (10.0 - 9.81 * h * h * N * (N + 1) / 2)) < 1e-9, "y follows the sum");
}

namespace {

    // A hinge about world Z: the joint frame's X must point along Z.
    const Quat kXtoZ = angleAxis(-kPi / 2.0, Vec3{0, 1, 0});
    Quat aboutZ(double angle) { return angleAxis(angle, Vec3{0, 0, 1}); }

} // namespace

TEST_CASE("a pendulum on a hinge to the world swings with its period") {
    World w;
    const double L = 1.0, th0 = 0.1, m = 1.0;
    const Vec3 pivot{0, 2, 0}, half{0.05, 0.05, 0.05};
    BodyDesc d;
    d.position = pivot + L * Vec3{std::sin(th0), -std::cos(th0), 0};
    d.orientation = aboutZ(th0);
    d.mass = m;
    d.inertia = boxInertia(m, half);
    const BodyId bob = w.addBody(d);
    JointDesc j;
    j.type = JointType::Hinge;
    j.a = bob;
    j.frameA = {Vec3{0, L, 0}, kXtoZ};
    j.frameB = {pivot, kXtoZ};
    const JointId hinge = w.addJoint(j);

    // Small swings: T = 2π √(I_pivot / (m g L)) (1 + θ0²/16), I_pivot = m (0.1² + 0.1²) / 12 + m L².
    const double Ip = m * 0.02 / 12.0 + m * L * L;
    const double T = 2 * kPi * std::sqrt(Ip / (m * 9.81 * L)) * (1 + th0 * th0 / 16);

    std::vector<double> crossings;
    double t = 0.0, xPrev = w.position(bob).x, worstArm = 0.0;
    while (crossings.size() < 3 && t < 10.0) {
        w.step(kDt);
        t += kDt;
        const double x = w.position(bob).x;
        if ((xPrev > 0) != (x > 0)) crossings.push_back(t - kDt + kDt * xPrev / (xPrev - x));
        xPrev = x;
        worstArm = std::max(worstArm, std::fabs(length(w.position(bob) - pivot) - L));
    }
    REQUIRE(crossings.size() == 3);
    CHECK_MESSAGE(std::fabs((crossings[2] - crossings[0]) - T) < 5e-3 * T,
                  "the period matches the small-swing formula");
    CHECK_MESSAGE(worstArm < 1e-9, "a tree hinge holds exactly");
    CHECK_MESSAGE(w.jointInfo(hinge).machine == 0, "it is a one-link machine");
}

TEST_CASE("a four-bar built from bodies and joints") {
    // Ground pivots at (0, 0) and (2.4, 0); crank 1 m, coupler 2.6 m, rocker 2.2 m: a Grashof crank-rocker.
    World w;
    const double lc = 1.0, lk = 2.6, lr = 2.2;
    const Vec3 B{2.4, 0, 0};
    const double th1 = 80.0 * kPi / 180.0;
    const Vec3 P{lc * std::cos(th1), lc * std::sin(th1), 0};
    const double dPB = length(B - P), alpha = std::atan2(B.y - P.y, B.x - P.x);
    const double th2 = alpha + std::acos((lk * lk + dPB * dPB - lr * lr) / (2 * lk * dPB));
    const Vec3 Q = P + lk * Vec3{std::cos(th2), std::sin(th2), 0};
    const double th3 = std::atan2(Q.y - B.y, Q.x - B.x);

    const auto rod = [&](const Vec3 &centre, double angle, double mass, double len) {
        BodyDesc d;
        d.position = centre;
        d.orientation = aboutZ(angle);
        d.mass = mass;
        d.inertia = Mat3(0.1, 0, 0, 0, mass * len * len / 12, 0, 0, 0, mass * len * len / 12);
        return w.addBody(d);
    };
    const BodyId crank = rod(0.5 * P, th1, 2.0, lc);
    const BodyId coupler = rod(P + 0.5 * lk * Vec3{std::cos(th2), std::sin(th2), 0}, th2, 3.0, lk);
    const BodyId rocker = rod(B + 0.5 * lr * Vec3{std::cos(th3), std::sin(th3), 0}, th3, 1.5, lr);

    const auto hinge = [&](JointType type, BodyId a, const Vec3 &onA, BodyId b, const Vec3 &onB) {
        JointDesc j;
        j.type = type;
        j.a = a;
        j.b = b;
        j.frameA = {onA, kXtoZ};
        j.frameB = {onB, kXtoZ};
        j.target = 3.0;
        j.maxForce = 400.0;
        return w.addJoint(j);
    };
    const JointId motor = hinge(JointType::Motor, crank, Vec3{-0.5 * lc, 0, 0}, BodyId{}, Vec3{});
    hinge(JointType::Hinge, coupler, Vec3{-0.5 * lk, 0, 0}, crank, Vec3{0.5 * lc, 0, 0});
    hinge(JointType::Hinge, rocker, Vec3{-0.5 * lr, 0, 0}, BodyId{}, B);
    hinge(JointType::Hinge, coupler, Vec3{0.5 * lk, 0, 0}, rocker, Vec3{0.5 * lr, 0, 0});

    const auto gap = [&] {
        return length((w.position(coupler) + w.orientation(coupler) * Vec3{0.5 * lk, 0, 0}) -
                      (w.position(rocker) + w.orientation(rocker) * Vec3{0.5 * lr, 0, 0}));
    };
    CHECK_MESSAGE(gap() < 1e-9, "it assembles closed");
    double worstGap = 0.0;
    for (int s = 0; s < 120; ++s) {
        w.step(kDt);
        worstGap = std::max(worstGap, gap());
    }

    int cuts = 0;
    for (std::uint32_t k = 0; k < 4; ++k) cuts += w.jointInfo({k, 0}).tree ? 0 : 1;
    CHECK_MESSAGE(cuts == 1, "one loop, one cut");
    CHECK_MESSAGE(w.jointInfo(motor).tree, "the motor stays in the tree");
    CHECK_MESSAGE(std::fabs(std::fabs(w.angularVelocity(crank).z) - 3.0) < 0.03,
                  "the motor turns the crank at 3 rad/s");
    CHECK_MESSAGE(worstGap < 1e-3, "and the loop holds");
}

TEST_CASE("a piston drives a body to its target") {
    World w;
    const BodyId rod = addBox(w, Vec3{}, 10.0, Vec3{0.2, 0.1, 0.1});
    JointDesc j;
    j.type = JointType::Piston;
    j.b = rod; // a is the world: the cylinder
    j.target = 0.5;
    j.speed = 1.0;
    j.maxForce = 1e4;
    const JointId piston = w.addJoint(j);
    run(w, 1.5);
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - 0.5) < 1e-3, "it reaches 0.5 m");
    CHECK_MESSAGE(std::fabs(w.position(rod).y) < 1e-6, "and holds the rod on its axis against gravity");
}

TEST_CASE("a box comes to rest on the ground") {
    World w;
    const BodyId ground = addGround(w);
    const Vec3 half{0.5, 0.5, 0.5};
    const BodyId box = addBox(w, Vec3{0, 0.6, 0}, 10.0, half);
    for (int s = 0; s < 120; ++s) {
        w.clearContacts();
        boxOnGround(w, ground, box, half);
        w.step(kDt);
    }
    CHECK_MESSAGE(std::fabs(w.position(box).y - 0.5) < 5e-3, "it rests on its face, within the slop");
    CHECK_MESSAGE(length(w.velocity(box)) < 1e-3, "and is still");
    CHECK_MESSAGE(length(w.angularVelocity(box)) < 1e-3, "and doesn't turn");
}

namespace {

    // A 20° slope, by tilting gravity instead of the ground: down the slope is +x.
    double slideAfter(double mu, double seconds) {
        const double th = 20.0 * kPi / 180.0;
        WorldSettings settings;
        settings.gravity = 9.81 * Vec3{std::sin(th), -std::cos(th), 0};
        World w(settings);
        const BodyId ground = addGround(w);
        const Vec3 half{0.5, 0.5, 0.5};
        const BodyId box = addBox(w, Vec3{0, 0.5, 0}, 10.0, half);
        for (int s = 0; s < int(std::lround(seconds / kDt)); ++s) {
            w.clearContacts();
            boxOnGround(w, ground, box, half, mu);
            w.step(kDt);
        }
        return w.velocity(box).x;
    }

} // namespace

TEST_CASE("friction on a slope") {
    // tan 20° = 0.364. Above it the box sticks; below, it slides at g (sin θ - μ cos θ).
    CHECK_MESSAGE(std::fabs(slideAfter(0.6, 2.0)) < 1e-2, "μ = 0.6 holds it");
    const double th = 20.0 * kPi / 180.0;
    const double expected = 9.81 * (std::sin(th) - 0.2 * std::cos(th)) * 1.0;
    CHECK_MESSAGE(std::fabs(slideAfter(0.2, 1.0) - expected) < 0.05 * expected, "μ = 0.2 slides at the right rate");
}

TEST_CASE("a dropped box bounces to e² of its drop") {
    World w;
    const BodyId ground = addGround(w);
    const Vec3 half{0.5, 0.5, 0.5};
    const BodyId box = addBox(w, Vec3{0, 1.5, 0}, 10.0, half);
    bool landed = false;
    double peak = 0.0;
    for (int s = 0; s < 180; ++s) {
        w.clearContacts();
        boxOnGround(w, ground, box, half, 0.6, 0.5);
        w.step(kDt);
        if (w.velocity(box).y > 0.0) landed = true;
        if (landed) peak = std::max(peak, w.position(box).y - 0.5);
    }
    CHECK_MESSAGE(landed, "it lands");
    CHECK_MESSAGE(std::fabs(peak - 0.25) < 0.08, "a 1 m drop at e = 0.5 rebounds about 0.25 m");
}

TEST_CASE("a fast-rolling ball keeps its height, whichever side of the contact it is on") {
    // Over a step the ball turns half a radian: a contact point carried round with it would lift off the ground.
    for (const bool ballIsA : {false, true}) {
        World w;
        const BodyId ground = addGround(w);
        const double r = 0.5, m = 10.0, spin = 30.0;
        BodyDesc d;
        d.position = Vec3{0, r, 0};
        d.mass = m;
        d.inertia = Mat3(0.4 * m * r * r);
        d.velocity = Vec3{spin * r, 0, 0};
        d.angularVelocity = Vec3{0, 0, -spin};
        const BodyId ball = w.addBody(d);
        double worstSink = 0.0, worstVy = 0.0;
        for (int s = 0; s < 120; ++s) {
            const Vec3 c = w.position(ball);
            Contact k;
            k.a = ballIsA ? ball : ground;
            k.b = ballIsA ? ground : ball;
            k.point = Vec3{c.x, 0, c.z};
            k.normal = ballIsA ? Vec3{0, -1, 0} : Vec3{0, 1, 0};
            k.depth = r - c.y;
            k.feature = 1;
            w.clearContacts();
            w.addContact(k);
            w.step(kDt);
            worstSink = std::max(worstSink, r - w.position(ball).y);
            worstVy = std::max(worstVy, std::fabs(w.velocity(ball).y));
        }
        CAPTURE(ballIsA);
        CHECK_MESSAGE(worstSink < 2e-3, "it rolls on the ground, not into it");
        CHECK_MESSAGE(worstVy < 0.02, "no vertical speed from the spin");
        CHECK_MESSAGE(std::fabs(w.velocity(ball).x - spin * r) < 0.05 * spin * r, "it keeps rolling");
    }
}

TEST_CASE("a hinged bar falls onto the ground and rests there") {
    World w;
    const BodyId ground = addGround(w);
    const Vec3 half{1.0, 0.05, 0.05}, pivot{0, 0.5, 0};
    const double up = 30.0 * kPi / 180.0;
    BodyDesc d;
    d.position = pivot + Vec3{std::cos(up), std::sin(up), 0};
    d.orientation = aboutZ(up);
    d.mass = 5.0;
    d.inertia = boxInertia(5.0, half);
    const BodyId bar = w.addBody(d);
    JointDesc j;
    j.type = JointType::Hinge;
    j.a = bar;
    j.frameA = {Vec3{-1.0, 0, 0}, kXtoZ};
    j.frameB = {pivot, kXtoZ};
    w.addJoint(j);

    double worstArm = 0.0;
    for (int s = 0; s < 180; ++s) {
        w.clearContacts();
        boxOnGround(w, ground, bar, half);
        w.step(kDt);
        const Vec3 end = w.position(bar) + w.orientation(bar) * Vec3{-1.0, 0, 0};
        worstArm = std::max(worstArm, length(end - pivot));
    }
    double lowest = 1e9;
    for (double sx : {-1.0, 1.0})
        for (double sy : {-1.0, 1.0})
            lowest = std::min(lowest, (w.position(bar) + w.orientation(bar) * Vec3{sx * 1.0, sy * 0.05, 0}).y);
    CHECK_MESSAGE((lowest > -0.01 && lowest < 0.005), "its tip rests on the ground");
    CHECK_MESSAGE(length(w.angularVelocity(bar)) < 1e-2, "and stays there");
    CHECK_MESSAGE(worstArm < 1e-9, "the hinge held exactly through the impact");
}

TEST_CASE("welded bodies move as one") {
    World w;
    const BodyId a = addBox(w, Vec3{0, 0, 0}, 2.0, Vec3{0.5, 0.5, 0.5});
    const BodyId b = addBox(w, Vec3{1, 0, 0}, 1.0, Vec3{0.5, 0.5, 0.5});
    JointDesc j;
    j.type = JointType::Fixed;
    j.a = a;
    j.b = b;
    const JointId weld = w.addJoint(j);
    for (int s = 0; s < 60; ++s) {
        w.applyTorque(a, Vec3{1, 2, 3});
        w.step(kDt);
    }
    const Vec3 rel = conjugate(w.orientation(a)) * (w.position(b) - w.position(a));
    CHECK_MESSAGE(length(rel - Vec3{1, 0, 0}) < 1e-12, "b keeps its place on a");
    CHECK_MESSAGE(length(w.angularVelocity(a) - w.angularVelocity(b)) < 1e-12, "and its spin");
    CHECK_MESSAGE(w.jointInfo(weld).machine == -1, "a weld isn't a machine joint");
}

TEST_CASE("a body hangs on a spring from the world") {
    World w;
    const BodyId bob = addBox(w, Vec3{0, 2, 0}, 3.0, Vec3{0.1, 0.1, 0.1});
    JointDesc j;
    j.type = JointType::Spring;
    j.a = bob;
    j.frameB = {Vec3{0, 3, 0}, Quat{}};
    j.restLength = 1.0;
    j.hertz = 2.0;
    j.damping = 1.0;
    const JointId spring = w.addJoint(j);
    run(w, 5.0);
    // k = m (2π f)², so the static extension is g / (2π f)², whatever the mass.
    const double stretch = 9.81 / std::pow(2 * kPi * 2.0, 2);
    CHECK_MESSAGE(std::fabs((2.0 - w.position(bob).y) - stretch) < 2e-3, "it settles at m g / k below its rest length");
    CHECK_MESSAGE(w.jointInfo(spring).machine == -1, "a lone body's spring is solved with the contacts");
}

TEST_CASE("handles of removed bodies stay invalid") {
    World w;
    const BodyId a = addBox(w, Vec3{}, 1.0, Vec3{0.5, 0.5, 0.5});
    w.removeBody(a);
    const BodyId b = addBox(w, Vec3{}, 1.0, Vec3{0.5, 0.5, 0.5});
    CHECK_MESSAGE(!w.contains(a), "the old handle is dead");
    CHECK_MESSAGE(w.contains(b), "the new one is live");
    CHECK_MESSAGE((a.index == b.index && a.generation != b.generation), "in the same slot, a new generation");
}

namespace {

    // The rotation taking +X onto axis.
    Quat xTo(const Vec3 &axis) {
        const Vec3 a = normalize(axis);
        return normalize(Quat{1.0 + a.x, cross(Vec3{1, 0, 0}, a)});
    }

    // A joint at a world point about a world axis, with the bodies where they are: their frames coincide.
    JointDesc jointAt(World &w, JointType type, BodyId a, BodyId b, const Vec3 &point, const Vec3 &axis) {
        const auto local = [&](BodyId body) -> Frame {
            if (!body.valid()) return {point, xTo(axis)};
            const Quat q = w.orientation(body);
            return {conjugate(q) * (point - w.position(body)), normalize(conjugate(q) * xTo(axis))};
        };
        JointDesc j;
        j.type = type;
        j.a = a;
        j.b = b;
        j.frameA = local(a);
        j.frameB = local(b);
        return j;
    }

    BodyId rodBetween(World &w, const Vec3 &p, const Vec3 &q, double mass) {
        const Vec3 d = q - p;
        const double half = 0.5 * length(d);
        BodyDesc b;
        b.position = 0.5 * (p + q);
        b.orientation = xTo(d);
        b.mass = mass;
        b.inertia = boxInertia(mass, Vec3{half, 0.1, 0.1});
        return w.addBody(b);
    }

} // namespace

TEST_CASE("a cylinder in a loop lifts a heavy boom to a jumped target at its speed") {
    World w;
    const Vec3 pivot{0.5, 1.6, 0}, tip{3.5, 3.6, 0}, base{1.0, 1.3, 0}, z{0, 0, 1};
    const Vec3 onBoom = pivot + 1.5 * normalize(tip - pivot), mid = 0.5 * (base + onBoom);
    const Vec3 u = normalize(onBoom - base);
    const BodyId boom = rodBetween(w, pivot, tip, 600.0);
    const BodyId barrel = rodBetween(w, base, mid + 0.2 * u, 40.0);
    const BodyId rod = rodBetween(w, mid - 0.2 * u, onBoom, 20.0);
    w.addJoint(jointAt(w, JointType::Hinge, boom, BodyId{}, pivot, z));
    w.addJoint(jointAt(w, JointType::Hinge, barrel, BodyId{}, base, z));
    w.addJoint(jointAt(w, JointType::Hinge, rod, boom, onBoom, z));
    JointDesc p = jointAt(w, JointType::Piston, barrel, rod, mid, u);
    p.maxForce = 2e5;
    p.speed = 0.5;
    const JointId piston = w.addJoint(p);

    run(w, 0.5);
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston)) < 1e-3, "it holds the boom where it was built");
    w.setTarget(piston, 0.2);
    run(w, 0.25);
    // 0.2 m at 0.5 m/s takes 0.4 s: after 0.25 s, about 0.125 m.
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - 0.125) < 0.02, "it moves at its speed, not at the servo's creep");
    run(w, 0.75);
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - 0.2) < 1e-3, "and arrives");
}
