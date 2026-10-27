#include <algorithm>
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

    constexpr double kG = 9.81;
    constexpr double kPi = std::numbers::pi;
    const Vec3 kZ{0, 0, 1};

    WorldSettings weightless() {
        WorldSettings s;
        s.gravity = Vec3{};
        return s;
    }

    BodyId ball(World &w, const Vec3 &at, double mass, double radius = 0.1) {
        BodyDesc d;
        d.position = at;
        d.mass = mass;
        d.inertia = Mat3(0.4 * mass * radius * radius);
        return w.addBody(d);
    }

    BodyId box(World &w, const Vec3 &at, double mass, const Vec3 &half) {
        BodyDesc d;
        d.position = at;
        d.mass = mass;
        d.inertia = boxInertia(mass, half);
        return w.addBody(d);
    }

    // A route point on a body, at a world point where the body is now.
    RopePoint on(World &w, BodyId b, const Vec3 &at) {
        RopePoint p;
        p.body = b;
        p.position = conjugate(w.orientation(b)) * (at - w.position(b));
        return p;
    }

    RopePoint fixedAt(const Vec3 &at) {
        RopePoint p;
        p.position = at;
        return p;
    }

    RopePoint sheaveAt(const Vec3 &at, double radius, const Vec3 &axis = kZ, int side = 0) {
        RopePoint p = fixedAt(at);
        p.radius = radius;
        p.axis = axis;
        p.side = side;
        return p;
    }

    RopeDesc rigid(double length) {
        RopeDesc d;
        d.length = length;
        return d;
    }

    RopeId rope(World &w, const RopeDesc &d, std::vector<RopePoint> route) { return w.addRope(d, route); }

} // namespace

TEST_CASE("a rope holds a hanging load at its length, pulling with its weight") {
    World w;
    const BodyId load = ball(w, Vec3{0, -2, 0}, 50.0);
    const RopeId r = rope(w, rigid(2.0), {fixedAt(Vec3{}), on(w, load, Vec3{0, -2, 0})});
    REQUIRE(r.valid());
    run(w, 10.0);
    CHECK_MESSAGE(std::fabs(w.ropeTension(r) / (50.0 * kG) - 1.0) < 1e-6, "its tension is the load's weight");
    CHECK_MESSAGE(std::fabs(w.position(load).y + 2.0) < 1e-5, "it hangs at the rope's length, not below");
    CHECK_MESSAGE(std::fabs(w.ropeSpan(r) - 2.0) < 1e-5, "the rope doesn't stretch");
    CHECK(w.ropeInfo(r).taut);
    CHECK_MESSAGE(w.asleep(load), "still, it sleeps, and keeps its tension");
}

TEST_CASE("a rope never pushes") {
    World w(weightless());
    const BodyId a = ball(w, Vec3{0, 0, 0}, 1.0), b = ball(w, Vec3{2, 0, 0}, 1.0);
    const RopeId r = rope(w, rigid(2.0), {on(w, a, Vec3{0, 0, 0}), on(w, b, Vec3{2, 0, 0})});
    w.setVelocity(a, Vec3{1, 0, 0}, Vec3{});
    w.setVelocity(b, Vec3{-1, 0, 0}, Vec3{});
    run(w, 0.5);
    CHECK_MESSAGE(w.ropeTension(r) == 0.0, "taut at first, then slack: it pulls nothing as its ends close");
    CHECK_MESSAGE(length(w.velocity(a) - Vec3{1, 0, 0}) < 1e-12, "nor slows them");
    CHECK(length(w.velocity(b) - Vec3{-1, 0, 0}) < 1e-12);
    CHECK_MESSAGE(!w.ropeInfo(r).taut, "slack");
}

TEST_CASE("a slack rope lets its load fall freely until it is taut, then catches it dead") {
    World w;
    const BodyId load = ball(w, Vec3{0, -1, 0}, 10.0);
    const RopeId r = rope(w, rigid(2.0), {fixedAt(Vec3{}), on(w, load, Vec3{0, -1, 0})});
    // 0.4 s of free fall: well before it has fallen the metre of slack. The world integrates gravity sub-step by
    // sub-step, v first: after n of them it has fallen g h² n(n+1)/2.
    run(w, 0.4);
    const double h = kDt / w.settings().substeps;
    const double n = std::round(0.4 / kDt) * w.settings().substeps;
    CHECK_MESSAGE(std::fabs(w.position(load).y - (-1.0 - kG * h * h * n * (n + 1) / 2.0)) < 1e-9,
                  "until taut it falls as if there were no rope");
    CHECK(w.ropeTension(r) == 0.0);
    double furthest = 0.0, after = 0.0;
    for (int s = 0; s < 120; ++s) {
        w.step(kDt);
        furthest = std::max(furthest, -w.position(load).y);
        if (s > 30) after = std::max(after, std::fabs(w.velocity(load).y));
    }
    CHECK_MESSAGE(furthest < 2.0 + 1e-4, "the rope stops it at its length, with nothing past it");
    CHECK_MESSAGE(after < 1e-3, "and it stays there: no bounce back up the rope");
}

TEST_CASE("a load on a rope swings as a pendulum, keeping its energy") {
    World w;
    const double L = 2.0, a0 = 5.0 * kPi / 180.0;
    const Vec3 at{L * std::sin(a0), -L * std::cos(a0), 0};
    const BodyId bob = ball(w, at, 5.0, 0.01);
    rope(w, rigid(L), {fixedAt(Vec3{}), on(w, bob, at)});
    const auto energy = [&] {
        const Vec3 v = w.velocity(bob);
        return 0.5 * 5.0 * dot(v, v) + 5.0 * kG * w.position(bob).y;
    };
    const double e0 = energy();
    // Times the bob crosses the bottom going one way: a period apart.
    std::vector<double> crossings;
    double prev = w.position(bob).x;
    for (int s = 1; s <= 20 * 60; ++s) {
        w.step(kDt);
        const double x = w.position(bob).x;
        if (prev > 0.0 && x <= 0.0) crossings.push_back(s * kDt - kDt * x / (x - prev));
        prev = x;
    }
    REQUIRE(crossings.size() >= 6);
    const double period = (crossings.back() - crossings.front()) / double(crossings.size() - 1);
    // Small swings: 2π√(L/g), and (1 + a0²/16) for the 5° amplitude.
    const double expected = 2.0 * kPi * std::sqrt(L / kG) * (1.0 + a0 * a0 / 16.0);
    CHECK_MESSAGE(std::fabs(period / expected - 1.0) < 1e-2, "its period is a pendulum's");
    CHECK_MESSAGE(std::fabs(energy() - e0) < 1e-2 * 5.0 * kG * L * (1.0 - std::cos(a0)),
                  "over 20 s it keeps its swing's energy to a percent");
}

TEST_CASE("an Atwood machine over two pulley points accelerates as (m1 - m2) g / (m1 + m2)") {
    World w;
    const BodyId heavy = ball(w, Vec3{-0.5, -1, 0}, 3.0), light = ball(w, Vec3{0.5, -1, 0}, 1.0);
    const RopeId r = rope(w, rigid(3.0),
                          {on(w, heavy, Vec3{-0.5, -1, 0}), fixedAt(Vec3{-0.5, 0, 0}), fixedAt(Vec3{0.5, 0, 0}),
                           on(w, light, Vec3{0.5, -1, 0})});
    run(w, 0.2);
    const double v0 = w.velocity(heavy).y;
    run(w, 0.3);
    const double a = -(w.velocity(heavy).y - v0) / 0.3;
    CHECK_MESSAGE(std::fabs(a / (0.5 * kG) - 1.0) < 1e-3, "the heavy side falls at half of g");
    CHECK_MESSAGE(std::fabs(w.velocity(light).y + w.velocity(heavy).y) < 1e-6, "the light side rises as fast");
    CHECK_MESSAGE(std::fabs(w.ropeTension(r) / (1.5 * kG) - 1.0) < 1e-3, "its tension is 2 m1 m2 g / (m1 + m2)");
    CHECK_MESSAGE(std::fabs(w.ropeSpan(r) - 3.0) < 1e-5, "the rope keeps its length round both pulleys");
}

namespace {

    // A 400 kg block hung on four strands of one rope, eyes on the block and on a beam above it, the rope's free end
    // over a last pulley to a counterweight. The block's eyes are spread along z about its centre, so the four strands
    // lift it without turning it.
    struct Tackle {
        BodyId block, weight;
        RopeId rope;
    };

    Tackle tackle(World &w, double counterweight) {
        Tackle t;
        t.block = box(w, Vec3{0, -2, 0.15}, 400.0, Vec3{0.1, 0.1, 0.2});
        t.weight = ball(w, Vec3{1, -1, 0.3}, counterweight);
        const auto eye = [&](double z) { return on(w, t.block, Vec3{0, -2, z}); };
        t.rope = rope(w, rigid(8.0 + 0.3 + 1.0 + 1.0),
                      {fixedAt(Vec3{0, 0, 0}), eye(0.0), eye(0.1), fixedAt(Vec3{0, 0, 0.1}), fixedAt(Vec3{0, 0, 0.2}),
                       eye(0.2), eye(0.3), fixedAt(Vec3{0, 0, 0.3}), fixedAt(Vec3{1, 0, 0.3}),
                       on(w, t.weight, Vec3{1, -1, 0.3})});
        return t;
    }

} // namespace

TEST_CASE("four strands give four times the pull: a quarter of the load balances it") {
    World w;
    const Tackle t = tackle(w, 100.0);
    run(w, 5.0);
    CHECK_MESSAGE(length(w.velocity(t.block)) < 1e-6, "100 kg holds 400 kg still");
    CHECK(length(w.velocity(t.weight)) < 1e-6);
    CHECK_MESSAGE(std::fabs(w.ropeTension(t.rope) / (100.0 * kG) - 1.0) < 1e-6,
                  "the rope pulls with the counterweight's weight, and the block hangs on four of it");

    World moving;
    const Tackle m = tackle(moving, 120.0);
    run(moving, 0.2);
    const Vec3 b0 = moving.position(m.block), c0 = moving.position(m.weight);
    const double vb0 = moving.velocity(m.block).y;
    run(moving, 0.5);
    const double rise = moving.position(m.block).y - b0.y, fall = c0.y - moving.position(m.weight).y;
    CHECK_MESSAGE(std::fabs(fall / rise - 4.0) < 1e-3, "the block rises a quarter of what the weight falls");
    // Block up x, weight down 4x: (M + 16 m) x'' = (4 m - M) g.
    const double a = (moving.velocity(m.block).y - vb0) / 0.5;
    CHECK_MESSAGE(std::fabs(a / ((4.0 * 120.0 - 400.0) * kG / (400.0 + 16.0 * 120.0)) - 1.0) < 1e-3,
                  "at the acceleration the lever of four strands gives");
}

TEST_CASE("a rope between two machines leaves them two machines, solved in one island") {
    World w;
    const BodyId a = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 10.0);
    const BodyId b = rodBetween(w, Vec3{3, 0, 0}, Vec3{5, 0, 0}, 10.0);
    const JointId ha = w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, a, Vec3{0, 0, 0}, kZ));
    const JointId hb = w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, b, Vec3{5, 0, 0}, kZ));
    const RopeId r = rope(w, rigid(1.0), {on(w, a, Vec3{2, 0, 0}), on(w, b, Vec3{3, 0, 0})});
    // Both swing down about their hinges: the rope between their tips holds them together.
    double widest = 0.0;
    for (int s = 0; s < 180; ++s) {
        w.step(kDt);
        const Vec3 tipA = w.position(a) + w.orientation(a) * Vec3{1, 0, 0};
        const Vec3 tipB = w.position(b) - w.orientation(b) * Vec3{1, 0, 0};
        widest = std::max(widest, length(tipA - tipB));
    }
    CHECK_MESSAGE(w.jointInfo(ha).machine != w.jointInfo(hb).machine, "each arm is still its own machine");
    CHECK(w.stats().machines == 2);
    CHECK_MESSAGE(w.ropeInfo(r).machine == -1, "the rope is swept between them");
    CHECK_MESSAGE(widest < 1.0 + 1e-4, "and it holds their tips together all the same");
}

TEST_CASE("a rope inside one machine is held with its joints: an arm held level by a rope at its tip") {
    World w;
    const BodyId arm = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 10.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, arm, Vec3{0, 0, 0}, kZ));
    const RopeId r = rope(w, rigid(1.0), {fixedAt(Vec3{2, 1, 0}), on(w, arm, Vec3{2, 0, 0})});
    run(w, 3.0);
    CHECK_MESSAGE(w.ropeInfo(r).machine >= 0, "the rope is the machine's");
    CHECK_MESSAGE(std::fabs(w.ropeTension(r) / (10.0 * kG / 2.0) - 1.0) < 1e-6,
                  "it carries half the arm's weight: the hinge the other half");
    const Vec3 tip = w.position(arm) + w.orientation(arm) * Vec3{1, 0, 0};
    CHECK_MESSAGE(std::fabs(tip.y) < 1e-5, "and holds it level");

    // The same arm unhinged, on two ropes, one at each end: each takes half its weight. Swept, not held.
    World free;
    const BodyId rod = rodBetween(free, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 10.0);
    const RopeId left = rope(free, rigid(1.0), {fixedAt(Vec3{0, 1, 0}), on(free, rod, Vec3{0, 0, 0})});
    const RopeId right = rope(free, rigid(1.0), {fixedAt(Vec3{2, 1, 0}), on(free, rod, Vec3{2, 0, 0})});
    run(free, 3.0);
    CHECK(free.ropeInfo(left).machine == -1);
    CHECK_MESSAGE(std::fabs(free.ropeTension(left) / (10.0 * kG / 2.0) - 1.0) < 1e-6, "swept, the same half");
    CHECK(std::fabs(free.ropeTension(right) / (10.0 * kG / 2.0) - 1.0) < 1e-6);
}

TEST_CASE("a 20 t load on a 10 kg hook: the ropes neither stretch nor jitter") {
    World w;
    const BodyId hook = ball(w, Vec3{0, -2, 0}, 10.0);
    const BodyId load = box(w, Vec3{0, -3.5, 0}, 20000.0, Vec3{0.5, 0.5, 0.5});
    const RopeId crane = rope(w, rigid(2.0), {fixedAt(Vec3{}), on(w, hook, Vec3{0, -2, 0})});
    const RopeId sling = rope(w, rigid(1.0), {on(w, hook, Vec3{0, -2, 0}), on(w, load, Vec3{0, -3, 0})});
    double stretch = 0.0, jitter = 0.0;
    for (int s = 0; s < 600; ++s) {
        w.step(kDt);
        stretch = std::max({stretch, w.ropeSpan(crane) - 2.0, w.ropeSpan(sling) - 1.0});
        if (s > 60) jitter = std::max(jitter, length(w.velocity(hook)));
    }
    CHECK_MESSAGE(stretch < 1e-4, "neither rope stretches by a tenth of a millimetre");
    CHECK_MESSAGE(jitter < 1e-4, "and the light hook between them stays still");
    CHECK_MESSAGE(std::fabs(w.ropeTension(crane) / (20010.0 * kG) - 1.0) < 1e-4, "the crane holds hook and load");
    CHECK(std::fabs(w.ropeTension(sling) / (20000.0 * kG) - 1.0) < 1e-4);
}

TEST_CASE("four slings from a light hook share a 20 t load exactly") {
    World w;
    const BodyId hook = ball(w, Vec3{0, -2, 0}, 10.0);
    const BodyId load = box(w, Vec3{0, -3.25, 0}, 20000.0, Vec3{0.5, 0.25, 0.5});
    const RopeId crane = rope(w, rigid(2.0), {fixedAt(Vec3{}), on(w, hook, Vec3{0, -2, 0})});
    std::vector<RopeId> slings;
    const double len = std::sqrt(1.5);
    for (const Vec3 corner : {Vec3{0.5, -3, 0.5}, Vec3{-0.5, -3, 0.5}, Vec3{0.5, -3, -0.5}, Vec3{-0.5, -3, -0.5}})
        slings.push_back(rope(w, rigid(len), {on(w, hook, Vec3{0, -2, 0}), on(w, load, corner)}));
    double stretch = 0.0;
    for (int s = 0; s < 600; ++s) {
        w.step(kDt);
        for (RopeId r : slings) stretch = std::max(stretch, w.ropeSpan(r) - len);
    }
    // Each takes a quarter of the weight along its slant: 1 m down in √1.5 of sling.
    const double each = 20000.0 * kG / 4.0 * len;
    for (RopeId r : slings) CHECK_MESSAGE(std::fabs(w.ropeTension(r) / each - 1.0) < 1e-3, "a quarter each");
    CHECK_MESSAGE(stretch < 1e-4, "none stretches");
    CHECK(std::fabs(w.ropeTension(crane) / (20010.0 * kG) - 1.0) < 1e-4);
}

TEST_CASE("a stretchy rope stretches mgL/EA and bounces at its spring's frequency") {
    World w;
    RopeDesc d = rigid(2.0);
    d.stiffness = 1e5;
    d.damping = 0.0;
    const BodyId load = ball(w, Vec3{0, -2, 0}, 100.0);
    const RopeId r = rope(w, d, {fixedAt(Vec3{}), on(w, load, Vec3{0, -2, 0})});
    std::vector<double> peaks; // the bottom of each bounce
    double prevV = 0.0;
    for (int s = 1; s <= 120; ++s) {
        w.step(kDt);
        const double v = w.velocity(load).y;
        if (prevV < 0.0 && v >= 0.0) peaks.push_back(s * kDt);
        prevV = v;
    }
    REQUIRE(peaks.size() >= 4);
    const double period = (peaks.back() - peaks.front()) / double(peaks.size() - 1);
    const double expected = 2.0 * kPi / std::sqrt(1e5 / 2.0 / 100.0);
    CHECK_MESSAGE(std::fabs(period / expected - 1.0) < 2e-2, "a spring of EA / L on its load");

    d.damping = 1.0;
    w.setRope(r, d);
    run(w, 5.0);
    CHECK_MESSAGE(std::fabs((w.ropeSpan(r) - 2.0) / (100.0 * kG * 2.0 / 1e5) - 1.0) < 1e-2,
                  "damped, it settles stretched by mgL/EA");
}

TEST_CASE("ropes the world can't use are refused, and a rope goes with a body removed") {
    World w;
    const BodyId b = ball(w, Vec3{0, -1, 0}, 1.0);
    const RopePoint top = fixedAt(Vec3{}), end = on(w, b, Vec3{0, -1, 0});
    RopeDesc bad = rigid(1.0);
    bad.length = std::nan("");
    CHECK(!w.addRope(bad, std::vector<RopePoint>{top, end}).valid());
    bad = rigid(-1.0);
    CHECK(!w.addRope(bad, std::vector<RopePoint>{top, end}).valid());
    bad = rigid(1.0);
    bad.minLength = 2.0;
    bad.maxLength = 1.0;
    CHECK(!w.addRope(bad, std::vector<RopePoint>{top, end}).valid());
    CHECK(!w.addRope(rigid(1.0), std::vector<RopePoint>{top}).valid());
    RopePoint flat = sheaveAt(Vec3{}, 0.1, Vec3{});
    CHECK(!w.addRope(rigid(1.0), std::vector<RopePoint>{flat, end}).valid());
    w.step(kDt);
    CHECK(w.stats().rejected == 5);

    const RopeId r = rope(w, rigid(1.0), {top, end});
    REQUIRE(w.contains(r));
    CHECK(!w.setSpoolLever(r, std::nan("")));
    run(w, 3.0);
    CHECK_MESSAGE(w.asleep(b), "a load still on its rope sleeps");
    w.removeBody(b);
    CHECK_MESSAGE(!w.contains(r), "removed with the body it held");
}

TEST_CASE("a hook hangs under a sheave's tangent, not its centre") {
    World w;
    const double r = 0.2;
    // In from the left along the sheave's top, round a quarter turn clockwise, and down its right side.
    const BodyId hook = ball(w, Vec3{r, -2, 0}, 20.0);
    const RopeId rp = rope(w, rigid(3.0 + 2.0 + 0.5 * kPi * r),
                           {fixedAt(Vec3{-3, r, 0}), sheaveAt(Vec3{}, r), on(w, hook, Vec3{r, -2, 0})});
    run(w, 5.0);
    CHECK_MESSAGE(std::fabs(w.position(hook).x - r) < 1e-6, "it hangs straight down from where the rope leaves");
    CHECK_MESSAGE(std::fabs(w.position(hook).y + 2.0) < 1e-5, "with the quarter turn round the sheave in the rope");
    RopeBend bends[3];
    REQUIRE(w.ropePath(rp, bends) == 3);
    CHECK(length(bends[1].in - Vec3{0, r, 0}) < 1e-6);
    CHECK(length(bends[1].out - Vec3{r, 0, 0}) < 1e-6);
    CHECK(bends[1].touching);
}

TEST_CASE("a route's span lays where it runs as a step would, and follows a body moved by hand at once") {
    World w;
    const double r = 0.2;
    const BodyId hook = ball(w, Vec3{r, -2, 0}, 20.0);
    const std::vector<RopePoint> route{fixedAt(Vec3{-3, r, 0}), sheaveAt(Vec3{}, r), on(w, hook, Vec3{r, -2, 0})};
    const RopeId rp = rope(w, rigid(3.0 + 2.0 + 0.5 * kPi * r), route);
    run(w, 2.0);
    RopeBend stepped[3], laid[3];
    REQUIRE(w.ropePath(rp, stepped) == 3);
    CHECK(std::fabs(w.routeSpan(route, laid) - w.ropeSpan(rp)) < 1e-6);
    for (int k = 0; k < 3; ++k) {
        CHECK(length(laid[k].in - stepped[k].in) < 1e-6);
        CHECK(length(laid[k].out - stepped[k].out) < 1e-6);
        CHECK(laid[k].touching == stepped[k].touching);
    }
    // Put a metre to the right, nothing stepped: the rope leaves the sheave toward it, on a tangent.
    const Vec3 to = w.position(hook) + Vec3{1, 0, 0};
    REQUIRE(w.setPose(hook, to, w.orientation(hook)));
    const double span = w.routeSpan(route, laid);
    CHECK(length(laid[2].in - to) < 1e-9);
    CHECK(std::fabs(length(laid[1].out) - r) < 1e-9);
    CHECK(std::fabs(dot(laid[1].out, to - laid[1].out)) < 1e-9);
    CHECK(std::fabs(span - (3.0 + length(to - laid[1].out) + r * std::acos(dot(laid[1].in, laid[1].out) / (r * r)))) < 1e-9);
}

TEST_CASE("an Atwood machine over a sheave: half a turn of rope round it, the same motion") {
    World w;
    const double r = 0.5;
    const BodyId heavy = ball(w, Vec3{-r, -1, 0}, 3.0), light = ball(w, Vec3{r, -1, 0}, 1.0);
    const RopeId rp = rope(w, rigid(2.0 + kPi * r),
                           {on(w, heavy, Vec3{-r, -1, 0}), sheaveAt(Vec3{}, r), on(w, light, Vec3{r, -1, 0})});
    w.step(kDt);
    CHECK_MESSAGE(std::fabs(w.ropeSpan(rp) - (2.0 + kPi * r)) < 1e-6, "two straight metres and half a turn");
    run(w, 0.2);
    const double v0 = w.velocity(heavy).y;
    run(w, 0.3);
    CHECK(std::fabs(-(w.velocity(heavy).y - v0) / 0.3 / (0.5 * kG) - 1.0) < 1e-3);
    CHECK(std::fabs(w.ropeTension(rp) / (1.5 * kG) - 1.0) < 1e-3);
    CHECK_MESSAGE(std::fabs(w.position(heavy).x + r) < 1e-6, "each side hangs straight down from its tangent");
}

namespace {

    // A 100 kg load 3 m under a spool on the beam above it.
    struct Winch {
        World w;
        BodyId load;
        RopeId rope;
    };

    void winch(Winch &k, double mass, RopeDesc d) {
        k.load = ball(k.w, Vec3{0, -3, 0}, mass);
        d.spool = RopeSpool::Start;
        k.rope = rope(k.w, d, {fixedAt(Vec3{}), on(k.w, k.load, Vec3{0, -3, 0})});
    }

} // namespace

TEST_CASE("a spool winds its load up at its speed") {
    Winch k;
    RopeDesc d = rigid(3.0);
    d.speed = 0.5;
    d.maxForce = 5000.0;
    winch(k, 100.0, d);
    run(k.w, 1.0);
    k.w.setSpoolLever(k.rope, -1.0);
    run(k.w, 1.0);
    CHECK_MESSAGE(std::fabs(k.w.velocity(k.load).y - 0.5) < 5e-3, "it rises at the spool's speed");
    CHECK_MESSAGE(std::fabs(k.w.ropeLength(k.rope) - 2.5) < 1e-2, "with half a metre less rope paid out");
    CHECK_MESSAGE(std::fabs(k.w.position(k.load).y + k.w.ropeLength(k.rope)) < 1e-4, "hanging on what is out");
}

TEST_CASE("a spool too weak for its load lets it down instead of stretching") {
    Winch k;
    RopeDesc d = rigid(3.0);
    d.maxForce = 500.0;
    d.brakeForce = 5000.0;
    winch(k, 100.0, d);
    run(k.w, 0.5);
    k.w.setSpoolLever(k.rope, -1.0);
    run(k.w, 0.3);
    const double v0 = k.w.velocity(k.load).y;
    run(k.w, 0.3);
    CHECK_MESSAGE(std::fabs(k.w.ropeTension(k.rope) / 500.0 - 1.0) < 1e-2, "the motor pulls its cap");
    CHECK_MESSAGE(std::fabs((k.w.velocity(k.load).y - v0) / 0.3 / -((981.0 - 500.0) / 100.0) - 1.0) < 2e-2,
                  "and the load comes down at what the rest of its weight gives");
    CHECK(k.w.ropeInfo(k.rope).slipping);
    CHECK_MESSAGE(std::fabs(k.w.ropeSpan(k.rope) - k.w.ropeLength(k.rope)) < 1e-3, "paying out, not stretching");
}

TEST_CASE("an idle spool's brake holds what it can and slips under more") {
    Winch held;
    RopeDesc d = rigid(3.0);
    d.brakeForce = 1500.0;
    winch(held, 100.0, d);
    run(held.w, 1.0);
    const double y = held.w.position(held.load).y;
    run(held.w, 5.0);
    CHECK_MESSAGE(std::fabs(held.w.position(held.load).y - y) < 1e-5, "981 N under a 1500 N brake: it holds");

    Winch slips;
    winch(slips, 200.0, d);
    run(slips.w, 0.3);
    const double v0 = slips.w.velocity(slips.load).y;
    run(slips.w, 0.3);
    CHECK_MESSAGE(std::fabs((slips.w.velocity(slips.load).y - v0) / 0.3 / -((1962.0 - 1500.0) / 200.0) - 1.0) < 2e-2,
                  "1962 N: the brake slips, and the load comes down against it");
}

TEST_CASE("a spool stops at the ends of its rope") {
    Winch in;
    RopeDesc d = rigid(3.0);
    d.speed = 1.0;
    d.minLength = 1.0;
    winch(in, 100.0, d);
    in.w.setSpoolLever(in.rope, -1.0);
    run(in.w, 4.0);
    CHECK_MESSAGE(std::fabs(in.w.ropeLength(in.rope) - 1.0) < 1e-9, "wound in, it stops at its shortest");
    CHECK_MESSAGE(length(in.w.velocity(in.load)) < 1e-3, "and holds the load there");

    Winch out;
    d = rigid(3.0);
    d.speed = 1.0;
    d.maxLength = 4.0;
    d.brakeForce = 100.0; // far too weak for the load: only the knot on the drum holds it
    winch(out, 100.0, d);
    out.w.setSpoolLever(out.rope, 1.0);
    run(out.w, 4.0);
    CHECK_MESSAGE(std::fabs(out.w.ropeLength(out.rope) - 4.0) < 1e-9, "paid out to its end, it is tied there");
    CHECK_MESSAGE(std::fabs(out.w.position(out.load).y + 4.0) < 1e-3, "and holds the load at full length");
}

TEST_CASE("a spool at work keeps its load awake; braked and still, it sleeps") {
    Winch k;
    RopeDesc d = rigid(3.0);
    d.speed = 0.1;
    winch(k, 10.0, d);
    k.w.setSpoolLever(k.rope, -1.0);
    run(k.w, 3.0);
    CHECK_MESSAGE(!k.w.asleep(k.load), "winding, it is awake");
    k.w.setSpoolLever(k.rope, 0.0);
    run(k.w, 3.0);
    CHECK_MESSAGE(k.w.asleep(k.load), "braked, it falls asleep");
    k.w.setSpoolLever(k.rope, 1.0);
    CHECK_MESSAGE(!k.w.asleep(k.load), "the lever wakes it");
}

TEST_CASE("a rope pulled harder than it is strong breaks, and is mended by setRope") {
    World w;
    const BodyId load = ball(w, Vec3{0, -1, 0}, 100.0);
    RopeDesc d = rigid(1.0);
    d.strength = 500.0;
    const RopeId r = rope(w, d, {fixedAt(Vec3{}), on(w, load, Vec3{0, -1, 0})});
    w.step(kDt);
    REQUIRE(w.brokenRopes().size() == 1);
    CHECK(w.brokenRopes()[0] == r);
    CHECK(w.ropeInfo(r).broken);
    run(w, 0.5);
    CHECK_MESSAGE(w.position(load).y < -1.5, "broken, it lets go");
    CHECK(w.brokenRopes().empty());
    CHECK(w.ropeTension(r) == 0.0);
    d.strength = 2000.0;
    w.setPose(load, Vec3{0, -1, 0}, Quat{});
    w.setVelocity(load, Vec3{}, Vec3{});
    REQUIRE(w.setRope(r, d));
    run(w, 1.0);
    CHECK_MESSAGE(!w.ropeInfo(r).broken, "mended, it holds again");
    CHECK(std::fabs(w.position(load).y + 1.0) < 1e-4);
}

TEST_CASE("a spool lifts a crate off the floor: the floor lets go of it as the rope takes its weight") {
    host::Host h;
    const BodyId crate = h.addBox(Vec3{0, 0.5, 0}, Quat{}, Vec3{0.5, 0.5, 0.5}, 500.0);
    RopeDesc d = rigid(4.4);
    d.spool = RopeSpool::Start;
    d.speed = 0.25;
    d.maxForce = 2e4;
    const RopeId r = h.world.addRope(d, std::vector<RopePoint>{fixedAt(Vec3{0, 5, 0}), on(h.world, crate, Vec3{0, 1, 0})});
    for (int s = 0; s < 60; ++s) h.step(kDt);
    CHECK_MESSAGE(h.world.ropeTension(r) == 0.0, "slack, the crate stands on the floor");
    h.world.setSpoolLever(r, -1.0);
    double maxTension = 0.0;
    for (int s = 0; s < 240; ++s) {
        h.step(kDt);
        maxTension = std::max(maxTension, h.world.ropeTension(r));
    }
    // 4.4 m out, 4 m to the crate's top: 0.4 m of slack at 0.25 m/s, then two and a half seconds of lifting.
    CHECK_MESSAGE(std::fabs(h.world.velocity(crate).y - 0.25) < 5e-3, "off the floor, it rises at the spool's speed");
    CHECK_MESSAGE(h.world.position(crate).y > 0.5 + 0.25, "well clear of it");
    CHECK_MESSAGE(std::fabs(h.world.ropeTension(r) / (500.0 * kG) - 1.0) < 1e-2, "on the rope alone");
    CHECK_MESSAGE(maxTension < 1.5 * 500.0 * kG, "taking the weight with no snatch worth the name");
    CHECK(std::fabs(h.world.ropeSpan(r) - h.world.ropeLength(r)) < 1e-4);
}

TEST_CASE("a crane jib on a capped cylinder hoists a load from its tip over a sheave") {
    host::Host h;
    const BodyId mast = h.addBox(Vec3{0, 2.5, 0}, Quat{}, Vec3{0.2, 2.5, 0.2}, 0.0);
    const BodyId jib = h.addBox(Vec3{3, 5, 0}, Quat{}, Vec3{3.0, 0.15, 0.15}, 300.0);
    h.addJoint(h.jointAt(JointType::Hinge, mast, jib, Vec3{0, 5, 0}, kZ));
    // The luffing cylinder: a capped piston from the mast to the jib's underside.
    const BodyId barrel = h.addBox(Vec3{0.5, 3.75, 0}, xTo(Vec3{1, 0.6, 0}), Vec3{0.6, 0.08, 0.08}, 40.0);
    const BodyId rod = h.addBox(Vec3{1.5, 4.35, 0}, xTo(Vec3{1, 0.6, 0}), Vec3{0.6, 0.04, 0.04}, 20.0);
    h.addJoint(h.jointAt(JointType::Hinge, barrel, mast, Vec3{0, 3.45, 0}, kZ));
    h.addJoint(h.jointAt(JointType::Hinge, rod, jib, Vec3{2.0, 4.65, 0}, kZ));
    JointDesc p = h.jointAt(JointType::Piston, barrel, rod, Vec3{1.0, 4.05, 0}, Vec3{1, 0.6, 0});
    p.maxForce = 2e5;
    p.min = -0.3;
    p.max = 0.3;
    const JointId luff = h.addJoint(p);
    h.ignore(barrel, rod);
    h.ignore(barrel, jib);
    h.ignore(rod, mast);
    const BodyId load = h.addBox(Vec3{6, 0.5, 0}, Quat{}, Vec3{0.5, 0.5, 0.5}, 800.0);
    RopeDesc d = rigid(5.0);
    d.spool = RopeSpool::Start;
    d.speed = 0.5;
    d.maxForce = 5e4;
    // From a winch at the mast's foot, up to a sheave at the jib's tip, down to the load.
    RopePoint tip = on(h.world, jib, Vec3{6, 5, 0});
    tip.radius = 0.2;
    tip.axis = conjugate(h.world.orientation(jib)) * kZ;
    const RopeId r = h.world.addRope(d, std::vector<RopePoint>{fixedAt(Vec3{0.3, 0.2, 0}), tip, on(h.world, load, Vec3{6, 1, 0})});
    h.world.setSpoolLever(r, -1.0);
    double stretch = 0.0;
    for (int s = 0; s < 360; ++s) {
        h.step(kDt);
        stretch = std::max(stretch, h.world.ropeSpan(r) - h.world.ropeLength(r));
        REQUIRE(h.world.faults().empty());
    }
    CHECK_MESSAGE(h.world.ropeInfo(r).machine == -1, "swept between the crane and the load");
    CHECK_MESSAGE(h.world.position(load).y > 1.0, "the load is up off the ground");
    CHECK_MESSAGE(stretch < 1e-3, "the hoist rope never stretches a millimetre");
    CHECK_MESSAGE(std::fabs(h.world.jointPosition(luff)) < 0.05, "and the cylinder holds the jib where it was");
    CHECK_MESSAGE(std::fabs(h.world.ropeTension(r) / (800.0 * kG) - 1.0) < 2e-2, "the rope carries the load");
}

TEST_CASE("a trolley moved by the host carries its load on a rope, which swings and never stretches") {
    World w;
    BodyDesc t;
    t.position = Vec3{0, 0, 0};
    const BodyId trolley = w.addBody(t);
    w.setKinematic(trolley, true);
    const BodyId load = ball(w, Vec3{0, -3, 0}, 200.0);
    const RopeId r = rope(w, rigid(3.0), {on(w, trolley, Vec3{}), on(w, load, Vec3{0, -3, 0})});
    double stretch = 0.0;
    for (int s = 1; s <= 180; ++s) {
        w.setPose(trolley, Vec3{std::min(1.0, s * kDt), 0, 0}, Quat{});
        w.step(kDt);
        stretch = std::max(stretch, w.ropeSpan(r) - 3.0);
    }
    CHECK_MESSAGE(stretch < 1e-4, "the rope holds its length with the trolley driving it");
    CHECK_MESSAGE(w.position(load).x > 0.3, "the load follows");
    CHECK_MESSAGE(w.ropeTension(r) > 0.5 * 200.0 * kG, "swinging under it");
}

TEST_CASE("a movable pulley: a sheave on a free block halves the pull and moves at half the rope's speed") {
    World w;
    const double r = 0.25;
    const BodyId block = box(w, Vec3{0, -2, 0}, 200.0, Vec3{0.1, 0.1, 0.1});
    const BodyId weight = ball(w, Vec3{2, -1, 0}, 100.0);
    RopePoint wheel = on(w, block, Vec3{0, -2, 0});
    wheel.radius = r;
    wheel.axis = kZ;
    // Down from the beam on the left of the wheel, round under it, up on the right over a beam pulley, down to 100 kg.
    const RopeId rp = rope(w, rigid(2.0 + kPi * r + 2.0 + (2.0 - r) + 1.0),
                           {fixedAt(Vec3{-r, 0, 0}), wheel, fixedAt(Vec3{r, 0, 0}), fixedAt(Vec3{2, 0, 0}),
                            on(w, weight, Vec3{2, -1, 0})});
    run(w, 3.0);
    CHECK_MESSAGE(length(w.velocity(block)) < 1e-6, "100 kg holds 200 kg on two strands");
    CHECK(length(w.velocity(weight)) < 1e-6);
    CHECK_MESSAGE(std::fabs(w.ropeTension(rp) / (100.0 * kG) - 1.0) < 1e-6, "with the rope pulling 100 kg's weight");
    CHECK_MESSAGE(std::fabs(w.position(block).x) < 1e-6, "and the block hangs centred under its two strands");
    w.applyForce(weight, Vec3{0, -500, 0}, w.centerOfMass(weight));
    w.step(kDt);
    const Vec3 b0 = w.position(block), c0 = w.position(weight);
    for (int s = 0; s < 30; ++s) {
        w.applyForce(weight, Vec3{0, -500, 0}, w.centerOfMass(weight));
        w.step(kDt);
    }
    CHECK_MESSAGE(std::fabs((c0.y - w.position(weight).y) / (w.position(block).y - b0.y) - 2.0) < 1e-3,
                  "pulled down, the weight goes twice as far as the block rises");
}

TEST_CASE("a rope re-routed through a new eye keeps what is paid out, and swings from the eye") {
    World w;
    const BodyId load = ball(w, Vec3{0, -3, 0}, 30.0);
    const RopeId r = rope(w, rigid(3.0), {fixedAt(Vec3{0, 0, 0}), on(w, load, Vec3{0, -3, 0})});
    run(w, 1.0);
    // An eye a metre down the rope, where it already runs: nothing moves, but the load now hangs 2 m below it.
    const std::vector<RopePoint> through{fixedAt(Vec3{0, 0, 0}), fixedAt(Vec3{0, -1, 0}), on(w, load, Vec3{0, -3, 0})};
    REQUIRE(w.setRopeRoute(r, through));
    CHECK(w.ropeLength(r) == 3.0);
    w.setVelocity(load, Vec3{2, 0, 0}, Vec3{});
    double furthest = 0.0;
    for (int s = 0; s < 120; ++s) {
        w.step(kDt);
        furthest = std::max(furthest, length(w.position(load) - Vec3{0, -1, 0}));
    }
    RopeBend bends[3];
    CHECK(w.ropePath(r, bends) == 3);
    CHECK_MESSAGE(w.position(load).x != 0.0, "pushed, it swings");
    CHECK_MESSAGE(furthest < 2.0 + 1e-4, "on a circle of 2 m about the eye, not 3 m about the top");
    CHECK_MESSAGE(std::fabs(w.ropeSpan(r) - 3.0) < 1e-4, "its span through the eye is still what is paid out");
}

TEST_CASE("spools at both ends wind in twice as fast") {
    Winch k;
    RopeDesc d = rigid(3.0);
    d.speed = 0.25;
    d.spool = RopeSpool::Both;
    k.load = ball(k.w, Vec3{0, -3, 0}, 50.0);
    k.rope = rope(k.w, d, {fixedAt(Vec3{}), on(k.w, k.load, Vec3{0, -3, 0})});
    k.w.setSpoolLever(k.rope, -1.0);
    run(k.w, 1.0);
    CHECK(std::fabs(k.w.velocity(k.load).y - 0.5) < 5e-3);
}

namespace {

    // A NaN in a swept rope's impulse, once its sweeps are done.
    void poisonRope(World::Impl &I, int substep, FaultPhase phase) {
        if (phase != FaultPhase::Contacts || substep != 4 || I.ropePreps.empty()) return;
        I.ropePreps[0].lambda = std::numeric_limits<double>::quiet_NaN();
    }

} // namespace

TEST_CASE("a rope gone wrong is put back with its island, and its spool's length with it") {
    Winch k;
    RopeDesc d = rigid(3.0);
    d.speed = 1.0;
    winch(k, 100.0, d);
    run(k.w, 1.0);
    k.w.setSpoolLever(k.rope, -1.0);
    const double length = k.w.ropeLength(k.rope);
    const Vec3 at = k.w.position(k.load);
    testing::phaseHook = poisonRope;
    k.w.step(kDt);
    testing::phaseHook = nullptr;
    REQUIRE(k.w.faults().size() == 1);
    CHECK(k.w.ropeLength(k.rope) == length);
    CHECK((k.w.position(k.load) == at));
    run(k.w, 1.0);
    CHECK_MESSAGE(k.w.faults().empty(), "and it carries on winding");
    CHECK(std::fabs(k.w.velocity(k.load).y - 1.0) < 1e-2);
}

TEST_CASE("a winch raises a hinged boom by a rope at its tip: a spool held in the boom's own solve") {
    World w;
    const BodyId boom = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 50.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, boom, Vec3{0, 0, 0}, kZ));
    RopeDesc d = rigid(std::sqrt(8.0));
    d.spool = RopeSpool::Start;
    d.speed = 0.2;
    d.maxForce = 2000.0;
    const RopeId r = rope(w, d, {fixedAt(Vec3{0, 2, 0}), on(w, boom, Vec3{2, 0, 0})});
    run(w, 0.5);
    w.setSpoolLever(r, -1.0);
    const auto tip = [&] { return w.position(boom) + w.orientation(boom) * Vec3{1, 0, 0}; };
    const double before = tip().y;
    double stretch = 0.0;
    for (int s = 0; s < 120; ++s) {
        w.step(kDt);
        stretch = std::max(stretch, w.ropeSpan(r) - w.ropeLength(r));
    }
    CHECK_MESSAGE(w.ropeInfo(r).machine >= 0, "the boom's solve holds the rope, spool and all");
    CHECK_MESSAGE(tip().y > before + 0.2, "winding in raises the boom");
    CHECK_MESSAGE(std::fabs(length(tip() - Vec3{0, 2, 0}) - w.ropeLength(r)) < 1e-4, "on the rope that is out");
    CHECK_MESSAGE(stretch < 1e-4, "which never stretches");

    // Too weak: 100 N at the tip can't hold 490 N·m about the hinge. It lets the boom down, slipping.
    World weak;
    const BodyId b2 = rodBetween(weak, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 50.0);
    weak.addJoint(jointAt(weak, JointType::Hinge, BodyId{}, b2, Vec3{0, 0, 0}, kZ));
    d.maxForce = 100.0;
    d.brakeForce = 100.0;
    const RopeId r2 = rope(weak, d, {fixedAt(Vec3{0, 2, 0}), on(weak, b2, Vec3{2, 0, 0})});
    run(weak, 0.5);
    CHECK_MESSAGE(weak.ropeInfo(r2).slipping, "the brake slips");
    CHECK_MESSAGE(std::fabs(weak.ropeTension(r2) / 100.0 - 1.0) < 1e-2, "holding back with all it has");
    CHECK(weak.ropeLength(r2) > std::sqrt(8.0) + 0.1);
}

TEST_CASE("a rope between two links of one machine is held inside it: it keeps an elbow from opening") {
    World w;
    const BodyId upper = rodBetween(w, Vec3{0, 0, 0}, Vec3{2, 0, 0}, 20.0);
    const BodyId lower = rodBetween(w, Vec3{2, 0, 0}, Vec3{2, -2, 0}, 20.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, upper, Vec3{0, 0, 0}, kZ));
    w.addJoint(jointAt(w, JointType::Hinge, upper, lower, Vec3{2, 0, 0}, kZ));
    // Between the links' middles across the elbow's right angle: it may close, never open.
    const RopeId r = rope(w, rigid(std::sqrt(2.0)), {on(w, upper, Vec3{1, 0, 0}), on(w, lower, Vec3{2, -1, 0})});
    double stretch = 0.0, pulled = 0.0;
    for (int s = 0; s < 300; ++s) {
        const Vec3 tip = w.position(lower) + w.orientation(lower) * Vec3{1, 0, 0};
        w.applyForce(lower, w.orientation(lower) * Vec3{0, 300, 0}, tip);
        w.step(kDt);
        stretch = std::max(stretch, w.ropeSpan(r) - std::sqrt(2.0));
        pulled = std::max(pulled, w.ropeTension(r));
    }
    CHECK_MESSAGE(w.ropeInfo(r).machine >= 0, "held in the machine's solve");
    CHECK_MESSAGE(w.stats().machines == 1, "one machine");
    CHECK_MESSAGE(pulled > 100.0, "pushed open, the rope takes the push");
    CHECK_MESSAGE(stretch < 1e-4, "and never stretches");
}
