#include <cmath>
#include <cstring>
#include <limits>
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
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    bool finite(const Vec3 &v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

    bool finite(World &w, BodyId b) {
        const Quat q = w.orientation(b);
        return finite(w.position(b)) && finite(w.velocity(b)) && finite(w.angularVelocity(b)) && std::isfinite(q.w) &&
               std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z);
    }

    // A body's whole state, to compare bit for bit.
    struct State {
        Vec3 c, v, w;
        Quat q;

        bool operator==(const State &o) const { return std::memcmp(this, &o, sizeof(State)) == 0; }
    };

    State stateOf(World &w, BodyId b) { return {w.position(b), w.velocity(b), w.angularVelocity(b), w.orientation(b)}; }

    // The poison: a NaN written into a body's state after the machines' own solve of one sub-step. A machine link's
    // state is its machine's rates, so the NaN goes there; a free body's is its group's.
    struct Poison {
        BodyId body;
        int substep = 3;
        bool everyStep = false;
        bool done = false;
    };
    Poison poison;

    void poisonHook(World::Impl &I, int substep, FaultPhase phase) {
        if (phase != FaultPhase::Internal || substep != poison.substep || poison.done) return;
        Group &g = I.groups[std::size_t(I.body(poison.body)->group)];
        if (g.machine >= 0) {
            MachineData &md = I.machines[std::size_t(g.machine)];
            md.m.art.qdot[std::size_t(md.m.art.links[std::size_t(g.link)].dof0)] = kNaN;
        } else {
            g.v.x = kNaN;
        }
        poison.done = !poison.everyStep;
    }

    // Steps once with the poison in, and takes it out again.
    void poisonedStep(World &w, BodyId body, int substep = 3) {
        poison = {body, substep};
        testing::phaseHook = poisonHook;
        w.step(kDt);
        testing::phaseHook = nullptr;
    }

    void handContacts(host::Host &h) {
        h.detect();
        h.world.clearContacts();
        for (const Contact &k : h.contacts) h.world.addContact(k);
    }

    // A three-link machine falling onto the floor, and five cubes falling in a stack well away from it: two islands
    // standing on one ground.
    struct TwoIslands {
        host::Host h;
        std::vector<BodyId> links, cubes;

        TwoIslands() {
            const Vec3 half{0.5, 0.1, 0.2};
            for (int i = 0; i < 3; ++i) links.push_back(h.addBox(Vec3{double(i), 0.3, 0}, Quat{}, half, 20.0));
            for (std::size_t i = 0; i < 2; ++i)
                h.addJoint(h.jointAt(JointType::Hinge, links[i], links[i + 1], Vec3{double(i) + 0.5, 0.3, 0}, kZ));
            for (int i = 0; i < 5; ++i)
                cubes.push_back(h.addBox(Vec3{8, 0.6 + 1.05 * i, 0}, Quat{}, Vec3{0.5, 0.5, 0.5}, 5.0));
        }
    };

} // namespace

TEST_CASE("a NaN is put back in the island it appeared in, before anything else sees it") {
    TwoIslands sick, well;
    for (int i = 0; i < 4; ++i) {
        sick.h.step(kDt);
        well.h.step(kDt);
    }
    std::vector<Vec3> before;
    for (BodyId link : sick.links) before.push_back(sick.h.world.position(link));

    handContacts(sick.h);
    poisonedStep(sick.h.world, sick.links[1]);
    handContacts(well.h);
    well.h.world.step(kDt);

    World &w = sick.h.world;
    for (std::size_t i = 0; i < sick.links.size(); ++i) {
        CHECK_MESSAGE(length(w.position(sick.links[i]) - before[i]) < 1e-9, "the machine is back where it began");
        CHECK_MESSAGE((w.velocity(sick.links[i]) == Vec3{} && w.angularVelocity(sick.links[i]) == Vec3{}),
                      "and at rest");
    }
    CHECK_MESSAGE(w.velocity(sick.h.ground()) == Vec3{}, "the ground both islands stand on was only read");
    for (std::size_t i = 0; i < sick.cubes.size(); ++i)
        CHECK_MESSAGE(stateOf(w, sick.cubes[i]) == stateOf(well.h.world, well.cubes[i]),
                      "the cubes stepped as if nothing had happened, bit for bit");

    REQUIRE(w.faults().size() == 1);
    const Fault &f = w.faults()[0];
    // Spread through the machine by the end of the step, the NaN is in every link: the first one found is named.
    CHECK_MESSAGE((f.body == sick.links[0] || f.body == sick.links[1] || f.body == sick.links[2]),
                  "the fault names a link of the poisoned machine");
    CHECK_MESSAGE(f.kind == FaultKind::NonFinite, "as a NaN");
    CHECK_MESSAGE(f.phase == FaultPhase::EndOfStep, "seen once the step was done");
    CHECK_MESSAGE(f.machine >= 0, "in its machine");
    CHECK_MESSAGE(f.bodies == 3, "and puts back its island of three");
    CHECK_MESSAGE(f.streak == 1, "for the first time");
    CHECK_MESSAGE(f.step == 5, "on the fifth step");
    CHECK_MESSAGE(w.stats().faults == 1, "counted in the stats");

    sick.h.step(kDt);
    for (BodyId link : sick.links) CHECK_MESSAGE(finite(w, link), "the next step is a number");
    CHECK_MESSAGE(w.faults().empty(), "and no fault");
}

TEST_CASE("looking after every phase, a fault names the phase it appeared in") {
    TwoIslands s;
    WorldSettings settings = s.h.world.settings();
    settings.checkPhases = true;
    s.h.world.setSettings(settings);
    s.h.step(kDt);
    handContacts(s.h);
    poisonedStep(s.h.world, s.links[1]);
    REQUIRE(s.h.world.faults().size() == 1);
    CHECK_MESSAGE(s.h.world.faults()[0].phase == FaultPhase::Internal, "the machines' own solve");
    CHECK_MESSAGE(s.h.world.faults()[0].body == s.links[1], "caught before it spread, it names the poisoned link");
}

TEST_CASE("an island that faults step after step counts its streak, and stays awake") {
    World w;
    BodyDesc d;
    d.position = Vec3{0, 5, 0};
    d.inertia = boxInertia(1.0, Vec3{0.5, 0.5, 0.5});
    const BodyId cube = w.addBody(d);
    poison = {cube, 0, true};
    testing::phaseHook = poisonHook;
    for (int i = 1; i <= 12; ++i) {
        w.step(kDt);
        REQUIRE(w.faults().size() == 1);
        CHECK_MESSAGE(w.faults()[0].streak == i, "one more each step");
    }
    testing::phaseHook = nullptr;
    CHECK_MESSAGE(!w.asleep(cube), "at rest twelve steps, but kept awake to try again");
    CHECK_MESSAGE((w.position(cube) == Vec3{0, 5, 0}), "where it began, every time");
    w.step(kDt);
    CHECK_MESSAGE(w.faults().empty(), "the poison gone, it falls");
    CHECK_MESSAGE(w.velocity(cube).y < 0.0, "under gravity again");
}

TEST_CASE("a body faster than the fault speed is put back like a NaN") {
    WorldSettings settings;
    settings.gravity = Vec3{};
    settings.faultSpeed = 20.0;
    World w(settings);
    BodyDesc d;
    d.inertia = boxInertia(1.0, Vec3{0.5, 0.5, 0.5});
    d.velocity = Vec3{10, 0, 0};
    const BodyId slow = w.addBody(d);
    d.position = Vec3{0, 10, 0};
    d.velocity = Vec3{30, 0, 0};
    const BodyId fast = w.addBody(d);
    w.step(kDt);
    REQUIRE(w.faults().size() == 1);
    CHECK_MESSAGE(w.faults()[0].kind == FaultKind::Speed, "a speed fault");
    CHECK_MESSAGE(w.faults()[0].body == fast, "the fast body's");
    CHECK_MESSAGE((w.position(fast) == Vec3{0, 10, 0}), "back where it was");
    CHECK_MESSAGE(w.velocity(fast) == Vec3{}, "at rest");
    CHECK_MESSAGE((w.velocity(slow) == Vec3{10, 0, 0}), "and the slow one carries on");
}

TEST_CASE("a lever let go on a step that faults still parks its drive where it is") {
    // The fault puts the drive back as the step found it, the lever still held, so the next step lets it go again.
    World w;
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
    p.maxForce = 2e5;
    p.speed = 0.5;
    p.min = -0.3;
    p.max = 0.2;
    const JointId piston = w.addJoint(p);

    run(w, 0.5);
    w.setLever(piston, 1.0);
    run(w, 0.3);
    w.setLever(piston, 0.0);
    poisonedStep(w, boom);
    REQUIRE(w.faults().size() == 1);
    w.step(kDt);
    const double parked = w.jointPosition(piston);
    CHECK_MESSAGE(parked > 0.08, "held, the lever moved it");
    // A ram under load trails its target: parked anywhere else, it would finish that travel by itself.
    CHECK_MESSAGE(std::fabs(w.target(piston) - parked) < 1e-6, "let go, its target parks where it is");
    run(w, 2.0);
    CHECK_MESSAGE(std::fabs(w.jointPosition(piston) - parked) < 1e-6, "and it stays there");
}

TEST_CASE("a thin link spun far past what a sub-step can follow is put back at rest") {
    // At 5000 rad/s a sub-step turns the link about five radians, and the midpoint passes return NaN.
    WorldSettings settings;
    settings.gravity = Vec3{};
    World w(settings);
    BodyDesc d;
    d.mass = 1.0;
    d.position = Vec3{0.5, 2, 0};
    d.inertia = boxInertia(1.0, Vec3{0.5, 0.01, 0.01});
    d.angularVelocity = Vec3{5000.0, 0.0, 1000.0};
    d.velocity = cross(d.angularVelocity, Vec3{0.5, 0, 0});
    const BodyId rod = w.addBody(d);
    JointDesc j;
    j.type = JointType::Ball;
    j.b = rod;
    j.frameA = Frame{Vec3{0, 2, 0}, Quat{}};
    j.frameB = Frame{Vec3{-0.5, 0, 0}, Quat{}};
    w.addJoint(j);
    w.step(kDt);
    REQUIRE(w.faults().size() == 1);
    CHECK_MESSAGE(w.faults()[0].kind == FaultKind::NonFinite, "a NaN");
    CHECK_MESSAGE((length(w.position(rod) - Vec3{0.5, 2, 0}) < 1e-9), "put back where it began");
    CHECK_MESSAGE(w.angularVelocity(rod) == Vec3{}, "at rest");
    w.step(kDt);
    CHECK_MESSAGE((w.faults().empty() && finite(w, rod)), "and it stays there");
}

TEST_CASE("values the world can't use are refused, changing nothing, and counted") {
    World w;
    BodyDesc g;
    g.isStatic = true;
    const BodyId ground = w.addBody(g);
    BodyDesc d;
    d.position = Vec3{0, 0.5, 0};
    d.inertia = boxInertia(1.0, Vec3{0.5, 0.5, 0.5});
    const BodyId cube = w.addBody(d);
    const BodyId arm = rodBetween(w, Vec3{3, 1, 0}, Vec3{4, 1, 0}, 5.0);
    JointDesc s = jointAt(w, JointType::Servo, BodyId{}, arm, Vec3{3, 1, 0}, kZ);
    s.maxForce = 1e3;
    const JointId servo = w.addJoint(s);

    Contact k;
    k.a = ground;
    k.b = cube;
    k.normal = Vec3{};
    CHECK_MESSAGE(!w.addContact(k), "a normal of no length");
    k.normal = Vec3{0, 1, 0};
    k.friction = -0.1;
    CHECK_MESSAGE(!w.addContact(k), "a negative friction");
    k.friction = 0.6;
    k.restitution = 3.0;
    CHECK_MESSAGE(w.addContact(k), "a restitution above 1 is clamped, not refused");

    CHECK_MESSAGE((!w.setPose(cube, Vec3{kNaN, 0, 0}, Quat{})), "a pose that isn't a number");
    CHECK_MESSAGE((w.position(cube) == Vec3{0, 0.5, 0}), "leaves the body where it was");
    CHECK_MESSAGE(!w.setTarget(servo, kNaN), "a target that isn't a number");
    CHECK_MESSAGE(w.target(servo) == 0.0, "leaves the target as it was");
    d.inertia = Mat3(1.0, 0, 0, 0, -1.0, 0, 0, 0, 1.0);
    CHECK_MESSAGE(!w.addBody(d).valid(), "an inertia that isn't positive definite");
    CHECK_MESSAGE(w.stats().bodies == 3, "adds no body");

    w.step(kDt);
    CHECK_MESSAGE(w.stats().rejected == 5, "five calls refused before that step");
    w.step(kDt);
    CHECK_MESSAGE(w.stats().rejected == 0, "and none before the next");
    CHECK_MESSAGE((!w.setPose(BodyId{}, Vec3{}, Quat{})), "a handle that's gone is refused too");
    w.step(kDt);
    CHECK_MESSAGE(w.stats().rejected == 0, "but isn't a bad value, so isn't counted");
}

TEST_CASE("each kind of value the world can't use is refused") {
    World w;
    const BodyId a = w.addBody(BodyDesc{}), b = w.addBody(BodyDesc{});
    const auto refusedBody = [&](auto edit) {
        BodyDesc x;
        edit(x);
        return !w.addBody(x).valid();
    };
    CHECK(refusedBody([](BodyDesc &x) { x.mass = 0.0; }));
    CHECK(refusedBody([](BodyDesc &x) { x.mass = kInf; }));
    CHECK(refusedBody([](BodyDesc &x) { x.inertia = Mat3(1, 0.5, 0, 0, 1, 0, 0, 0, 1); })); // not symmetric
    CHECK(refusedBody([](BodyDesc &x) { x.inertia = Mat3(1, 0, 0, 0, 0, 0, 0, 0, 1); }));   // not definite
    CHECK(refusedBody([](BodyDesc &x) { x.orientation = Quat{0, 0, 0, 0}; }));
    CHECK(refusedBody([](BodyDesc &x) { x.velocity = Vec3{0, kNaN, 0}; }));
    CHECK_MESSAGE(!refusedBody([](BodyDesc &x) {
        x.isStatic = true;
        x.mass = 0.0;
        x.inertia = Mat3{};
    }),
                  "a static body needs no mass");
    CHECK(!w.setMassProperties(a, -1.0, Mat3(1.0)));
    CHECK(!w.setVelocity(a, Vec3{}, Vec3{kNaN, 0, 0}));
    CHECK(!w.applyForce(a, Vec3{kNaN, 0, 0}, Vec3{}));
    CHECK(!w.applyForce(a, Vec3{1, 0, 0}, Vec3{0, kInf, 0}));
    CHECK(!w.applyTorque(a, Vec3{0, 0, kNaN}));

    const JointDesc j = jointAt(w, JointType::Piston, a, b, Vec3{}, Vec3{1, 0, 0});
    const auto refusedJoint = [&](auto edit) {
        JointDesc x = j;
        edit(x);
        return !w.addJoint(x).valid();
    };
    CHECK(refusedJoint([](JointDesc &x) { x.frameA.position.x = kNaN; }));
    CHECK(refusedJoint([](JointDesc &x) { x.frameB.rotation = Quat{0, 0, 0, 0}; }));
    CHECK(refusedJoint([](JointDesc &x) { x.maxForce = -1.0; }));
    CHECK(refusedJoint([](JointDesc &x) { x.forceSlew = kNaN; }));
    CHECK(refusedJoint([](JointDesc &x) { x.speed = -1.0; }));
    CHECK(refusedJoint([](JointDesc &x) { x.hertz = kNaN; }));
    CHECK(refusedJoint([](JointDesc &x) { x.damping = -0.5; }));
    CHECK(refusedJoint([](JointDesc &x) { x.restLength = -1.0; }));
    CHECK(refusedJoint([](JointDesc &x) {
        x.min = 1.0;
        x.max = 0.0;
    }));
    CHECK(refusedJoint([](JointDesc &x) { x.target = kNaN; }));
    const JointId piston = w.addJoint(j);
    REQUIRE(piston.valid());
    JointDesc bad = j;
    bad.min = kNaN;
    CHECK_MESSAGE(!w.setJoint(piston, bad), "an edit is checked like a new joint");
    CHECK(!w.setLever(piston, kNaN));

    Contact k;
    k.a = a;
    k.b = b;
    k.normal = Vec3{0, 1, 0};
    const auto refusedContact = [&](auto edit) {
        Contact x = k;
        edit(x);
        return !w.addContact(x);
    };
    CHECK(refusedContact([](Contact &x) { x.point.z = kNaN; }));
    CHECK(refusedContact([](Contact &x) { x.depth = kInf; }));
    CHECK(refusedContact([](Contact &x) { x.normal = Vec3{0, 3, 0}; }));
    CHECK(refusedContact([](Contact &x) { x.rolling = -0.01; }));
    CHECK(refusedContact([](Contact &x) { x.spinning = kNaN; }));
    CHECK(refusedContact([](Contact &x) { x.restitution = kNaN; }));
    CHECK_MESSAGE((!refusedContact([](Contact &x) { x.normal = Vec3{0, 1.5, 0}; })),
                  "a normal off unit length is fine");
}

TEST_CASE("a spring between two machines anchored on their hinge axes does nothing") {
    // Neither anchor can move, so the spring's row has no mobility at all: it was 0/0.
    World w;
    const BodyId a = rodBetween(w, Vec3{0, 1, 0}, Vec3{1, 1, 0}, 10.0);
    const BodyId b = rodBetween(w, Vec3{2, 1, 0}, Vec3{3, 1, 0}, 10.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, a, Vec3{0, 1, 0}, kZ));
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, b, Vec3{2, 1, 0}, kZ));
    JointDesc s = jointAt(w, JointType::Spring, a, b, Vec3{0, 1, 0}, kZ);
    s.frameB = jointAt(w, JointType::Spring, a, b, Vec3{2, 1, 0}, kZ).frameB;
    s.hertz = 2.0;
    s.restLength = 1.0;
    w.addJoint(s);
    int faults = 0;
    for (int i = 0; i < 60; ++i) {
        w.step(kDt);
        faults += w.stats().faults;
    }
    CHECK_MESSAGE((finite(w, a) && finite(w, b)), "both arms stay finite");
    CHECK_MESSAGE(faults == 0, "with nothing to put back");
    // Each arm swings down under gravity alone, about its own hinge: its centre ends up below its hinge.
    CHECK_MESSAGE(w.position(a).y < 1.0, "the first arm swings as if the spring weren't there");
    CHECK_MESSAGE(w.position(b).y < 1.0, "and so does the second");
}

TEST_CASE("a thin rotor spun up by a torque spins at the rate the torque gives it") {
    // ω × Iω, solved at the midpoint, has no answer for a Newton step whose Jacobian is singular; then the spin
    // keeps its value for the sub-step rather than dividing by zero. Spun up about its axis it never gets there.
    WorldSettings settings;
    settings.gravity = Vec3{};
    World w(settings);
    BodyDesc d;
    d.mass = 1.0;
    d.inertia = Mat3(0.002, 0, 0, 0, 0.5, 0, 0, 0, 0.5);
    const BodyId rotor = w.addBody(d);
    const double torque = 1.0, seconds = 10.0;
    for (int i = 0; i < int(std::lround(seconds / kDt)); ++i) {
        w.applyTorque(rotor, w.orientation(rotor) * Vec3{torque, 0, 0});
        w.step(kDt);
    }
    const double expected = torque * seconds / 0.002;
    CHECK_MESSAGE(finite(w, rotor), "it stays finite");
    CHECK_MESSAGE(std::fabs(length(w.angularVelocity(rotor)) - expected) < 1e-6 * expected,
                  "τ t / I: 5000 rad/s after 10 s");
}

TEST_CASE("a capped motor in a truss that locks it stalls at its cap, and nothing moves") {
    // Its row has no mobility once the loop holds it: the sweeps must not divide by that.
    World w;
    const Vec3 p0{0, 1, 0}, p1{2, 1, 0}, p2{1, 2.5, 0};
    const BodyId a = rodBetween(w, p0, p2, 5.0), b = rodBetween(w, p1, p2, 5.0);
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, a, p0, kZ));
    w.addJoint(jointAt(w, JointType::Hinge, BodyId{}, b, p1, kZ));
    JointDesc m = jointAt(w, JointType::Motor, a, b, p2, kZ);
    m.maxForce = 100.0;
    m.target = 3.0;
    const JointId motor = w.addJoint(m);
    const Vec3 at = w.position(a);
    run(w, 2.0);
    CHECK_MESSAGE((finite(w, a) && finite(w, b)), "the truss stays finite");
    CHECK_MESSAGE(length(w.position(a) - at) < 1e-3, "and stands where it was");
    CHECK_MESSAGE(std::fabs(std::fabs(w.driveForce(motor)) - 100.0) < 1e-9, "the motor pushes at its cap, stalled");
}
