#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kZ{0, 0, 1};

    bool finite(const Vec3 &v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

    bool finite(World &w, BodyId b) {
        const Quat q = w.orientation(b);
        return finite(w.position(b)) && finite(w.velocity(b)) && finite(w.angularVelocity(b)) && std::isfinite(q.w) &&
               std::isfinite(q.x) && std::isfinite(q.y) && std::isfinite(q.z);
    }

    BodyId addCube(World &w, const Vec3 &at) {
        BodyDesc d;
        d.position = at;
        d.inertia = boxInertia(1.0, Vec3{0.5, 0.5, 0.5});
        return w.addBody(d);
    }

    // The host's detection for unit cubes over the plane y = 0: their four lowest corners.
    void cubesOnGround(World &w, BodyId ground, std::initializer_list<BodyId> cubes) {
        w.clearContacts();
        std::uint64_t feature = 0;
        for (BodyId cube : cubes) {
            const Vec3 c = w.position(cube);
            const Quat q = w.orientation(cube);
            for (double sx : {-0.5, 0.5})
                for (double sy : {-0.5, 0.5})
                    for (double sz : {-0.5, 0.5}) {
                        const Vec3 p = c + q * Vec3{sx, sy, sz};
                        ++feature;
                        if (!(p.y < 0.05)) continue;
                        Contact k;
                        k.a = ground;
                        k.b = cube;
                        k.point = p;
                        k.normal = Vec3{0, 1, 0};
                        k.depth = -p.y;
                        k.feature = feature;
                        w.addContact(k);
                    }
        }
    }

} // namespace

TEST_CASE("a NaN on one body never reaches another through the ground they share") {
    World w;
    BodyDesc g;
    g.isStatic = true;
    const BodyId ground = w.addBody(g);
    const BodyId sick = addCube(w, Vec3{0, 0.5, 0}), well = addCube(w, Vec3{3, 0.5, 0});
    for (int i = 0; i < 3; ++i) { // still awake: a sleeper is never stepped, and would carry nothing anywhere
        cubesOnGround(w, ground, {sick, well});
        w.step(kDt);
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    w.setVelocity(sick, Vec3{nan, 0, 0}, Vec3{});
    for (int i = 0; i < 10; ++i) {
        cubesOnGround(w, ground, {sick, well});
        w.step(kDt);
    }
    CHECK_MESSAGE(!finite(w, sick), "the poisoned cube is poisoned: nothing repairs it yet");
    CHECK_MESSAGE(finite(w, well), "the cube beside it never saw it");
    CHECK_MESSAGE(w.velocity(ground) == Vec3{}, "and the ground they stand on was only read");
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
    run(w, 1.0);
    CHECK_MESSAGE((finite(w, a) && finite(w, b)), "both arms stay finite");
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
