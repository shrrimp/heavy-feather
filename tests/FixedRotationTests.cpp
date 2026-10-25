#include <cmath>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"
#include "host/Host.hpp"

using namespace hf;
using namespace build;

namespace {

    WorldSettings weightless() {
        WorldSettings s;
        s.gravity = Vec3{};
        return s;
    }

    // How far a body has turned from where it started, in radians.
    double turned(const Quat &from, const Quat &to) {
        return 2.0 * std::acos(std::fmin(1.0, std::fabs(dot(from, to))));
    }

    // A tall box standing with its centre of mass past a ledge's edge: free, it tips off.
    double tipAfter(bool fixed) {
        host::Host h;
        h.addBox(Vec3{0, 0.5, 0}, Quat{}, Vec3{1, 0.5, 1}, 0.0);
        const BodyId post = h.addBox(Vec3{1.1, 2.0, 0}, Quat{}, Vec3{0.3, 1.0, 0.3}, 10.0);
        h.world.setFixedRotation(post, fixed);
        for (int s = 0; s < 120; ++s) h.step(kDt);
        return turned(Quat{}, h.world.orientation(post));
    }

} // namespace

TEST_CASE("a box past a ledge tips off it, and with its rotation fixed stands upright on the edge") {
    CHECK_MESSAGE(tipAfter(false) > 0.5, "free, it tips");
    CHECK_MESSAGE(tipAfter(true) < 1e-9, "fixed, it does not turn at all");
}

TEST_CASE("a body with its rotation fixed is pushed off its centre and moves without turning") {
    World w(weightless());
    BodyDesc d;
    d.mass = 2.0;
    d.inertia = boxInertia(2.0, Vec3{0.5, 0.5, 0.5});
    d.angularVelocity = Vec3{0, 3, 0};
    d.fixedRotation = true;
    const BodyId b = w.addBody(d);
    CHECK_MESSAGE(length(w.angularVelocity(b)) == 0.0, "the spin it was given is dropped");
    for (int s = 0; s < 60; ++s) {
        w.applyForce(b, Vec3{4, 0, 0}, w.centerOfMass(b) + Vec3{0, 0.5, 0.5});
        w.applyTorque(b, Vec3{0, 0, 10});
        w.step(kDt);
    }
    CHECK_MESSAGE(turned(Quat{}, w.orientation(b)) == 0.0, "no force or torque turns it");
    CHECK_MESSAGE(std::fabs(w.velocity(b).x - 2.0) < 1e-9, "the push moves it as it moves any body: F t / m");
}

TEST_CASE("its rotation let go, a body turns again from rest") {
    World w(weightless());
    BodyDesc d;
    d.inertia = Mat3(1.0);
    d.fixedRotation = true;
    const BodyId b = w.addBody(d);
    w.applyTorque(b, Vec3{0, 1, 0});
    w.step(kDt);
    CHECK(length(w.angularVelocity(b)) == 0.0);
    w.setFixedRotation(b, false);
    for (int s = 0; s < 60; ++s) {
        w.applyTorque(b, Vec3{0, 1, 0});
        w.step(kDt);
    }
    CHECK_MESSAGE(std::fabs(w.angularVelocity(b).y - 1.0) < 1e-9, "a second of a unit torque on a unit inertia");
}
