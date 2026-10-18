#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "host/Host.hpp"

using namespace hf;

namespace {

    constexpr double kDt = 1.0 / 60.0;
    constexpr double kPi = std::numbers::pi;

    // A direct-acting valve train with gravity off, as in an engine where the spring dwarfs the valve's weight: an
    // eccentric round cam on a motor lifts a flat follower on a slider, and a spring holds the follower down on it.
    // Everything hangs off a heavy free head, so the cam, the follower and their contact are one machine.
    struct ValveTrain {
        static constexpr double e = 0.05;       // m: the cam's offset, half the lift
        static constexpr double R = 0.10;       // m: the cam's radius
        static constexpr double half = 0.02;    // m: the follower's half height
        static constexpr double preload = 0.05; // m: the spring's compression on the base circle
        static constexpr double hertz = 2.0;    // the spring's frequency on the follower

        host::Host h;
        BodyId cam, follower;

        // The follower stays on the lobe while the spring can decelerate it at the top, e ω² ≤ (2πf)² (δ0 + 2e).
        static double floatSpeed() { return 2.0 * kPi * hertz * std::sqrt((preload + 2.0 * e) / e); }

        explicit ValveTrain(double omega) : h(settings()) {
            const Vec3 axle{0, 1, 0};
            BodyDesc head;
            head.position = Vec3{0, 1.2, 0};
            head.mass = 1e4;
            head.inertia = host::boxInertia(1e4, Vec3{0.5, 0.5, 0.5});
            const BodyId frame = h.world.addBody(head);

            BodyDesc c;
            c.position = axle;
            c.mass = 4.0;
            c.inertia = Mat3(0.5 * 4.0 * R * R);
            c.angularVelocity = Vec3{0, 0, omega};
            cam = h.world.addBody(c);
            host::Shape lobe;
            lobe.body = cam;
            lobe.kind = host::ShapeKind::Sphere;
            lobe.radius = R;
            lobe.offset = Vec3{e, 0, 0};
            lobe.mechanism = true;
            h.shapes.push_back(lobe);
            JointDesc motor = h.jointAt(JointType::Motor, frame, cam, axle, Vec3{0, 0, 1});
            motor.target = omega;
            h.addJoint(motor);

            const Vec3 at = axle + Vec3{0, R + half, 0};
            BodyDesc f;
            f.position = at;
            f.mass = 0.2;
            f.inertia = host::boxInertia(0.2, Vec3{0.08, half, 0.03});
            f.velocity = Vec3{0, e * omega, 0};
            follower = h.world.addBody(f);
            host::Shape tappet;
            tappet.body = follower;
            tappet.half = Vec3{0.08, half, 0.03};
            tappet.mechanism = true;
            h.shapes.push_back(tappet);
            h.addJoint(h.jointAt(JointType::Slider, frame, follower, at, Vec3{0, 1, 0}));

            // Lift runs from -e (the base circle) to +e: the spring's anchor sits 0.3 m above the follower.
            JointDesc spring = h.jointAt(JointType::Spring, frame, follower, at, Vec3{0, 1, 0});
            spring.frameA.position += Vec3{0, 0.3, 0};
            spring.hertz = hertz;
            spring.damping = 0.0;
            spring.restLength = 0.3 + e + preload;
            h.addJoint(spring);
        }

        // A lobe turning several degrees a frame outruns contacts detected once a frame: the host detects at the
        // sub-step rate instead, as a real valve train's simulation must.
        static WorldSettings settings() {
            WorldSettings s;
            s.gravity = Vec3{};
            s.substeps = 1;
            return s;
        }

        // How far the follower's face is off the lobe's top: positive is daylight between them.
        double gap() const {
            const Vec3 lobeCentre = h.world.position(cam) + h.world.orientation(cam) * Vec3{e, 0, 0};
            return (h.world.position(follower).y - half) - (lobeCentre.y + R);
        }
    };

} // namespace


TEST_CASE("a valve floats at the cam speed its spring predicts") {
    // The widest daylight between the follower and the lobe over a few turns, once the start has settled.
    const auto widestGap = [](double speed) {
        ValveTrain v(speed);
        constexpr double dt = kDt / 16.0;
        double most = 0.0;
        for (int i = 0; i < 960; ++i) {
            v.h.step(dt);
            if (i >= 240) most = std::max(most, v.gap());
        }
        return most;
    };
    const double wf = ValveTrain::floatSpeed();
    CHECK_MESSAGE(widestGap(0.95 * wf) < 1e-4, "5% under, the spring keeps it on the lobe");
    CHECK_MESSAGE(widestGap(1.05 * wf) > 1e-3, "5% over, it leaves the lobe at the top");
}

namespace {

    // A wheel on a motor and a brake pad weighing down on it from a slider, both on one frame standing on the
    // ground: the pad's weight is the clamp, and the motor's cap the torque it fights.
    struct Brake {
        static constexpr double R = 0.2, pad = 20.0, mu = 0.5;
        host::Host h;
        BodyId wheel;
        JointId motor;

        static double holdingTorque() { return mu * pad * 9.81 * R; }

        explicit Brake(double torque) {
            const BodyId frame = h.addBox(Vec3{0, 0.1, 0}, Quat{}, Vec3{1.0, 0.1, 1.0}, 1000.0);
            const Vec3 axle{0, 0.6, 0};
            wheel = h.addSphere(axle, R, 10.0);
            h.shapes.back().friction = mu;
            h.shapes.back().mechanism = true;
            JointDesc m = h.jointAt(JointType::Motor, frame, wheel, axle, Vec3{0, 0, 1});
            m.speed = 1.0;
            m.maxForce = torque;
            motor = h.addJoint(m);
            const Vec3 at = axle + Vec3{0, R + 0.05, 0};
            const BodyId shoe = h.addBox(at, Quat{}, Vec3{0.1, 0.05, 0.1}, pad, mu);
            h.shapes.back().mechanism = true;
            h.addJoint(h.jointAt(JointType::Slider, frame, shoe, at, Vec3{0, 1, 0}));
        }
    };

} // namespace

TEST_CASE("a brake holds its wheel until the torque passes mu times its clamp") {
    const auto spin = [](double torque) {
        Brake b(torque);
        for (int i = 0; i < 30; ++i) b.h.step(kDt);
        b.h.world.setLever(b.motor, 1.0);
        for (int i = 0; i < 60; ++i) b.h.step(kDt);
        return b.h.world.angularVelocity(b.wheel).z;
    };
    const double hold = Brake::holdingTorque();
    CHECK_MESSAGE(std::fabs(spin(0.97 * hold)) < 1e-3, "3% under, the pad holds it");
    CHECK_MESSAGE(spin(1.03 * hold) > 0.9, "3% over, it slips and the motor takes it to speed");
}

TEST_CASE("a cart whose driven wheels spin pulls away at mu g, straight") {
    // The wheels' contacts with the ground solved with the cart's joints. Sliding, friction is mu times the push and
    // points against the slip: no more, and nothing sideways.
    host::Host h;
    h.shapes[0].mechanism = true;
    const BodyId cart = h.addBox(Vec3{0, 0.6, 0}, Quat{}, Vec3{1.0, 0.15, 0.6}, 100.0);
    std::vector<JointId> motors;
    for (double x : {-0.8, 0.8})
        for (double z : {-0.75, 0.75}) {
            const BodyId wheel = h.addSphere(Vec3{x, 0.3, z}, 0.3, 5.0);
            h.shapes.back().mechanism = true;
            JointDesc m = h.jointAt(JointType::Motor, cart, wheel, Vec3{x, 0.3, z}, Vec3{0, 0, 1});
            m.speed = 10.0;
            m.maxForce = 200.0; // four times what the grip can pass on
            motors.push_back(h.addJoint(m));
        }
    for (int i = 0; i < 30; ++i) h.step(kDt);
    for (JointId m : motors) h.world.setLever(m, -1.0);
    h.step(kDt);
    const double v0 = h.world.velocity(cart).x;
    h.step(kDt);
    h.step(kDt);
    const double a = (h.world.velocity(cart).x - v0) / (2.0 * kDt), mug = std::sqrt(0.6 * 0.8) * 9.81;
    CHECK_MESSAGE(std::fabs(a - mug) < 0.01 * mug, "Coulomb's law, to the percent: " << a << " against " << mug);
    for (int i = 0; i < 120; ++i) h.step(kDt);
    CHECK_MESSAGE(std::fabs(h.world.velocity(cart).x - 3.0) < 1e-3, "then rolls at the wheels' speed");
    CHECK_MESSAGE(std::fabs(h.world.position(cart).z) < 0.01, "on the line it started on");
}

TEST_CASE("a mechanism contact bounces like any other") {
    // A ball on a vertical slider dropped 1 m onto the ground, the contact solved with the slider's machine.
    host::Host h;
    h.shapes[0].mechanism = true;
    const BodyId ball = h.addSphere(Vec3{0, 1.2, 0}, 0.2, 5.0, 0.5);
    h.shapes.back().mechanism = true;
    h.addJoint(h.jointAt(JointType::Slider, BodyId{}, ball, Vec3{0, 1.2, 0}, Vec3{0, 1, 0}));
    bool landed = false;
    double peak = 0.0;
    for (int i = 0; i < 180; ++i) {
        h.step(kDt);
        if (h.world.velocity(ball).y > 0.0) landed = true;
        if (landed) peak = std::max(peak, h.world.position(ball).y - 0.2);
    }
    CHECK_MESSAGE(landed, "it lands");
    CHECK_MESSAGE(std::fabs(peak - 0.25) < 0.08, "a 1 m drop at e = 0.5 rebounds about 0.25 m: " << peak);
}
