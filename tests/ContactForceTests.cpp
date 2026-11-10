#include <cmath>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "host/Host.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kHalf{0.5, 0.5, 0.5};

    // What the contacts on a body push it with, summed: a takes the opposite of each.
    Vec3 pushOn(const World &w, BodyId body) {
        Vec3 sum{};
        for (const ContactForce &f : w.contactForces()) {
            if (f.b == body) sum += f.force;
            if (f.a == body) sum -= f.force;
        }
        return sum;
    }

    void settle(host::Host &h, double seconds) {
        for (int s = 0; s < int(std::lround(seconds / kDt)); ++s) h.step(kDt);
    }

} // namespace

TEST_CASE("a box at rest stands on its weight, awake and asleep") {
    host::Host h;
    const BodyId box = h.addBox(Vec3{0, 0.5, 0}, Quat{}, kHalf, 40.0);
    settle(h, 0.1);
    REQUIRE_FALSE(h.world.asleep(box));
    const double weight = 40.0 * 9.81;
    CHECK(pushOn(h.world, box).y == doctest::Approx(weight).epsilon(0.02));
    CHECK(std::fabs(pushOn(h.world, box).x) < 0.01 * weight);

    for (int i = 0; i < 600 && !h.world.asleep(box); ++i) h.step(kDt);
    REQUIRE(h.world.asleep(box));
    h.step(kDt);
    CHECK_MESSAGE(pushOn(h.world, box).y == doctest::Approx(weight).epsilon(0.02), "the kept contacts still hold it");
}

TEST_CASE("a hinged pair on the ground carries both weights through its contacts") {
    host::Host h;
    const BodyId a = h.addBox(Vec3{0, 0.5, 0}, Quat{}, kHalf, 30.0);
    const BodyId b = h.addBox(Vec3{1.02, 0.5, 0}, Quat{}, kHalf, 50.0);
    h.addJoint(h.jointAt(JointType::Hinge, a, b, Vec3{0.51, 0.5, 0}, Vec3{0, 0, 1}));
    settle(h, 1.0);
    const Vec3 total = pushOn(h.world, a) + pushOn(h.world, b);
    CHECK(total.y == doctest::Approx(80.0 * 9.81).epsilon(0.02));
}

TEST_CASE("a box resting on another passes its weight down") {
    host::Host h;
    const BodyId low = h.addBox(Vec3{0, 0.5, 0}, Quat{}, kHalf, 20.0);
    const BodyId top = h.addBox(Vec3{0, 1.5, 0}, Quat{}, kHalf, 10.0);
    settle(h, 1.0);
    CHECK(pushOn(h.world, top).y == doctest::Approx(10.0 * 9.81).epsilon(0.03));
    CHECK(pushOn(h.world, low).y == doctest::Approx(20.0 * 9.81).epsilon(0.03));
}
