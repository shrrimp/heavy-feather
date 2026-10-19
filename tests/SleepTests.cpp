#include <algorithm>
#include <cmath>
#include <vector>

#include <doctest/doctest.h>

#include "host/Host.hpp"
#include "host/Scenes.hpp"

using namespace hf;

namespace {

    constexpr double kDt = 1.0 / 60.0;

    void run(host::Host &h, double seconds) {
        for (int i = 0; i < int(std::lround(seconds / kDt)); ++i) h.step(kDt);
    }

    BodyId box(host::Host &h, const Vec3 &at, double half = 0.25, double mass = 5.0) {
        return h.addBox(at, Quat{}, Vec3{half, half, half}, mass);
    }

} // namespace

TEST_CASE("a box at rest falls asleep, and a push wakes it") {
    host::Host h;
    const BodyId b = box(h, Vec3{0, 0.3, 0});
    run(h, 2.0);
    CHECK_MESSAGE(h.world.asleep(b), "still for a few steps, it sleeps");
    CHECK_MESSAGE(h.world.velocity(b) == Vec3{}, "with no velocity left at all");
    CHECK_MESSAGE(h.world.stats().stepped == 0, "and costs nothing");
    const Vec3 at = h.world.position(b);
    run(h, 1.0);
    CHECK_MESSAGE(h.world.position(b) == at, "asleep, it doesn't drift by a bit");

    h.world.applyForce(b, Vec3{3000, 0, 0}, h.world.position(b));
    h.step(kDt);
    CHECK_MESSAGE(!h.world.asleep(b), "a push wakes it");
    CHECK_MESSAGE(h.world.velocity(b).x > 1.0, "and moves it");
}

TEST_CASE("a host writing back the pose it read doesn't wake anything") {
    host::Host h;
    const BodyId b = box(h, Vec3{0, 0.3, 0});
    run(h, 2.0);
    REQUIRE(h.world.asleep(b));
    for (int i = 0; i < 30; ++i) {
        h.world.setPose(b, h.world.position(b), h.world.orientation(b));
        h.world.setVelocity(b, h.world.velocity(b), h.world.angularVelocity(b));
        h.world.applyForce(b, Vec3{}, h.world.position(b));
        h.step(kDt);
    }
    CHECK_MESSAGE(h.world.asleep(b), "the same pose, no velocity and no force are no push");
    h.world.setPose(b, h.world.position(b) + Vec3{0, 1.0, 0}, h.world.orientation(b));
    CHECK_MESSAGE(!h.world.asleep(b), "a teleport is");
}

TEST_CASE("a ball landing on a sleeping stack steps the top of it, not the whole stack") {
    host::Host h;
    std::vector<BodyId> stack;
    for (int i = 0; i < 8; ++i) stack.push_back(box(h, Vec3{0, 0.25 + 0.5 * i, 0}));
    run(h, 3.0);
    for (BodyId b : stack) REQUIRE(h.world.asleep(b));

    const BodyId ball = h.addSphere(Vec3{0, 4.4, 0}, 0.1, 0.5);
    int most = 0;
    bool bottomWoke = false;
    for (int i = 0; i < 90; ++i) {
        h.step(kDt);
        most = std::max(most, h.world.stats().stepped);
        bottomWoke = bottomWoke || !h.world.asleep(stack[0]) || !h.world.asleep(stack[1]);
    }
    CHECK_MESSAGE(most >= 3, "the ball, the box it lands on and the one under that are stepped: " << most);
    CHECK_MESSAGE(most < 9, "a couple of contacts deep, not the whole stack");
    CHECK_MESSAGE(!bottomWoke, "the bottom of the stack sleeps through it");
    CHECK_MESSAGE(std::fabs(h.world.position(ball).y - 4.1) < 0.01, "and the ball rests on top");
    run(h, 3.0);
    CHECK_MESSAGE(h.world.stats().stepped == 0, "then everything sleeps again");
}

TEST_CASE("a machine sleeps whole, and its lever wakes it") {
    host::Host h;
    const host::Excavator e = host::addExcavator(h, Vec3{});
    run(h, 5.0);
    CHECK_MESSAGE(h.world.stats().stepped == 0, "parked, every link of it sleeps");
    const double boom = h.world.jointPosition(e.boom);
    h.world.setLever(e.boom, 1.0);
    run(h, 0.5);
    CHECK_MESSAGE(h.world.stats().stepped > 0, "the lever wakes it");
    CHECK_MESSAGE(h.world.jointPosition(e.boom) > boom + 0.05, "and the boom lifts");
    h.world.setLever(e.boom, 0.0);
    run(h, 5.0);
    CHECK_MESSAGE(h.world.stats().stepped == 0, "let go, it parks and sleeps again");
}

TEST_CASE("taking the support away drops what slept on it") {
    for (bool remove : {true, false}) {
        CAPTURE(remove);
        host::Host h;
        const BodyId under = box(h, Vec3{0, 0.25, 0}, 0.25, 50.0), top = box(h, Vec3{0, 0.75, 0});
        run(h, 2.0);
        REQUIRE(h.world.asleep(top));
        if (remove) h.removeBody(under);
        else h.world.setVelocity(under, Vec3{20.0, 0, 0}, Vec3{}); // yanked out from under it
        run(h, 1.5);
        CHECK_MESSAGE(std::fabs(h.world.position(top).y - 0.25) < 0.01, "it falls to the ground");
    }
}

TEST_CASE("sleeping can be turned off") {
    WorldSettings settings;
    settings.sleeping = false;
    host::Host h(settings);
    const BodyId b = box(h, Vec3{0, 0.3, 0});
    run(h, 2.0);
    CHECK_MESSAGE(!h.world.asleep(b), "stepped every step");
    CHECK(h.world.stats().stepped == 1);
}

TEST_CASE("bodies apart are islands apart: a joint makes one, the ground doesn't") {
    host::Host h;
    box(h, Vec3{0, 1.0, 0});
    box(h, Vec3{3, 1.0, 0});
    const BodyId c = box(h, Vec3{6, 1.0, 0}), d = box(h, Vec3{6.6, 1.0, 0});
    h.addJoint(h.jointAt(JointType::Ball, c, d, Vec3{6.3, 1.0, 0}, Vec3{1, 0, 0}));
    h.step(kDt);
    CHECK_MESSAGE(h.world.stats().islands == 3, "two boxes alone and a jointed pair");
    run(h, 0.5);
    CHECK_MESSAGE(h.world.stats().islands == 3, "landed on the same ground, still three");
    run(h, 3.0);
    CHECK_MESSAGE(h.world.stats().islands == 0, "asleep, none");
}
