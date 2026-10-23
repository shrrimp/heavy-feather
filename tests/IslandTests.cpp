#include <cmath>
#include <vector>

#include <doctest/doctest.h>

#include "Builders.hpp"
#include "hf/World.hpp"
#include "host/Host.hpp"

using namespace hf;
using namespace build;

namespace {

    const Vec3 kHalf{0.5, 0.5, 0.5};

    bool allAsleep(host::Host &h, const std::vector<BodyId> &bodies) {
        for (BodyId b : bodies)
            if (!h.world.asleep(b)) return false;
        return true;
    }

    // Steps until every body sleeps, or the time is up; the seconds it took, or more than the time if it never did.
    double secondsToSleep(host::Host &h, const std::vector<BodyId> &bodies, double limit) {
        for (int i = 1; i <= int(std::lround(limit / kDt)); ++i) {
            h.step(kDt);
            if (allAsleep(h, bodies)) return i * kDt;
        }
        return 2.0 * limit;
    }

} // namespace

TEST_CASE("a stack of ten boxes on the rigid budget settles and sleeps within three seconds") {
    host::Host h;
    std::vector<BodyId> boxes;
    for (int i = 0; i < 10; ++i) boxes.push_back(h.addBox(Vec3{0, 0.5 + 1.0 * i, 0}, Quat{}, kHalf, 10.0));
    const double t = secondsToSleep(h, boxes, 3.0);
    CHECK_MESSAGE(t <= 3.0, "asleep within three seconds");
    const Vec3 top = h.world.position(boxes.back());
    CHECK_MESSAGE((std::fabs(top.x) < 0.01 && std::fabs(top.z) < 0.01), "standing straight");
    CHECK_MESSAGE(std::fabs(top.y - 9.5) < 0.05, "and as tall as it was built");
}

TEST_CASE("a pyramid of fifty-five boxes on the rigid budget settles and sleeps within three seconds") {
    host::Host h;
    std::vector<BodyId> boxes;
    for (int row = 0; row < 10; ++row)
        for (int i = 0; i < 10 - row; ++i)
            boxes.push_back(h.addBox(Vec3{1.05 * (i - 0.5 * (9 - row)), 0.5 + 1.0 * row, 0}, Quat{}, kHalf, 10.0));
    REQUIRE(boxes.size() == 55);
    const double t = secondsToSleep(h, boxes, 3.0);
    CHECK_MESSAGE(t <= 3.0, "asleep within three seconds");
    CHECK_MESSAGE(std::fabs(h.world.position(boxes.back()).y - 9.5) < 0.1, "with its top box still on top");
}

TEST_CASE("an island holding a machine keeps the machine budget, and the rows say which is which") {
    host::Host h;
    // A pendulum hanging onto a box, and a pile of boxes well away: one island of each kind.
    const BodyId arm = h.addBox(Vec3{0, 2.0, 0}, Quat{}, Vec3{0.1, 1.0, 0.1}, 5.0);
    h.addJoint(h.jointAt(JointType::Hinge, BodyId{}, arm, Vec3{0, 3.0, 0}, Vec3{0, 0, 1}));
    h.addBox(Vec3{0, 0.5, 0}, Quat{}, kHalf, 10.0);
    for (int i = 0; i < 3; ++i) h.addBox(Vec3{8, 0.5 + 1.0 * i, 0}, Quat{}, kHalf, 10.0);
    h.step(kDt);
    const WorldStats s = h.world.stats();
    CHECK_MESSAGE(s.islands == 2, "two islands");
    CHECK_MESSAGE(s.machineRows > 0, "the machine's swept rows");
    CHECK_MESSAGE(s.rigidRows > 0, "and the pile's");
}

TEST_CASE("a step's profile is off unless asked for, and adds up when it is") {
    host::Host h;
    for (int i = 0; i < 5; ++i) h.addBox(Vec3{0, 0.5 + 1.0 * i, 0}, Quat{}, kHalf, 10.0);
    h.step(kDt);
    CHECK_MESSAGE(h.world.stats().profile.total == 0.0, "nothing timed by default");
    WorldSettings settings = h.world.settings();
    settings.profile = true;
    h.world.setSettings(settings);
    h.step(kDt);
    const StepProfile p = h.world.stats().profile;
    CHECK_MESSAGE(p.total > 0.0, "a step takes some time");
    CHECK_MESSAGE(std::fabs(p.begin + p.machines + p.contacts + p.positions + p.integrate + p.relax + p.end - p.total) <
                      1e-9,
                  "and its phases add up to it");
}
