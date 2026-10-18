#include <doctest/doctest.h>

#include "host/Allocations.hpp"
#include "host/Scenes.hpp"

using namespace host;

namespace {

    constexpr double kDt = 1.0 / 60.0;

#ifdef NDEBUG
    constexpr int kSteps = 4;
#else
    constexpr int kSteps = 1; // see StressTests
#endif

    // The host's own detection runs outside the count: only the library's step is measured.
    std::size_t stepAllocations(Scene &s, int i) {
        if (s.operate) s.operate(*s.host, i);
        s.host->detect();
        s.host->world.clearContacts();
        for (const hf::Contact &c : s.host->contacts) s.host->world.addContact(c);
        const std::size_t before = allocations().count;
        s.host->world.step(kDt);
        return allocations().count - before;
    }

} // namespace

TEST_CASE("a warm step allocates nothing") {
    for (auto build : {+[] { return pile(kSteps > 1 ? 120 : 40); }, +[] { return chain(kSteps > 1 ? 80 : 30); },
                       +[] { return scissorLift(kSteps > 1 ? 8 : 4); }}) {
        Scene s = build();
        CAPTURE(s.name);
        for (int i = 0; i < 60 * kSteps; ++i) s.step(i, kDt);
        std::size_t count = 0;
        for (int i = 60 * kSteps; i < 75 * kSteps; ++i) count += stepAllocations(s, i);
        CHECK_MESSAGE(count == 0, "storage made once is reused every step");
    }
}

TEST_CASE("a busy scene only allocates to grow") {
    Scene s = yard();
    for (int i = 0; i < 30 * kSteps; ++i) s.step(i, kDt);
    std::size_t count = 0;
    for (int i = 30 * kSteps; i < 105 * kSteps; ++i) count += stepAllocations(s, i);
    CHECK_MESSAGE(count < 40, "only growth to new peak sizes, a few dozen at most, never per step");
}

TEST_CASE("bodies and joints made and destroyed leave nothing behind") {
    Scene s = chain(10);
    Host &h = *s.host;
    const auto cycle = [&](int i) {
        const hf::BodyId box = h.addBox(hf::Vec3{0.5, 3.0, 0.5}, hf::Quat{}, hf::Vec3{0.2, 0.2, 0.2}, 5.0);
        const hf::BodyId ball = h.addSphere(hf::Vec3{-0.5, 3.0, 0.5}, 0.2, 5.0);
        h.addJoint(h.jointAt(hf::JointType::Hinge, box, ball, hf::Vec3{0, 3.0, 0.5}, hf::Vec3{0, 0, 1}));
        for (int k = 0; k < 4; ++k) s.step(i * 4 + k, kDt);
        h.removeBody(box);
        h.removeBody(ball);
        s.step(i * 4 + 4, kDt);
    };
    // Containers reach their working sizes first; after that, nothing may grow.
    for (int i = 0; i < 75 * kSteps; ++i) cycle(i);
    const std::size_t warm = allocations().live;
    for (int i = 75 * kSteps; i < 150 * kSteps; ++i) cycle(i);
    CHECK_MESSAGE(allocations().live <= warm, "churn doesn't grow the world");
}

TEST_CASE("a world gives back everything it took") {
    const std::size_t before = allocations().live;
    {
        Scene s = yard();
        for (int i = 0; i < 15 * kSteps; ++i) s.step(i, kDt);
    }
    CHECK_MESSAGE(allocations().live == before, "destroying a world frees all of it");
}
