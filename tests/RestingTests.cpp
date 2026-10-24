#include <cmath>
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

    const Vec3 kHalf{0.5, 0.5, 0.5};

    // The host that only hands a pair with an awake end: the ground is static, so it never counts as one.
    void handAwake(host::Host &h) {
        h.detect();
        h.world.clearContacts();
        const auto awake = [&](BodyId b) { return b != h.ground() && !h.world.asleep(b); };
        for (const Contact &k : h.contacts)
            if (awake(k.a) || awake(k.b)) h.world.addContact(k);
    }

    struct Pile {
        host::Host h;
        std::vector<BodyId> boxes;

        // Three columns of three boxes, settled until every one sleeps, every pair handed all along.
        Pile() {
            for (int col = 0; col < 3; ++col)
                for (int row = 0; row < 3; ++row)
                    boxes.push_back(h.addBox(Vec3{1.01 * col, 0.5 + 1.0 * row, 0}, Quat{}, kHalf, 10.0));
            for (int i = 0; i < 600 && !asleep(); ++i) h.step(kDt);
        }
        bool asleep() {
            for (BodyId b : boxes)
                if (!h.world.asleep(b)) return false;
            return true;
        }
    };

    std::vector<double> warmStarts;
    void readWarmStarts(World::Impl &I, int substep, FaultPhase phase) {
        if (substep != 0 || phase != FaultPhase::Internal) return;
        warmStarts.clear();
        for (const ContactData &c : I.contacts) warmStarts.push_back(c.ln);
    }

    // What the world keeps this step: what it kept from the last, and what it takes in now, merged once it's done.
    std::size_t keptCount = 0;
    double keptImpulse = 0.0;
    void readKept(World::Impl &I, int, FaultPhase phase) {
        if (phase != FaultPhase::EndOfStep) return;
        keptCount = I.kept.size() + I.newlyKept.size();
        keptImpulse = 0.0;
        for (const ContactData &c : I.kept) keptImpulse += c.ln;
        for (const ContactData &c : I.newlyKept) keptImpulse += c.ln;
    }

} // namespace

TEST_CASE("a sleeping pile stands on the contacts the world kept once the host stops handing them") {
    Pile p;
    REQUIRE(p.asleep());
    std::vector<Vec3> at;
    for (BodyId b : p.boxes) at.push_back(p.h.world.position(b));
    for (int i = 0; i < 300; ++i) {
        handAwake(p.h);
        p.h.world.step(kDt);
    }
    CHECK_MESSAGE(p.asleep(), "five seconds on, handed nothing, it still sleeps");
    for (std::size_t i = 0; i < p.boxes.size(); ++i)
        CHECK_MESSAGE((p.h.world.position(p.boxes[i]) == at[i]), "exactly where it was");

    // A ball dropped on the middle column: what it reaches wakes, and stands on what the world kept for it.
    const BodyId ball = p.h.addSphere(Vec3{1.01, 3.5, 0}, 0.3, 5.0);
    for (int i = 0; i < 120; ++i) {
        handAwake(p.h);
        p.h.world.step(kDt);
    }
    CHECK_MESSAGE(p.h.world.position(ball).y > 3.0 + 0.3 - 0.05, "the ball rests on the pile");
    for (std::size_t i = 0; i < p.boxes.size(); ++i)
        CHECK_MESSAGE(length(p.h.world.position(p.boxes[i]) - at[i]) < 0.02, "and nothing sank into what it stood on");
}

TEST_CASE("a body woken among sleepers solves its contacts warm") {
    Pile p;
    REQUIRE(p.asleep());
    for (int i = 0; i < 10; ++i) {
        handAwake(p.h);
        p.h.world.step(kDt);
    }
    // The middle box of the middle column, woken: its contacts above and below were kept between sleepers.
    p.h.world.wake(p.boxes[4]);
    testing::phaseHook = readWarmStarts;
    handAwake(p.h);
    p.h.world.step(kDt);
    testing::phaseHook = nullptr;
    REQUIRE(!warmStarts.empty());
    int warmed = 0;
    for (double l : warmStarts) warmed += l > 0.0;
    CHECK_MESSAGE(warmed >= 4, "the contacts under it and over it start with what held them");
}

TEST_CASE("a pair the host hands again replaces what the world kept, and starts from it") {
    Pile p;
    REQUIRE(p.asleep());
    testing::phaseHook = readKept;
    for (int i = 0; i < 3; ++i) {
        handAwake(p.h);
        p.h.world.step(kDt);
    }
    const std::size_t keptAlone = keptCount;
    const double heldAlone = keptImpulse;
    CHECK_MESSAGE(keptAlone > 0, "the world keeps the sleeping pile's contacts");
    CHECK_MESSAGE(heldAlone > 0.0, "with what held them");
    // A host that hands every pair still: each replaces its kept self, and keeps what it held.
    for (int i = 0; i < 3; ++i) p.h.step(kDt);
    testing::phaseHook = nullptr;
    CHECK_MESSAGE(keptCount == p.h.contacts.size(), "each pair handed kept once, none twice");
    CHECK_MESSAGE(std::fabs(keptImpulse - heldAlone) < 1e-9 * heldAlone, "nor lost what held it");
    CHECK_MESSAGE(p.asleep(), "and the pile sleeps on");
}
