#include <cmath>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "host/Scenes.hpp"

using namespace host;

namespace {

    constexpr double kDt = 1.0 / 60.0;

    // Optimised builds run the scenes at full length; debug and sanitizer builds a shorter run of smaller scenes,
    // still every scene and every check.
#ifdef NDEBUG
    constexpr int kSteps = 4;
#else
    constexpr int kSteps = 1;
#endif

    struct Case {
        const char *name;
        std::function<Scene()> build;
    };

    // Sized to run in a debug build; the benchmark runs the large ones.
    const std::vector<Case> &cases() {
        static const std::vector<Case> list = {
            {"pile", [] { return pile(kSteps > 1 ? 150 : 60); }},
            {"chain", [] { return chain(kSteps > 1 ? 120 : 40); }},
            {"ring", [] { return ring(kSteps > 1 ? 48 : 24); }},
            {"branches", [] { return branches(kSteps > 1 ? 6 : 4); }},
            {"scissor lift", [] { return scissorLift(kSteps > 1 ? 10 : 5); }},
            {"yard", [] { return yard(); }},
            {"ropes", [] { return ropes(kSteps > 1 ? 6 : 2); }},
        };
        return list;
    }

    // Every body's pose and velocity, in shape order.
    std::vector<double> state(const Host &h) {
        std::vector<double> s;
        for (const Shape &shape : h.shapes) {
            if (!h.world.contains(shape.body)) continue;
            const hf::Vec3 p = h.world.position(shape.body), v = h.world.velocity(shape.body);
            const hf::Vec3 w = h.world.angularVelocity(shape.body);
            const hf::Quat q = h.world.orientation(shape.body);
            s.insert(s.end(), {p.x, p.y, p.z, q.w, q.x, q.y, q.z, v.x, v.y, v.z, w.x, w.y, w.z});
        }
        return s;
    }

    bool sameBits(const std::vector<double> &a, const std::vector<double> &b) {
        return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
    }

} // namespace

TEST_CASE("the same scene stepped twice is the same, bit for bit") {
    for (const Case &c : cases()) {
        CAPTURE(std::string(c.name));
        Scene a = c.build(), b = c.build();
        // Stepped in turns, so one world can't lean on anything the other left behind.
        bool same = true;
        int firstDiff = -1;
        for (int i = 0; i < 45 * kSteps && same; ++i) {
            a.step(i, kDt);
            b.step(i, kDt);
            same = sameBits(state(*a.host), state(*b.host));
            if (!same) firstDiff = i;
        }
        CAPTURE(firstDiff);
        CHECK_MESSAGE(same, "two runs of one scene never part");
    }
}

TEST_CASE("no scene blows up") {
    for (const Case &c : cases()) {
        CAPTURE(std::string(c.name));
        Scene s = c.build();
        double worstSpeed = 0.0, farthest = 0.0;
        bool finite = true;
        for (int i = 0; i < 150 * kSteps; ++i) {
            s.step(i, kDt);
            const std::vector<double> st = state(*s.host);
            for (std::size_t k = 0; k < st.size(); k += 13) {
                for (std::size_t j = 0; j < 13; ++j) finite &= std::isfinite(st[k + j]);
                farthest = std::max(farthest, std::hypot(st[k], st[k + 1], st[k + 2]));
                worstSpeed = std::max(worstSpeed, std::hypot(st[k + 7], st[k + 8], st[k + 9]));
            }
        }
        CAPTURE(worstSpeed);
        CAPTURE(farthest);
        CHECK_MESSAGE(finite, "every pose and velocity stays a number");
        CHECK_MESSAGE(farthest < 200.0, "nothing is flung away");
        CHECK_MESSAGE(worstSpeed < 50.0, "nothing moves faster than anything in the scene could make it");
    }
}
