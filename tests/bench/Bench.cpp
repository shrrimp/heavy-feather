// The cost of a step on the large stress scenes: the world's step alone, the host's detection left out. Budgets are a
// few times what a desktop CPU takes, unthrottled: a laptop saving power can run several times slower.
//   hf-bench            prints the table
//   hf-bench --check    also fails when a scene's mean is over its budget
//   hf-bench <name>     only the scenes whose name contains it
//   --steps N           times N steps instead of each scene's own count
//   --trace             prints every step that takes over twice the mean so far

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <vector>

#include "host/Scenes.hpp"

using namespace host;

namespace {

    struct Case {
        const char *name;
        std::function<Scene()> build;
        int steps;
        double budgetMs; // the mean a step may take
    };

    const std::vector<Case> &cases() {
        static const std::vector<Case> list = {
            {"pile 600", [] { return pile(600); }, 600, 35.0},
            {"chain 500 piled", [] { return chain(500); }, 300, 1000.0},
            {"ring 200", [] { return ring(200); }, 300, 160.0},
            {"branches 511", [] { return branches(9); }, 300, 12.0},
            {"scissor lift 30", [] { return scissorLift(30); }, 300, 35.0},
            {"yard", [] { return yard(); }, 600, 40.0},
        };
        return list;
    }

    double ms(std::chrono::steady_clock::duration d) { return std::chrono::duration<double, std::milli>(d).count(); }

} // namespace

int main(int argc, char **argv) {
    bool check = false;
    const char *only = nullptr;
    int stepsOverride = 0;
    bool trace = false;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--check")) check = true;
        else if (!std::strcmp(argv[i], "--steps") && i + 1 < argc) stepsOverride = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--trace")) trace = true;
        else only = argv[i];
    }

    constexpr double kDt = 1.0 / 60.0;
    int failed = 0;
    std::printf("%-18s %8s %8s %8s %8s %8s %9s\n", "scene", "bodies", "contacts", "mean ms", "p95 ms", "max ms",
                "budget");
    for (const Case &c : cases()) {
        if (only && !std::strstr(c.name, only)) continue;
        Scene s = c.build();
        Host &h = *s.host;
        for (int i = 0; i < 60; ++i) s.step(i, kDt);

        std::vector<double> times;
        const int steps = stepsOverride > 0 ? stepsOverride : c.steps;
        for (int i = 60; i < 60 + steps; ++i) {
            if (s.operate) s.operate(h, i);
            h.detect();
            h.world.clearContacts();
            for (const hf::Contact &k : h.contacts) h.world.addContact(k);
            const auto t0 = std::chrono::steady_clock::now();
            h.world.step(kDt);
            times.push_back(ms(std::chrono::steady_clock::now() - t0));
            if (trace && times.size() > 5) {
                double sum = 0.0;
                for (double t : times) sum += t;
                const hf::WorldStats st = h.world.stats();
                if (times.back() > 2.0 * sum / double(times.size()))
                    std::printf("  step %d: %.3f ms, contacts %d, rows %d, pivots %d, fallbacks %d\n", i, times.back(),
                                st.contacts, st.rows, st.lcpPivots, st.lcpFallbacks);
            }
        }
        double mean = 0.0;
        for (double t : times) mean += t;
        mean /= double(times.size());
        std::sort(times.begin(), times.end());
        const double p95 = times[std::size_t(0.95 * double(times.size() - 1))], worst = times.back();
        const hf::WorldStats st = h.world.stats();
        const bool over = mean > c.budgetMs;
        failed += over;
        std::printf("%-18s %8d %8d %8.3f %8.3f %8.3f %8.1f%s\n", c.name, st.bodies, st.contacts, mean, p95, worst,
                    c.budgetMs, over ? "  OVER" : "");
    }
    return check && failed ? 1 : 0;
}
