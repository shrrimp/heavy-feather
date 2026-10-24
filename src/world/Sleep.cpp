#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <utility>
#include <vector>

#include "world/WorldImpl.hpp"

// Sleeping, as the prototype does it: a body still for a few steps falls asleep unless something it touches still
// moves, and each step the awake bodies flood out to whatever they could disturb; the rest stand still, as if static.
// A machine sleeps and wakes whole, and a drive at work keeps it awake.
namespace hf {

    namespace {

        constexpr int kStillSteps = 8;
        constexpr double kStillSpeed2 = 5e-4;      // m²/s²
        constexpr double kStillSpin2 = 1e-4;       // (ω × radius)², m²/s²
        constexpr double kStillCorrection = 1e-5;  // m the position pass may still be moving it
        constexpr int kMaxDepth = 2;               // contacts a wake travels through; joints don't count
        constexpr std::uint64_t kConductAge = 2;   // steps a pair still carries force after it last touched
        constexpr std::uint64_t kLostSupport = 2;  // steps apart before a sleeper's support is gone: one is flicker
        constexpr std::uint64_t kForgetAge = 4;

        std::uint64_t keyOf(BodyId id) { return (std::uint64_t(id.index) << 32) | id.generation; }

        std::uint64_t hashOf(BodyId a, BodyId b, std::uint64_t feature) {
            std::uint64_t h = keyOf(a) * 0x9E3779B97F4A7C15ull;
            h ^= (keyOf(b) + 0x632BE59BD9B4E019ull + (h << 6) + (h >> 2)) * 0xC2B2AE3D27D4EB4Full;
            h ^= (feature + 0x165667B19E3779F9ull + (h << 6) + (h >> 2)) * 0x94D049BB133111EBull;
            return h ^ (h >> 31);
        }

        // An open-addressing table of contact indices, sized for count entries at most half full.
        void clearTable(std::vector<int> &table, std::size_t count) {
            std::size_t size = 16;
            while (size < 2 * count) size *= 2;
            table.assign(size, -1);
        }

        // The slot holding the entry same() matches, or the empty slot it would go in.
        template <typename Same>
        std::size_t probe(const std::vector<int> &table, std::uint64_t hash, Same same) {
            const std::size_t mask = table.size() - 1;
            std::size_t i = std::size_t(hash) & mask;
            while (table[i] >= 0 && !same(table[i])) i = (i + 1) & mask;
            return i;
        }
        BodyId idOf(std::uint64_t key) { return {std::uint32_t(key >> 32), std::uint32_t(key)}; }

        bool sweeping(const Group &g) { return g.kinematic && (g.v != Vec3{} || g.w != Vec3{}); }

        Touch touchOf(BodyId x, BodyId y, std::uint64_t step) {
            const std::uint64_t a = keyOf(x), b = keyOf(y);
            return {std::min(a, b), std::max(a, b), step};
        }

        // A drive moving its machine, or pushing to: an operator's lever held, a setpoint on its way, a drive still
        // bringing itself back to it, a motor on.
        bool atWork(const JointSlot &j) {
            constexpr double kOff = 1e-4; // m or rad
            const JointDesc &d = j.desc;
            if (d.type == JointType::Motor) return d.target != 0.0;
            if (!positionDrive(d.type)) return false;
            const bool held = std::isfinite(d.maxForce) && j.locked;
            const double want = !held ? j.targetEff : j.commanded ? j.holdPos : d.target;
            const double off = d.type == JointType::Servo ? std::remainder(j.position - want, 2.0 * std::numbers::pi)
                                                          : j.position - want;
            return j.lever != 0.0 || j.targetFrom != d.target || std::fabs(off) > kOff;
        }

        void wakeNext(World::Impl &I, std::uint64_t key) {
            if (BodySlot *b = I.body(idOf(key))) b->wakeNext = true;
        }

        // Stepped this step: the awake; one hop more to a sleeper an awake body is closing on, as it's about to be
        // hit; then along every link that carries force, a few contacts deep but through any number of joints.
        void flood(World::Impl &I) {
            const std::size_t G = I.groups.size();
            I.reached = I.seeded;
            I.floodDepth.assign(G, 0);
            I.floodQueue.clear();
            for (std::size_t g = 0; g < G; ++g)
                if (I.seeded[g]) I.floodQueue.push_back(int(g));
            for (const SleepLink &l : I.sleepLinks) {
                if (!l.approaching || I.seeded[std::size_t(l.a)] == I.seeded[std::size_t(l.b)]) continue;
                const int sleeper = I.seeded[std::size_t(l.a)] ? l.b : l.a;
                if (I.reached[std::size_t(sleeper)]) continue;
                I.reached[std::size_t(sleeper)] = 1;
                I.floodQueue.push_back(sleeper);
            }
            for (std::size_t head = 0; head < I.floodQueue.size(); ++head) {
                const int g = I.floodQueue[head];
                const bool deep = I.floodDepth[std::size_t(g)] >= kMaxDepth;
                for (int e = I.linkStart[std::size_t(g)]; e < I.linkStart[std::size_t(g) + 1]; ++e) {
                    const SleepLink &l = I.sleepLinks[std::size_t(I.linkOf[std::size_t(e)])];
                    if (!l.conducts || (deep && !l.jointed)) continue;
                    const int o = l.a == g ? l.b : l.a;
                    if (I.reached[std::size_t(o)]) continue;
                    I.reached[std::size_t(o)] = 1;
                    I.floodDepth[std::size_t(o)] = I.floodDepth[std::size_t(g)] + (l.jointed ? 0 : 1);
                    I.floodQueue.push_back(o);
                }
            }
            I.steppedGroups.clear();
            for (std::size_t g = 0; g < G; ++g) {
                Group &gr = I.groups[g];
                gr.frozen = !gr.isStatic && !I.reached[g];
                if (moves(gr) || gr.kinematic) I.steppedGroups.push_back(int(g));
            }
            I.steppedMachines.clear();
            for (std::size_t mi = 0; mi < I.machines.size(); ++mi) {
                MachineData &md = I.machines[mi];
                md.frozen = I.groups[std::size_t(md.groupOfLink[0])].frozen;
                if (!md.frozen) I.steppedMachines.push_back(int(mi));
            }
        }

        // The stepped groups split into islands: what touches or is jointed to what. Each one is a solve of its own.
        void findIslands(World::Impl &I) {
            std::vector<int> &root = I.islandOf;
            root.resize(I.groups.size());
            for (std::size_t g = 0; g < root.size(); ++g) root[g] = int(g);
            const auto find = [&](int x) {
                while (root[std::size_t(x)] != x) x = root[std::size_t(x)] = root[std::size_t(root[std::size_t(x)])];
                return x;
            };
            for (const SleepLink &l : I.sleepLinks) {
                if (!(l.touching || l.jointed) || !I.reached[std::size_t(l.a)] || !I.reached[std::size_t(l.b)]) continue;
                // A kinematic body joins no island: nothing it touches can push it.
                if (I.groups[std::size_t(l.a)].isStatic || I.groups[std::size_t(l.b)].isStatic) continue;
                const int a = find(l.a), b = find(l.b);
                if (a != b) root[std::size_t(std::max(a, b))] = std::min(a, b);
            }
            I.islands = 0;
            for (std::size_t g = 0; g < root.size(); ++g) {
                Group &gr = I.groups[g];
                gr.island = -1;
                if (gr.isStatic || gr.frozen) continue;
                const int r = find(int(g));
                if (r == int(g)) gr.island = I.islands++;
                else gr.island = I.groups[std::size_t(r)].island;
            }
        }

    } // namespace

    void World::Impl::beginSleep() {
        ++stepCount;
        // What the world kept between sleepers is gone with a body that's gone, the host moved, or reshaped, whose
        // touches wake. A contact the host hands again is the host's fresh look at one it keeps: that takes it in and
        // keeps its impulses. Any other it kept for a pair the host handed is gone.
        if (!kept.empty()) {
            clearTable(handedByKey, contacts.size());
            clearTable(handedByPair, contacts.size());
            for (std::size_t i = 0; i < contacts.size(); ++i) {
                const Contact &in = contacts[i].in;
                const auto sameKey = [&](int j) {
                    const Contact &o = contacts[std::size_t(j)].in;
                    return o.a == in.a && o.b == in.b && o.feature == in.feature;
                };
                const auto samePair = [&](int j) {
                    const Contact &o = contacts[std::size_t(j)].in;
                    return o.a == in.a && o.b == in.b;
                };
                int &byKey = handedByKey[probe(handedByKey, hashOf(in.a, in.b, in.feature), sameKey)];
                if (byKey < 0) byKey = int(i);
                int &byPair = handedByPair[probe(handedByPair, hashOf(in.a, in.b, 0), samePair)];
                if (byPair < 0) byPair = int(i);
            }
        }
        std::size_t k = 0;
        for (std::size_t i = 0; i < kept.size(); ++i) {
            ContactData &c = kept[i];
            BodySlot *a = body(c.in.a), *b = body(c.in.b);
            if (!a || !b || a->teleported || b->teleported) continue;
            if (a->reshaped || b->reshaped) {
                a->wakeNext = b->wakeNext = true;
                continue;
            }
            const int byKey = handedByKey[probe(handedByKey, hashOf(c.in.a, c.in.b, c.in.feature), [&](int j) {
                const Contact &o = contacts[std::size_t(j)].in;
                return o.a == c.in.a && o.b == c.in.b && o.feature == c.in.feature;
            })];
            if (byKey >= 0) {
                ContactData &handed = contacts[std::size_t(byKey)];
                c.in = handed.in;
                handed.absorbed = true;
            } else if (handedByPair[probe(handedByPair, hashOf(c.in.a, c.in.b, 0), [&](int j) {
                           const Contact &o = contacts[std::size_t(j)].in;
                           return o.a == c.in.a && o.b == c.in.b;
                       })] >= 0) {
                continue;
            }
            if (k != i) kept[k] = c;
            ++k;
        }
        kept.resize(k);
        if (reshapedAny) {
            for (const ContactData &c : contacts) {
                BodySlot *a = body(c.in.a), *b = body(c.in.b);
                if (a && b && (a->reshaped || b->reshaped)) a->wakeNext = b->wakeNext = true;
            }
            for (BodySlot &b : bodies) b.reshaped = false;
            reshapedAny = false;
        }
        const std::size_t G = groups.size();
        seeded.assign(G, 0);
        for (BodySlot &b : bodies) {
            if (!b.alive) continue;
            if (b.wakeNext || !settings.sleeping) {
                b.asleep = b.wakeNext = false;
                b.stillSteps = 0;
            }
            if (!b.asleep && !groups[std::size_t(b.group)].isStatic) seeded[std::size_t(b.group)] = 1;
        }
        // A kinematic body on the move wakes what it reaches, as an awake body does.
        for (std::size_t g = 0; g < G; ++g)
            if (sweeping(groups[g])) seeded[g] = 1;
        const auto awake = [&](std::uint64_t key) {
            const BodySlot *b = body(idOf(key));
            return b && seeded[std::size_t(b->group)];
        };

        // The pairs touching now. One that touched lately, has an awake end so the host looked, and touches no
        // more has lost its support: both wake next step, which is soon enough for something that was resting.
        touchScratch.clear();
        // A contact the world keeps, handed again, is between sleepers: no support is lost there.
        for (const ContactData &c : contacts)
            if (!c.absorbed) touchScratch.push_back(touchOf(c.in.a, c.in.b, stepCount));
        std::sort(touchScratch.begin(), touchScratch.end());
        touchScratch.erase(std::unique(touchScratch.begin(), touchScratch.end(),
                                       [](const Touch &x, const Touch &y) { return !(x < y) && !(y < x); }),
                           touchScratch.end());
        const std::size_t now = touchScratch.size();
        bool keptPairsBuilt = false;
        for (const Touch &t : touches) {
            const std::uint64_t age = stepCount - t.step;
            if (age > kForgetAge || std::binary_search(touchScratch.begin(), touchScratch.begin() + long(now), t))
                continue;
            if (age >= kLostSupport && (awake(t.a) || awake(t.b))) {
                // A pair the world keeps still touches, though the host no longer hands it.
                if (!keptPairsBuilt) {
                    keptPairs.clear();
                    for (const ContactData &c : kept) {
                        const std::uint64_t a = keyOf(c.in.a), b = keyOf(c.in.b);
                        keptPairs.push_back({std::min(a, b), std::max(a, b)});
                    }
                    std::sort(keptPairs.begin(), keptPairs.end());
                    keptPairsBuilt = true;
                }
                if (std::binary_search(keptPairs.begin(), keptPairs.end(), std::pair{t.a, t.b})) {
                    touchScratch.push_back({t.a, t.b, stepCount});
                    continue;
                }
                wakeNext(*this, t.a);
                wakeNext(*this, t.b);
            }
            touchScratch.push_back(t);
        }
        touches.swap(touchScratch);
        std::sort(touches.begin(), touches.end());

        // Links between dynamic groups: the contacts, the pairs that touched a step or two ago, and the joints.
        sleepLinks.clear();
        const auto link = [&](const BodySlot *a, const BodySlot *b, SleepLink l) {
            if (!a || !b || a->group == b->group) return;
            const Group &ga = groups[std::size_t(a->group)], &gb = groups[std::size_t(b->group)];
            if ((ga.isStatic && !sweeping(ga)) || (gb.isStatic && !sweeping(gb))) return;
            l.a = a->group;
            l.b = b->group;
            sleepLinks.push_back(l);
        };
        for (const ContactData &c : contacts) {
            const BodySlot *a = body(c.in.a), *b = body(c.in.b);
            if (!a || !b || c.absorbed) continue; // its kept twin links them
            Vec3 ca, va, wa, cb, vb, wb;
            Quat q;
            bodyState(*a, ca, q, va, wa);
            bodyState(*b, cb, q, vb, wb);
            const Vec3 rel = vb + cross(wb, c.in.point - cb) - va - cross(wa, c.in.point - ca);
            link(a, b, {.touching = true, .approaching = dot(rel, c.in.normal) < 0.0, .conducts = true});
        }
        for (const Touch &t : touches)
            if (stepCount - t.step <= kConductAge && t.step != stepCount)
                link(body(idOf(t.a)), body(idOf(t.b)), {.conducts = true});
        for (const ContactData &c : kept) link(body(c.in.a), body(c.in.b), {.touching = true, .conducts = true});
        for (const JointSlot &j : joints)
            if (j.alive && j.desc.a.valid() && j.desc.b.valid())
                link(body(j.desc.a), body(j.desc.b), {.conducts = true, .jointed = true});

        linkStart.assign(G + 1, 0);
        for (const SleepLink &l : sleepLinks) {
            ++linkStart[std::size_t(l.a) + 1];
            ++linkStart[std::size_t(l.b) + 1];
        }
        for (std::size_t g = 0; g < G; ++g) linkStart[g + 1] += linkStart[g];
        linkOf.resize(std::size_t(linkStart[G]));
        floodQueue.assign(linkStart.begin(), linkStart.end() - 1); // write heads
        for (std::size_t i = 0; i < sleepLinks.size(); ++i) {
            linkOf[std::size_t(floodQueue[std::size_t(sleepLinks[i].a)]++)] = int(i);
            linkOf[std::size_t(floodQueue[std::size_t(sleepLinks[i].b)]++)] = int(i);
        }
        flood(*this);
        findIslands(*this);
    }

    // Per unit, a group or a whole machine: still this step counts on, moving wakes it. Still long enough, it falls
    // asleep unless something it touches or is jointed to is still moving; asleep, it keeps no velocity at all.
    void World::Impl::endSleep() {
        if (!settings.sleeping) return;
        const std::size_t G = groups.size(), U = G + machines.size();
        unitOf.assign(G, -1);
        for (std::size_t g = 1; g < G; ++g)
            if (!groups[g].isStatic) unitOf[g] = groups[g].machine >= 0 ? int(G) + groups[g].machine : int(g);
        unitStepped.assign(U, 0);
        unitStill.assign(U, 1);
        unitSteps.assign(U, INT_MAX);
        unitAsleep.assign(U, 1);
        unitVetoed.assign(U, 0);

        for (std::size_t g = 1; g < G; ++g) {
            const Group &gr = groups[g];
            if (unitOf[g] < 0 || gr.frozen) continue;
            const std::size_t u = std::size_t(unitOf[g]);
            unitStepped[u] = 1;
            // Put back after a fault, it is at rest, but not still: it stays awake to try again.
            const double spin = length(gr.w) * gr.radius;
            if (dot(gr.v, gr.v) >= kStillSpeed2 || spin * spin >= kStillSpin2 || gr.correction >= kStillCorrection ||
                faulted(gr))
                unitStill[u] = 0;
        }
        // A grabbed body is held by the host's hand: never still enough to sleep.
        for (const GrabSlot &s : grabs) {
            const BodySlot *b = s.alive ? body(s.body) : nullptr;
            if (b && b->group >= 0 && unitOf[std::size_t(b->group)] >= 0)
                unitStill[std::size_t(unitOf[std::size_t(b->group)])] = 0;
        }
        for (std::size_t mi = 0; mi < machines.size(); ++mi) {
            const MachineData &md = machines[mi];
            bool working = md.correction >= kStillCorrection;
            for (const JointRef &r : md.joints) working = working || atWork(joints[r.joint]);
            if (working) unitStill[G + mi] = 0;
        }
        for (const BodySlot &b : bodies) {
            if (!b.alive || unitOf[std::size_t(b.group)] < 0) continue;
            const std::size_t u = std::size_t(unitOf[std::size_t(b.group)]);
            unitSteps[u] = std::min(unitSteps[u], b.stillSteps);
            unitAsleep[u] = unitAsleep[u] && b.asleep;
        }
        for (std::size_t u = 0; u < U; ++u) {
            if (!unitStepped[u]) continue;
            if (unitStill[u]) {
                unitSteps[u] = std::min(unitSteps[u] + 1, 1 << 20);
            } else {
                unitSteps[u] = 0;
                unitAsleep[u] = 0;
            }
        }

        // The veto travels one link a step. Falling asleep late costs nothing; a step early is a visible freeze.
        for (const SleepLink &l : sleepLinks) {
            if (!l.touching && !l.jointed) continue;
            if (unitOf[std::size_t(l.a)] < 0 || unitOf[std::size_t(l.b)] < 0) continue;
            const std::size_t a = std::size_t(unitOf[std::size_t(l.a)]), b = std::size_t(unitOf[std::size_t(l.b)]);
            if (a == b) continue;
            if (unitSteps[a] == 0) unitVetoed[b] = 1;
            if (unitSteps[b] == 0) unitVetoed[a] = 1;
        }
        for (std::size_t u = 0; u < U; ++u)
            if (unitStepped[u] && unitSteps[u] >= kStillSteps && (unitAsleep[u] || !unitVetoed[u])) unitAsleep[u] = 1;
            else if (unitStepped[u]) unitAsleep[u] = 0;

        for (std::size_t g = 1; g < G; ++g) {
            Group &gr = groups[g];
            if (unitOf[g] < 0 || !unitStepped[std::size_t(unitOf[g])] || !unitAsleep[std::size_t(unitOf[g])]) continue;
            gr.v = gr.w = Vec3{};
        }
        for (BodySlot &b : bodies) {
            if (!b.alive || unitOf[std::size_t(b.group)] < 0) continue;
            const std::size_t u = std::size_t(unitOf[std::size_t(b.group)]);
            if (!unitStepped[u]) continue;
            b.stillSteps = unitSteps[u];
            b.asleep = unitAsleep[u];
        }
    }

    void World::Impl::wakePartners(BodyId id) {
        const std::uint64_t key = keyOf(id);
        for (const Touch &t : touches)
            if (t.a == key || t.b == key) {
                wakeNext(*this, t.a);
                wakeNext(*this, t.b);
            }
        for (const ContactData &c : kept)
            if (c.in.a == id || c.in.b == id) {
                wakeNext(*this, keyOf(c.in.a));
                wakeNext(*this, keyOf(c.in.b));
            }
    }

} // namespace hf
