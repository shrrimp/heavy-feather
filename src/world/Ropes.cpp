#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include "math/Dense.hpp"
#include "world/Checks.hpp"
#include "world/WorldImpl.hpp"

// Ropes: one row each, on the rope's slack, s = L - span >= 0, where L is what is paid out and span its length along
// the route as it lies. It only ever pulls, and it pulls every point of its route toward its neighbours at once: a
// sheave feels both strands that run round it, and the tension is the same all along the rope. Nothing in it grows
// with the rope's length.
//
// A rope never joins two machines into one. When every moving point of it is on one machine it is a row of that
// machine's solve, held exactly with the joints (MachineRows.cpp). Otherwise it is swept with the contacts, as a loose
// spring is, and the ropes of an island are solved together and exactly, so four slings on a light hook share a heavy
// load whatever the masses. A stretched rope is pulled back in the position pass; one with a spring is soft instead.
namespace hf {

    namespace {

        constexpr double kInf = std::numeric_limits<double>::infinity();
        constexpr double kPi = std::numbers::pi;
        constexpr double kRopeBeta = 0.8;         // fraction of a stretch corrected per sub-step, as a joint's error
        constexpr double kMaxPosBias = 1000.0;    // m/s: a pose-only ceiling, it never becomes momentum
        constexpr double kMinSpringLength = 1e-3; // m: a stretchy rope shorter than this is as stiff as this one
        constexpr double kSlipping = 1.0 - 1e-6;  // a spool pulling this much of its cap is overloaded
        constexpr double kRelEps = 1e-9;          // a rigid rope's compliance relative to its own mobility

        void wakeUp(BodySlot *b) {
            if (!b) return;
            b->asleep = false;
            b->stillSteps = 0;
        }

        bool nonNegative(double x, bool infinite = false) { return x >= 0.0 && (infinite || std::isfinite(x)); }

        bool usable(const RopeDesc &d) {
            return nonNegative(d.length) && nonNegative(d.stiffness, true) && d.stiffness > 0.0 &&
                   nonNegative(d.damping) && nonNegative(d.strength, true) && nonNegative(d.speed) &&
                   nonNegative(d.maxForce, true) && nonNegative(d.brakeForce, true) && nonNegative(d.minLength) &&
                   nonNegative(d.maxLength, true) && d.minLength <= d.maxLength && d.spool <= RopeSpool::Both;
        }

        bool usable(const RopePoint &p) {
            return finite(p.position) && finite(p.axis) && nonNegative(p.radius) &&
                   (p.radius == 0.0 || length(p.axis) > 1e-9) && p.side >= -1 && p.side <= 1;
        }

        int spools(const RopeDesc &d) { return d.spool == RopeSpool::Both ? 2 : d.spool == RopeSpool::None ? 0 : 1; }

        double withinSpools(const RopeDesc &d, double length) {
            return spools(d) ? std::clamp(length, d.minLength, d.maxLength) : length;
        }

        // Where a rope running round a circle of radius r about c, in the plane square to n, touches it: arriving from
        // a, or leaving toward a. Side +1 goes round anticlockwise about n, -1 clockwise. A point inside the circle
        // has no tangent: the rim's nearest point stands in.
        Vec3 tangentPoint(const Vec3 &a, const Vec3 &c, const Vec3 &n, double r, double side, bool arriving) {
            const Vec3 d = (a - c) - dot(a - c, n) * n;
            const double D = length(d);
            if (!(D > r)) return D > 0.0 ? c + (r / D) * d : c;
            const Vec3 e = d / D, f = cross(n, e);
            const double cosPhi = r / D, sinPhi = std::sqrt(std::max(0.0, 1.0 - cosPhi * cosPhi));
            return c + r * (cosPhi * e + (arriving ? side : -side) * sinPhi * f);
        }

        // How far round from where the rope arrives to where it leaves, the way it goes, in (-π, π].
        double arcAngle(const Vec3 &in, const Vec3 &out, const Vec3 &c, const Vec3 &n, double side) {
            const Vec3 a = in - c, b = out - c;
            return side * std::atan2(dot(n, cross(a, b)), dot(a, b));
        }

        bool sheave(const RopePoint &p) { return p.radius > 0.0; }

        // The rope through its route, each bend's p and n where its point is: each sheave's tangent points and angle,
        // which way the span grows as each point moves, and the span. A sheave starts from its last answer, and its
        // angle stays continuous; one unwound past nothing is off, and the rope runs straight past it. The first time,
        // a sheave's angle is taken in [-π/2, 3π/2): a half turn round it, as a block's, is a half turn.
        double laySpan(const std::vector<RopePoint> &route, std::vector<RopeBendState> &bends, bool first) {
            const int count = int(route.size());
            for (int k = 0; k < count; ++k) {
                RopeBendState &b = bends[std::size_t(k)];
                b.dir = Vec3{};
                if (!sheave(route[std::size_t(k)])) {
                    b.in = b.out = b.p;
                    b.touching = true;
                } else if (first) {
                    b.in = b.out = b.p;
                }
            }
            const auto prevOut = [&](int k) {
                for (int j = k - 1; j >= 0; --j)
                    if (bends[std::size_t(j)].touching) return bends[std::size_t(j)].out;
                return bends[std::size_t(k)].p;
            };
            const auto nextIn = [&](int k) {
                for (int j = k + 1; j < count; ++j)
                    if (bends[std::size_t(j)].touching) return bends[std::size_t(j)].in;
                return bends[std::size_t(k)].p;
            };
            // Neighbouring sheaves depend on each other's tangent points: two passes, each from the last.
            for (int pass = 0; pass < 2; ++pass)
                for (int k = 0; k < count; ++k) {
                    const RopePoint &rp = route[std::size_t(k)];
                    if (!sheave(rp)) continue;
                    RopeBendState &b = bends[std::size_t(k)];
                    const double side = double(rp.side);
                    // A drum at an end: the rope leaves it, or arrives at it, along one tangent.
                    if (k == 0 || k == count - 1) {
                        b.in = b.out = k == 0 ? tangentPoint(nextIn(k), b.p, b.n, rp.radius, side, false)
                                              : tangentPoint(prevOut(k), b.p, b.n, rp.radius, side, true);
                        b.touching = true;
                        continue;
                    }
                    const Vec3 in = tangentPoint(prevOut(k), b.p, b.n, rp.radius, side, true);
                    const Vec3 out = tangentPoint(nextIn(k), b.p, b.n, rp.radius, side, false);
                    const double raw = arcAngle(in, out, b.p, b.n, side);
                    double theta = raw;
                    if (first && pass == 0) theta = raw < -0.5 * kPi ? raw + 2.0 * kPi : raw;
                    else theta = raw + 2.0 * kPi * std::round((b.theta - raw) / (2.0 * kPi));
                    b.in = in;
                    b.out = out;
                    b.theta = theta;
                    b.touching = theta > 0.0;
                }

            double span = 0.0;
            int prev = 0;
            Vec3 u{};
            for (int k = 1; k < count; ++k) {
                RopeBendState &b = bends[std::size_t(k)];
                if (!b.touching) continue;
                const Vec3 seg = b.in - bends[std::size_t(prev)].out;
                const double len = length(seg);
                if (len > 1e-12) u = seg / len; // a segment of no length runs on the way the last one did
                span += len;
                bends[std::size_t(prev)].dir -= u;
                b.dir += u;
                prev = k;
            }
            for (int k = 1; k + 1 < count; ++k) {
                const RopeBendState &b = bends[std::size_t(k)];
                if (sheave(route[std::size_t(k)]) && b.touching) span += route[std::size_t(k)].radius * b.theta;
            }
            return span;
        }

        // A sheave left to choose its side goes round the way the route bends there; a drum at an end, anticlockwise.
        void settleSides(std::vector<RopePoint> &route, const std::vector<RopeBendState> &bends) {
            const int count = int(route.size());
            for (int k = 0; k < count; ++k) {
                RopePoint &rp = route[std::size_t(k)];
                if (!sheave(rp) || rp.side != 0) continue;
                rp.side = 1;
                if (k == 0 || k == count - 1) continue;
                const RopeBendState &a = bends[std::size_t(k - 1)], &b = bends[std::size_t(k)],
                                    &c = bends[std::size_t(k + 1)];
                if (dot(b.n, cross(b.p - a.p, c.p - b.p)) < 0.0) rp.side = -1;
            }
        }

    } // namespace

    // A rope's velocity row for slack s, the spools paying out at their rate: a rigid rope lets the slack close within
    // the sub-step and no further, a stretchy one is a spring on its stretch, as stiff as EA over its length.
    void ropeSlackRow(const RopeSlot &r, double k, double h, double &gamma, double &target) {
        const double s = r.length - r.span;
        if (std::isfinite(r.desc.stiffness)) {
            const double ks = r.desc.stiffness / std::max(r.length, kMinSpringLength);
            const double hertz = k > 0.0 ? std::sqrt(ks * k) / (2.0 * kPi) : 0.0;
            double stretched = 0.0;
            springRow(k, hertz, r.desc.damping, std::min(s, 0.0), h, gamma, stretched);
            target = (s < 0.0 ? stretched : -s / h) - r.lengthRate;
            return;
        }
        gamma = kRelEps * k;
        target = -std::max(s, 0.0) / h - r.lengthRate;
    }

    // The position pass's target for a rigid rope: a stretch pulled back a share at a time, a slack one free to close.
    // An overloaded spool is paying out, not stretching.
    double positionTarget(const RopeSlot &r, double h) {
        const double s = r.length - r.span;
        if (s >= 0.0) return -s / h;
        return r.slipping ? 0.0 : std::min(kRopeBeta * -s / h, kMaxPosBias);
    }

    namespace {

        double velocityOf(const World::Impl &I, const RopePrep &p, bool pseudo) {
            double v = 0.0;
            for (int e = p.endBegin; e < p.endBegin + p.endCount; ++e) {
                const RopeEnd &end = I.ropeEnds[std::size_t(e)];
                const Group &g = I.groups[std::size_t(end.group)];
                v += pseudo ? dot(end.lin, g.pv) + dot(end.ang, g.pw) : dot(end.lin, g.v) + dot(end.ang, g.w);
            }
            for (int k = p.blockBegin; k < p.blockBegin + p.blockCount; ++k) {
                const ContactBlock &blk = I.ropeBlocks[std::size_t(k)];
                const Machine &m = I.machines[std::size_t(blk.machine)].m;
                v += dot(I.arena.data() + blk.off, pseudo ? m.pseudo().data() : m.art.qdot.data(), std::size_t(blk.n));
            }
            return v;
        }

        void applyRope(World::Impl &I, const RopePrep &p, double dl, bool pseudo) {
            if (dl == 0.0) return;
            for (int e = p.endBegin; e < p.endBegin + p.endCount; ++e) {
                const RopeEnd &end = I.ropeEnds[std::size_t(e)];
                if (!end.write) continue;
                Group &g = I.groups[std::size_t(end.group)];
                (pseudo ? g.pv : g.v) += dl * end.imLin;
                (pseudo ? g.pw : g.w) += dl * end.imAng;
            }
            for (int k = p.blockBegin; k < p.blockBegin + p.blockCount; ++k) {
                const ContactBlock &blk = I.ropeBlocks[std::size_t(k)];
                I.machines[std::size_t(blk.machine)].m.applyContact(
                    I.arena.data() + blk.off + std::size_t(blk.n), dl, pseudo);
            }
        }

        // How an impulse on rope b moves rope a's row: through the groups and machines they share.
        double coupling(const World::Impl &I, const RopePrep &a, const RopePrep &b) {
            double s = 0.0;
            for (int i = a.endBegin; i < a.endBegin + a.endCount; ++i)
                for (int j = b.endBegin; j < b.endBegin + b.endCount; ++j) {
                    const RopeEnd &x = I.ropeEnds[std::size_t(i)], &y = I.ropeEnds[std::size_t(j)];
                    if (x.group == y.group && y.write) s += dot(x.lin, y.imLin) + dot(x.ang, y.imAng);
                }
            for (int i = a.blockBegin; i < a.blockBegin + a.blockCount; ++i)
                for (int j = b.blockBegin; j < b.blockBegin + b.blockCount; ++j) {
                    const ContactBlock &x = I.ropeBlocks[std::size_t(i)], &y = I.ropeBlocks[std::size_t(j)];
                    if (x.machine == y.machine)
                        s += dot(I.arena.data() + x.off, I.arena.data() + y.off + std::size_t(y.n), std::size_t(x.n));
                }
            return s;
        }

    } // namespace

    RopeSlot *World::Impl::rope(RopeId id) {
        if (id.index >= ropes.size()) return nullptr;
        RopeSlot &r = ropes[id.index];
        return r.alive && r.generation == id.generation ? &r : nullptr;
    }

    const RopeSlot *World::Impl::rope(RopeId id) const { return const_cast<Impl *>(this)->rope(id); }

    void World::Impl::wakeRope(const RopeSlot &r) {
        for (const RopePoint &p : r.route) wakeUp(body(p.body));
    }

    // Once a step: which ropes are stepped, each one's points on their groups, and how it is solved. A rope with
    // nothing moving on it waits, keeping its tension; one inside one machine is that machine's; the rest are swept,
    // and make their islands machine islands.
    void World::Impl::beginRopes() {
        steppedRopes.clear();
        brokenRopes.clear();
        for (MachineData &md : machines) md.ropes.clear();
        for (std::size_t ri = 0; ri < ropes.size(); ++ri) {
            RopeSlot &r = ropes[ri];
            if (!r.alive) continue;
            r.machine = -2;
            r.row = -1;
            if (r.broken) continue;
            int machine = -2, island = -1;
            bool swept = false;
            for (std::size_t k = 0; k < r.route.size(); ++k) {
                const RopePoint &rp = r.route[k];
                RopeBendState &b = r.bends[k];
                const BodySlot *s = body(rp.body);
                if (s && s->teleported) r.placed = false; // a sheave's angle can't follow a jump
                b.group = s ? s->group : 0;
                b.local = s ? s->offset + s->rot * (rp.position - s->desc.centerOfMass) : rp.position;
                b.axis = s ? s->rot * rp.axis : rp.axis;
                const Group &g = groups[std::size_t(b.group)];
                if (!moves(g)) continue;
                if (island < 0) island = g.island;
                const int m = machineOf(g);
                if (m < 0) swept = true;
                else if (machine == -2) machine = m;
                else if (machine != m) swept = true;
            }
            if (island < 0) continue;
            r.machine = r.held = swept ? -1 : machine;
            r.island = island;
            r.length0 = r.length;
            r.peak = 0.0;
            if (r.machine >= 0) machines[std::size_t(r.machine)].ropes.push_back(int(ri));
            else machineIsland[std::size_t(island)] = 1;

            // The spools: a rate within their stock, the motor's cap while they turn, the brake's while they don't.
            // Paid out to the end, the rope is tied to the drum.
            const RopeDesc &d = r.desc;
            r.lengthRate = 0.0;
            r.cap = kInf;
            if (spools(d)) {
                double rate = r.lever * d.speed * spools(d);
                if ((rate > 0.0 && r.length >= d.maxLength) || (rate < 0.0 && r.length <= d.minLength)) rate = 0.0;
                r.lengthRate = rate;
                const bool tied = r.length >= d.maxLength;
                r.cap = rate < 0.0 ? d.maxForce : tied ? kInf : rate > 0.0 ? d.maxForce : d.brakeForce;
            }
            steppedRopes.push_back(int(ri));
        }
    }

    // A rope as it lies, from its groups' poses. A spool that slipped has paid out what the load took.
    void World::Impl::layRope(RopeSlot &r) {
        for (RopeBendState &b : r.bends) {
            const Group &g = groups[std::size_t(b.group)];
            b.p = g.c + g.q * b.local;
            b.n = normalize(g.q * b.axis);
        }
        const bool first = !r.placed;
        if (first) settleSides(r.route, r.bends);
        r.span = laySpan(r.route, r.bends, first);
        r.placed = true;
        if (r.slipping) r.length = withinSpools(r.desc, std::max(r.length, r.span));
    }

    // Each sub-step, before anything is solved: every stepped rope as it lies.
    void World::Impl::measureRopes() {
        for (int ri : steppedRopes) layRope(ropes[std::size_t(ri)]);
    }

    // The span of a rope with md's links where the line search tried them: the rest stay where they are.
    double World::Impl::ropeSpanIn(const RopeSlot &r, const MachineData &md) {
        ropeScratch = r.bends;
        for (RopeBendState &b : ropeScratch) {
            const Group &g = groups[std::size_t(b.group)];
            const bool link = g.machine >= 0 && &machines[std::size_t(g.machine)] == &md;
            const Vec3 c = link ? md.m.art.links[std::size_t(g.link)].c : g.c;
            const Quat q = link ? md.m.art.links[std::size_t(g.link)].R : g.q;
            b.p = c + q * b.local;
            b.n = normalize(q * b.axis);
        }
        return laySpan(r.route, ropeScratch, false);
    }

    // Each sub-step, after the springs: the rows of the ropes no one machine holds, warm from the last sub-step's
    // tension, island by island.
    void World::Impl::prepareRopes() {
        ropePreps.clear();
        ropeEnds.clear();
        ropeBlocks.clear();
        for (int ri : steppedRopes) {
            RopeSlot &r = ropes[std::size_t(ri)];
            if (r.machine != -1) continue;
            RopePrep p;
            p.rope = ri;
            p.island = r.island;
            p.endBegin = int(ropeEnds.size());
            p.blockBegin = int(ropeBlocks.size());
            for (const RopeBendState &b : r.bends) {
                if (!b.touching || b.dir == Vec3{}) continue;
                const Group &g = groups[std::size_t(b.group)];
                const int mi = machineOf(g);
                if (mi >= 0) {
                    ContactBlock *blk = nullptr;
                    for (int k = p.blockBegin; k < int(ropeBlocks.size()); ++k)
                        if (ropeBlocks[std::size_t(k)].machine == mi) blk = &ropeBlocks[std::size_t(k)];
                    Machine &m = machines[std::size_t(mi)].m;
                    if (!blk) {
                        ropeBlocks.push_back({mi, arena.size(), m.dofs()});
                        blk = &ropeBlocks.back();
                        arena.resize(blk->off + 2 * std::size_t(blk->n), 0.0);
                    }
                    m.art.addPointRow(g.link, b.p, b.dir, -1.0, arena.data() + blk->off);
                    continue;
                }
                // A group that doesn't move is only read: a kinematic one's motion counts, the ground's is nothing.
                if (!moves(g) && !g.kinematic) continue;
                RopeEnd *end = nullptr;
                for (int e = p.endBegin; e < int(ropeEnds.size()); ++e)
                    if (ropeEnds[std::size_t(e)].group == b.group) end = &ropeEnds[std::size_t(e)];
                if (!end) {
                    ropeEnds.push_back({});
                    end = &ropeEnds.back();
                    end->group = b.group;
                    end->write = moves(g);
                }
                end->lin -= b.dir;
                end->ang -= cross(b.p - g.c, b.dir);
            }
            p.endCount = int(ropeEnds.size()) - p.endBegin;
            p.blockCount = int(ropeBlocks.size()) - p.blockBegin;
            for (int e = p.endBegin; e < p.endBegin + p.endCount; ++e) {
                RopeEnd &end = ropeEnds[std::size_t(e)];
                if (!end.write) continue;
                const Group &g = groups[std::size_t(end.group)];
                end.imLin = g.invMass * end.lin;
                end.imAng = g.invInertiaWorld() * end.ang;
                p.k += dot(end.lin, end.imLin) + dot(end.ang, end.imAng);
            }
            for (int k = p.blockBegin; k < p.blockBegin + p.blockCount; ++k) {
                const ContactBlock &blk = ropeBlocks[std::size_t(k)];
                double *g = arena.data() + blk.off, *w = g + blk.n;
                const double unheld = machines[std::size_t(blk.machine)].m.mobility(g, w);
                p.k += std::max(dot(g, w, std::size_t(blk.n)), kLockedMobility * unheld);
            }
            ropeSlackRow(r, p.k, h, p.gamma, p.target);
            // Nothing on it can move along it: a rope between two anchored hinges' axles, say.
            if (!(p.k + p.gamma > 0.0)) {
                ropeEnds.resize(std::size_t(p.endBegin));
                ropeBlocks.resize(std::size_t(p.blockBegin));
                continue;
            }
            p.hi = r.cap * h;
            p.posTarget = positionTarget(r, h);
            p.lambda = std::clamp(r.force * h, 0.0, p.hi);
            applyRope(*this, p, p.lambda, false);
            ropePreps.push_back(p);
        }

        // Island by island, and each island's block of couplings, made once a sub-step.
        const std::size_t L = std::size_t(islands);
        ropeStart.assign(L + 1, 0);
        for (const RopePrep &p : ropePreps) ++ropeStart[std::size_t(p.island) + 1];
        for (std::size_t i = 1; i <= L; ++i) ropeStart[i] += ropeStart[i - 1];
        ropeSorted.resize(ropePreps.size());
        heads.assign(ropeStart.begin(), ropeStart.end() - 1);
        for (const RopePrep &p : ropePreps) ropeSorted[std::size_t(heads[std::size_t(p.island)]++)] = p;
        ropePreps.swap(ropeSorted);
        ropeBlockStart.assign(L + 1, 0);
        for (std::size_t i = 0; i < L; ++i) {
            const std::size_t n = std::size_t(ropeStart[i + 1] - ropeStart[i]);
            ropeBlockStart[i + 1] = ropeBlockStart[i] + (n > 1 ? n * n : 0);
        }
        ropeS.resize(ropeBlockStart[L]);
        for (std::size_t i = 0; i < L; ++i) {
            const int b = ropeStart[i], n = ropeStart[i + 1] - b;
            if (n < 2) continue;
            double *S = ropeS.data() + ropeBlockStart[i];
            for (int x = 0; x < n; ++x)
                for (int y = 0; y < n; ++y)
                    S[x * n + y] = x == y ? ropePreps[std::size_t(b + x)].k
                                          : coupling(*this, ropePreps[std::size_t(b + x)], ropePreps[std::size_t(b + y)]);
        }
    }

    // One sweep of an island's swept ropes, on the velocity or on the position pass's pseudo-velocity. One rope is a
    // row; several are one block, solved exactly given what everything else in the island is doing.
    void World::Impl::sweepRopes(int island, bool pseudo) {
        if (ropeStart.empty()) return;
        const int begin = ropeStart[std::size_t(island)], n = ropeStart[std::size_t(island) + 1] - begin;
        if (n == 0) return;
        // The position pass pulls back a stretched rigid rope; a stretchy one's spring does that itself.
        const auto rigid = [&](const RopePrep &p) { return !std::isfinite(ropes[std::size_t(p.rope)].desc.stiffness); };
        if (n == 1) {
            RopePrep &p = ropePreps[std::size_t(begin)];
            if (pseudo && !rigid(p)) return;
            const double v = velocityOf(*this, p, pseudo);
            double &l = pseudo ? p.pseudo : p.lambda;
            const double nl = pseudo ? std::max(0.0, l + (p.posTarget - v) / (p.k + p.gamma))
                                     : std::clamp(l + (p.target - v - p.gamma * l) / (p.k + p.gamma), 0.0, p.hi);
            const double dl = nl - l;
            l = nl;
            applyRope(*this, p, dl, pseudo);
            return;
        }

        // w = (S + Γ) x + (v - S λ - target): free rows meet their target, rows at 0 open, rows at their cap pull it.
        const double *S = ropeS.data() + ropeBlockStart[std::size_t(island)];
        const std::size_t m = std::size_t(n);
        ropeA.assign(S, S + m * m);
        ropeB.resize(m);
        ropeLo.assign(m, 0.0);
        ropeHi.resize(m);
        ropeX.resize(m);
        ropeV.resize(m);
        ropeState.resize(m);
        for (std::size_t i = 0; i < m; ++i) ropeV[i] = velocityOf(*this, ropePreps[std::size_t(begin) + i], pseudo);
        for (std::size_t i = 0; i < m; ++i) {
            const RopePrep &p = ropePreps[std::size_t(begin) + i];
            const double l = pseudo ? p.pseudo : p.lambda;
            double Sl = 0.0;
            for (std::size_t j = 0; j < m; ++j)
                Sl += S[i * m + j] * (pseudo ? ropePreps[std::size_t(begin) + j].pseudo
                                             : ropePreps[std::size_t(begin) + j].lambda);
            ropeA[i * m + i] += p.gamma;
            ropeB[i] = ropeV[i] - Sl - (pseudo ? p.posTarget : p.target);
            ropeHi[i] = pseudo ? (rigid(p) ? kInf : 0.0) : p.hi;
            ropeX[i] = l;
            ropeState[i] = l <= 0.0 ? RowState::AtLo : l >= ropeHi[i] ? RowState::AtHi : RowState::Free;
        }
        solveBoxLcp(ropeA, ropeB, ropeLo, ropeHi, n, ropeX, ropeState, ropeLcp);
        for (std::size_t i = 0; i < m; ++i) {
            RopePrep &p = ropePreps[std::size_t(begin) + i];
            double &l = pseudo ? p.pseudo : p.lambda;
            const double nl = std::isfinite(ropeX[i]) ? ropeX[i] : l;
            const double dl = nl - l;
            l = nl;
            applyRope(*this, p, dl, pseudo);
        }
    }

    // A machine's ropes' tensions, from its solve.
    void World::Impl::recordRopes(MachineData &md) {
        const std::vector<double> &lambda = md.m.lambda();
        for (int ri : md.ropes) {
            RopeSlot &r = ropes[std::size_t(ri)];
            if (r.row >= 0 && std::size_t(r.row) < lambda.size()) r.force = lambda[std::size_t(r.row)] / h;
        }
    }

    // Each sub-step, its velocities solved: every stepped rope's tension, and what its spools did. A spool at its cap
    // is slipping: it pays out what the load takes, measured once the rope has moved. One that reaches the end of its
    // stock stops there: wound in, its brake holds; paid out, the rope is tied to the drum.
    void World::Impl::afterRopes() {
        for (const RopePrep &p : ropePreps) ropes[std::size_t(p.rope)].force = p.lambda / h;
        for (int ri : steppedRopes) {
            RopeSlot &r = ropes[std::size_t(ri)];
            r.peak = std::max(r.peak, r.force);
            const RopeDesc &d = r.desc;
            if (!spools(d)) continue;
            r.slipping = std::isfinite(r.cap) && r.force >= kSlipping * r.cap;
            r.length = withinSpools(d, r.length + r.lengthRate * h);
            if (r.lengthRate > 0.0 && r.length >= d.maxLength) {
                r.lengthRate = 0.0;
                r.cap = kInf;
            } else if (r.lengthRate < 0.0 && r.length <= d.minLength) {
                r.lengthRate = 0.0;
                r.cap = d.brakeForce;
            }
        }
    }

    // Once the step is done: every rope as the step left it, for the host to read, a rope with nothing
    // moving on it too. And a rope pulled harder than it is strong breaks, and pulls nothing from the
    // next step on.
    void World::Impl::endRopes() {
        for (RopeSlot &r : ropes)
            if (r.alive && !r.broken) layRope(r);
        for (int ri : steppedRopes) {
            RopeSlot &r = ropes[std::size_t(ri)];
            if (r.broken || !(r.peak > r.desc.strength)) continue;
            r.broken = true;
            r.force = 0.0;
            r.slipping = false;
            brokenRopes.push_back({std::uint32_t(ri), r.generation});
            wakeRope(r);
        }
    }

    RopeId World::addRope(const RopeDesc &desc, std::span<const RopePoint> route) {
        Impl &I = *impl_;
        for (const RopePoint &p : route)
            if (p.body.valid() && !I.body(p.body)) return {};
        if (route.size() < 2 || !usable(desc)) return I.reject(RopeId{});
        for (const RopePoint &p : route)
            if (!usable(p)) return I.reject(RopeId{});
        std::uint32_t index;
        if (!I.freeRopes.empty()) {
            index = I.freeRopes.back();
            I.freeRopes.pop_back();
        } else {
            index = std::uint32_t(I.ropes.size());
            I.ropes.emplace_back();
        }
        RopeSlot &r = I.ropes[index];
        RopeSlot fresh;
        fresh.generation = r.generation;
        r = std::move(fresh);
        r.alive = true;
        r.desc = desc;
        r.route.assign(route.begin(), route.end());
        for (RopePoint &p : r.route)
            if (sheave(p)) p.axis = normalize(p.axis);
        r.bends.assign(route.size(), RopeBendState{});
        r.length = withinSpools(desc, desc.length);
        I.wakeRope(r);
        return {index, r.generation};
    }

    double World::routeSpan(std::span<const RopePoint> route, std::span<RopeBend> path) const {
        const Impl &I = *impl_;
        if (route.size() < 2) return 0.0;
        std::vector<RopePoint> points(route.begin(), route.end());
        std::vector<RopeBendState> bends(route.size());
        for (std::size_t k = 0; k < points.size(); ++k) {
            RopePoint &p = points[k];
            if (!usable(p) || (p.body.valid() && !I.body(p.body))) return 0.0;
            if (sheave(p)) p.axis = normalize(p.axis);
            Vec3 c{}, v, w;
            Quat q{};
            if (const BodySlot *b = I.body(p.body)) {
                I.bodyState(*b, c, q, v, w);
                bends[k].p = c + q * (p.position - b->desc.centerOfMass);
            } else {
                bends[k].p = p.position;
            }
            bends[k].n = q * p.axis;
        }
        settleSides(points, bends);
        const double span = laySpan(points, bends, true);
        for (std::size_t k = 0; k < std::min(path.size(), bends.size()); ++k)
            path[k] = {bends[k].in, bends[k].out, bends[k].touching};
        return span;
    }

    void World::removeRope(RopeId id) {
        Impl &I = *impl_;
        RopeSlot *r = I.rope(id);
        if (!r) return;
        I.wakeRope(*r);
        r->alive = false;
        ++r->generation;
        r->route.clear();
        r->bends.clear();
        I.freeRopes.push_back(id.index);
    }

    bool World::contains(RopeId id) const { return impl_->rope(id) != nullptr; }

    bool World::setRope(RopeId id, const RopeDesc &desc) {
        Impl &I = *impl_;
        RopeSlot *r = I.rope(id);
        if (!r) return false;
        if (!usable(desc)) return I.reject(false);
        if (desc.length != r->desc.length) r->length = desc.length;
        r->desc = desc;
        r->length = withinSpools(desc, r->length);
        r->broken = false;
        I.wakeRope(*r);
        return true;
    }

    RopeDesc World::ropeDesc(RopeId id) const {
        const RopeSlot *r = impl_->rope(id);
        return r ? r->desc : RopeDesc{};
    }

    bool World::setRopeRoute(RopeId id, std::span<const RopePoint> route) {
        Impl &I = *impl_;
        RopeSlot *r = I.rope(id);
        if (!r) return false;
        for (const RopePoint &p : route)
            if (p.body.valid() && !I.body(p.body)) return false;
        if (route.size() < 2) return I.reject(false);
        for (const RopePoint &p : route)
            if (!usable(p)) return I.reject(false);
        // A point that stays keeps what it was: a sheave its tangents and its angle, which a host
        // adding and dropping the points a rope catches on, step after step, must not reset.
        const auto same = [](const RopePoint &a, const RopePoint &b) {
            return a.body == b.body && a.position == b.position && a.radius == b.radius &&
                   (a.radius == 0.0 || (a.axis == normalize(b.axis) && (b.side == 0 || a.side == b.side)));
        };
        std::vector<RopeBendState> bends(route.size());
        std::vector<RopePoint> points(route.begin(), route.end());
        bool changed = route.size() != r->route.size(), kept = true;
        std::size_t from = 0;
        for (std::size_t k = 0; k < points.size(); ++k) {
            if (sheave(points[k])) points[k].axis = normalize(points[k].axis);
            std::size_t found = r->route.size();
            for (std::size_t j = from; j < r->route.size(); ++j)
                if (same(r->route[j], route[k])) {
                    found = j;
                    break;
                }
            if (found < r->route.size()) {
                bends[k] = r->bends[found];
                points[k].side = r->route[found].side;
                changed = changed || found != k;
                from = found + 1;
            } else {
                changed = true;
                kept = kept && !sheave(points[k]);
            }
        }
        if (!changed) return true;
        I.wakeRope(*r);
        r->route = std::move(points);
        r->bends = std::move(bends);
        r->placed = r->placed && kept;
        I.wakeRope(*r);
        return true;
    }

    bool World::setRopeLength(RopeId id, double length) {
        Impl &I = *impl_;
        RopeSlot *r = I.rope(id);
        if (!r) return false;
        if (!nonNegative(length)) return I.reject(false);
        r->length = withinSpools(r->desc, length);
        r->slipping = false;
        I.wakeRope(*r);
        return true;
    }

    bool World::setSpoolLever(RopeId id, double lever) {
        Impl &I = *impl_;
        RopeSlot *r = I.rope(id);
        if (!r) return false;
        if (!std::isfinite(lever)) return I.reject(false);
        const double l = std::clamp(lever, -1.0, 1.0);
        if (l != r->lever) I.wakeRope(*r);
        r->lever = l;
        return true;
    }

    double World::ropeLength(RopeId id) const {
        const RopeSlot *r = impl_->rope(id);
        return r ? r->length : 0.0;
    }

    double World::ropeSpan(RopeId id) const {
        const RopeSlot *r = impl_->rope(id);
        return r ? r->span : 0.0;
    }

    double World::ropeTension(RopeId id) const {
        const RopeSlot *r = impl_->rope(id);
        return r && !r->broken ? r->force : 0.0;
    }

    RopeInfo World::ropeInfo(RopeId id) const {
        const RopeSlot *r = impl_->rope(id);
        if (!r) return {};
        RopeInfo info;
        info.machine = r->held;
        info.taut = !r->broken && r->force > 0.0;
        info.slipping = r->slipping;
        info.broken = r->broken;
        return info;
    }

    int World::ropePath(RopeId id, std::span<RopeBend> out) const {
        const RopeSlot *r = impl_->rope(id);
        if (!r) return 0;
        const std::size_t count = std::min(out.size(), r->bends.size());
        for (std::size_t k = 0; k < count; ++k) out[k] = {r->bends[k].in, r->bends[k].out, r->bends[k].touching};
        return int(count);
    }

    std::span<const RopeId> World::brokenRopes() const { return impl_->brokenRopes; }

} // namespace hf
