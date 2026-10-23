#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <tuple>
#include <vector>

#include "contact/ContactParams.hpp"
#include "math/Dense.hpp"
#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        constexpr int kRestitutionIterations = 4;
        constexpr int kSolverIterations = 16; // per step, split across the sub-steps
        constexpr int kProximalOuter = 4;
        constexpr int kRelaxIterations = 2;
        constexpr int kPositionIterations = 16;

        // The machine that solves a mechanism contact with its joints: the one both ends are on, or one end with the
        // other static. Any other contact is swept.
        int directMachine(const Contact &k, const Group &a, const Group &b) {
            if (!k.mechanism) return -1;
            const int ma = machineOf(a), mb = machineOf(b);
            if (ma >= 0 && (mb == ma || !moves(b))) return ma;
            if (mb >= 0 && !moves(a)) return mb;
            return -1;
        }

    } // namespace

    SoftParams SoftParams::fromHertz(double hertz, double zeta, double h) {
        const double omega = 2.0 * std::numbers::pi * std::min(hertz, 0.25 / h);
        const double a1 = 2.0 * zeta + h * omega;
        const double a2 = h * omega * a1;
        const double a3 = 1.0 / (1.0 + a2);
        return {omega / a1, a2 * a3, a3};
    }

    // Once per step: the host's contacts in group frames, their warm start, and the approach speed a bounce
    // reflects, sampled before anything has touched a velocity.
    void World::Impl::beginContacts() {
        for (Group &g : groups) {
            g.c0 = g.c;
            g.q0 = g.q;
        }
        std::size_t kept = 0;
        for (std::size_t i = 0; i < contacts.size(); ++i) {
            ContactData &c = contacts[i];
            const BodySlot *a = body(c.in.a), *b = body(c.in.b);
            if (!a || !b || a->group == b->group) continue;
            const Group &ga = groups[std::size_t(a->group)], &gb = groups[std::size_t(b->group)];
            if (!moves(ga) && !moves(gb)) continue;
            c.ga = a->group;
            c.gb = b->group;
            c.rA = c.in.point - ga.c;
            c.rB = c.in.point - gb.c;
            c.depth0 = c.in.depth;
            c.n = c.in.normal;
            c.depth = c.in.depth;
            // A gap starts cold: replaying an impulse across it would pull the pair together.
            const WarmStart key{{c.in.a, c.in.b, c.in.feature}};
            const auto it = std::lower_bound(warm.begin(), warm.end(), key);
            const bool found = c.in.depth >= 0.0 && it != warm.end() && it->key == key.key;
            c.ln = found ? it->ln : 0.0;
            c.lt1 = found ? it->lt1 : 0.0;
            c.lt2 = found ? it->lt2 : 0.0;
            c.lr1 = found ? it->lr1 : 0.0;
            c.lr2 = found ? it->lr2 : 0.0;
            c.ls = found ? it->ls : 0.0;
            c.lnMax = 0.0;
            const Vec3 va = ga.v + cross(ga.w, c.rA), vb = gb.v + cross(gb.w, c.rB);
            c.vnImpact = dot(vb - va, c.n);
            c.machine = directMachine(c.in, ga, gb);
            c.row = -1;
            contacts[kept++] = c;
        }
        contacts.resize(kept);

        // Support travels one contact per iteration, so contacts are solved from the ground up: by how many contacts
        // separate them from something static. A counting sort, so equal depths keep the host's order.
        groupDepth.assign(groups.size(), 1 << 30);
        for (std::size_t g = 0; g < groups.size(); ++g)
            if (!moves(groups[g])) groupDepth[g] = 0;
        for (bool changed = true; changed;) {
            changed = false;
            for (const ContactData &c : contacts) {
                int &a = groupDepth[std::size_t(c.ga)], &b = groupDepth[std::size_t(c.gb)];
                if (a > b + 1) a = b + 1, changed = true;
                if (b > a + 1) b = a + 1, changed = true;
            }
        }
        const auto depthOf = [&](const ContactData &c) {
            return std::size_t(
                std::min({groupDepth[std::size_t(c.ga)], groupDepth[std::size_t(c.gb)], int(contacts.size())}));
        };
        depthStart.resize(contacts.size() + 2); // resized, not assigned: it grows geometrically
        std::fill(depthStart.begin(), depthStart.end(), 0);
        for (const ContactData &c : contacts) ++depthStart[depthOf(c) + 1];
        for (std::size_t d = 1; d < depthStart.size(); ++d) depthStart[d] += depthStart[d - 1];
        contactScratch.resize(contacts.size());
        for (const ContactData &c : contacts) contactScratch[depthStart[depthOf(c)]++] = c;
        contacts.swap(contactScratch);

        for (MachineData &md : machines) md.contacts.clear();
        for (std::size_t i = 0; i < contacts.size(); ++i)
            if (contacts[i].machine >= 0) machines[std::size_t(contacts[i].machine)].contacts.push_back(int(i));
    }

    // The warm start for the next step. Only this step's contacts are kept, so the cache can't grow.
    void World::Impl::storeContacts() {
        warm.clear();
        for (const ContactData &c : contacts)
            warm.push_back({{c.in.a, c.in.b, c.in.feature}, c.ln, c.lt1, c.lt2, c.lr1, c.lr2, c.ls});
        std::sort(warm.begin(), warm.end());
    }

    namespace {

        void apply(World::Impl &I, const ContactPrep &p, const ContactData &c, double ln, double lt1, double lt2);

        // Angular impulses about t1, t2 (rolling) and n (spin) on both ends.
        void applyTurn(World::Impl &I, const ContactPrep &p, double l1, double l2, double ls) {
            Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
            if (p.writeA) a.w -= l1 * p.imAr1 + l2 * p.imAr2 + ls * p.imAs;
            if (p.writeB) b.w += l1 * p.imBr1 + l2 * p.imBr2 + ls * p.imBs;
            const double l[6] = {0.0, 0.0, 0.0, l1, l2, ls};
            machineImpulse(I, p, false, 6, l);
        }

        // ωB - ωA about t1, t2 and n.
        Vec3 turnVelocity(World::Impl &I, const ContactPrep &p, const ContactData &c) {
            const Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
            Vec3 w{};
            if (machineOf(b) < 0) w += b.w;
            if (machineOf(a) < 0) w -= a.w;
            double m[6] = {};
            machineVelocity(I, p, false, 6, m);
            return {dot(w, p.t1) + m[3], dot(w, p.t2) + m[4], dot(w, c.n) + m[5]};
        }

    } // namespace

    // The poses have moved since the host detected: the depth is re-derived, not re-detected. Each end moves with its
    // body, turning to first order, and the normal holds, so a rolling contact stays on the ground.
    double World::Impl::contactDepth(const ContactData &c) const {
        const Group &a = groups[std::size_t(c.ga)], &b = groups[std::size_t(c.gb)];
        const Vec3 da = a.c - a.c0 + cross(a.turned, c.rA);
        const Vec3 db = b.c - b.c0 + cross(b.turned, c.rB);
        return c.depth0 - dot(db - da, c.n);
    }

    void World::Impl::updateTurned() {
        for (Group &g : groups) g.turned = rotationVectorOf(g.q * conjugate(g.q0));
    }

    void World::Impl::prepareContacts() {
        prep.clear();
        arena.clear();
        for (const LinkResponse &R : responses) machines[std::size_t(R.machine)].response[std::size_t(R.link)] = -1;
        responses.clear();
        responseArena.clear();
        for (std::vector<LinkRun> &r : runs) r.assign(machines.size(), LinkRun{});
        for (std::size_t ci = 0; ci < contacts.size(); ++ci) {
            ContactData &c = contacts[ci];
            if (c.machine >= 0) continue;
            Group &a = groups[std::size_t(c.ga)], &b = groups[std::size_t(c.gb)];
            c.depth = contactDepth(c);

            // A gap is speculative: a rigid row that lets the pair close it within the sub-step and no faster, so a
            // fast body stops on a surface instead of passing through it. One the pair isn't closing fast enough to
            // reach this sub-step, with room to spare, could only solve to zero: the next sub-step looks again.
            const bool gap = c.depth < 0.0;
            if (gap && c.lnMax <= 0.0) {
                const Vec3 va = a.v + cross(a.w, c.rA), vb = b.v + cross(b.w, c.rB);
                if (-dot(vb - va, c.n) * h * kGapReach < -c.depth) continue;
            }

            ContactPrep p;
            p.contact = int(ci);
            p.a = c.ga;
            p.b = c.gb;

            // A machine end's mass is its whole linkage's, through W: its rigid terms are zero.
            const bool rigidA = moves(a) && machineOf(a) < 0, rigidB = moves(b) && machineOf(b) < 0;
            p.writeA = rigidA;
            p.writeB = rigidB;
            if (rigidA) a.pv = a.pw = Vec3{};
            if (rigidB) b.pv = b.pw = Vec3{};
            if (!rigidA && !rigidB && machineOf(a) < 0 && machineOf(b) < 0) continue;
            p.imA = rigidA ? a.invMass : 0.0;
            p.imB = rigidB ? b.invMass : 0.0;
            p.invIwA = rigidA ? a.invInertiaWorld() : Mat3{};
            p.invIwB = rigidB ? b.invInertiaWorld() : Mat3{};

            // The lever arms hold for the step: a contact stays under a body that moves across it.
            p.r1 = c.rA;
            p.r2 = c.rB;
            contactBasis(c.n, p.t1, p.t2);
            p.r1cn = cross(p.r1, c.n);
            p.r2cn = cross(p.r2, c.n);
            p.r1ct1 = cross(p.r1, p.t1);
            p.r1ct2 = cross(p.r1, p.t2);
            p.r2ct1 = cross(p.r2, p.t1);
            p.r2ct2 = cross(p.r2, p.t2);
            p.imR1cn = p.invIwA * p.r1cn;
            p.imR2cn = p.invIwB * p.r2cn;
            p.imR1ct1 = p.invIwA * p.r1ct1;
            p.imR2ct1 = p.invIwB * p.r2ct1;
            p.imR1ct2 = p.invIwA * p.r1ct2;
            p.imR2ct2 = p.invIwB * p.r2ct2;
            p.kn = p.imA + p.imB + dot(p.r1cn, p.imR1cn) + dot(p.r2cn, p.imR2cn);
            p.kt1 = p.imA + p.imB + dot(p.r1ct1, p.imR1ct1) + dot(p.r2ct1, p.imR2ct1);
            p.kt2 = p.imA + p.imB + dot(p.r1ct2, p.imR1ct2) + dot(p.r2ct2, p.imR2ct2);
            p.kt12 = dot(p.r1ct1, p.imR1ct2) + dot(p.r2ct1, p.imR2ct2);
            p.turns = c.in.rolling > 0.0 || c.in.spinning > 0.0;
            if (p.turns) {
                p.imAr1 = p.invIwA * p.t1;
                p.imAr2 = p.invIwA * p.t2;
                p.imAs = p.invIwA * c.n;
                p.imBr1 = p.invIwB * p.t1;
                p.imBr2 = p.invIwB * p.t2;
                p.imBs = p.invIwB * c.n;
                p.kr1 = dot(p.t1, p.imAr1 + p.imBr1);
                p.kr2 = dot(p.t2, p.imAr2 + p.imBr2);
                p.ks = dot(c.n, p.imAs + p.imBs);
            }
            p.soft = gap ? SoftParams::rigid() : SoftParams::fromHertz(kContactHertz, kContactDamping, h);

            addMachineSides(*this, p, c);
            // Neither end can move along the normal, so the contact can do nothing, and its rows would divide by zero.
            if (!(p.kn > 0.0)) continue;
            apply(*this, p, c, c.ln, c.lt1, c.lt2);
            if (p.turns) applyTurn(*this, p, c.lr1, c.lr2, c.ls);

            // Penetration is recovered in the position pass, where it moves the pose without leaving momentum.
            p.vtarget = gap ? c.depth / h : 0.0;
            p.posBias = gap ? 0.0 : std::min(kMaxBias, p.soft.biasRate * std::max(0.0, c.depth - kSlop));
            prep.push_back(p);
        }
        flushRuns(*this, false);
    }

    namespace {

        struct RowVelocity {
            double n = 0.0, t1 = 0.0, t2 = 0.0;
        };

        // The relative velocity vB - vA at the contact: a rigid end from its group, a machine end from qdot (or its
        // pseudo-velocity, in the position pass) through its rows.
        RowVelocity rowVelocity(World::Impl &I, const ContactPrep &p, const ContactData &c, bool pseudo) {
            const Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
            Vec3 vrel{};
            if (machineOf(b) < 0) vrel += pseudo ? b.pv + cross(b.pw, p.r2) : b.v + cross(b.w, p.r2);
            if (machineOf(a) < 0) vrel -= pseudo ? a.pv + cross(a.pw, p.r1) : a.v + cross(a.w, p.r1);
            double m[3] = {};
            machineVelocity(I, p, pseudo, pseudo ? 1 : 3, m);
            return {dot(vrel, c.n) + m[0], dot(vrel, p.t1) + m[1], dot(vrel, p.t2) + m[2]};
        }

        void apply(World::Impl &I, const ContactPrep &p, const ContactData &c, double ln, double lt1, double lt2) {
            Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
            const Vec3 Pt = lt1 * p.t1 + lt2 * p.t2;
            if (p.writeA) {
                a.v -= p.imA * ln * c.n;
                a.w -= ln * p.imR1cn;
                a.v -= p.imA * Pt;
                a.w -= lt1 * p.imR1ct1 + lt2 * p.imR1ct2;
            }
            if (p.writeB) {
                b.v += p.imB * ln * c.n;
                b.w += ln * p.imR2cn;
                b.v += p.imB * Pt;
                b.w += lt1 * p.imR2ct1 + lt2 * p.imR2ct2;
            }
            const double l[3] = {ln, lt1, lt2};
            machineImpulse(I, p, false, 3, l);
        }

        // The friction pair's new accumulated impulse, clamped to the cone. With a machine end the 2×2 block is solved
        // exactly: a wheel's mobility lies along its rolling direction, and a row-by-row step diverges on it.
        void frictionStep(const ContactPrep &p, const RowVelocity &v, double maxT, double &lt1, double &lt2) {
            // A direction neither end can move along takes no impulse: a hinge's axis on a face square to it.
            double d1 = p.kt1 > 0.0 ? -v.t1 / p.kt1 : 0.0, d2 = p.kt2 > 0.0 ? -v.t2 / p.kt2 : 0.0;
            if (p.ns > 0) {
                const double a = p.kt1, b = p.kt12, d = p.kt2;
                const double det = a * d - b * b;
                if (det > 1e-12 * a * d) {
                    d1 = -(d * v.t1 - b * v.t2) / det;
                    d2 = -(a * v.t2 - b * v.t1) / det;
                }
            }
            lt1 += d1;
            lt2 += d2;
            // Outside the cone, back onto its rim; not a number at all, nothing.
            const double len = std::sqrt(lt1 * lt1 + lt2 * lt2);
            if (!(len <= maxT)) {
                const double s = std::isfinite(len) && len > 0.0 ? maxT / len : 0.0;
                lt1 *= s;
                lt2 *= s;
            }
        }

        // One pass over a contact: soft when biased, the pure non-penetration row when relaxing.
        // Rolling and spin, resisted by up to their lever arm times the normal force, once the linear rows are done.
        void turnStep(World::Impl &I, const ContactPrep &p, ContactData &c, double ln) {
            const Vec3 w = turnVelocity(I, p, c);
            double l1 = p.kr1 > 0.0 ? c.lr1 - w.x / p.kr1 : 0.0, l2 = p.kr2 > 0.0 ? c.lr2 - w.y / p.kr2 : 0.0;
            const double maxR = c.in.rolling * ln, len = std::sqrt(l1 * l1 + l2 * l2);
            if (!(len <= maxR)) {
                const bool onRim = std::isfinite(len) && len > 0.0;
                l1 = onRim ? l1 * maxR / len : 0.0;
                l2 = onRim ? l2 * maxR / len : 0.0;
            }
            const double maxS = c.in.spinning * ln;
            const double ls = p.ks > 0.0 ? std::clamp(c.ls - w.z / p.ks, -maxS, maxS) : 0.0;
            applyTurn(I, p, l1 - c.lr1, l2 - c.lr2, ls - c.ls);
            c.lr1 = l1;
            c.lr2 = l2;
            c.ls = ls;
        }

        void sweepContact(World::Impl &I, ContactPrep &p, bool relax) {
            ContactData &c = I.contacts[std::size_t(p.contact)];
            const RowVelocity v = rowVelocity(I, p, c, false);
            const double lnn = relax ? std::max(0.0, c.ln + (p.vtarget - v.n) / p.kn)
                                     : std::max(0.0, c.ln + p.soft.massCoeff * (p.vtarget - v.n) / p.kn -
                                                         p.soft.impulseCoeff * (c.ln - p.refLn));
            double lt1 = c.lt1, lt2 = c.lt2;
            frictionStep(p, v, c.in.friction * lnn, lt1, lt2);
            const double dn = lnn - c.ln, d1 = lt1 - c.lt1, d2 = lt2 - c.lt2;
            c.ln = lnn;
            c.lnMax = std::max(c.lnMax, lnn);
            c.lt1 = lt1;
            c.lt2 = lt2;
            apply(I, p, c, dn, d1, d2);
            if (p.turns) turnStep(I, p, c, lnn);
        }

        void sweepBounded(World::Impl &I) {
            for (MachineData &md : I.machines)
                if (!md.frozen && md.m.boundedRows()) md.m.sweepBounded();
        }

    } // namespace

    void World::Impl::prepareSprings() {
        springs.clear();
        for (const JointRef &r : looseSprings) {
            const JointDesc &d = joints[r.joint].desc;
            const Group &a = groups[std::size_t(r.A.group)], &b = groups[std::size_t(r.B.group)];
            const Vec3 pA = a.c + a.q * r.A.anchor, pB = b.c + b.q * r.B.anchor;
            const double len = length(pA - pB);
            if (!(len > 1e-9)) continue;

            SpringPrep sp;
            sp.a = r.A.group;
            sp.b = r.B.group;
            sp.n = (pA - pB) / len;
            sp.r1 = pA - a.c;
            sp.r2 = pB - b.c;
            const bool rigidA = moves(a) && machineOf(a) < 0, rigidB = moves(b) && machineOf(b) < 0;
            sp.imA = rigidA ? a.invMass : 0.0;
            sp.imB = rigidB ? b.invMass : 0.0;
            sp.imR1n = (rigidA ? a.invInertiaWorld() : Mat3{}) * cross(sp.r1, sp.n);
            sp.imR2n = (rigidB ? b.invInertiaWorld() : Mat3{}) * cross(sp.r2, sp.n);
            sp.k = sp.imA + sp.imB + dot(cross(sp.r1, sp.n), sp.imR1n) + dot(cross(sp.r2, sp.n), sp.imR2n);
            sp.writeA = rigidA;
            sp.writeB = rigidB;

            if (!moves(a) && !moves(b)) continue;
            for (int mi : {machineOf(a), machineOf(b)}) {
                if (mi < 0 || (sp.nb == 1 && sp.blocks[0].machine == mi)) continue;
                Machine &m = machines[std::size_t(mi)].m;
                const int n = m.dofs();
                ContactBlock blk{mi, arena.size(), n};
                arena.resize(blk.off + 2 * std::size_t(n), 0.0);
                double *g = arena.data() + blk.off, *w = g + n;
                if (machineOf(a) == mi) m.art.addPointRow(a.link, pA, sp.n, 1.0, g);
                if (machineOf(b) == mi) m.art.addPointRow(b.link, pB, sp.n, -1.0, g);
                const double unheld = m.mobility(g, w);
                sp.k += std::max(dot(g, w, std::size_t(n)), kLockedMobility * unheld);
                sp.blocks[sp.nb++] = blk;
            }
            springRow(sp.k, d.hertz, d.damping, len - d.restLength, h, sp.gamma, sp.target);
            // Neither end can move along it: between two machines anchored on their hinge axes, say. It can do
            // nothing, and its row would be 0/0.
            if (!(sp.k + sp.gamma > 0.0)) continue;
            springs.push_back(sp);
        }
    }

    namespace {

        void sweepSpring(World::Impl &I, SpringPrep &sp) {
            Group &a = I.groups[std::size_t(sp.a)], &b = I.groups[std::size_t(sp.b)];
            double v = 0.0;
            if (machineOf(a) < 0) v += dot(a.v + cross(a.w, sp.r1), sp.n);
            if (machineOf(b) < 0) v -= dot(b.v + cross(b.w, sp.r2), sp.n);
            for (int k = 0; k < sp.nb; ++k) {
                const ContactBlock &blk = sp.blocks[k];
                const std::vector<double> &qd = I.machines[std::size_t(blk.machine)].m.art.qdot;
                v += dot(I.arena.data() + blk.off, qd.data(), std::size_t(blk.n));
            }
            const double dl = (sp.target - v - sp.gamma * sp.lambda) / (sp.k + sp.gamma);
            sp.lambda += dl;
            if (sp.writeA) {
                a.v += sp.imA * dl * sp.n;
                a.w += dl * sp.imR1n;
            }
            if (sp.writeB) {
                b.v -= sp.imB * dl * sp.n;
                b.w -= dl * sp.imR2n;
            }
            for (int k = 0; k < sp.nb; ++k) {
                const ContactBlock &blk = sp.blocks[k];
                I.machines[std::size_t(blk.machine)].m.applyContact(I.arena.data() + blk.off + std::size_t(blk.n), dl);
            }
        }

    } // namespace

    // The biased sweeps. The proximal anchor is re-taken each outer pass, so the impulse bleed cancels at the fixed
    // point. A machine's capped drives and limits are swept in the same loop as the contacts that load them.
    void World::Impl::solveContacts() {
        const int iterations =
            std::max(std::max(1, settings.contactIterations), kSolverIterations / std::max(1, settings.substeps));
        const int outer = std::min(kProximalOuter, iterations);
        const int inner = std::max(1, iterations / outer);
        for (int o = 0; o < outer; ++o) {
            for (ContactPrep &p : prep) p.refLn = contacts[std::size_t(p.contact)].ln;
            for (int k = 0; k < inner; ++k) {
                for (SpringPrep &sp : springs) sweepSpring(*this, sp);
                for (ContactPrep &p : prep) sweepContact(*this, p, false);
                flushRuns(*this, false);
                sweepBounded(*this);
            }
        }
    }

    // After the integration, and nothing is integrated from it: drains the separating velocity the bias added, and
    // leaves impulses the next warm start can trust.
    void World::Impl::relaxContacts() {
        // A gap keeps what is left of it after the integration; anything touching is held to no approach at all.
        updateTurned();
        for (ContactPrep &p : prep) {
            const ContactData &c = contacts[std::size_t(p.contact)];
            p.vtarget = std::min(0.0, contactDepth(c)) / h;
        }
        for (int k = 0; k < kRelaxIterations; ++k) {
            for (ContactPrep &p : prep) sweepContact(*this, p, true);
            flushRuns(*this, false);
            sweepBounded(*this);
        }
    }

    void World::Impl::solveContactPositions() {
        for (int k = 0; k < kPositionIterations; ++k)
            for (ContactPrep &p : prep) {
                if (p.posBias == 0.0) continue;
                const ContactData &c = contacts[std::size_t(p.contact)];
                const double dl = (p.posBias - rowVelocity(*this, p, c, true).n) / p.kn;
                const double nl = std::max(0.0, p.pseudo + dl); // push only, like the velocity row
                const double app = nl - p.pseudo;
                p.pseudo = nl;
                Group &a = groups[std::size_t(p.a)], &b = groups[std::size_t(p.b)];
                if (p.writeA) {
                    a.pv -= p.imA * app * c.n;
                    a.pw -= app * p.imR1cn;
                }
                if (p.writeB) {
                    b.pv += p.imB * app * c.n;
                    b.pw += app * p.imR2cn;
                }
                machineImpulse(*this, p, true, 1, &app);
            }
        flushRuns(*this, true);
    }

    // Once per step, after every sub-step: a bounce is a property of the impact, solved against the speed the pair
    // actually met at.
    void World::Impl::restitution() {
        for (int k = 0; k < kRestitutionIterations; ++k)
            for (ContactPrep &p : prep) {
                ContactData &c = contacts[std::size_t(p.contact)];
                if (c.in.restitution <= 0.0 || c.vnImpact > -kRestitutionMinSpeed || c.lnMax <= 0.0) continue;
                const double vn = rowVelocity(*this, p, c, false).n;
                const double lnn = std::max(0.0, c.ln + (-c.in.restitution * c.vnImpact - vn) / p.kn);
                const double dn = lnn - c.ln;
                c.ln = lnn;
                apply(*this, p, c, dn, 0.0, 0.0);
            }
        flushRuns(*this, false);
        for (int k = 0; k < kRestitutionIterations; ++k)
            for (MachineData &md : machines)
                for (int ci : md.contacts) {
                    const ContactData &c = contacts[std::size_t(ci)];
                    if (c.row < 0 || c.in.restitution <= 0.0 || c.vnImpact > -kRestitutionMinSpeed || c.lnMax <= 0.0)
                        continue;
                    md.m.pushRow(c.row, -c.in.restitution * c.vnImpact);
                }
    }

} // namespace hf
