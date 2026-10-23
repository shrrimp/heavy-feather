#include <algorithm>
#include <cmath>

#include "contact/ContactParams.hpp"
#include "world/WorldImpl.hpp"

// A rigid island's contacts, lean. A row is built once a step: its basis, its lever arms and its effective masses
// hold for the step, as the body inertias it was built with do, and only its depth is refreshed each sub-step. The
// bodies it pushes sit together in one array, a moving body's slot written, any other only read.
namespace hf {

    namespace {

        constexpr int kProximalOuter = 4;
        constexpr int kRestitutionIterations = 4;
        constexpr int kRelaxIterations = 2; // a machine island's: what an impact's rows get

        struct Velocity {
            double n = 0.0, t1 = 0.0, t2 = 0.0;
        };

        Velocity velocityOf(const std::vector<SolverBody> &B, const RigidRow &r, bool pseudo) {
            const SolverBody &a = B[std::size_t(r.a)], &b = B[std::size_t(r.b)];
            const Vec3 v = pseudo ? (b.pv + cross(b.pw, r.r2)) - (a.pv + cross(a.pw, r.r1))
                                  : (b.v + cross(b.w, r.r2)) - (a.v + cross(a.w, r.r1));
            return {dot(v, r.n), dot(v, r.t1), dot(v, r.t2)};
        }

        // An impulse P at the contact: -P on a, +P on b.
        void push(std::vector<SolverBody> &B, const RigidRow &r, const Vec3 &P) {
            SolverBody &a = B[std::size_t(r.a)], &b = B[std::size_t(r.b)];
            if (a.moves) {
                a.v -= a.invMass * P;
                a.w -= a.invI * cross(r.r1, P);
            }
            if (b.moves) {
                b.v += b.invMass * P;
                b.w += b.invI * cross(r.r2, P);
            }
        }

        void turn(std::vector<SolverBody> &B, const RigidRow &r, const Vec3 &L) {
            SolverBody &a = B[std::size_t(r.a)], &b = B[std::size_t(r.b)];
            if (a.moves) a.w -= a.invI * L;
            if (b.moves) b.w += b.invI * L;
        }

        void turnStep(std::vector<SolverBody> &B, const RigidRow &r, RigidTurn &t, double ln) {
            const SolverBody &a = B[std::size_t(r.a)], &b = B[std::size_t(r.b)];
            const Vec3 w = b.w - a.w;
            double l1 = t.kr1 > 0.0 ? t.lr1 - dot(w, r.t1) / t.kr1 : 0.0;
            double l2 = t.kr2 > 0.0 ? t.lr2 - dot(w, r.t2) / t.kr2 : 0.0;
            const double maxR = t.rolling * ln, len = std::sqrt(l1 * l1 + l2 * l2);
            if (!(len <= maxR)) {
                const bool onRim = std::isfinite(len) && len > 0.0;
                l1 = onRim ? l1 * maxR / len : 0.0;
                l2 = onRim ? l2 * maxR / len : 0.0;
            }
            const double maxS = t.spinning * ln;
            const double ls = t.ks > 0.0 ? std::clamp(t.ls - dot(w, r.n) / t.ks, -maxS, maxS) : 0.0;
            turn(B, r, (l1 - t.lr1) * r.t1 + (l2 - t.lr2) * r.t2 + (ls - t.ls) * r.n);
            t.lr1 = l1;
            t.lr2 = l2;
            t.ls = ls;
        }

        // One pass over a row: soft while solving, the bare non-penetration row when relaxing.
        void sweep(std::vector<SolverBody> &B, std::vector<RigidTurn> &turns, RigidRow &r, bool relax) {
            const Velocity v = velocityOf(B, r, false);
            const double lnn =
                relax
                    ? std::max(0.0, r.ln + (r.vtarget - v.n) / r.kn)
                    : std::max(0.0, r.ln + r.massCoeff * (r.vtarget - v.n) / r.kn - r.impulseCoeff * (r.ln - r.refLn));
            double lt1 = r.lt1 + (r.kt1 > 0.0 ? -v.t1 / r.kt1 : 0.0);
            double lt2 = r.lt2 + (r.kt2 > 0.0 ? -v.t2 / r.kt2 : 0.0);
            const double maxT = r.friction * lnn, len = std::sqrt(lt1 * lt1 + lt2 * lt2);
            if (!(len <= maxT)) {
                const double s = std::isfinite(len) && len > 0.0 ? maxT / len : 0.0;
                lt1 *= s;
                lt2 *= s;
            }
            push(B, r, (lnn - r.ln) * r.n + (lt1 - r.lt1) * r.t1 + (lt2 - r.lt2) * r.t2);
            r.ln = lnn;
            r.lnMax = std::max(r.lnMax, lnn);
            r.lt1 = lt1;
            r.lt2 = lt2;
            if (r.turn >= 0) turnStep(B, r, turns[std::size_t(r.turn)], lnn);
        }

        // The proximal loop: the anchor the soft rows bleed toward is re-taken each outer pass, so the bleed cancels
        // at the fixed point.
        template <typename Rows>
        void solve(World::Impl &I, const Rows &rows, int count, int iterations) {
            if (count == 0 || iterations <= 0) return;
            const int outer = std::min(kProximalOuter, iterations);
            const int inner = std::max(1, iterations / outer);
            for (int o = 0; o < outer; ++o) {
                for (int k = 0; k < count; ++k) {
                    RigidRow &r = rows(k);
                    r.refLn = r.ln;
                }
                for (int i = 0; i < inner; ++i)
                    for (int k = 0; k < count; ++k) {
                        RigidRow &r = rows(k);
                        if (r.active) sweep(I.solverBodies, I.rigidTurns, r, false);
                    }
            }
        }

        template <typename Rows>
        void relaxRows(World::Impl &I, const Rows &rows, int count, int iterations) {
            for (int i = 0; i < iterations; ++i)
                for (int k = 0; k < count; ++k) {
                    RigidRow &r = rows(k);
                    if (r.active) sweep(I.solverBodies, I.rigidTurns, r, true);
                }
        }

    } // namespace

    void World::Impl::beginRigidRows() {
        rigidRows.clear();
        rigidTurns.clear();
        solverBodies.clear();
        slotOfGroup.assign(groups.size(), -1);
        rigidStart.assign(std::size_t(islands) + 1, 0);
        const auto slot = [&](int gi) {
            int &s = slotOfGroup[std::size_t(gi)];
            if (s >= 0) return s;
            s = int(solverBodies.size());
            const Group &g = groups[std::size_t(gi)];
            SolverBody b;
            b.group = gi;
            b.moves = moves(g);
            b.invMass = b.moves ? g.invMass : 0.0;
            b.invI = b.moves ? g.invInertiaWorld() : Mat3{};
            solverBodies.push_back(b);
            return s;
        };
        for (std::size_t ci = 0; ci < contacts.size(); ++ci) {
            const ContactData &c = contacts[ci];
            if (machineIsland[std::size_t(c.island)]) continue;
            RigidRow r;
            r.contact = int(ci);
            r.a = slot(c.ga);
            r.b = slot(c.gb);
            r.n = c.n;
            contactBasis(c.n, r.t1, r.t2);
            r.r1 = c.rA;
            r.r2 = c.rB;
            const SolverBody &a = solverBodies[std::size_t(r.a)], &b = solverBodies[std::size_t(r.b)];
            const auto mass = [&](const Vec3 &d) {
                const Vec3 ca = cross(r.r1, d), cb = cross(r.r2, d);
                return a.invMass + b.invMass + dot(ca, a.invI * ca) + dot(cb, b.invI * cb);
            };
            r.kn = mass(r.n);
            r.kt1 = mass(r.t1);
            r.kt2 = mass(r.t2);
            // Neither end can move along the normal: it can do nothing.
            if (!(r.kn > 0.0)) continue;
            r.friction = c.in.friction;
            r.ln = c.ln;
            r.lt1 = c.lt1;
            r.lt2 = c.lt2;
            if (c.in.rolling > 0.0 || c.in.spinning > 0.0) {
                RigidTurn t;
                const Mat3 I = a.invI + b.invI;
                t.kr1 = dot(r.t1, I * r.t1);
                t.kr2 = dot(r.t2, I * r.t2);
                t.ks = dot(r.n, I * r.n);
                t.lr1 = c.lr1;
                t.lr2 = c.lr2;
                t.ls = c.ls;
                t.rolling = c.in.rolling;
                t.spinning = c.in.spinning;
                r.turn = int(rigidTurns.size());
                rigidTurns.push_back(t);
            }
            rigidRows.push_back(r);
            ++rigidStart[std::size_t(c.island) + 1];
        }
        for (std::size_t i = 1; i < rigidStart.size(); ++i) rigidStart[i] += rigidStart[i - 1];
    }

    void World::Impl::prepareRigidRows() {
        const std::size_t L = std::size_t(islands);
        impactRows.clear();
        impactStart.assign(L + 1, 0);
        rowsOfKind[0] = 0;
        for (SolverBody &b : solverBodies) {
            const Group &g = groups[std::size_t(b.group)];
            b.v = g.v;
            b.w = g.w;
            b.pv = b.pw = Vec3{};
        }
        const SoftParams soft = SoftParams::fromHertz(kContactHertz, kContactDamping, h);
        for (std::size_t island = 0; island < L; ++island) {
            impactStart[island] = int(impactRows.size());
            if (machineIsland[island]) continue;
            for (int k = rigidStart[island]; k < rigidStart[island + 1]; ++k) {
                RigidRow &r = rigidRows[std::size_t(k)];
                ContactData &c = contacts[std::size_t(r.contact)];
                c.depth = contactDepth(c);
                if (c.impact > 0) --c.impact;
                // A gap is speculative: the pair may close it within the sub-step and no faster. One it isn't closing
                // fast enough to reach could only solve to zero: the next sub-step looks again.
                const bool gap = c.depth < 0.0;
                const double approach = -velocityOf(solverBodies, r, false).n;
                r.active = !(gap && r.lnMax <= 0.0 && approach * h * kGapReach < -c.depth);
                if (!r.active) continue;
                ++rowsOfKind[0];
                if (approach > kRestitutionMinSpeed) c.impact = 2;
                if (c.impact > 0) impactRows.push_back(k);
                const SoftParams s = gap ? SoftParams::rigid() : soft;
                r.massCoeff = s.massCoeff;
                r.impulseCoeff = s.impulseCoeff;
                r.vtarget = gap ? c.depth / h : 0.0;
                r.posBias = gap ? 0.0 : std::min(kMaxBias, s.biasRate * std::max(0.0, c.depth - kSlop));
                r.pseudo = 0.0;
                push(solverBodies, r, r.ln * r.n + r.lt1 * r.t1 + r.lt2 * r.t2);
                if (r.turn >= 0) {
                    const RigidTurn &t = rigidTurns[std::size_t(r.turn)];
                    turn(solverBodies, r, t.lr1 * r.t1 + t.lr2 * r.t2 + t.ls * r.n);
                }
            }
        }
        impactStart[L] = int(impactRows.size());
    }

    void World::Impl::solveRigidIsland(int island) {
        const std::size_t i = std::size_t(island);
        const int begin = rigidStart[i], count = rigidStart[i + 1] - begin;
        solve(
            *this, [&](int k) -> RigidRow & { return rigidRows[std::size_t(begin + k)]; }, count,
            std::max(1, settings.rigidIterations));
        // An impact's rows, on a machine island's budget of their own.
        const int *impacts = impactRows.data() + impactStart[i];
        solve(
            *this, [&](int k) -> RigidRow & { return rigidRows[std::size_t(impacts[k])]; },
            impactStart[i + 1] - impactStart[i], std::max(1, settings.contactIterations));
    }

    void World::Impl::relaxRigidIsland(int island) {
        const std::size_t i = std::size_t(island);
        const int begin = rigidStart[i], count = rigidStart[i + 1] - begin;
        for (int k = begin; k < begin + count; ++k) {
            RigidRow &r = rigidRows[std::size_t(k)];
            if (r.active) r.vtarget = std::min(0.0, contactDepth(contacts[std::size_t(r.contact)])) / h;
        }
        relaxRows(
            *this, [&](int k) -> RigidRow & { return rigidRows[std::size_t(begin + k)]; }, count,
            std::max(0, settings.rigidRelaxIterations));
        const int *impacts = impactRows.data() + impactStart[i];
        relaxRows(
            *this, [&](int k) -> RigidRow & { return rigidRows[std::size_t(impacts[k])]; },
            impactStart[i + 1] - impactStart[i], kRelaxIterations);
    }

    // Penetration past the slop, pushed out with a pseudo-velocity the integration spends and forgets.
    void World::Impl::rigidPositions(int island) {
        const std::size_t i = std::size_t(island);
        const int iterations = std::max(0, settings.rigidPositionIterations);
        for (int it = 0; it < iterations; ++it)
            for (int k = rigidStart[i]; k < rigidStart[i + 1]; ++k) {
                RigidRow &r = rigidRows[std::size_t(k)];
                if (!r.active || r.posBias == 0.0) continue;
                const double nl = std::max(0.0, r.pseudo + (r.posBias - velocityOf(solverBodies, r, true).n) / r.kn);
                const double d = nl - r.pseudo;
                r.pseudo = nl;
                SolverBody &a = solverBodies[std::size_t(r.a)], &b = solverBodies[std::size_t(r.b)];
                const Vec3 P = d * r.n;
                if (a.moves) {
                    a.pv -= a.invMass * P;
                    a.pw -= a.invI * cross(r.r1, P);
                }
                if (b.moves) {
                    b.pv += b.invMass * P;
                    b.pw += b.invI * cross(r.r2, P);
                }
            }
    }

    void World::Impl::rigidVelocitiesOut() {
        for (const SolverBody &b : solverBodies) {
            if (!b.moves) continue;
            Group &g = groups[std::size_t(b.group)];
            g.v = b.v;
            g.w = b.w;
        }
    }

    void World::Impl::rigidPositionsOut() {
        for (const SolverBody &b : solverBodies) {
            if (!b.moves) continue;
            Group &g = groups[std::size_t(b.group)];
            g.pv = b.pv;
            g.pw = b.pw;
        }
    }

    // A bounce is a property of the impact: solved once the sub-steps are done, against the speed the pair met at.
    void World::Impl::rigidRestitution() {
        for (int it = 0; it < kRestitutionIterations; ++it)
            for (RigidRow &r : rigidRows) {
                const ContactData &c = contacts[std::size_t(r.contact)];
                if (c.in.restitution <= 0.0 || c.vnImpact > -kRestitutionMinSpeed || r.lnMax <= 0.0) continue;
                const double vn = velocityOf(solverBodies, r, false).n;
                const double lnn = std::max(0.0, r.ln + (-c.in.restitution * c.vnImpact - vn) / r.kn);
                push(solverBodies, r, (lnn - r.ln) * r.n);
                r.ln = lnn;
            }
        rigidVelocitiesOut();
    }

    void World::Impl::endRigidRows() {
        for (const RigidRow &r : rigidRows) {
            ContactData &c = contacts[std::size_t(r.contact)];
            c.ln = r.ln;
            c.lt1 = r.lt1;
            c.lt2 = r.lt2;
            c.lnMax = r.lnMax;
            if (r.turn >= 0) {
                const RigidTurn &t = rigidTurns[std::size_t(r.turn)];
                c.lr1 = t.lr1;
                c.lr2 = t.lr2;
                c.ls = t.ls;
            }
        }
    }

} // namespace hf
