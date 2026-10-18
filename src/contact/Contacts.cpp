#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        constexpr double kContactHertz = 120.0;
        constexpr double kContactDamping = 4.0; // a larger ratio makes the row stiffer: 4 is the quiet point
        constexpr double kSlop = 0.002;         // m of penetration left alone
        constexpr double kMaxBias = 3.0;        // m/s: the fastest penetration is pushed out
        constexpr double kRestitutionMinSpeed = 0.5;
        constexpr int kRestitutionIterations = 4;
        constexpr int kSolverIterations = 16; // per step, split across the sub-steps
        constexpr int kMinSolverIterations = 8;
        constexpr int kProximalOuter = 4;
        constexpr int kRelaxIterations = 2;
        constexpr int kPositionIterations = 16;
        // Two links of one machine its held rows lock together have ~0 relative mobility through W. Floored at this
        // much of their unconstrained mobility, far below any real mass ratio, their contact stays finite.
        constexpr double kLockedMobility = 1e-7;

        // A stable tangent pair for any unit normal; (n, t1, t2) is right-handed.
        void contactBasis(const Vec3 &n, Vec3 &t1, Vec3 &t2) {
            const double s = std::copysign(1.0, n.z);
            const double a = -1.0 / (s + n.z);
            const double b = n.x * n.y * a;
            t1 = {1.0 + s * n.x * n.x * a, s * b, -s * n.x};
            t2 = {b, s + n.y * n.y * a, -n.y};
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
        std::vector<ContactData> kept;
        kept.reserve(contacts.size());
        for (ContactData &c : contacts) {
            const BodySlot *a = body(c.in.a), *b = body(c.in.b);
            if (!a || !b || a->group == b->group) continue;
            const Group &ga = groups[std::size_t(a->group)], &gb = groups[std::size_t(b->group)];
            if (ga.isStatic && gb.isStatic) continue;
            c.ga = a->group;
            c.gb = b->group;
            c.rA = c.in.point - ga.c;
            c.rB = c.in.point - gb.c;
            c.depth0 = c.in.depth;
            c.n = c.in.normal;
            c.depth = c.in.depth;
            const auto it = warm.find({c.in.a, c.in.b, c.in.feature});
            c.ln = it != warm.end() ? it->second.first : 0.0;
            c.lt1 = it != warm.end() ? it->second.second.x : 0.0;
            c.lt2 = it != warm.end() ? it->second.second.y : 0.0;
            c.lnMax = 0.0;
            const Vec3 va = ga.v + cross(ga.w, c.rA), vb = gb.v + cross(gb.w, c.rB);
            c.vnImpact = dot(vb - va, c.n);
            kept.push_back(c);
        }
        contacts.swap(kept);

        // Support travels one contact per iteration, so contacts are solved from the ground up: by how many contacts
        // separate them from something static. Stable, so equal depths keep the host's order.
        std::vector<int> depth(groups.size(), 1 << 30);
        for (std::size_t g = 0; g < groups.size(); ++g)
            if (groups[g].isStatic) depth[g] = 0;
        for (bool changed = true; changed;) {
            changed = false;
            for (const ContactData &c : contacts) {
                const std::size_t a = std::size_t(c.ga), b = std::size_t(c.gb);
                if (depth[a] > depth[b] + 1) {
                    depth[a] = depth[b] + 1;
                    changed = true;
                }
                if (depth[b] > depth[a] + 1) {
                    depth[b] = depth[a] + 1;
                    changed = true;
                }
            }
        }
        std::stable_sort(contacts.begin(), contacts.end(), [&](const ContactData &x, const ContactData &y) {
            return std::min(depth[std::size_t(x.ga)], depth[std::size_t(x.gb)]) <
                   std::min(depth[std::size_t(y.ga)], depth[std::size_t(y.gb)]);
        });
    }

    // The warm start for the next step. Only this step's contacts are kept, so the cache can't grow.
    void World::Impl::storeContacts() {
        warm.clear();
        for (const ContactData &c : contacts) warm[{c.in.a, c.in.b, c.in.feature}] = {c.ln, Vec3{c.lt1, c.lt2, 0.0}};
    }

    namespace {

        // A machine end of a contact: its normal and tangent rows through the machine's mobility, and their effective
        // inverse masses added to the row's. When both ends are links of the same machine, one block carries both.
        void addMachineBlock(World::Impl &I, ContactPrep &p, const ContactData &c, int mi) {
            MachineData &md = I.machines[std::size_t(mi)];
            Machine &m = md.m;
            const int n = m.dofs();
            ContactBlock blk{mi, I.arena.size(), n};
            I.arena.resize(blk.off + 6 * std::size_t(n), 0.0);
            double *gn = I.arena.data() + blk.off, *gt1 = gn + n, *gt2 = gt1 + n;
            double *wn = gt2 + n, *wt1 = wn + n, *wt2 = wt1 + n;

            const Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
            // B's end with +, A's with -: the row's relative velocity is vB - vA.
            if (b.machine == mi) {
                const Vec3 at = b.c + p.r2;
                m.art.addPointRow(b.link, at, c.n, 1.0, gn);
                m.art.addPointRow(b.link, at, p.t1, 1.0, gt1);
                m.art.addPointRow(b.link, at, p.t2, 1.0, gt2);
            }
            if (a.machine == mi) {
                const Vec3 at = a.c + p.r1;
                m.art.addPointRow(a.link, at, c.n, -1.0, gn);
                m.art.addPointRow(a.link, at, p.t1, -1.0, gt1);
                m.art.addPointRow(a.link, at, p.t2, -1.0, gt2);
            }
            m.mobility(gn, wn);
            m.mobility(gt1, wt1);
            m.mobility(gt2, wt2);
            double kn = 0.0, k1 = 0.0, k2 = 0.0, k12 = 0.0;
            for (int i = 0; i < n; ++i) {
                kn += gn[i] * wn[i];
                k1 += gt1[i] * wt1[i];
                k2 += gt2[i] * wt2[i];
                k12 += gt1[i] * wt2[i];
            }
            p.kn += std::max(kn, kLockedMobility * m.unconstrainedMobility(gn));
            p.kt1 += std::max(k1, kLockedMobility * m.unconstrainedMobility(gt1));
            p.kt2 += std::max(k2, kLockedMobility * m.unconstrainedMobility(gt2));
            p.kt12 += k12;
            p.blocks[p.nb++] = blk;
        }

    } // namespace

    void World::Impl::prepareContacts() {
        prep.clear();
        arena.clear();
        for (std::size_t ci = 0; ci < contacts.size(); ++ci) {
            ContactData &c = contacts[ci];
            Group &a = groups[std::size_t(c.ga)], &b = groups[std::size_t(c.gb)];

            // The poses have moved since the host detected: the depth is re-derived, not re-detected. Each end moves
            // with its body, turning to first order, and the normal holds, so a rolling contact stays on the ground.
            const Vec3 da = a.c - a.c0 + cross(rotationVectorOf(a.q * conjugate(a.q0)), c.rA);
            const Vec3 db = b.c - b.c0 + cross(rotationVectorOf(b.q * conjugate(b.q0)), c.rB);
            c.depth = c.depth0 - dot(db - da, c.n);
            if (c.depth < 0.0) continue; // a gap applies no force

            ContactPrep p;
            p.contact = int(ci);
            p.a = c.ga;
            p.b = c.gb;
            a.pv = a.pw = b.pv = b.pw = Vec3{};

            // A machine end's mass is its whole linkage's, through W: its rigid terms are zero.
            const bool rigidA = !a.isStatic && a.machine < 0, rigidB = !b.isStatic && b.machine < 0;
            if (!rigidA && !rigidB && a.machine < 0 && b.machine < 0) continue;
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
            p.soft = SoftParams::fromHertz(kContactHertz, kContactDamping, h);

            const Vec3 P = c.ln * c.n + c.lt1 * p.t1 + c.lt2 * p.t2;
            a.v -= p.imA * P;
            a.w -= p.invIwA * cross(p.r1, P);
            b.v += p.imB * P;
            b.w += p.invIwB * cross(p.r2, P);

            if (a.machine >= 0) addMachineBlock(*this, p, c, a.machine);
            if (b.machine >= 0 && b.machine != a.machine) addMachineBlock(*this, p, c, b.machine);
            for (int k = 0; k < p.nb; ++k) {
                const ContactBlock &blk = p.blocks[k];
                Machine &m = machines[std::size_t(blk.machine)].m;
                const double *gn = arena.data() + blk.off;
                const std::size_t n = std::size_t(blk.n);
                m.applyContact(gn + 3 * n, c.ln);
                m.applyContact(gn + 4 * n, c.lt1);
                m.applyContact(gn + 5 * n, c.lt2);
            }

            // Penetration is recovered in the position pass, where it moves the pose without leaving momentum.
            p.vtarget = 0.0;
            p.posBias = std::min(kMaxBias, p.soft.biasRate * std::max(0.0, c.depth - kSlop));
            prep.push_back(p);
        }
    }

    namespace {

        struct RowVelocity {
            double n = 0.0, t1 = 0.0, t2 = 0.0;
        };

        // The relative velocity vB - vA at the contact: a rigid end from its group, a machine end from qdot (or its
        // pseudo-velocity, in the position pass) through its rows.
        RowVelocity rowVelocity(const World::Impl &I, const ContactPrep &p, const ContactData &c, bool pseudo) {
            const Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
            Vec3 vrel{};
            if (b.machine < 0) vrel += pseudo ? b.pv + cross(b.pw, p.r2) : b.v + cross(b.w, p.r2);
            if (a.machine < 0) vrel -= pseudo ? a.pv + cross(a.pw, p.r1) : a.v + cross(a.w, p.r1);
            RowVelocity v{dot(vrel, c.n), dot(vrel, p.t1), dot(vrel, p.t2)};
            for (int k = 0; k < p.nb; ++k) {
                const ContactBlock &blk = p.blocks[k];
                const Machine &m = I.machines[std::size_t(blk.machine)].m;
                const double *gn = I.arena.data() + blk.off;
                const std::size_t n = std::size_t(blk.n);
                const std::vector<double> &qd = pseudo ? m.pseudo() : m.art.qdot;
                for (std::size_t i = 0; i < n; ++i) {
                    v.n += gn[i] * qd[i];
                    if (!pseudo) {
                        v.t1 += gn[n + i] * qd[i];
                        v.t2 += gn[2 * n + i] * qd[i];
                    }
                }
            }
            return v;
        }

        void apply(World::Impl &I, const ContactPrep &p, const ContactData &c, double ln, double lt1, double lt2) {
            Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
            a.v -= p.imA * ln * c.n;
            a.w -= ln * p.imR1cn;
            b.v += p.imB * ln * c.n;
            b.w += ln * p.imR2cn;
            const Vec3 Pt = lt1 * p.t1 + lt2 * p.t2;
            a.v -= p.imA * Pt;
            a.w -= lt1 * p.imR1ct1 + lt2 * p.imR1ct2;
            b.v += p.imB * Pt;
            b.w += lt1 * p.imR2ct1 + lt2 * p.imR2ct2;
            for (int k = 0; k < p.nb; ++k) {
                const ContactBlock &blk = p.blocks[k];
                Machine &m = I.machines[std::size_t(blk.machine)].m;
                const double *gn = I.arena.data() + blk.off;
                const std::size_t n = std::size_t(blk.n);
                if (ln != 0.0) m.applyContact(gn + 3 * n, ln);
                if (lt1 != 0.0) m.applyContact(gn + 4 * n, lt1);
                if (lt2 != 0.0) m.applyContact(gn + 5 * n, lt2);
            }
        }

        // The friction pair's new accumulated impulse, clamped to the cone. With a machine end the 2×2 block is solved
        // exactly: a wheel's mobility lies along its rolling direction, and a row-by-row step diverges on it.
        void frictionStep(const ContactPrep &p, const RowVelocity &v, double maxT, double &lt1, double &lt2) {
            double d1 = -v.t1 / p.kt1, d2 = -v.t2 / p.kt2;
            if (p.nb > 0) {
                const double a = p.kt1, b = p.kt12, d = p.kt2;
                const double det = a * d - b * b;
                if (det > 1e-12 * a * d) {
                    d1 = -(d * v.t1 - b * v.t2) / det;
                    d2 = -(a * v.t2 - b * v.t1) / det;
                }
            }
            lt1 += d1;
            lt2 += d2;
            const double len = std::sqrt(lt1 * lt1 + lt2 * lt2);
            if (len > maxT) {
                lt1 *= maxT / len;
                lt2 *= maxT / len;
            }
        }

        // One pass over a contact: soft when biased, the pure non-penetration row when relaxing.
        void sweepContact(World::Impl &I, ContactPrep &p, bool relax) {
            ContactData &c = I.contacts[std::size_t(p.contact)];
            const RowVelocity v = rowVelocity(I, p, c, false);
            const double lnn = relax ? std::max(0.0, c.ln - v.n / p.kn)
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
        }

        void sweepBounded(World::Impl &I) {
            for (MachineData &md : I.machines)
                if (md.m.boundedRows()) md.m.sweepBounded();
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
            const bool rigidA = !a.isStatic && a.machine < 0, rigidB = !b.isStatic && b.machine < 0;
            sp.imA = rigidA ? a.invMass : 0.0;
            sp.imB = rigidB ? b.invMass : 0.0;
            sp.imR1n = (rigidA ? a.invInertiaWorld() : Mat3{}) * cross(sp.r1, sp.n);
            sp.imR2n = (rigidB ? b.invInertiaWorld() : Mat3{}) * cross(sp.r2, sp.n);
            sp.k = sp.imA + sp.imB + dot(cross(sp.r1, sp.n), sp.imR1n) + dot(cross(sp.r2, sp.n), sp.imR2n);

            for (int mi : {a.machine, b.machine}) {
                if (mi < 0 || (sp.nb == 1 && sp.blocks[0].machine == mi)) continue;
                Machine &m = machines[std::size_t(mi)].m;
                const int n = m.dofs();
                ContactBlock blk{mi, arena.size(), n};
                arena.resize(blk.off + 2 * std::size_t(n), 0.0);
                double *g = arena.data() + blk.off, *w = g + n;
                if (a.machine == mi) m.art.addPointRow(a.link, pA, sp.n, 1.0, g);
                if (b.machine == mi) m.art.addPointRow(b.link, pB, sp.n, -1.0, g);
                m.mobility(g, w);
                double k = 0.0;
                for (int i = 0; i < n; ++i) k += g[i] * w[i];
                sp.k += std::max(k, kLockedMobility * m.unconstrainedMobility(g));
                sp.blocks[sp.nb++] = blk;
            }
            springRow(sp.k, d.hertz, d.damping, len - d.restLength, h, sp.gamma, sp.target);
            springs.push_back(sp);
        }
    }

    namespace {

        void sweepSpring(World::Impl &I, SpringPrep &sp) {
            Group &a = I.groups[std::size_t(sp.a)], &b = I.groups[std::size_t(sp.b)];
            double v = 0.0;
            if (a.machine < 0) v += dot(a.v + cross(a.w, sp.r1), sp.n);
            if (b.machine < 0) v -= dot(b.v + cross(b.w, sp.r2), sp.n);
            for (int k = 0; k < sp.nb; ++k) {
                const ContactBlock &blk = sp.blocks[k];
                const std::vector<double> &qd = I.machines[std::size_t(blk.machine)].m.art.qdot;
                const double *g = I.arena.data() + blk.off;
                for (std::size_t i = 0; i < std::size_t(blk.n); ++i) v += g[i] * qd[i];
            }
            const double dl = (sp.target - v - sp.gamma * sp.lambda) / (sp.k + sp.gamma);
            sp.lambda += dl;
            a.v += sp.imA * dl * sp.n;
            a.w += dl * sp.imR1n;
            b.v -= sp.imB * dl * sp.n;
            b.w -= dl * sp.imR2n;
            for (int k = 0; k < sp.nb; ++k) {
                const ContactBlock &blk = sp.blocks[k];
                I.machines[std::size_t(blk.machine)].m.applyContact(I.arena.data() + blk.off + std::size_t(blk.n), dl);
            }
        }

    } // namespace

    // The biased sweeps. The proximal anchor is re-taken each outer pass, so the impulse bleed cancels at the fixed
    // point. A machine's capped drives and limits are swept in the same loop as the contacts that load them.
    void World::Impl::solveContacts() {
        const int iterations = std::max(kMinSolverIterations, kSolverIterations / std::max(1, settings.substeps));
        const int outer = std::min(kProximalOuter, iterations);
        const int inner = std::max(1, iterations / outer);
        for (int o = 0; o < outer; ++o) {
            for (ContactPrep &p : prep) p.refLn = contacts[std::size_t(p.contact)].ln;
            for (int k = 0; k < inner; ++k) {
                for (SpringPrep &sp : springs) sweepSpring(*this, sp);
                for (ContactPrep &p : prep) sweepContact(*this, p, false);
                sweepBounded(*this);
            }
        }
    }

    // After the integration, and nothing is integrated from it: drains the separating velocity the bias added, and
    // leaves impulses the next warm start can trust.
    void World::Impl::relaxContacts() {
        for (int k = 0; k < kRelaxIterations; ++k) {
            for (ContactPrep &p : prep) sweepContact(*this, p, true);
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
                a.pv -= p.imA * app * c.n;
                a.pw -= app * p.imR1cn;
                b.pv += p.imB * app * c.n;
                b.pw += app * p.imR2cn;
                for (int i = 0; i < p.nb; ++i) {
                    const ContactBlock &blk = p.blocks[i];
                    machines[std::size_t(blk.machine)].m.applyPseudo(arena.data() + blk.off + 3 * std::size_t(blk.n),
                                                                     app);
                }
            }
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
    }

} // namespace hf
