#include <algorithm>
#include <cstdint>
#include <vector>

#include "contact/ContactParams.hpp"
#include "world/WorldImpl.hpp"

// Mechanism contacts: rows in their machine's direct solve, beside its joints'.
namespace hf {

    namespace {

        // A contact row's tag: its key hashed, top bit set so it never meets a joint's, and the row in the low bits.
        std::uint64_t tagOf(const Contact &k, std::uint64_t slot) {
            std::uint64_t h = k.feature * 0x9E3779B97F4A7C15ull;
            h ^= ((std::uint64_t(k.a.index) << 32) | k.b.index) * 0xC2B2AE3D27D4EB4Full;
            h ^= ((std::uint64_t(k.a.generation) << 32) | k.b.generation) * 0x165667B19E3779F9ull;
            return (std::uint64_t(1) << 63) | (h >> 5 << 4) | slot;
        }

        // Each end's link (-1: static) and its point, moving with it as it did when the step began.
        struct Ends {
            int la = -1, lb = -1;
            Vec3 pA{}, pB{};
        };

        Ends endsOf(const World::Impl &I, const ContactData &c) {
            const Group &a = I.groups[std::size_t(c.ga)], &b = I.groups[std::size_t(c.gb)];
            return {a.machine >= 0 ? a.link : -1, b.machine >= 0 ? b.link : -1, a.c + c.rA, b.c + c.rB};
        }

        // B's point less A's, along d; or B's turn less A's, about d.
        void pointRow(const Articulation &art, const Ends &e, const Vec3 &d, double *J) {
            art.addPointRow(e.lb, e.pB, d, 1.0, J);
            art.addPointRow(e.la, e.pA, d, -1.0, J);
        }

        void turnRow(const Articulation &art, const Ends &e, const Vec3 &d, double *J) {
            art.addAngularRow(e.lb, d, 1.0, J);
            art.addAngularRow(e.la, d, -1.0, J);
        }

    } // namespace

    // A push-only normal, soft like a swept contact's but anchored at its last answer, so it holds rigidly once its
    // load settles; then friction rows bounded by mu times whatever the normal's answer turns out to be.
    void World::Impl::contactRows(MachineData &md) {
        constexpr double kFrictionYield = 1e-2; // rigid, a contact's redundant friction rows set the pivoting cycling
        constexpr double kSlideSpeed = 1e-2;    // m/s
        Machine &m = md.m;
        const SoftParams soft = SoftParams::fromHertz(kContactHertz, kContactDamping, h);
        const double softness = soft.impulseCoeff / soft.massCoeff; // compliance over the row's own mobility
        for (int ci : md.contacts) {
            ContactData &c = contacts[std::size_t(ci)];
            const Ends e = endsOf(*this, c);
            const auto velocityAt = [&](int link, const Vec3 &x) {
                return link < 0 ? Vec3{} : pointVelocity(m.art.twistOf(link, m.art.qdot.data()), x - m.art.O);
            };
            const Vec3 v = velocityAt(e.lb, e.pB) - velocityAt(e.la, e.pA);
            c.depth = contactDepth(c);
            c.row = -1;
            const bool gap = c.depth < 0.0;
            if (gap && c.lnMax <= 0.0 && -dot(v, c.n) * h * kGapReach < -c.depth) continue;

            // A box of two friction rows reaches √2 mu off its axes. Sliding, it's one row against the slip; holding,
            // a pair the sweeps keep in its circle.
            Vec3 T1, T2;
            contactBasis(c.n, T1, T2);
            const Vec3 F = c.lt1 * T1 + c.lt2 * T2, R = c.lr1 * T1 + c.lr2 * T2;
            const Vec3 slide = v - dot(v, c.n) * c.n;
            const bool sliding = length(slide) > kSlideSpeed;
            c.t1 = sliding ? normalize(slide) : T1;
            c.t2 = cross(c.n, c.t1);

            const auto add = [&](MachineRow row, std::uint64_t slot, bool turn, const Vec3 &d) {
                md.velTags.push_back(tagOf(c.in, slot));
                double *J = m.addRow(row);
                if (turn) turnRow(m.art, e, d, J);
                else pointRow(m.art, e, d, J);
            };
            const auto following = [&](double ratio, double anchor, bool paired = false) {
                MachineRow row;
                row.follows = c.row;
                row.ratio = ratio;
                row.pairedWithNext = paired;
                row.relEps = kFrictionYield;
                row.anchor = gap ? 0.0 : anchor;
                return row;
            };
            MachineRow normal;
            normal.lo = 0.0;
            normal.target = gap ? c.depth / h : 0.0;
            if (!gap) {
                normal.relEps = softness;
                normal.anchor = c.ln;
            }
            c.row = m.rowCount();
            add(normal, 0, false, c.n);
            add(following(c.in.friction, dot(F, c.t1), !sliding), 1, false, c.t1);
            add(following(sliding ? 0.0 : c.in.friction, dot(F, c.t2)), 2, false, c.t2);
            if (c.in.rolling > 0.0 || c.in.spinning > 0.0) {
                add(following(c.in.rolling, dot(R, c.t1), true), 3, true, c.t1);
                add(following(c.in.rolling, dot(R, c.t2)), 4, true, c.t2);
                add(following(c.in.spinning, c.ls), 5, true, c.n);
            }
        }
    }

    // Penetration past the slop, pushed out in the position pass: the pose moves and no momentum is left.
    void World::Impl::contactPositionRows(MachineData &md) {
        const double rate = SoftParams::fromHertz(kContactHertz, kContactDamping, h).biasRate;
        for (int ci : md.contacts) {
            const ContactData &c = contacts[std::size_t(ci)];
            if (c.row < 0 || c.depth <= kSlop) continue;
            MachineRow row;
            row.lo = 0.0;
            row.target = std::min(kMaxBias, rate * (c.depth - kSlop));
            md.posTags.push_back(tagOf(c.in, 8));
            pointRow(md.m.art, endsOf(*this, c), c.n, md.m.addPosRow(row));
        }
    }

    // The answers, kept in the contact's own basis: the next sub-step's anchors, and the next step's warm start.
    void World::Impl::recordContacts(MachineData &md) {
        const std::vector<double> &l = md.m.lambda();
        for (int ci : md.contacts) {
            ContactData &c = contacts[std::size_t(ci)];
            if (c.row < 0) continue;
            const std::size_t r = std::size_t(c.row);
            Vec3 T1, T2;
            contactBasis(c.n, T1, T2);
            const Vec3 F = l[r + 1] * c.t1 + l[r + 2] * c.t2;
            c.ln = l[r];
            c.lt1 = dot(F, T1);
            c.lt2 = dot(F, T2);
            if (c.in.rolling > 0.0 || c.in.spinning > 0.0) {
                const Vec3 R = l[r + 3] * c.t1 + l[r + 4] * c.t2;
                c.lr1 = dot(R, T1);
                c.lr2 = dot(R, T2);
                c.ls = l[r + 5];
            }
            c.lnMax = std::max(c.lnMax, c.ln);
        }
    }

} // namespace hf
