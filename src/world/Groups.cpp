#include <cmath>
#include <numeric>
#include <vector>

#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        // Half the diagonal of the box with this inertia: a fair bound on a body's extent about its centre.
        double radiusFromInertia(double mass, const Mat3 &I) {
            if (!(mass > 0.0)) return 0.0;
            return std::sqrt(1.5 * (I[0][0] + I[1][1] + I[2][2]) / mass);
        }

    } // namespace

    void World::Impl::bodyState(const BodySlot &b, Vec3 &c, Quat &q, Vec3 &v, Vec3 &w) const {
        if (b.group < 0) {
            c = b.desc.position;
            q = b.desc.orientation;
            v = b.desc.velocity;
            w = b.desc.angularVelocity;
            return;
        }
        const Group &g = groups[std::size_t(b.group)];
        q = normalize(g.q * b.rot);
        c = g.c + g.q * b.offset;
        w = g.w;
        v = g.v + cross(g.w, c - g.c);
    }

    // A weld holds its frames together: each member is placed from the one it is welded to, out from what can't
    // move (the world, a static body) or else from the group's first member.
    template <typename Find>
    void World::Impl::placeWelds(Find &find) {
        const std::size_t n = bodies.size();
        std::vector<std::uint8_t> placed(n, 0), rooted(n, 0);
        const auto place = [&](BodySlot &b, const Frame &onB, const Vec3 &at, const Quat &q) {
            b.desc.orientation = normalize(q * conjugate(onB.rotation));
            b.desc.position = at - b.desc.orientation * onB.position;
        };
        for (const JointSlot &j : joints) {
            if (!j.alive || j.desc.type != JointType::Fixed || (j.desc.a.valid() && j.desc.b.valid())) continue;
            const bool aBody = j.desc.a.valid();
            const std::uint32_t i = (aBody ? j.desc.a : j.desc.b).index;
            const Frame &world = aBody ? j.desc.frameB : j.desc.frameA;
            if (!bodies[i].desc.isStatic)
                place(bodies[i], aBody ? j.desc.frameA : j.desc.frameB, world.position, world.rotation);
            placed[i] = rooted[find(i)] = 1;
        }
        for (std::size_t i = 0; i < n; ++i)
            if (bodies[i].alive && bodies[i].desc.isStatic) placed[i] = rooted[find(i)] = 1;
        for (std::size_t i = 0; i < n; ++i)
            if (bodies[i].alive && !rooted[find(i)]) placed[i] = rooted[find(i)] = 1;

        for (bool changed = true; changed;) {
            changed = false;
            for (const JointSlot &j : joints) {
                if (!j.alive || j.desc.type != JointType::Fixed || !j.desc.a.valid() || !j.desc.b.valid()) continue;
                const std::uint32_t a = j.desc.a.index, b = j.desc.b.index;
                if (placed[a] == placed[b]) continue;
                const bool fromA = placed[a] != 0;
                const BodySlot &from = bodies[fromA ? a : b];
                const Frame &f = fromA ? j.desc.frameA : j.desc.frameB;
                place(bodies[fromA ? b : a], fromA ? j.desc.frameB : j.desc.frameA,
                      from.desc.position + from.desc.orientation * f.position, from.desc.orientation * f.rotation);
                placed[fromA ? b : a] = 1;
                changed = true;
            }
        }
    }

    void World::Impl::rebuildGroups() {
        ++groupRebuilds;
        // Where the step left each body, unless the host has written it since.
        for (BodySlot &b : bodies) {
            if (!b.alive || b.group < 0) continue;
            Vec3 c, v, w;
            Quat q;
            bodyState(b, c, q, v, w);
            if (!b.poseWritten) {
                b.desc.position = c;
                b.desc.orientation = q;
            }
            if (!b.velocityWritten) {
                b.desc.velocity = v;
                b.desc.angularVelocity = w;
            }
        }

        // Welds: union-find over Fixed joints. A body welded to the world is static.
        const std::size_t n = bodies.size();
        std::vector<std::size_t> uf(n);
        std::iota(uf.begin(), uf.end(), std::size_t(0));
        const auto find = [&](std::size_t x) {
            while (uf[x] != x) x = uf[x] = uf[uf[x]];
            return x;
        };
        std::vector<std::uint8_t> toWorld(n, 0);
        for (const JointSlot &j : joints) {
            if (!j.alive || j.desc.type != JointType::Fixed) continue;
            if (!j.desc.a.valid() || !j.desc.b.valid()) {
                toWorld[(j.desc.a.valid() ? j.desc.a : j.desc.b).index] = 1;
                continue;
            }
            const std::size_t a = find(j.desc.a.index), b = find(j.desc.b.index);
            if (a != b) uf[std::max(a, b)] = std::min(a, b);
        }
        for (std::size_t i = 0; i < n; ++i)
            if (toWorld[i]) toWorld[find(i)] = 1;
        placeWelds(find);

        groups.assign(1, Group{});
        std::vector<int> groupOfRoot(n, -1);
        for (std::size_t i = 0; i < n; ++i) {
            BodySlot &b = bodies[i];
            if (!b.alive) continue;
            const std::size_t r = find(i);
            if (groupOfRoot[r] < 0) {
                groupOfRoot[r] = int(groups.size());
                Group g;
                g.q = b.desc.orientation; // a group's frame is its first member's
                g.c = b.desc.position;
                g.isStatic = toWorld[r] != 0;
                groups.push_back(g);
            }
            b.group = groupOfRoot[r];
        }

        // Mass and centre of mass.
        std::vector<double> M(groups.size(), 0.0);
        std::vector<Vec3> mc(groups.size());
        for (const BodySlot &b : bodies) {
            if (!b.alive) continue;
            Group &g = groups[std::size_t(b.group)];
            if (b.desc.isStatic || !(b.desc.mass > 0.0)) g.isStatic = true;
            if (b.kinematic) g.kinematic = true;
            if (b.desc.mass > 0.0) {
                M[std::size_t(b.group)] += b.desc.mass;
                mc[std::size_t(b.group)] += b.desc.mass * b.desc.position;
            }
        }
        for (std::size_t gi = 1; gi < groups.size(); ++gi)
            if (M[gi] > 0.0) groups[gi].c = mc[gi] / M[gi];

        // Inertia, momentum and extent about that centre.
        std::vector<Mat3> Iw(groups.size());
        std::vector<Vec3> P(groups.size()), L(groups.size());
        for (const BodySlot &b : bodies) {
            if (!b.alive) continue;
            const std::size_t gi = std::size_t(b.group);
            Group &g = groups[gi];
            const double m = std::max(b.desc.mass, 0.0);
            const Vec3 d = b.desc.position - g.c;
            const Mat3 R = toMat3(b.desc.orientation);
            const Mat3 Ib = R * b.desc.inertia * transpose(R);
            Iw[gi] = Iw[gi] + Ib + m * (dot(d, d) * Mat3(1.0) - outerProduct(d, d));
            P[gi] += m * b.desc.velocity;
            L[gi] += Ib * b.desc.angularVelocity + m * cross(d, b.desc.velocity);
            const double r = b.desc.radius > 0.0 ? b.desc.radius : radiusFromInertia(m, b.desc.inertia);
            g.radius = std::max(g.radius, length(d) + r);
        }

        for (std::size_t gi = 1; gi < groups.size(); ++gi) {
            Group &g = groups[gi];
            if (g.isStatic) continue;
            g.invMass = 1.0 / M[gi];
            const Mat3 R = toMat3(g.q);
            const Mat3 Ibody = transpose(R) * Iw[gi] * R;
            const bool solid = std::fabs(determinant(Ibody)) > 1e-300;
            g.inertiaBody = solid ? Ibody : Mat3(1e-3 * M[gi]);
            g.invInertiaBody = inverse(g.inertiaBody);
            g.v = P[gi] / M[gi];
            g.w = std::fabs(determinant(Iw[gi])) > 1e-300 ? inverse(Iw[gi]) * L[gi] : Vec3{};
        }
        // A kinematic group keeps its mass for when it is let go, but is static to the solver until then.
        for (Group &g : groups) {
            if (!g.kinematic) continue;
            g.isStatic = true;
            g.v = g.w = Vec3{};
            g.targetC = g.c;
            g.targetQ = g.q;
        }

        for (BodySlot &b : bodies) {
            if (!b.alive) continue;
            const Group &g = groups[std::size_t(b.group)];
            b.offset = conjugate(g.q) * (b.desc.position - g.c);
            b.rot = normalize(conjugate(g.q) * b.desc.orientation);
        }
    }

} // namespace hf
