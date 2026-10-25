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
            b.desc.position = at - b.desc.orientation * (onB.position - b.desc.centerOfMass);
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
                      from.desc.position + from.desc.orientation * (f.position - from.desc.centerOfMass),
                      from.desc.orientation * f.rotation);
                placed[fromA ? b : a] = 1;
                changed = true;
            }
        }
    }

    int World::Impl::addGroup(BodySlot &b) {
        int gi;
        if (!freeGroups.empty()) {
            gi = freeGroups.back();
            freeGroups.pop_back();
        } else {
            gi = int(groups.size());
            groups.emplace_back();
        }
        Group &g = groups[std::size_t(gi)];
        g = Group{};
        g.c = g.c0 = b.desc.position;
        g.q = g.q0 = b.desc.orientation;
        g.members = 1;
        b.group = gi;
        b.offset = Vec3{};
        b.rot = Quat{};
        setGroupMass(gi, b);
        if (!g.isStatic) {
            g.v = b.desc.velocity;
            g.w = b.desc.angularVelocity;
        }
        return gi;
    }

    void World::Impl::freeGroup(int gi) {
        groups[std::size_t(gi)] = Group{};
        freeGroups.push_back(gi);
    }

    // A group of one body: its mass, inertia and extent are the body's, in the body's own axes.
    void World::Impl::setGroupMass(int gi, const BodySlot &b) {
        Group &g = groups[std::size_t(gi)];
        g.isStatic = b.desc.isStatic || !(b.desc.mass > 0.0) || b.kinematic;
        g.kinematic = b.kinematic;
        const double m = std::max(b.desc.mass, 0.0);
        g.radius = b.desc.radius > 0.0 ? b.desc.radius : radiusFromInertia(m, b.desc.inertia);
        if (b.desc.isStatic || !(b.desc.mass > 0.0)) {
            g.invMass = 0.0;
            g.v = g.w = Vec3{};
            return;
        }
        g.invMass = 1.0 / b.desc.mass;
        g.inertiaBody = b.desc.inertia;
        g.fixedRotation = b.desc.fixedRotation;
        g.invInertiaBody = g.fixedRotation ? Mat3{} : inverse(g.inertiaBody);
        if (g.fixedRotation) g.w = Vec3{};
    }

    void World::Impl::touch(int gi) {
        if (gi <= 0) return;
        touched.push_back(gi);
    }

    void World::Impl::touchJoint(const JointDesc &d) {
        if (d.type == JointType::Fixed) {
            dirty = true;
            return;
        }
        for (BodyId id : {d.a, d.b}) {
            const BodySlot *b = body(id);
            if (!b) continue;
            if (b->group < 0) dirty = true;
            else touch(b->group);
        }
    }

    void World::Impl::touchJointsOf(BodyId id) {
        for (const JointSlot &j : joints)
            if (j.alive && (j.desc.a == id || j.desc.b == id)) touchJoint(j.desc);
        const BodySlot *b = body(id);
        if (b && b->group > 0 && groups[std::size_t(b->group)].machine >= 0)
            touchedMachines.push_back(groups[std::size_t(b->group)].machine);
    }

    // The centre of mass moved within the body, its origin where it was. Inside the world the body is placed by its
    // centre, so the centre moves, and everything measured from it, joints' ends and grabs, shifts back to stay put.
    void World::Impl::moveCenterOfMass(BodySlot &b, const Vec3 &centerOfMass) {
        Vec3 c, v, w;
        Quat q;
        bodyState(b, c, q, v, w);
        const Vec3 shift = centerOfMass - b.desc.centerOfMass, moved = q * shift;
        b.desc.centerOfMass = centerOfMass;
        b.desc.position = c + moved;
        b.desc.velocity = v + cross(w, moved);
        const BodyId id{std::uint32_t(&b - bodies.data()), b.generation};
        for (GrabSlot &g : grabs)
            if (g.alive && g.body == id) g.local -= shift;
        if (dirty || b.group < 0 || groups[std::size_t(b.group)].members > 1) {
            // A weld group's centre is all its bodies': rebuilt, from the body as it is now.
            b.poseWritten = b.velocityWritten = true;
            noteWritten(b);
            dirty = true;
            return;
        }
        const int gi = b.group;
        Group &g = groups[std::size_t(gi)];
        g.c += moved;
        g.c0 = g.c;
        if (!g.isStatic) g.v += cross(g.w, moved);
        const Vec3 inGroup = b.rot * shift;
        const auto shiftEnd = [&](End &e) {
            if (e.group == gi) e.anchor -= inGroup;
        };
        for (JointRef &r : looseSprings) {
            shiftEnd(r.A);
            shiftEnd(r.B);
        }
        for (MachineData &md : machines) {
            for (JointRef &r : md.joints) {
                shiftEnd(r.A);
                shiftEnd(r.B);
            }
            Articulation &art = md.m.art;
            for (std::size_t l = 0; l < art.links.size(); ++l) {
                Link &L = art.links[l];
                // This link's own end of its joint, and its children's ends on it, both measured from its centre.
                if (md.groupOfLink[l] == gi) {
                    L.pc -= inGroup;
                    L.c = g.c;
                }
                if (L.parent >= 0 && md.groupOfLink[std::size_t(L.parent)] == gi) L.pp -= inGroup;
            }
            // Rooted on it, as on the ground: its root ends shift, and its pose as an anchor is where it now is.
            for (Anchor &a : md.anchors) {
                if (a.group != gi) continue;
                if (a.link >= 0) {
                    Link &L = art.links[std::size_t(a.link)];
                    L.pp -= inGroup;
                    L.anchorPos = g.c;
                } else {
                    a.c = g.c;
                }
            }
        }
    }

    void World::Impl::bodyMassChanged(BodySlot &b) {
        if (dirty || b.group < 0 || groups[std::size_t(b.group)].members > 1) {
            dirty = true;
            return;
        }
        Group &g = groups[std::size_t(b.group)];
        const bool wasStatic = g.isStatic;
        const Vec3 v = g.v, w = g.w;
        setGroupMass(b.group, b);
        if (g.isStatic != wasStatic) {
            // Static to moving or back: its joints make it a part of a machine, or no longer.
            if (!g.isStatic) {
                g.v = b.velocityWritten ? b.desc.velocity : Vec3{};
                g.w = b.velocityWritten ? b.desc.angularVelocity : Vec3{};
            }
            touchJointsOf({std::uint32_t(&b - bodies.data()), b.generation});
            return;
        }
        if (!g.isStatic) {
            g.v = v;
            g.w = w;
        }
        // A machine's link carries on with its new mass, the machine unbuilt.
        if (g.machine >= 0) {
            Link &L = machines[std::size_t(g.machine)].m.art.links[std::size_t(g.link)];
            L.mass = b.desc.mass;
            L.inertiaBody = g.inertiaBody;
        }
    }

    void World::Impl::rebuildGroups() {
        ++groupRebuilds;
        freeGroups.clear();
        touched.clear();
        touchedMachines.clear();
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
            ++groups[std::size_t(b.group)].members;
        }

        // Mass and centre of mass.
        std::vector<double> M(groups.size(), 0.0);
        std::vector<Vec3> mc(groups.size());
        for (const BodySlot &b : bodies) {
            if (!b.alive) continue;
            Group &g = groups[std::size_t(b.group)];
            if (b.desc.isStatic || !(b.desc.mass > 0.0)) g.isStatic = true;
            if (b.kinematic) g.kinematic = true;
            if (b.desc.fixedRotation) g.fixedRotation = true;
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
            g.invInertiaBody = g.fixedRotation ? Mat3{} : inverse(g.inertiaBody);
            g.v = P[gi] / M[gi];
            g.w = !g.fixedRotation && std::fabs(determinant(Iw[gi])) > 1e-300 ? inverse(Iw[gi]) * L[gi] : Vec3{};
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
