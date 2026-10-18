#include <algorithm>
#include <cmath>

#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        bool rigid(const Group &g) { return !g.isStatic && g.machine < 0; }

        void integrateRigid(Group &g, double h) {
            g.c += (g.v + g.pv) * h;
            const Vec3 w = g.w + g.pw;
            g.q = normalize(g.q + Quat{0.0, 0.5 * h * w} * g.q);
            g.pv = g.pw = Vec3{};
        }

        // The ceiling on how far a position correction may move a body in one sub-step: scaled, not zeroed, so a
        // correction that is merely too fast arrives late rather than not at all.
        void clampCorrection(Group &g, double h, double maxMove) {
            const double move = (length(g.pv) + length(g.pw) * g.radius) * h;
            if (!std::isfinite(move)) {
                g.pv = g.pw = Vec3{};
            } else if (move > maxMove) {
                g.pv *= maxMove / move;
                g.pw *= maxMove / move;
            }
        }

    } // namespace

    void World::Impl::step(double dt) {
        if (!(dt > 0.0)) return;
        const int n = std::max(1, settings.substeps);
        h = dt / n;

        if (dirty) {
            rebuildGroups();
            buildMachines();
            dirty = false;
        }

        for (Group &g : groups) g.force = g.torque = Vec3{};
        for (BodySlot &b : bodies) {
            if (!b.alive) continue;
            Group &g = groups[std::size_t(b.group)];
            const Vec3 d = g.q * b.offset;
            g.force += b.force;
            g.torque += b.torque + cross(d, b.force);
            b.force = b.torque = Vec3{};
        }

        for (MachineData &md : machines) readState(md);
        beginDrives(dt);
        beginContacts();

        for (int i = 0; i < n; ++i) {
            // At the top of sub-step i a moving target is i/n of the way along.
            const double alpha = double(i) / double(n);
            for (JointSlot &j : joints)
                if (j.alive) j.targetEff = j.targetFrom + (j.targetTo - j.targetFrom) * alpha;

            for (Group &g : groups) {
                if (!rigid(g)) continue;
                g.v += h * (g.invMass * g.force + settings.gravity);
                g.w += h * (g.invInertiaWorld() * g.torque);
            }
            for (MachineData &md : machines) solveInternal(md);

            prepareContacts();
            prepareSprings();
            solveContacts();
            for (MachineData &md : machines) afterContacts(md);

            for (MachineData &md : machines) solvePositions(md);
            solveContactPositions();
            for (Group &g : groups)
                if (rigid(g)) clampCorrection(g, h, settings.maxCorrection);

            for (MachineData &md : machines) integrate(md);
            for (Group &g : groups)
                if (rigid(g)) integrateRigid(g, h);

            relaxContacts();
            for (MachineData &md : machines) writeVelocities(md);
        }

        for (JointSlot &j : joints)
            if (j.alive) j.targetFrom = j.targetTo;
        restitution();
        for (MachineData &md : machines) writeVelocities(md);
        storeContacts();

        stats = WorldStats{};
        stats.contacts = int(prep.size());
        stats.machines = int(machines.size());
        for (const MachineData &md : machines) {
            stats.links += int(md.groupOfLink.size());
            stats.dofs += md.m.dofs();
            stats.rows += md.m.rowCount();
            stats.lcpPivots += md.m.stats.lcpPivots;
            stats.lcpFallbacks += md.m.stats.lcpFallbacks;
        }
    }

} // namespace hf
