#include <algorithm>
#include <cmath>

#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        bool rigid(const Group &g) { return !g.isStatic && g.machine < 0; }

        // The spin's own torque, ω × Iω, at the step's midpoint: I(ω' − ω) + h ω_m × Iω_m = 0 in the body frame,
        // solved by Newton. The midpoint keeps a free body's energy and angular momentum; taken at the start, a thin
        // body spinning near its long axis gains energy every step until it explodes, and left out it never precesses.
        void gyroscopic(Group &g, double h) {
            const Mat3 R = toMat3(g.q);
            const Mat3 &I = g.inertiaBody;
            const Vec3 w = transpose(R) * g.w;
            Vec3 next = w;
            for (int k = 0; k < 2; ++k) {
                const Vec3 m = 0.5 * (w + next), Im = I * m;
                const Vec3 f = I * (next - w) + h * cross(m, Im);
                const Mat3 J = I + 0.5 * h * (skew(m) * I - skew(Im));
                next = next - inverse(J) * f;
            }
            g.w = R * next;
        }

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
                gyroscopic(g, h);
            }
            updateTurned();
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
        carriedRows.clear();

        stats = WorldStats{};
        stats.contacts = int(prep.size());
        stats.machines = int(machines.size());
        for (const MachineData &md : machines) {
            for (int ci : md.contacts) stats.contacts += contacts[std::size_t(ci)].row >= 0;
            stats.links += int(md.groupOfLink.size());
            stats.dofs += md.m.dofs();
            stats.rows += md.m.rowCount();
            stats.lcpPivots += md.m.stats.lcpPivots;
            stats.lcpFallbacks += md.m.stats.lcpFallbacks;
        }
    }

} // namespace hf
