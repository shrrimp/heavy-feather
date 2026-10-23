#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        bool rigid(const Group &g) { return !g.isStatic && !g.frozen && g.machine < 0; }

        // The spin's own torque, ω × Iω, at the step's midpoint: I(ω' − ω) + h ω_m × Iω_m = 0 in the body frame,
        // solved by Newton. The midpoint keeps a free body's energy and angular momentum; taken at the start, a thin
        // body spinning near its long axis gains energy every step until it explodes, and left out it never precesses.
        //
        // A Newton step with no inverse, or none worth trusting, leaves the spin as it was for the sub-step.
        void gyroscopic(Group &g, double h) {
            const Mat3 R = toMat3(g.q);
            const Mat3 &I = g.inertiaBody;
            const Vec3 w = transpose(R) * g.w;
            const double tiny = 1e-12 * std::fabs(determinant(I));
            Vec3 next = w;
            for (int k = 0; k < 2; ++k) {
                const Vec3 m = 0.5 * (w + next), Im = I * m;
                const Vec3 f = I * (next - w) + h * cross(m, Im);
                const Mat3 J = I + 0.5 * h * (skew(m) * I - skew(Im));
                const double det = determinant(J);
                if (!std::isfinite(det) || !(std::fabs(det) > tiny)) return;
                next = next - inverse(J) * f;
            }
            g.w = R * next;
        }

#ifndef NDEBUG
        // What doesn't move this step is only ever read: a row that wrote into it would carry its numbers, a NaN
        // among them, to every island that touches it.
        void checkStill(const std::vector<Group> &groups) {
            for (const Group &g : groups) {
                if (moves(g) || g.kinematic) continue;
                assert(g.v == Vec3{} && g.w == Vec3{} && g.pv == Vec3{} && g.pw == Vec3{});
            }
        }
#endif

        void integrateRigid(Group &g, double h) {
            g.correction = (length(g.pv) + length(g.pw) * g.radius) * h;
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

        build();
        // A machine the host moved starts its rows cold.
        for (std::uint32_t i : written) {
            const BodySlot &b = bodies[i];
            if (!b.alive || !b.teleported || b.group < 0) continue;
            const int m = groups[std::size_t(b.group)].machine;
            if (m < 0) continue;
            machines[std::size_t(m)].lastVel.clear();
            machines[std::size_t(m)].lastPos.clear();
        }
        beginKinematic(dt);
        beginSleep();
        const auto stepped = [](const MachineData &md) { return !md.frozen; };

        for (Group &g : groups) g.force = g.torque = Vec3{};
        for (BodySlot &b : bodies) {
            if (!b.alive) continue;
            Group &g = groups[std::size_t(b.group)];
            const Vec3 d = g.q * b.offset;
            g.force += b.force;
            g.torque += b.torque + cross(d, b.force);
            b.force = b.torque = Vec3{};
        }

        for (MachineData &md : machines)
            if (stepped(md)) readState(md);
        beginFaults();
        beginDrives(dt);
        beginContacts();

        for (int i = 0; i < n; ++i) {
#ifndef NDEBUG
            checkStill(groups);
#endif
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
            for (MachineData &md : machines)
                if (stepped(md)) solveInternal(md);
            afterPhase(FaultPhase::Internal, i);

            prepareContacts();
            prepareSprings();
            solveContacts();
            afterPhase(FaultPhase::Contacts, i);
            for (MachineData &md : machines)
                if (stepped(md)) afterContacts(md);
            afterPhase(FaultPhase::AfterContacts, i);

            for (MachineData &md : machines)
                if (stepped(md)) solvePositions(md);
            afterPhase(FaultPhase::Positions, i);
            solveContactPositions();
            for (Group &g : groups)
                if (rigid(g)) clampCorrection(g, h, settings.maxCorrection);
            afterPhase(FaultPhase::ContactPositions, i);

            for (MachineData &md : machines)
                if (stepped(md)) integrate(md);
            for (Group &g : groups) {
                if (rigid(g)) integrateRigid(g, h);
                else if (g.kinematic) moveKinematic(g, double(i + 1) * h, i == n - 1);
            }
            afterPhase(FaultPhase::Integrate, i);

            relaxContacts();
            for (MachineData &md : machines)
                if (stepped(md)) writeVelocities(md);
            afterPhase(FaultPhase::Relax, i);
        }

        for (JointSlot &j : joints)
            if (j.alive) j.targetFrom = j.targetTo;
        restitution();
        for (MachineData &md : machines)
            if (stepped(md)) writeVelocities(md);
        afterPhase(FaultPhase::Restitution, n);

        // Nothing has read this step's state yet: an island gone wrong is put back before anything does.
        checkFaults(FaultPhase::EndOfStep);
        if (!faults.empty()) rollBack();
        storeContacts();
        carriedRows.clear();
        endSleep();
        endStepBodies();

        stats = WorldStats{};
        stats.faults = int(faults.size());
        stats.rejected = rejected;
        stats.groupRebuilds = groupRebuilds;
        stats.machineRebuilds = machineRebuilds;
        rejected = groupRebuilds = machineRebuilds = 0;
        stats.contacts = int(prep.size());
        stats.machines = int(machines.size());
        for (const MachineData &md : machines) {
            stats.links += int(md.groupOfLink.size());
            stats.dofs += md.m.dofs();
            if (!stepped(md)) continue;
            for (int ci : md.contacts) stats.contacts += contacts[std::size_t(ci)].row >= 0;
            stats.rows += md.m.rowCount();
            stats.lcpPivots += md.m.stats.lcpPivots;
            stats.lcpFallbacks += md.m.stats.lcpFallbacks;
        }
        for (BodySlot &b : bodies) {
            if (!b.alive) continue;
            const Group &g = groups[std::size_t(b.group)];
            stats.stepped += moves(g);
            stats.asleep += b.asleep;
            if (moves(g) && !faulted(g)) b.faultStreak = 0;
        }
        stats.islands = islands;
    }

} // namespace hf
