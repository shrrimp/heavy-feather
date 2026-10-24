#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <vector>

#include "world/TestHooks.hpp"
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

        // A step's phases, timed into a StepProfile when the host asks.
        struct Laps {
            using Clock = std::chrono::steady_clock;
            bool on = false;
            Clock::time_point t{};

            explicit Laps(bool profile) : on(profile) {
                if (on) t = Clock::now();
            }
            void lap(double &into) {
                if (!on) return;
                const Clock::time_point now = Clock::now();
                into += std::chrono::duration<double, std::milli>(now - t).count();
                t = now;
            }
        };

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

    // One sub-step of a pass: the machine islands', with their machines, springs and grabs; or the rigid islands'.
    template <typename Laps>
    void World::Impl::subStep(int i, int n, Laps &laps, StepProfile &profile) {
#ifndef NDEBUG
        checkStill(groups);
#endif
        const bool machinesPass = pass == 0;
        const std::vector<int> &mine = passGroups[pass];
        if (machinesPass) {
            // At the top of sub-step i a moving target is i/n of the way along.
            const double alpha = double(i) / double(n);
            for (int mi : steppedMachines)
                for (const JointRef &r : machines[std::size_t(mi)].joints) {
                    JointSlot &j = joints[r.joint];
                    j.targetEff = j.targetFrom + (j.targetTo - j.targetFrom) * alpha;
                }
        }
        for (int gi : mine) {
            Group &g = groups[std::size_t(gi)];
            if (!rigid(g)) continue;
            g.v += h * (g.invMass * g.force + settings.gravity);
            g.w += h * (g.invInertiaWorld() * g.torque);
            gyroscopic(g, h);
        }
        updateTurned();
        if (machinesPass)
            for (int mi : steppedMachines) solveInternal(machines[std::size_t(mi)]);
        afterPhase(FaultPhase::Internal, i);
        laps.lap(profile.machines);

        if (machinesPass) {
            prepareContacts();
            prepareSprings();
            prepareGrabs();
        } else {
            prepareRigidRows();
        }
        solveContacts();
        afterPhase(FaultPhase::Contacts, i);
        if (machinesPass)
            for (int mi : steppedMachines) afterContacts(machines[std::size_t(mi)]);
        afterPhase(FaultPhase::AfterContacts, i);
        laps.lap(profile.contacts);

        if (machinesPass)
            for (int mi : steppedMachines) solvePositions(machines[std::size_t(mi)]);
        afterPhase(FaultPhase::Positions, i);
        solveContactPositions();
        for (int gi : mine)
            if (rigid(groups[std::size_t(gi)])) clampCorrection(groups[std::size_t(gi)], h, settings.maxCorrection);
        afterPhase(FaultPhase::ContactPositions, i);
        laps.lap(profile.positions);

        if (machinesPass)
            for (int mi : steppedMachines) integrate(machines[std::size_t(mi)]);
        for (int gi : mine) {
            Group &g = groups[std::size_t(gi)];
            if (rigid(g)) integrateRigid(g, h);
            else if (g.kinematic) moveKinematic(g, double(i + 1) * h, i == n - 1);
        }
        afterPhase(FaultPhase::Integrate, i);
        laps.lap(profile.integrate);

        relaxContacts();
        if (machinesPass)
            for (int mi : steppedMachines) writeVelocities(machines[std::size_t(mi)]);
        afterPhase(FaultPhase::Relax, i);
        laps.lap(profile.relax);
    }

    void World::Impl::step(double dt) {
        if (!(dt > 0.0)) return;
        const int n = std::max(1, settings.substeps);
        h = dt / n;
        StepProfile profile;
        Laps laps(settings.profile);

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

        // The host's forces, on whatever is stepped: a force on a sleeper woke it.
        for (int gi : steppedGroups) groups[std::size_t(gi)].force = groups[std::size_t(gi)].torque = Vec3{};
        std::sort(forced.begin(), forced.end());
        for (std::uint32_t i : forced) {
            BodySlot &b = bodies[i];
            b.forced = false;
            if (!b.alive || b.group < 0) continue;
            Group &g = groups[std::size_t(b.group)];
            const Vec3 d = g.q * b.offset;
            g.force += b.force;
            g.torque += b.torque + cross(d, b.force);
            b.force = b.torque = Vec3{};
        }
        forced.clear();

        for (MachineData &md : machines)
            if (stepped(md)) readState(md);
        beginFaults();
        beginDrives(dt);
        beginContacts();
        // The two passes' groups: the machine islands', the rigid islands', the kinematic ones in both.
        passGroups[0].clear();
        passGroups[1].clear();
        bool kind[2] = {false, false}; // a machine island, a rigid one
        for (std::size_t island = 0; island < machineIsland.size(); ++island)
            kind[machineIsland[island] ? 0 : 1] = true;
        for (int gi : steppedGroups) {
            const Group &g = groups[std::size_t(gi)];
            if (g.kinematic) {
                passGroups[0].push_back(gi);
                passGroups[1].push_back(gi);
            } else {
                passGroups[machineIsland[std::size_t(g.island)] ? 0 : 1].push_back(gi);
            }
        }
        laps.lap(profile.begin);

        const int substeps[2] = {n, settings.rigidSubsteps > 0 ? settings.rigidSubsteps : n};
        for (pass = 0; pass < 2; ++pass) {
            if (!kind[pass]) continue;
            h = dt / substeps[pass];
            // Each pass moves what is kinematic from where the step found it, at its own sub-steps.
            for (int gi : passGroups[pass]) {
                Group &g = groups[std::size_t(gi)];
                if (!g.kinematic) continue;
                g.c = g.c0;
                g.q = g.q0;
            }
            for (int i = 0; i < substeps[pass]; ++i) subStep(i, substeps[pass], laps, profile);
        }
        pass = 0;
        h = dt / n;
        for (int gi : steppedGroups) {
            Group &g = groups[std::size_t(gi)];
            if (g.kinematic) moveKinematic(g, dt, true);
        }

        for (int mi : steppedMachines)
            for (const JointRef &r : machines[std::size_t(mi)].joints)
                joints[r.joint].targetFrom = joints[r.joint].targetTo;
        restitution();
        for (MachineData &md : machines)
            if (stepped(md)) writeVelocities(md);
        afterPhase(FaultPhase::Restitution, n);

        if (testing::phaseHook) testing::phaseHook(*this, n, FaultPhase::EndOfStep);
        // Nothing has read this step's state yet: an island gone wrong is put back before anything does.
        checkFaults(FaultPhase::EndOfStep);
        if (!faults.empty()) rollBack();
        storeContacts();
        carriedRows.clear();
        endSleep();
        keepResting();
        endStepBodies();

        stats = WorldStats{};
        stats.faults = int(faults.size());
        stats.rejected = rejected;
        stats.groupRebuilds = groupRebuilds;
        stats.machineRebuilds = machineRebuilds;
        rejected = groupRebuilds = machineRebuilds = 0;
        stats.rigidRows = rowsOfKind[0];
        stats.machineRows = rowsOfKind[1];
        stats.contacts = int(prep.size()) + rowsOfKind[0];
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
        laps.lap(profile.end);
        profile.total = profile.begin + profile.machines + profile.contacts + profile.positions + profile.integrate +
                        profile.relax + profile.end;
        stats.profile = profile;
    }

} // namespace hf
