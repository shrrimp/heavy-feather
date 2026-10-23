#include <cmath>
#include <vector>

#include "world/TestHooks.hpp"
#include "world/WorldImpl.hpp"

// The tripwire. A row only ever writes into the moving groups of its own island, so a NaN born in an island stays in
// it until the step is done. Then, before the host or any other island reads it, the island is put back as the step
// found it: its poses, at rest, its drives as they were, its contacts cold. It stays awake, and the host is told.
namespace hf {

    namespace testing {
        PhaseHook phaseHook = nullptr;
    } // namespace testing

    namespace {

        bool finite(double x) { return std::isfinite(x); }
        bool finite(const Vec3 &v) { return finite(v.x) && finite(v.y) && finite(v.z); }
        bool finite(const Quat &q) { return finite(q.w) && finite(q.x) && finite(q.y) && finite(q.z); }

        bool finite(const Group &g) {
            return finite(g.c) && finite(g.q) && finite(g.v) && finite(g.w) && finite(g.pv) && finite(g.pw);
        }

        // What a step changes in a drive.
        bool finite(const JointSlot &j) {
            return finite(j.lambda) && finite(j.desc.target) && finite(j.targetFrom) && finite(j.targetTo) &&
                   finite(j.targetEff) && finite(j.holdPos);
        }

        bool finite(const ContactData &c) {
            return finite(c.ln) && finite(c.lt1) && finite(c.lt2) && finite(c.lr1) && finite(c.lr2) && finite(c.ls);
        }

    } // namespace

    void World::Impl::beginFaults() {
        faults.clear();
        faultSites.clear();
        islandFaulted.assign(std::size_t(islands), 0);
        driveSaves.clear();
        for (const MachineData &md : machines) {
            if (md.frozen) continue;
            for (const JointRef &r : md.joints) driveSaves.push_back({r.joint, joints[r.joint]});
        }
    }

    void World::Impl::afterPhase(FaultPhase phase, int substep) {
        if (testing::phaseHook) testing::phaseHook(*this, substep, phase);
        if (settings.checkPhases) checkFaults(phase);
    }

    void World::Impl::checkFaults(FaultPhase phase) {
        const bool end = phase == FaultPhase::EndOfStep;
        // One fault an island: the first thing found wrong in it.
        const auto fault = [&](int island, FaultKind kind, int group, BodyId body, int machine) {
            if (island < 0 || islandFaulted[std::size_t(island)]) return;
            islandFaulted[std::size_t(island)] = 1;
            Fault f;
            f.step = stepCount;
            f.kind = kind;
            f.phase = phase;
            f.body = body;
            f.machine = machine;
            faults.push_back(f);
            faultSites.push_back({island, group});
        };

        for (std::size_t gi = 0; gi < groups.size(); ++gi) {
            const Group &g = groups[gi];
            if (!moves(g) || faulted(g)) continue;
            if (!finite(g)) fault(g.island, FaultKind::NonFinite, int(gi), BodyId{}, machineOf(g));
            else if (end && length(g.v) + length(g.w) * g.radius > settings.faultSpeed)
                fault(g.island, FaultKind::Speed, int(gi), BodyId{}, machineOf(g));
        }
        for (std::size_t mi = 0; mi < machines.size(); ++mi) {
            const MachineData &md = machines[mi];
            if (md.frozen || md.groupOfLink.empty()) continue;
            const int island = groups[std::size_t(md.groupOfLink[0])].island;
            if (islandFaulted[std::size_t(island)]) continue;
            // The link whose rates went wrong first: mid-step, before they spread through the machine, the one that
            // went wrong.
            const Articulation &art = md.m.art;
            for (std::size_t l = 0; l < art.links.size(); ++l)
                for (int k = 0; k < art.links[l].ndof; ++k)
                    if (!finite(art.qdot[std::size_t(art.links[l].dof0 + k)]))
                        fault(island, FaultKind::NonFinite, md.groupOfLink[l], BodyId{}, int(mi));
            for (const JointRef &r : md.joints) {
                const JointDesc &d = joints[r.joint].desc;
                if (!finite(joints[r.joint])) fault(island, FaultKind::NonFinite, -1, d.b.valid() ? d.b : d.a, int(mi));
            }
        }
        if (!end) return;
        for (const ContactData &c : contacts) {
            if (finite(c)) continue;
            const bool onA = moves(groups[std::size_t(c.ga)]);
            const Group &g = groups[std::size_t(onA ? c.ga : c.gb)];
            fault(g.island, FaultKind::NonFinite, -1, onA ? c.in.a : c.in.b, machineOf(g));
        }
    }

    void World::Impl::rollBack() {
        for (Group &g : groups) {
            if (!moves(g) || !faulted(g)) continue;
            g.c = g.c0;
            g.q = g.q0;
            g.v = g.w = g.pv = g.pw = g.turned = Vec3{};
            g.correction = 0.0;
        }
        // Its machines' rows start cold; their articulations are read from the poses again next step.
        for (MachineData &md : machines) {
            if (md.frozen || md.groupOfLink.empty() || !faulted(groups[std::size_t(md.groupOfLink[0])])) continue;
            md.lastVel.clear();
            md.lastPos.clear();
            md.correction = 0.0;
        }
        for (const DriveSave &s : driveSaves) {
            JointSlot &j = joints[s.joint];
            const MachineData &md = machines[std::size_t(s.slot.info.machine)];
            if (!faulted(groups[std::size_t(md.groupOfLink[0])])) continue;
            j = s.slot;
            j.lambda = 0.0;
        }
        // Stored with no impulse, a contact starts cold.
        for (ContactData &c : contacts) {
            const Group &a = groups[std::size_t(c.ga)];
            if (!faulted(moves(a) ? a : groups[std::size_t(c.gb)])) continue;
            c.ln = c.lt1 = c.lt2 = c.lr1 = c.lr2 = c.ls = 0.0;
        }

        // Who to name: a group's first body, and how many bodies each island holds.
        bodyOfGroup.assign(groups.size(), -1);
        islandBodies.assign(islandFaulted.size(), 0);
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            BodySlot &b = bodies[i];
            if (!b.alive) continue;
            int &first = bodyOfGroup[std::size_t(b.group)];
            if (first < 0) first = int(i);
            const Group &g = groups[std::size_t(b.group)];
            if (!moves(g) || !faulted(g)) continue;
            ++b.faultStreak;
            ++islandBodies[std::size_t(g.island)];
        }
        for (std::size_t k = 0; k < faults.size(); ++k) {
            Fault &f = faults[k];
            const FaultSite &site = faultSites[k];
            if (site.group >= 0) {
                const int i = bodyOfGroup[std::size_t(site.group)];
                if (i >= 0) f.body = {std::uint32_t(i), bodies[std::size_t(i)].generation};
            }
            f.bodies = islandBodies[std::size_t(site.island)];
            if (const BodySlot *b = body(f.body)) f.streak = b->faultStreak;
        }
    }

} // namespace hf
