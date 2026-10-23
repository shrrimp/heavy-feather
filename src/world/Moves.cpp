#include <cmath>
#include <vector>

#include "math/Spatial.hpp"
#include "world/Checks.hpp"
#include "world/WorldImpl.hpp"

// The host's moves: whatever puts a body somewhere the step didn't take it. A move is rigid wherever it can be, so a
// machine moved whole keeps its joints and its drives exactly as they were; and whatever is moved by the host starts
// its contacts and its machine's rows cold, since the pushes it was taking no longer apply where it is now.
namespace hf {

    namespace {

        void wakeUp(BodySlot &b) {
            b.asleep = false;
            b.stillSteps = 0;
        }

        bool identity(const Vec3 &fromC, const Quat &fromQ, const Vec3 &toC, const Quat &toQ) {
            return length(toC - fromC) <= 1e-9 && std::fabs(dot(fromQ, toQ)) >= 1.0 - 1e-12;
        }

    } // namespace

    void World::Impl::build() {
        if (!dirty) return;
        rebuildGroups();
        buildMachines();
        dirty = false;
    }

    void World::Impl::noteWritten(BodySlot &b) {
        if (b.listed) return;
        b.listed = true;
        written.push_back(std::uint32_t(&b - bodies.data()));
    }

    void World::Impl::markMachine(int machine, std::vector<std::uint8_t> &mark) const {
        mark.assign(groups.size(), 0);
        for (int g : machines[std::size_t(machine)].groupOfLink) mark[std::size_t(g)] = 1;
    }

    void World::Impl::moveGroups(const std::vector<std::uint8_t> &mark, const Vec3 &fromC, const Quat &fromQ,
                                 const Vec3 &toC, const Quat &toQ) {
        const Quat d = normalize(toQ * conjugate(fromQ));
        for (std::size_t gi = 1; gi < groups.size(); ++gi) {
            if (!mark[gi]) continue;
            Group &g = groups[gi];
            g.c = toC + d * (g.c - fromC);
            g.q = normalize(d * g.q);
            g.v = d * g.v;
            g.w = d * g.w;
            g.targetC = toC + d * (g.targetC - fromC);
            g.targetQ = normalize(d * g.targetQ);
        }
        for (BodySlot &b : bodies) {
            if (!b.alive || b.group < 0 || !mark[std::size_t(b.group)]) continue;
            bodyState(b, b.desc.position, b.desc.orientation, b.desc.velocity, b.desc.angularVelocity);
            b.poseWritten = b.velocityWritten = b.teleported = true;
            noteWritten(b);
            wakeUp(b);
        }
        // Whatever they touched wakes: it may have lost what it rested on.
        for (const Touch &t : touches) {
            BodySlot *a = body({std::uint32_t(t.a >> 32), std::uint32_t(t.a)});
            BodySlot *b = body({std::uint32_t(t.b >> 32), std::uint32_t(t.b)});
            const bool moved = (a && a->group >= 0 && mark[std::size_t(a->group)]) ||
                               (b && b->group >= 0 && mark[std::size_t(b->group)]);
            if (!moved) continue;
            if (a) a->wakeNext = true;
            if (b) b->wakeNext = true;
        }
    }

    bool World::Impl::followAnchors(MachineData &md) {
        Articulation &art = md.m.art;
        const auto held = [&](Anchor &a, Vec3 *&c, Quat *&q) {
            c = a.link >= 0 ? &art.links[std::size_t(a.link)].anchorPos : &a.c;
            q = a.link >= 0 ? &art.links[std::size_t(a.link)].anchorRot : &a.q;
        };
        for (Anchor &a : md.anchors) {
            Vec3 *c;
            Quat *q;
            held(a, c, q);
            const Group &g = groups[std::size_t(a.group)];
            if (g.c == *c && g.q == *q) continue;
            // The first anchor that moved carries the machine; any other closes on it through the position pass.
            markMachine(int(&md - machines.data()), groupMark);
            moveGroups(groupMark, *c, *q, g.c, g.q);
            for (Anchor &b : md.anchors) {
                held(b, c, q);
                *c = groups[std::size_t(b.group)].c;
                *q = groups[std::size_t(b.group)].q;
            }
            md.lastVel.clear();
            md.lastPos.clear();
            return true;
        }
        return false;
    }

    void World::Impl::carryAnchoredOn(int group) {
        for (MachineData &md : machines)
            for (const Anchor &a : md.anchors)
                if (a.group == group) {
                    followAnchors(md);
                    break;
                }
    }

    // Where each kinematic group has to be once the step is done, and the velocity that takes it there.
    void World::Impl::beginKinematic(double dt) {
        for (Group &g : groups) {
            if (!g.kinematic) continue;
            g.targetC = g.c;
            g.targetQ = g.q;
        }
        for (std::uint32_t i : written) {
            const BodySlot &b = bodies[i];
            if (!b.alive || !b.kinematic || b.group < 0) continue;
            Group &g = groups[std::size_t(b.group)];
            if (b.hasTarget) {
                g.targetQ = normalize(b.targetOrientation * conjugate(b.rot));
                g.targetC = b.targetPosition - g.targetQ * b.offset;
            } else if (b.velocityWritten) {
                const Vec3 w = b.desc.angularVelocity, v = b.desc.velocity - cross(w, g.q * b.offset);
                g.targetC = g.c + dt * v;
                g.targetQ = normalize(rotationVector(dt * w) * g.q);
            }
        }
        for (Group &g : groups) {
            if (!g.kinematic) continue;
            g.v = (g.targetC - g.c) / dt;
            g.w = rotationVectorOf(g.targetQ * conjugate(g.q)) / dt;
        }
    }

    // Along its way at time t into the step, from where the step began; on its target, exactly, at the end.
    void World::Impl::moveKinematic(Group &g, double t, bool last) const {
        if (last) {
            g.c = g.targetC;
            g.q = g.targetQ;
            return;
        }
        g.c = g.c0 + t * g.v;
        g.q = normalize(rotationVector(t * g.w) * g.q0);
    }

    void World::Impl::endStepBodies() {
        movedBodies.clear();
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            const BodySlot &b = bodies[i];
            if (!b.alive) continue;
            const Group &g = groups[std::size_t(b.group)];
            if ((moves(g) || g.kinematic) && (g.c != g.c0 || g.q != g.q0))
                movedBodies.push_back({std::uint32_t(i), b.generation});
        }
        for (std::uint32_t i : written) {
            BodySlot &b = bodies[i];
            b.poseWritten = b.velocityWritten = b.teleported = b.listed = b.hasTarget = false;
        }
        written.clear();
    }

    BodyState World::state(BodyId id) const {
        BodyState s;
        if (const BodySlot *b = impl_->body(id))
            impl_->bodyState(*b, s.position, s.orientation, s.velocity, s.angularVelocity);
        return s;
    }

    std::span<const BodyId> World::moved() const { return impl_->movedBodies; }

    bool World::moveMachine(BodyId id, const Vec3 &position, const Quat &orientation) {
        Impl &I = *impl_;
        BodySlot *b = I.body(id);
        if (!b) return false;
        if (!finite(position) || !rotation(orientation)) return I.reject(false);
        I.build();
        const Group &g = I.groups[std::size_t(b->group)];
        if (g.isStatic && !g.kinematic) return setPose(id, position, orientation);
        if (g.kinematic) return false;

        Vec3 c, v, w;
        Quat q;
        I.bodyState(*b, c, q, v, w);
        const Quat to = normalize(orientation);
        if (g.machine >= 0) {
            MachineData &md = I.machines[std::size_t(g.machine)];
            // Anchored, it can only go where its anchor went, and an anchor moved carries it there already.
            if (!md.anchors.empty() || md.anchored) return identity(c, q, position, to);
            I.markMachine(g.machine, I.groupMark);
        } else {
            I.groupMark.assign(I.groups.size(), 0);
            I.groupMark[std::size_t(b->group)] = 1;
        }
        I.moveGroups(I.groupMark, c, q, position, to);
        return true;
    }

    bool World::setMachineVelocity(BodyId id, const Vec3 &velocity, const Vec3 &angularVelocity) {
        Impl &I = *impl_;
        BodySlot *b = I.body(id);
        if (!b) return false;
        if (!finite(velocity) || !finite(angularVelocity)) return I.reject(false);
        I.build();
        const Group &g = I.groups[std::size_t(b->group)];
        if (g.isStatic) return false;
        Vec3 c, v, w;
        Quat q;
        I.bodyState(*b, c, q, v, w);
        if (g.machine >= 0) {
            const MachineData &md = I.machines[std::size_t(g.machine)];
            if (!md.anchors.empty() || md.anchored) return false;
            I.markMachine(g.machine, I.groupMark);
        } else {
            I.groupMark.assign(I.groups.size(), 0);
            I.groupMark[std::size_t(b->group)] = 1;
        }
        for (std::size_t gi = 1; gi < I.groups.size(); ++gi) {
            if (!I.groupMark[gi]) continue;
            Group &m = I.groups[gi];
            m.v = velocity + cross(angularVelocity, m.c - c);
            m.w = angularVelocity;
        }
        for (BodySlot &m : I.bodies) {
            if (!m.alive || m.group < 0 || !I.groupMark[std::size_t(m.group)]) continue;
            I.bodyState(m, c, q, m.desc.velocity, m.desc.angularVelocity);
            m.velocityWritten = true;
            I.noteWritten(m);
            wakeUp(m);
        }
        return true;
    }

    bool World::setKinematic(BodyId id, bool on) {
        Impl &I = *impl_;
        BodySlot *b = I.body(id);
        if (!b) return false;
        if (b->kinematic == on) return true;
        I.build();
        Group &g = I.groups[std::size_t(b->group)];
        if (g.machine >= 0) return false;

        // A weld group moves as one, so it is kinematic as one.
        bool isStatic = false, jointed = false;
        for (BodySlot &m : I.bodies) {
            if (!m.alive || m.group != b->group) continue;
            m.kinematic = on;
            isStatic = isStatic || m.desc.isStatic || !(m.desc.mass > 0.0);
            wakeUp(m);
        }
        I.wakePartners(id);
        g.kinematic = on;
        if (on) {
            g.isStatic = true;
            g.v = g.w = Vec3{};
            g.targetC = g.c;
            g.targetQ = g.q;
            return true;
        }
        g.isStatic = isStatic || g.invMass == 0.0;
        if (g.isStatic) {
            g.v = g.w = Vec3{};
        } else if (b->velocityWritten) {
            // Let go with a throw of the host's own.
            g.w = b->desc.angularVelocity;
            g.v = b->desc.velocity - cross(g.w, g.q * b->offset);
        }
        // Let go, a body jointed to anything becomes a part of a machine, or the end of a loose spring.
        for (const JointSlot &j : I.joints) {
            if (!j.alive || j.desc.type == JointType::Fixed) continue;
            const BodySlot *a = I.body(j.desc.a), *c = I.body(j.desc.b);
            jointed = jointed || (a && a->group == b->group) || (c && c->group == b->group);
        }
        if (jointed && !g.isStatic) I.dirty = true;
        return true;
    }

    bool World::kinematic(BodyId id) const {
        const BodySlot *b = impl_->body(id);
        return b && b->kinematic;
    }

} // namespace hf
