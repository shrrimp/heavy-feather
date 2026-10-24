#include <cmath>

#include "math/Spatial.hpp"
#include "world/Checks.hpp"
#include "world/WorldImpl.hpp"

// The host's hand on a body: its point pulled softly toward a target, never harder than the cap, as a mouse drags it.
// On a free body the rows are swept with the contacts, as a loose spring is; on a machine link they are bounded rows
// of the machine's own solve, as a capped drive is (MachineRows.cpp), so the whole linkage follows.
namespace hf {

    namespace {

        void wakeUp(BodySlot &b) {
            b.asleep = false;
            b.stillSteps = 0;
        }

        bool usable(const GrabDesc &d) {
            const auto cap = [](double x) { return x >= 0.0; }; // infinite: no cap
            return cap(d.maxForce) && cap(d.maxTorque) && std::isfinite(d.hertz) && d.hertz >= 0.0 &&
                   std::isfinite(d.damping) && d.damping >= 0.0;
        }

        constexpr Vec3 kAxes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

        // Three rows, held within a cap together: a pull of more than the cap in any direction is scaled back
        // onto it, so pulling diagonally is no stronger than pulling straight.
        Vec3 sweepRows(const double *k, const double *gamma, const double *target, double *lambda, const Vec3 &v,
                       double cap) {
            Vec3 next;
            for (int i = 0; i < 3; ++i) {
                const double den = k[i] + gamma[i];
                next[i] = den > 0.0 ? lambda[i] + (target[i] - v[i] - gamma[i] * lambda[i]) / den : lambda[i];
            }
            const double len = length(next);
            if (!(len <= cap)) next *= std::isfinite(len) && len > 0.0 ? cap / len : 0.0;
            const Vec3 d = next - Vec3{lambda[0], lambda[1], lambda[2]};
            for (int i = 0; i < 3; ++i) lambda[i] = next[i];
            return d;
        }

    } // namespace

    GrabId World::grab(BodyId id, const Vec3 &local, const GrabDesc &desc) {
        Impl &I = *impl_;
        BodySlot *b = I.body(id);
        if (!b) return {};
        if (!finite(local) || !usable(desc)) return I.reject(GrabId{});
        std::uint32_t index;
        if (!I.freeGrabs.empty()) {
            index = I.freeGrabs.back();
            I.freeGrabs.pop_back();
        } else {
            index = std::uint32_t(I.grabs.size());
            I.grabs.emplace_back();
        }
        GrabSlot &g = I.grabs[index];
        const std::uint32_t generation = g.generation;
        g = GrabSlot{};
        g.generation = generation;
        g.alive = true;
        g.desc = desc;
        g.body = id;
        g.local = local - b->desc.centerOfMass; // kept from the centre of mass, as the world works
        // Held where it is until the host says where to.
        Vec3 c, v, w;
        Quat q;
        I.bodyState(*b, c, q, v, w);
        g.target = c + q * g.local;
        g.targetRot = q;
        wakeUp(*b);
        return {index, generation};
    }

    bool World::setGrabTarget(GrabId id, const Vec3 &point, const Quat *orientation) {
        Impl &I = *impl_;
        if (id.index >= I.grabs.size() || !I.grabs[id.index].alive || I.grabs[id.index].generation != id.generation)
            return false;
        if (!finite(point) || (orientation && !rotation(*orientation))) return I.reject(false);
        GrabSlot &g = I.grabs[id.index];
        g.target = point;
        if (orientation) {
            g.targetRot = normalize(*orientation);
            g.turns = true;
        }
        if (BodySlot *b = I.body(g.body)) wakeUp(*b);
        return true;
    }

    void World::release(GrabId id) {
        Impl &I = *impl_;
        if (id.index >= I.grabs.size() || !I.grabs[id.index].alive || I.grabs[id.index].generation != id.generation)
            return;
        GrabSlot &g = I.grabs[id.index];
        if (BodySlot *b = I.body(g.body)) wakeUp(*b);
        g.alive = false;
        ++g.generation;
        I.freeGrabs.push_back(id.index);
    }

    // Each sub-step, the free bodies' grabs: soft rows along the world's axes toward the target, and about them toward
    // its orientation, their stiffness from the body's mass along each.
    void World::Impl::prepareGrabs() {
        grabPreps.clear();
        for (std::size_t gi = 0; gi < grabs.size(); ++gi) {
            const GrabSlot &s = grabs[gi];
            if (!s.alive) continue;
            const BodySlot *b = body(s.body);
            if (!b || b->group < 0) continue;
            const Group &g = groups[std::size_t(b->group)];
            if (!moves(g) || machineOf(g) >= 0) continue;
            GrabPrep p;
            p.grab = int(gi);
            p.group = b->group;
            p.island = g.island;
            p.r = g.q * (b->offset + b->rot * s.local);
            p.invIw = g.invInertiaWorld();
            const Vec3 error = g.c + p.r - s.target;
            for (int i = 0; i < 3; ++i) {
                const Vec3 rc = cross(p.r, kAxes[i]);
                p.k[i] = g.invMass + dot(rc, p.invIw * rc);
                springRow(p.k[i], s.desc.hertz, s.desc.damping, error[i], h, p.gamma[i], p.target[i]);
            }
            p.turns = s.turns && s.desc.maxTorque > 0.0;
            if (p.turns) {
                const Vec3 turn = rotationVectorOf(s.targetRot * conjugate(normalize(g.q * b->rot)));
                for (int i = 0; i < 3; ++i) {
                    p.k[3 + i] = dot(kAxes[i], p.invIw * kAxes[i]);
                    springRow(p.k[3 + i], s.desc.hertz, s.desc.damping, -turn[i], h, p.gamma[3 + i], p.target[3 + i]);
                }
            }
            grabPreps.push_back(p);
        }
    }

    void World::Impl::sweepGrabs(int island) {
        for (GrabPrep &p : grabPreps) {
            if (p.island != island) continue;
            Group &g = groups[std::size_t(p.group)];
            const GrabDesc &d = grabs[std::size_t(p.grab)].desc;
            const Vec3 P = sweepRows(p.k, p.gamma, p.target, p.lambda, g.v + cross(g.w, p.r), d.maxForce * h);
            g.v += g.invMass * P;
            g.w += p.invIw * cross(p.r, P);
            if (!p.turns) continue;
            const Vec3 L = sweepRows(p.k + 3, p.gamma + 3, p.target + 3, p.lambda + 3, g.w, d.maxTorque * h);
            g.w += p.invIw * L;
        }
    }

} // namespace hf
