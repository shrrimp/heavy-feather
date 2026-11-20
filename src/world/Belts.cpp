#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include "world/Checks.hpp"
#include "world/Routes.hpp"
#include "world/WorldImpl.hpp"

// A belt's speed is a link of the machine its guides make (Machines.cpp), a wheel on its carrier: everything that ties
// it to a guide then sits in that machine's solve beside the guide's own drive (MachineRows.cpp).
namespace hf {

    namespace {

        constexpr double kPi = std::numbers::pi;
        constexpr Vec3 kX{1.0, 0.0, 0.0};
        constexpr double kParallel = 0.99984769515639; // cos 1°: guide axes further apart than this run on no one loop
        constexpr double kArcSamples = 16.0;           // per half turn: the area only sets the link's inertia

        bool usable(const BeltDesc &d) {
            return positive(d.length) && positive(d.mass) && nonNegative(d.pitch) && nonNegative(d.stiffness, true) &&
                   d.stiffness > 0.0 && nonNegative(d.damping) && nonNegative(d.strength, true) && d.strength > 0.0;
        }

        bool usable(const BeltGuide &g) {
            return positive(g.radius) && (g.side == 1 || g.side == -1) && nonNegative(g.friction);
        }

        Quat xTo(const Vec3 &a) {
            const Vec3 u = normalize(a);
            if (u.x < -1.0 + 1e-12) return Quat{0.0, 0.0, 1.0, 0.0};
            return normalize(Quat{1.0 + u.x, cross(kX, u)});
        }

        double enclosedArea(const std::vector<RopeBendState> &bends, const std::vector<RopePoint> &route,
                            const Vec3 &n) {
            const Vec3 e1 = normalize(std::fabs(n.x) < 0.9 ? cross(n, kX) : cross(n, Vec3{0.0, 1.0, 0.0}));
            const Vec3 e2 = cross(n, e1);
            std::vector<Vec3> ring;
            for (std::size_t k = 0; k < bends.size(); ++k) {
                const RopeBendState &b = bends[k];
                if (!b.touching) continue;
                const int steps = std::max(1, int(std::ceil(std::fabs(b.theta) / kPi * kArcSamples)));
                const Vec3 from = b.in - b.p;
                const double side = double(route[k].side);
                for (int s = 0; s <= steps; ++s)
                    ring.push_back(b.p + angleAxis(side * b.theta * double(s) / double(steps), b.n) * from);
            }
            double twice = 0.0;
            for (std::size_t k = 0; k < ring.size(); ++k) {
                const Vec3 &a = ring[k], &c = ring[(k + 1) % ring.size()];
                twice += dot(a, e1) * dot(c, e2) - dot(c, e1) * dot(a, e2);
            }
            return 0.5 * std::fabs(twice);
        }

    } // namespace

    BeltSlot *World::Impl::belt(BeltId id) {
        if (id.index >= belts.size()) return nullptr;
        BeltSlot &b = belts[id.index];
        return b.alive && b.generation == id.generation ? &b : nullptr;
    }

    const BeltSlot *World::Impl::belt(BeltId id) const { return const_cast<Impl *>(this)->belt(id); }

    int World::Impl::carrierGroup(const BeltSlot &b) const {
        const BodySlot *s = body(b.desc.carrier);
        return s ? s->group : 0;
    }

    void World::Impl::beltGroups(const BeltSlot &b, std::vector<int> &out) const {
        out.clear();
        out.push_back(carrierGroup(b));
        for (const RopePoint &p : b.route)
            if (const BodySlot *s = body(p.body)) out.push_back(s->group);
    }

    void World::Impl::touchBelt(const BeltSlot &b) {
        if (dirty) return;
        std::vector<int> mine;
        beltGroups(b, mine);
        for (int g : mine) {
            if (g < 0) {
                dirty = true;
                return;
            }
            touch(g);
            if (groups[std::size_t(g)].machine >= 0) touchedMachines.push_back(groups[std::size_t(g)].machine);
        }
    }

    void World::Impl::wakeBelt(const BeltSlot &b) {
        wakeUp(body(b.desc.carrier));
        for (const RopePoint &p : b.route) wakeUp(body(p.body));
    }

    void World::Impl::removeBeltsOn(JointId hinge) {
        for (std::uint32_t k = 0; k < belts.size(); ++k) {
            const BeltSlot &b = belts[k];
            if (!b.alive) continue;
            for (const BeltGuide &g : b.guides)
                if (g.hinge == hinge) {
                    dropBelt(k);
                    break;
                }
        }
    }

    void World::Impl::dropBelt(std::uint32_t index) {
        BeltSlot &b = belts[index];
        wakeBelt(b);
        touchBelt(b);
        b.alive = false;
        ++b.generation;
        freeBelts.push_back(index);
    }

    void World::Impl::beltLinkState(const BeltSlot &b, Vec3 &c, Quat &R, Vec3 &v, Vec3 &w) const {
        const End e = beltEnd(b);
        const Group &g = groups[std::size_t(e.group)];
        c = g.c + g.q * e.anchor;
        const Quat frame = g.q * e.frame;
        R = frame * angleAxis(std::remainder(b.position / b.radius, 2.0 * kPi), kX);
        v = g.v + cross(g.w, c - g.c);
        w = g.w + (frame * kX) * (b.speed / b.radius);
    }

    End World::Impl::beltEnd(const BeltSlot &b) const {
        End e;
        e.group = carrierGroup(b);
        e.isStatic = groups[std::size_t(e.group)].isStatic;
        if (const BodySlot *s = body(b.desc.carrier)) {
            e.anchor = s->offset + s->rot * (b.middle - s->desc.centerOfMass);
            e.frame = normalize(s->rot * xTo(b.axis));
        } else {
            e.anchor = b.middle;
            e.frame = xTo(b.axis);
        }
        return e;
    }

    void World::Impl::beginBelts() {
        brokenBelts.clear();
        for (BeltSlot &b : belts) {
            b.machine = b.link = -1;
            b.peak = 0.0;
        }
        for (std::size_t mi = 0; mi < machines.size(); ++mi) {
            const MachineData &md = machines[mi];
            if (md.frozen) continue;
            for (const MachineData::BeltLink &bl : md.belts) {
                BeltSlot &b = belts[std::size_t(bl.belt)];
                b.machine = int(mi);
                b.link = bl.link;
                b.angle0 = md.m.art.q[std::size_t(md.m.art.links[std::size_t(bl.link)].coord)];
            }
        }
        for (BeltSlot &b : belts) {
            if (!b.alive || b.machine < 0) continue;
            for (std::size_t k = 0; k < b.route.size(); ++k) {
                const RopePoint &rp = b.route[k];
                RopeBendState &bend = b.bends[k];
                const BodySlot &s = *body(rp.body);
                if (s.teleported) b.placed = false; // a wrap angle can't follow a jump
                bend.group = s.group;
                bend.local = s.offset + s.rot * (rp.position - s.desc.centerOfMass);
                bend.axis = s.rot * rp.axis;
            }
            if (b.phase.empty()) placeTeeth(b);
        }
    }

    // Teeth that come out whole keep their step on the links: where they sit the first step is where they stay.
    void World::Impl::placeTeeth(BeltSlot &b) {
        b.phase.assign(b.guides.size(), std::numeric_limits<double>::quiet_NaN());
        if (!(b.desc.pitch > 0.0)) return;
        const Quat carrier = groups[std::size_t(carrierGroup(b))].q;
        const Vec3 axis = normalize(carrier * beltEnd(b).frame * kX);
        for (std::size_t k = 0; k < b.guides.size(); ++k) {
            const BeltGuide &g = b.guides[k];
            const double teeth = 2.0 * kPi * g.radius / b.desc.pitch;
            if (!g.meshes || std::fabs(teeth - std::round(teeth)) > 1e-6 * teeth) continue;
            const Quat guide = groups[std::size_t(body(b.route[k].body)->group)].q;
            const double turn = turnAbout(conjugate(carrier) * guide, conjugate(carrier) * axis);
            b.phase[k] = wrapped(b.position - double(g.side) * g.radius * turn, b.desc.pitch);
        }
    }

    void World::Impl::layBelt(BeltSlot &b) {
        for (RopeBendState &bend : b.bends) {
            const Group &g = groups[std::size_t(bend.group)];
            bend.p = g.c + g.q * bend.local;
            bend.n = normalize(g.q * bend.axis);
        }
        b.span = laySpan(b.route, b.bends, !b.placed, true);
        b.placed = true;
    }

    void World::Impl::recordBelts(MachineData &md) {
        const std::vector<double> &lambda = md.m.lambda();
        const auto at = [&](int row) {
            return row >= 0 && std::size_t(row) < lambda.size() ? lambda[std::size_t(row)] : 0.0;
        };
        for (const MachineData::BeltLink &bl : md.belts) {
            BeltSlot &b = belts[std::size_t(bl.belt)];
            b.force = std::max(at(b.row) / h, b.stretched());
            b.peak = std::max(b.peak, b.force);
            b.slipping = 0;
            for (std::size_t k = 0; k < b.spinRows.size(); ++k)
                if (b.spinRows[k] >= 0 && std::fabs(at(b.spinRows[k])) >= kSlipping * b.spinCaps[k]) ++b.slipping;
        }
    }

    ContactData::Grip World::Impl::gripOf(const Contact &k, const Group &guide) const {
        ContactData::Grip g;
        const BeltSlot *b = belt(k.belt);
        if (!b || b->broken || b->machine < 0 || guide.machine != b->machine) return g;
        for (std::size_t i = 0; i < b->route.size(); ++i)
            if (b->route[i].body == k.b) {
                g.belt = int(k.belt.index);
                g.guide = int(i);
                g.axle = guide.q * b->bends[i].local;
                const double way = length(k.beltDirection);
                g.way = way > 0.0 ? k.beltDirection / way : Vec3{};
                return g;
            }
        return g;
    }

    Vec3 World::Impl::gripVelocity(const ContactData &c) const {
        const BeltSlot &b = belts[std::size_t(c.grip.belt)];
        const Group &guide = groups[std::size_t(c.gb)], &carrier = groups[std::size_t(carrierGroup(b))];
        return guide.v + cross(guide.w, c.grip.axle) + cross(carrier.w, c.rB - c.grip.axle) + b.speed * c.grip.way;
    }

    void World::Impl::recordGrips() {
        for (BeltSlot &b : belts)
            if (b.alive) b.grip.assign(b.guides.size(), 0.0);
        for (const ContactPrep &p : prep) {
            const ContactData &c = contacts[std::size_t(p.contact)];
            if (c.grip.belt >= 0) belts[std::size_t(c.grip.belt)].grip[std::size_t(c.grip.guide)] += c.ln / h;
        }
    }

    void World::Impl::measureBelts() {
        for (BeltSlot &b : belts)
            if (b.alive && b.machine >= 0) layBelt(b);
    }

    // The link's angle wraps when the machine is read again, so the belt's own position is what runs on.
    void World::Impl::endBelts() {
        for (std::uint32_t k = 0; k < belts.size(); ++k) {
            BeltSlot &b = belts[k];
            if (!b.alive || b.machine < 0) continue;
            const MachineData &md = machines[std::size_t(b.machine)];
            if (faulted(groups[std::size_t(md.groupOfLink[0])])) {
                b.force = b.peak = 0.0;
                continue;
            }
            if (!b.broken && b.peak > b.desc.strength) {
                b.broken = true;
                b.force = 0.0;
                brokenBelts.push_back({k, b.generation});
                wakeBelt(b);
            }
            const Link &L = md.m.art.links[std::size_t(b.link)];
            b.position += b.radius * (md.m.art.q[std::size_t(L.coord)] - b.angle0);
            b.speed = b.radius * md.m.art.qdot[std::size_t(L.dof0)];
        }
    }

    // The loop as given, laid where its guides are now: false for one no belt can run round.
    bool World::Impl::loopOf(const BeltDesc &desc, std::span<const BeltGuide> loop, BeltSlot &b) const {
        for (const BeltGuide &g : loop) {
            const JointSlot *j = joint(g.hinge);
            if (!j || !body(j->desc.b) || !usable(g) || !turns(j->desc.type)) return false;
        }
        if (loop.size() < 2) return false;

        std::vector<RopePoint> route;
        std::vector<RopeBendState> bends(loop.size());
        Vec3 axis{};
        for (std::size_t k = 0; k < loop.size(); ++k) {
            const JointSlot &j = *joint(loop[k].hinge);
            const BodySlot &guide = *body(j.desc.b);
            Vec3 c, v, w;
            Quat q;
            bodyState(guide, c, q, v, w);
            const Vec3 localAxis = normalize(j.desc.frameB.rotation * kX);
            const Vec3 worldAxis = q * localAxis;
            if (k == 0) axis = worldAxis;
            if (std::fabs(dot(axis, worldAxis)) < kParallel) return false;
            RopePoint p;
            p.body = j.desc.b;
            p.position = j.desc.frameB.position;
            p.axis = dot(axis, worldAxis) > 0.0 ? localAxis : -localAxis;
            p.radius = loop[k].radius;
            p.side = loop[k].side;
            route.push_back(p);
            bends[k].p = c + q * (p.position - guide.desc.centerOfMass);
        }

        // Guides listed clockwise run the belt clockwise: the axis it runs anticlockwise about is the one its sides
        // are taken about.
        double turning = 0.0;
        for (std::size_t k = 0; k < bends.size(); ++k)
            turning += dot(axis, cross(bends[k].p, bends[(k + 1) % bends.size()].p));
        if (turning < 0.0) {
            axis = -axis;
            for (RopePoint &p : route) p.axis = -p.axis;
        }
        for (RopeBendState &bend : bends) bend.n = axis;
        const double span = laySpan(route, bends, true, true);
        const double area = enclosedArea(bends, route, axis);
        if (!positive(span) || !positive(area)) return false;

        Vec3 middle{};
        for (const RopeBendState &bend : bends) middle += bend.p;
        middle = middle / double(bends.size());
        Vec3 c{}, v, w;
        Quat q{};
        if (const BodySlot *s = body(desc.carrier)) {
            bodyState(*s, c, q, v, w);
            c = c - q * s->desc.centerOfMass;
        }
        b.guides.assign(loop.begin(), loop.end());
        b.route = std::move(route);
        b.bends = std::move(bends);
        b.phase.clear();
        b.placed = true;
        b.span = span;
        // 2A/L, but never below its smallest guide's radius: a loop crossed on itself encloses next to nothing, and its
        // link would turn without bound.
        double smallest = loop[0].radius;
        for (const BeltGuide &g : loop) smallest = std::min(smallest, g.radius);
        b.radius = std::max(2.0 * area / span, smallest);
        b.middle = conjugate(q) * (middle - c);
        b.axis = conjugate(q) * axis;
        return true;
    }

    // A guide whose hinge or body is gone is a handle gone, not a value refused.
    bool World::Impl::guidesHere(std::span<const BeltGuide> loop) const {
        for (const BeltGuide &g : loop) {
            const JointSlot *j = joint(g.hinge);
            if (!j || !body(j->desc.b)) return false;
        }
        return true;
    }

    BeltId World::addBelt(const BeltDesc &desc, std::span<const BeltGuide> loop) {
        Impl &I = *impl_;
        if (desc.carrier.valid() && !I.body(desc.carrier)) return {};
        if (!I.guidesHere(loop)) return {};
        BeltSlot fresh;
        fresh.desc = desc;
        if (!usable(desc) || !I.loopOf(desc, loop, fresh)) return I.reject(BeltId{});

        std::uint32_t index;
        if (!I.freeBelts.empty()) {
            index = I.freeBelts.back();
            I.freeBelts.pop_back();
        } else {
            index = std::uint32_t(I.belts.size());
            I.belts.emplace_back();
        }
        BeltSlot &b = I.belts[index];
        fresh.generation = b.generation;
        fresh.alive = true;
        b = std::move(fresh);
        I.touchBelt(b);
        I.wakeBelt(b);
        return {index, b.generation};
    }

    // The belt keeps its place and speed: only the guides it runs round change.
    bool World::setBeltLoop(BeltId id, std::span<const BeltGuide> loop) {
        Impl &I = *impl_;
        BeltSlot *b = I.belt(id);
        if (!b || !I.guidesHere(loop)) return false;
        BeltSlot laid;
        if (!I.loopOf(b->desc, loop, laid)) return I.reject(false);
        // The machines of the guides it leaves are rebuilt as well as those it joins.
        I.wakeBelt(*b);
        I.touchBelt(*b);
        b->guides = std::move(laid.guides);
        b->route = std::move(laid.route);
        b->bends = std::move(laid.bends);
        b->phase.clear();
        b->placed = true;
        b->span = laid.span;
        b->radius = laid.radius;
        b->middle = laid.middle;
        b->axis = laid.axis;
        I.wakeBelt(*b);
        I.touchBelt(*b);
        return true;
    }

    void World::removeBelt(BeltId id) {
        if (impl_->belt(id)) impl_->dropBelt(id.index);
    }

    bool World::contains(BeltId id) const { return impl_->belt(id) != nullptr; }

    bool World::setBelt(BeltId id, const BeltDesc &desc) {
        Impl &I = *impl_;
        BeltSlot *b = I.belt(id);
        if (!b) return false;
        if (!usable(desc) || desc.carrier != b->desc.carrier) return I.reject(false);
        if (desc.pitch != b->desc.pitch) b->phase.clear();
        b->desc = desc;
        b->broken = false;
        I.wakeBelt(*b);
        I.touchBelt(*b);
        return true;
    }

    BeltDesc World::beltDesc(BeltId id) const {
        const BeltSlot *b = impl_->belt(id);
        return b ? b->desc : BeltDesc{};
    }

    double World::beltPosition(BeltId id) const {
        const BeltSlot *b = impl_->belt(id);
        return b ? b->position : 0.0;
    }

    bool World::setBeltMotion(BeltId id, double position, double speed) {
        BeltSlot *b = impl_->belt(id);
        if (!b) return false;
        if (!std::isfinite(position) || !std::isfinite(speed)) return impl_->reject(false);
        b->position = position;
        b->speed = speed;
        impl_->wakeBelt(*b);
        return true;
    }

    double World::beltSpeed(BeltId id) const {
        const BeltSlot *b = impl_->belt(id);
        return b ? b->speed : 0.0;
    }

    double World::beltTension(BeltId id) const {
        const BeltSlot *b = impl_->belt(id);
        return b ? b->force : 0.0;
    }

    BeltInfo World::beltInfo(BeltId id) const {
        BeltInfo info;
        if (const BeltSlot *b = impl_->belt(id)) {
            info.machine = b->machine;
            info.taut = b->force > 0.0;
            info.slipping = b->slipping;
            info.broken = b->broken;
        }
        return info;
    }

    std::span<const BeltId> World::brokenBelts() const { return impl_->brokenBelts; }

    int World::beltPath(BeltId id, std::span<RopeBend> out) const {
        const BeltSlot *b = impl_->belt(id);
        if (!b) return 0;
        const std::size_t n = std::min(out.size(), b->bends.size());
        for (std::size_t k = 0; k < n; ++k) out[k] = {b->bends[k].in, b->bends[k].out, b->bends[k].touching};
        return int(n);
    }

} // namespace hf
