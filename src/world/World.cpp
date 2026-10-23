#include "hf/World.hpp"

#include <algorithm>
#include <cmath>

#include "world/Checks.hpp"
#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        void wakeUp(BodySlot *b) {
            if (!b) return;
            b->asleep = false;
            b->stillSteps = 0;
        }

        void wakeEnds(World::Impl &I, const JointDesc &d) {
            wakeUp(I.body(d.a));
            wakeUp(I.body(d.b));
        }

        bool attachable(const World::Impl &I, const JointDesc &d) {
            const bool aOk = d.a.valid() ? I.body(d.a) != nullptr : d.b.valid();
            const bool bOk = d.b.valid() ? I.body(d.b) != nullptr : d.a.valid();
            return aOk && bOk && !(d.a == d.b);
        }

        // Symmetric to rounding, and positive definite: every leading minor above zero.
        bool inertia(const Mat3 &I) {
            double size = 0.0;
            for (int c = 0; c < 3; ++c)
                for (int r = 0; r < 3; ++r) {
                    if (!std::isfinite(I[c][r])) return false;
                    size = std::max(size, std::fabs(I[c][r]));
                }
            for (int c = 0; c < 3; ++c)
                for (int r = c + 1; r < 3; ++r)
                    if (std::fabs(I[c][r] - I[r][c]) > 1e-9 * size) return false;
            return I[0][0] > 0.0 && I[0][0] * I[1][1] - I[0][1] * I[1][0] > 0.0 && determinant(I) > 0.0;
        }

        // A dynamic body moves by its mass and inertia; a static one never uses them, so they only have to be numbers.
        bool massProperties(bool isStatic, double mass, const Mat3 &I) {
            if (isStatic) {
                for (int c = 0; c < 3; ++c)
                    if (!finite(I[c])) return false;
                return std::isfinite(mass);
            }
            return std::isfinite(mass) && mass > 0.0 && inertia(I);
        }

        bool usable(const BodyDesc &d) {
            return finite(d.position) && rotation(d.orientation) && finite(d.velocity) && finite(d.angularVelocity) &&
                   std::isfinite(d.radius) && massProperties(d.isStatic, d.mass, d.inertia);
        }

        // At or above zero, and a number; infinite only where infinite means no limit.
        bool nonNegative(double x, bool infinite = false) { return x >= 0.0 && (infinite || std::isfinite(x)); }

        bool usable(const JointDesc &d) {
            return finite(d.frameA.position) && rotation(d.frameA.rotation) && finite(d.frameB.position) &&
                   rotation(d.frameB.rotation) && std::isfinite(d.target) && nonNegative(d.maxForce, true) &&
                   nonNegative(d.forceSlew, true) && nonNegative(d.speed) && nonNegative(d.hertz) &&
                   nonNegative(d.damping) && nonNegative(d.restLength) && !std::isnan(d.min) && !std::isnan(d.max) &&
                   d.min <= d.max;
        }

        bool usable(const Contact &k) {
            const double n = length(k.normal);
            return finite(k.point) && std::isfinite(k.depth) && n >= 0.5 && n <= 2.0 && nonNegative(k.friction) &&
                   nonNegative(k.rolling) && nonNegative(k.spinning) && std::isfinite(k.restitution);
        }

        bool sameShape(const JointDesc &x, const JointDesc &y) {
            return x.type == y.type && x.a == y.a && x.b == y.b && x.frameA == y.frameA && x.frameB == y.frameB;
        }

        // A drive as a new joint has it. A piston's or servo's setpoint starts where it is, and travels to its target
        // at the drive's speed. One born a hair off starts on its target: twin cylinders authored a fraction of a
        // millimetre apart would otherwise ask the solve for two speeds on one freedom.
        void startDrive(World::Impl &I, JointSlot &j) {
            JointSlot fresh;
            fresh.desc = j.desc;
            fresh.generation = j.generation;
            fresh.alive = true;
            fresh.targetFrom = fresh.targetTo = fresh.targetEff = j.desc.target;
            j = fresh;
            if (positionDrive(j.desc.type)) {
                j.position = j.positionPrev = I.jointPosition(j);
                const bool there = std::fabs(j.position - j.desc.target) <= kServoDeadband;
                j.targetFrom = j.targetTo = j.targetEff = there ? j.desc.target : j.position;
            }
        }

    } // namespace

    World::World(const WorldSettings &settings) : impl_(std::make_unique<Impl>()) {
        impl_->settings = settings;
        impl_->groups.emplace_back();
    }

    World::~World() = default;
    World::World(World &&) noexcept = default;
    World &World::operator=(World &&) noexcept = default;

    const WorldSettings &World::settings() const { return impl_->settings; }
    void World::setSettings(const WorldSettings &settings) {
        impl_->settings = settings;
        for (BodySlot &b : impl_->bodies) wakeUp(&b);
    }

    BodySlot *World::Impl::body(BodyId id) {
        if (id.index >= bodies.size()) return nullptr;
        BodySlot &b = bodies[id.index];
        return b.alive && b.generation == id.generation ? &b : nullptr;
    }

    const BodySlot *World::Impl::body(BodyId id) const { return const_cast<Impl *>(this)->body(id); }

    JointSlot *World::Impl::joint(JointId id) {
        if (id.index >= joints.size()) return nullptr;
        JointSlot &j = joints[id.index];
        return j.alive && j.generation == id.generation ? &j : nullptr;
    }

    const JointSlot *World::Impl::joint(JointId id) const { return const_cast<Impl *>(this)->joint(id); }

    BodyId World::addBody(const BodyDesc &desc) {
        Impl &I = *impl_;
        if (!usable(desc)) return I.reject(BodyId{});
        std::uint32_t index;
        if (!I.freeBodies.empty()) {
            index = I.freeBodies.back();
            I.freeBodies.pop_back();
        } else {
            index = std::uint32_t(I.bodies.size());
            I.bodies.emplace_back();
        }
        // A slot reused starts as new: nothing of the body that had it is left, but its generation.
        BodySlot &b = I.bodies[index];
        BodySlot fresh;
        fresh.generation = b.generation;
        fresh.listed = b.listed;
        fresh.forced = b.forced;
        b = fresh;
        b.desc = desc;
        b.desc.orientation = normalize(desc.orientation);
        b.alive = true;
        I.dirty = true;
        return {index, b.generation};
    }

    void World::removeBody(BodyId id) {
        Impl &I = *impl_;
        BodySlot *b = I.body(id);
        if (!b) return;
        I.wakePartners(id);
        for (std::uint32_t k = 0; k < I.joints.size(); ++k) {
            const JointSlot &j = I.joints[k];
            if (j.alive && (j.desc.a == id || j.desc.b == id)) removeJoint({k, j.generation});
        }
        b->alive = false;
        ++b->generation;
        b->group = -1;
        I.freeBodies.push_back(id.index);
        I.dirty = true;
    }

    bool World::contains(BodyId id) const { return impl_->body(id) != nullptr; }

    bool World::asleep(BodyId id) const {
        const BodySlot *b = impl_->body(id);
        return b && b->asleep;
    }

    void World::wake(BodyId id) { wakeUp(impl_->body(id)); }

    void World::reshaped(BodyId id) {
        BodySlot *b = impl_->body(id);
        if (!b) return;
        wakeUp(b);
        impl_->wakePartners(id);
        b->reshaped = true;
        impl_->reshapedAny = true;
    }

    Vec3 World::position(BodyId id) const {
        const BodySlot *b = impl_->body(id);
        if (!b) return {};
        Vec3 c, v, w;
        Quat q;
        impl_->bodyState(*b, c, q, v, w);
        return c;
    }

    Quat World::orientation(BodyId id) const {
        const BodySlot *b = impl_->body(id);
        if (!b) return {};
        Vec3 c, v, w;
        Quat q;
        impl_->bodyState(*b, c, q, v, w);
        return q;
    }

    Vec3 World::velocity(BodyId id) const {
        const BodySlot *b = impl_->body(id);
        if (!b) return {};
        Vec3 c, v, w;
        Quat q;
        impl_->bodyState(*b, c, q, v, w);
        return v;
    }

    Vec3 World::angularVelocity(BodyId id) const {
        const BodySlot *b = impl_->body(id);
        if (!b) return {};
        Vec3 c, v, w;
        Quat q;
        impl_->bodyState(*b, c, q, v, w);
        return w;
    }

    bool World::setMassProperties(BodyId id, double mass, const Mat3 &inertia) {
        BodySlot *b = impl_->body(id);
        if (!b) return false;
        if (!massProperties(b->desc.isStatic, mass, inertia)) return impl_->reject(false);
        b->desc.mass = mass;
        b->desc.inertia = inertia;
        impl_->dirty = true;
        wakeUp(b);
        impl_->wakePartners(id);
        return true;
    }

    void World::setStatic(BodyId id, bool isStatic) {
        BodySlot *b = impl_->body(id);
        if (!b || b->desc.isStatic == isStatic) return;
        b->desc.isStatic = isStatic;
        impl_->dirty = true;
        wakeUp(b);
    }

    bool World::setPose(BodyId id, const Vec3 &position, const Quat &orientation) {
        Impl &I = *impl_;
        BodySlot *b = I.body(id);
        if (!b) return false;
        if (!finite(position) || !rotation(orientation)) return I.reject(false);
        const Quat qb = normalize(orientation);
        // A host writing back the pose it read isn't a push: only a real move wakes it.
        const bool moves =
            length(position - this->position(id)) > 1e-9 || std::fabs(dot(qb, this->orientation(id))) < 1.0 - 1e-12;
        if (moves) {
            wakeUp(b);
            I.wakePartners(id);
        }
        I.noteWritten(*b);
        if (b->kinematic) {
            b->targetPosition = position;
            b->targetOrientation = qb;
            b->hasTarget = true;
            return true;
        }
        b->desc.position = position;
        b->desc.orientation = qb;
        b->poseWritten = true;
        b->teleported = b->teleported || moves;
        if (b->group < 0) return true;
        Group &g = I.groups[std::size_t(b->group)];
        g.q = normalize(qb * conjugate(b->rot));
        g.c = position - g.q * b->offset;
        if (g.isStatic && moves) I.carryAnchoredOn(b->group);
        return true;
    }

    bool World::setVelocity(BodyId id, const Vec3 &velocity, const Vec3 &angularVelocity) {
        BodySlot *b = impl_->body(id);
        if (!b) return false;
        if (!finite(velocity) || !finite(angularVelocity)) return impl_->reject(false);
        if (length(velocity - this->velocity(id)) > 1e-9 || length(angularVelocity - this->angularVelocity(id)) > 1e-9)
            wakeUp(b);
        b->desc.velocity = velocity;
        b->desc.angularVelocity = angularVelocity;
        b->velocityWritten = true;
        impl_->noteWritten(*b);
        if (b->group < 0) return true;
        Group &g = impl_->groups[std::size_t(b->group)];
        // A sleeper keeps no velocity: a host writing back what it read leaves it still.
        if (g.isStatic || b->asleep) return true;
        g.w = angularVelocity;
        g.v = velocity - cross(angularVelocity, g.q * b->offset);
        return true;
    }

    bool World::applyForce(BodyId id, const Vec3 &force, const Vec3 &point) {
        BodySlot *b = impl_->body(id);
        if (!b) return false;
        if (!finite(force) || !finite(point)) return impl_->reject(false);
        if (force != Vec3{}) wakeUp(b);
        b->force += force;
        b->torque += cross(point - position(id), force);
        if (!b->forced) {
            b->forced = true;
            impl_->forced.push_back(id.index);
        }
        return true;
    }

    bool World::applyTorque(BodyId id, const Vec3 &torque) {
        BodySlot *b = impl_->body(id);
        if (!b) return false;
        if (!finite(torque)) return impl_->reject(false);
        if (torque != Vec3{}) wakeUp(b);
        b->torque += torque;
        if (!b->forced) {
            b->forced = true;
            impl_->forced.push_back(id.index);
        }
        return true;
    }

    JointId World::addJoint(const JointDesc &desc) {
        Impl &I = *impl_;
        if (!attachable(I, desc)) return {};
        if (!usable(desc)) return I.reject(JointId{});
        std::uint32_t index;
        if (!I.freeJoints.empty()) {
            index = I.freeJoints.back();
            I.freeJoints.pop_back();
        } else {
            index = std::uint32_t(I.joints.size());
            I.joints.emplace_back();
        }
        JointSlot &j = I.joints[index];
        j.desc = desc;
        startDrive(I, j);
        I.dirty = true;
        wakeEnds(I, desc);
        return {index, j.generation};
    }

    void World::removeJoint(JointId id) {
        Impl &I = *impl_;
        JointSlot *j = I.joint(id);
        if (!j) return;
        wakeEnds(I, j->desc);
        j->alive = false;
        ++j->generation;
        I.freeJoints.push_back(id.index);
        I.dirty = true;
    }

    bool World::contains(JointId id) const { return impl_->joint(id) != nullptr; }

    bool World::setJoint(JointId id, const JointDesc &desc) {
        Impl &I = *impl_;
        JointSlot *j = I.joint(id);
        if (!j || !attachable(I, desc)) return false;
        if (!usable(desc)) return I.reject(false);
        wakeEnds(I, j->desc);
        wakeEnds(I, desc);
        const bool reshape = !sameShape(j->desc, desc);
        const double target = j->desc.target;
        j->desc = desc;
        if (reshape) {
            startDrive(I, *j);
            I.dirty = true;
            return true;
        }
        if (positionDrive(desc.type)) {
            const auto inside = [&](double &x) { x = std::clamp(x, desc.min, desc.max); };
            for (double *x : {&j->desc.target, &j->targetFrom, &j->targetTo, &j->targetEff, &j->holdPos}) inside(*x);
        }
        if (desc.type == JointType::Motor && j->lever != 0.0) j->desc.target = j->lever * desc.speed;
        if (j->desc.target != target) j->commanded = true;
        return true;
    }

    JointDesc World::jointDesc(JointId id) const {
        const JointSlot *j = impl_->joint(id);
        return j ? j->desc : JointDesc{};
    }

    bool World::setTarget(JointId id, double target) {
        JointSlot *j = impl_->joint(id);
        if (!j) return false;
        if (!std::isfinite(target)) return impl_->reject(false);
        const double t = positionDrive(j->desc.type) ? std::clamp(target, j->desc.min, j->desc.max) : target;
        if (t != j->desc.target) wakeEnds(*impl_, j->desc);
        j->desc.target = t;
        j->commanded = true;
        return true;
    }

    bool World::setLever(JointId id, double lever) {
        JointSlot *j = impl_->joint(id);
        if (!j) return false;
        if (!std::isfinite(lever)) return impl_->reject(false);
        if (std::clamp(lever, -1.0, 1.0) != j->lever) wakeEnds(*impl_, j->desc);
        j->lever = std::clamp(lever, -1.0, 1.0);
        if (j->desc.type == JointType::Motor) j->desc.target = j->lever * j->desc.speed;
        return true;
    }

    double World::target(JointId id) const {
        const JointSlot *j = impl_->joint(id);
        return j ? j->desc.target : 0.0;
    }

    double World::jointPosition(JointId id) const {
        const JointSlot *j = impl_->joint(id);
        return j ? impl_->jointPosition(*j) : 0.0;
    }

    double World::Impl::jointPosition(const JointSlot &j) const {
        const auto frame = [&](BodyId id, const Frame &f, Vec3 &p, Quat &q) {
            const BodySlot *b = body(id);
            if (!b) {
                p = f.position;
                q = f.rotation;
                return;
            }
            Vec3 c, v, w;
            Quat qb;
            bodyState(*b, c, qb, v, w);
            p = c + qb * f.position;
            q = qb * f.rotation;
        };
        Vec3 pA, pB;
        Quat qA, qB;
        frame(j.desc.a, j.desc.frameA, pA, qA);
        frame(j.desc.b, j.desc.frameB, pB, qB);
        switch (j.desc.type) {
            case JointType::Hinge:
            case JointType::Motor:
            case JointType::AngularSpring:
            case JointType::Servo: {
                const Vec3 yA = qA * Vec3{0, 1, 0}, zA = qA * Vec3{0, 0, 1}, yB = qB * Vec3{0, 1, 0};
                return std::atan2(dot(yB, zA), dot(yB, yA));
            }
            case JointType::Slider:
            case JointType::Piston: return dot(pB - pA, qA * Vec3{1, 0, 0});
            case JointType::Spring: return length(pA - pB);
            case JointType::Fixed:
            case JointType::Ball: return 0.0;
        }
        return 0.0;
    }

    double World::driveForce(JointId id) const {
        const JointSlot *j = impl_->joint(id);
        return j && impl_->h > 0.0 ? j->lambda / impl_->h : 0.0;
    }

    JointInfo World::jointInfo(JointId id) const {
        const JointSlot *j = impl_->joint(id);
        return j ? j->info : JointInfo{};
    }

    void World::clearContacts() { impl_->contacts.clear(); }

    bool World::addContact(const Contact &contact) {
        if (!impl_->body(contact.a) || !impl_->body(contact.b) || contact.a == contact.b) return false;
        if (!usable(contact)) return impl_->reject(false);
        ContactData c;
        c.in = contact;
        c.in.normal = normalize(contact.normal);
        c.in.restitution = std::clamp(contact.restitution, 0.0, 1.0);
        impl_->contacts.push_back(c);
        return true;
    }

    void World::step(double dt) { impl_->step(dt); }

    WorldStats World::stats() const {
        WorldStats s = impl_->stats;
        s.bodies = int(impl_->bodies.size() - impl_->freeBodies.size());
        s.joints = int(impl_->joints.size() - impl_->freeJoints.size());
        return s;
    }

    std::span<const Fault> World::faults() const { return impl_->faults; }

} // namespace hf
