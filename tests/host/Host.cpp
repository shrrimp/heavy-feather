#include "host/Host.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

using hf::BodyId;
using hf::Quat;
using hf::Vec3;

namespace host {

    namespace {

        std::uint64_t pairKey(BodyId a, BodyId b) {
            const std::uint32_t lo = std::min(a.index, b.index), hi = std::max(a.index, b.index);
            return (std::uint64_t(lo) << 32) | hi;
        }

    } // namespace

    Quat xTo(const Vec3 &axis) {
        const Vec3 a = normalize(axis);
        if (a.x < -0.999999) return hf::angleAxis(std::numbers::pi, Vec3{0, 1, 0});
        return normalize(Quat{1.0 + a.x, cross(Vec3{1, 0, 0}, a)});
    }

    hf::Mat3 boxInertia(double m, const Vec3 &h) {
        return {m / 3.0 * (h.y * h.y + h.z * h.z), 0, 0, 0, m / 3.0 * (h.x * h.x + h.z * h.z), 0, 0, 0,
                m / 3.0 * (h.x * h.x + h.y * h.y)};
    }

    Host::Host(const hf::WorldSettings &settings) : world(settings) {
        hf::BodyDesc g;
        g.isStatic = true;
        ground_ = world.addBody(g);
        Shape plane;
        plane.body = ground_;
        plane.kind = ShapeKind::Plane;
        plane.friction = 0.8;
        shapes.push_back(plane);
    }

    BodyId Host::addBox(const Vec3 &position, const Quat &orientation, const Vec3 &half, double mass, double friction) {
        hf::BodyDesc d;
        d.position = position;
        d.orientation = orientation;
        d.mass = mass;
        d.inertia = boxInertia(mass, half);
        d.isStatic = mass <= 0.0;
        Shape s;
        s.body = world.addBody(d);
        s.half = half;
        s.friction = friction;
        shapes.push_back(s);
        return s.body;
    }

    BodyId Host::addSphere(const Vec3 &position, double radius, double mass, double restitution) {
        hf::BodyDesc d;
        d.position = position;
        d.mass = mass;
        d.inertia = hf::Mat3(0.4 * mass * radius * radius);
        Shape s;
        s.body = world.addBody(d);
        s.kind = ShapeKind::Sphere;
        s.radius = radius;
        s.restitution = restitution;
        shapes.push_back(s);
        return s.body;
    }

    void Host::removeBody(BodyId body) {
        world.removeBody(body);
        std::erase_if(shapes, [&](const Shape &s) { return s.body == body; });
    }

    hf::JointDesc Host::jointAt(hf::JointType type, BodyId a, BodyId b, const Vec3 &point, const Vec3 &axis) const {
        const Quat frame = xTo(axis);
        const auto local = [&](BodyId body) -> hf::Frame {
            if (!body.valid()) return {point, frame};
            const Quat q = world.orientation(body);
            return {conjugate(q) * (point - world.position(body)), normalize(conjugate(q) * frame)};
        };
        hf::JointDesc j;
        j.type = type;
        j.a = a;
        j.b = b;
        j.frameA = local(a);
        j.frameB = local(b);
        return j;
    }

    hf::JointId Host::addJoint(const hf::JointDesc &desc, bool collide) {
        if (!collide && desc.a.valid() && desc.b.valid()) ignore(desc.a, desc.b);
        return world.addJoint(desc);
    }

    void Host::ignore(BodyId a, BodyId b) { ignored_.insert(pairKey(a, b)); }

} // namespace host

namespace host {

    // Sort and sweep along x over bounding spheres; planes against everything.
    void Host::detect() {
        contacts.clear();
        poses_.resize(shapes.size());
        spans_.clear();
        planes_.clear();
        for (std::size_t i = 0; i < shapes.size(); ++i) {
            const Shape &s = shapes[i];
            const Quat q = world.orientation(s.body);
            poses_[i] = {world.position(s.body) + q * s.offset, toMat3(q)};
            if (s.kind == ShapeKind::Plane) {
                planes_.push_back(int(i));
                continue;
            }
            const double r = boundingRadius(s) + margin;
            spans_.push_back({poses_[i].c.x - r, poses_[i].c.x + r, int(i)});
        }
        std::sort(spans_.begin(), spans_.end(),
                  [](const Span &x, const Span &y) { return x.lo < y.lo || (x.lo == y.lo && x.shape < y.shape); });

        const auto test = [&](int i, int j) {
            const Shape &a = shapes[std::size_t(i)], &b = shapes[std::size_t(j)];
            if (a.body == b.body || ignored_.count(pairKey(a.body, b.body))) return;
            if (!world.contains(a.body) || !world.contains(b.body)) return;
            const std::uint64_t base = (std::uint64_t(i) << 40) | (std::uint64_t(j) << 16);
            collide(a, poses_[std::size_t(i)], b, poses_[std::size_t(j)], margin, base, contacts);
        };
        for (std::size_t k = 0; k < spans_.size(); ++k)
            for (std::size_t m = k + 1; m < spans_.size() && spans_[m].lo <= spans_[k].hi; ++m) {
                const int i = std::min(spans_[k].shape, spans_[m].shape),
                          j = std::max(spans_[k].shape, spans_[m].shape);
                const double reach =
                    boundingRadius(shapes[std::size_t(i)]) + boundingRadius(shapes[std::size_t(j)]) + margin;
                if (length(poses_[std::size_t(i)].c - poses_[std::size_t(j)].c) <= reach) test(i, j);
            }
        for (int p : planes_)
            for (const Span &s : spans_) test(std::min(p, s.shape), std::max(p, s.shape));
    }

    void Host::step(double dt) {
        detect();
        world.clearContacts();
        for (const hf::Contact &c : contacts) world.addContact(c);
        world.step(dt);
    }

} // namespace host
