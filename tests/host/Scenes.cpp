#include "host/Scenes.hpp"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

using hf::BodyId;
using hf::JointId;
using hf::JointType;
using hf::Quat;
using hf::Vec3;

namespace host {

    namespace {

        constexpr double kPi = std::numbers::pi;
        const Vec3 kX{1, 0, 0}, kY{0, 1, 0}, kZ{0, 0, 1};

        // A fixed sequence, so every build of a scene is the same scene.
        struct Rng {
            std::uint64_t s;
            double next() {
                s = s * 6364136223846793005ull + 1442695040888963407ull;
                return double(s >> 11) * (1.0 / 9007199254740992.0);
            }
            double in(double lo, double hi) { return lo + (hi - lo) * next(); }
        };

        // A bar from p to q, its X along the bar.
        BodyId bar(Host &h, const Vec3 &p, const Vec3 &q, double thickness, double mass) {
            const Vec3 d = q - p;
            return h.addBox(0.5 * (p + q), xTo(d), Vec3{0.5 * length(d), thickness, thickness}, mass);
        }

        void ignoreAll(Host &h, const std::vector<BodyId> &parts) {
            for (std::size_t i = 0; i < parts.size(); ++i)
                for (std::size_t j = i + 1; j < parts.size(); ++j) h.ignore(parts[i], parts[j]);
        }

        // Four static walls round a square pit, and the bodies above it in layers: boxes and balls of mixed sizes.
        void addPile(Host &h, const Vec3 &at, int count, double size, Rng &rng) {
            const double half = 0.5 * size, wall = 0.2;
            for (const auto &[c, e] : {std::pair{Vec3{half + wall, 1, 0}, Vec3{wall, 1, half + 2 * wall}},
                                       {Vec3{-half - wall, 1, 0}, Vec3{wall, 1, half + 2 * wall}},
                                       {Vec3{0, 1, half + wall}, Vec3{half, 1, wall}},
                                       {Vec3{0, 1, -half - wall}, Vec3{half, 1, wall}}})
                h.addBox(at + c, Quat{}, e, 0.0);
            const int perRow = std::max(1, int(size / 0.9));
            for (int i = 0; i < count; ++i) {
                const int layer = i / (perRow * perRow), cell = i % (perRow * perRow);
                const Vec3 p = at +
                               Vec3{-half + 0.45 + 0.9 * (cell % perRow), 0.6 + 0.9 * layer,
                                    -half + 0.45 + 0.9 * (cell / perRow)} +
                               Vec3{rng.in(-0.05, 0.05), 0.0, rng.in(-0.05, 0.05)};
                if (i % 3 == 2) {
                    const double r = rng.in(0.18, 0.3);
                    h.addSphere(p, r, 900.0 * r * r * r * 4.2, rng.in(0.0, 0.4));
                } else {
                    const Vec3 e{rng.in(0.15, 0.35), rng.in(0.15, 0.35), rng.in(0.15, 0.35)};
                    const Quat q = hf::angleAxis(rng.in(-0.3, 0.3), normalize(Vec3{rng.in(-1, 1), 1, rng.in(-1, 1)}));
                    h.addBox(p, q, e, 600.0 * 8.0 * e.x * e.y * e.z, rng.in(0.3, 0.9));
                }
            }
        }

    } // namespace

    Scene pile(int bodies) {
        Scene s{"pile", std::make_unique<Host>(), nullptr};
        Rng rng{7};
        addPile(*s.host, Vec3{}, bodies, 6.0, rng);
        return s;
    }

} // namespace host

namespace host {

    namespace {

        // Links laid out along +x from a static post, hinged about z and y in turn with a ball every tenth joint,
        // falling onto a row of crates.
        std::vector<BodyId> addChain(Host &h, const Vec3 &at, int links, double height) {
            const double len = 0.3;
            const BodyId post = h.addBox(at + Vec3{-0.1, 0.5 * height, 0}, Quat{}, Vec3{0.1, 0.5 * height, 0.1}, 0.0);
            std::vector<BodyId> parts{post};
            BodyId prev = post;
            for (int i = 0; i < links; ++i) {
                const Vec3 p = at + Vec3{len * i, height, 0}, q = p + Vec3{len, 0, 0};
                const BodyId link = bar(h, p, q, 0.04, 2.0);
                const JointType type = i % 10 == 9 ? JointType::Ball : JointType::Hinge;
                h.addJoint(h.jointAt(type, prev, link, p, i % 2 ? kY : kZ));
                parts.push_back(link);
                prev = link;
            }
            for (double x = 1.0; x < len * links; x += 2.5)
                h.addBox(at + Vec3{x, 0.4, 0.1}, Quat{}, Vec3{0.4, 0.4, 0.4}, 30.0);
            return parts;
        }

    } // namespace

    Scene chain(int links) {
        Scene s{"chain", std::make_unique<Host>(), nullptr};
        addChain(*s.host, Vec3{}, links, 1.5);
        return s;
    }

    // Standing on edge like a loose track, it folds down onto the crates and balls around it.
    Scene ring(int links) {
        Scene s{"ring", std::make_unique<Host>(), nullptr};
        Host &h = *s.host;
        const double len = 0.3, radius = links * len / (2.0 * kPi);
        const Vec3 centre{0, radius + 0.3, 0};
        const auto at = [&](int i) {
            const double a = 2.0 * kPi * i / links;
            return centre + radius * Vec3{std::cos(a), std::sin(a), 0};
        };
        std::vector<BodyId> parts;
        for (int i = 0; i < links; ++i) parts.push_back(bar(h, at(i), at(i + 1), 0.05, 3.0));
        for (int i = 0; i < links; ++i)
            h.addJoint(
                h.jointAt(JointType::Hinge, parts[std::size_t(i)], parts[std::size_t((i + 1) % links)], at(i + 1), kZ));
        ignoreAll(h, parts);
        Rng rng{11};
        for (int i = 0; i < 12; ++i) {
            const Vec3 p{rng.in(-radius, radius), 0.4, rng.in(-1.5, 1.5)};
            if (i % 2) h.addSphere(p, 0.3, 40.0);
            else h.addBox(p, Quat{}, Vec3{0.35, 0.35, 0.35}, 40.0);
        }
        return s;
    }

    // Each link carries two children at its tip, splayed in planes that turn a quarter at every level.
    Scene branches(int depth) {
        Scene s{"branches", std::make_unique<Host>(), nullptr};
        Host &h = *s.host;
        const double top = 2.0 + 0.9 * depth;
        const BodyId root = h.addBox(Vec3{0, top + 0.2, 0}, Quat{}, Vec3{0.3, 0.2, 0.3}, 0.0);
        std::vector<BodyId> parts{root};
        struct Node {
            BodyId parent;
            Vec3 tip, dir;
            int level;
        };
        std::vector<Node> open{{root, Vec3{0, top, 0}, Vec3{0, -1, 0}, 0}};
        while (!open.empty()) {
            const Node n = open.back();
            open.pop_back();
            const double len = 0.8 * std::pow(0.85, n.level);
            const BodyId link = bar(h, n.tip, n.tip + len * n.dir, 0.04, 3.0 * std::pow(0.8, n.level));
            const Vec3 hingeAxis = n.level % 2 ? kX : kZ;
            const JointType type = n.level % 3 == 2 ? JointType::Ball : JointType::Hinge;
            h.addJoint(h.jointAt(type, n.parent, link, n.tip, hingeAxis));
            parts.push_back(link);
            if (n.level + 1 >= depth) continue;
            const Vec3 splay = n.level % 2 ? kZ : kX;
            for (double side : {-1.0, 1.0})
                open.push_back({link, n.tip + len * n.dir, normalize(n.dir + 0.7 * side * splay), n.level + 1});
        }
        ignoreAll(h, parts);
        for (int i = 0; i < 8; ++i) h.addSphere(Vec3{-2.0 + 0.6 * i, 0.3, 0.2 * (i % 3)}, 0.3, 20.0);
        return s;
    }

} // namespace host

namespace host {

    namespace {

        struct Lift {
            JointId piston;
            std::vector<BodyId> parts;
        };

        // Crossed bars pinned at their middles, each stage's tops pinned to the next stage's bottoms. One bottom end
        // is pinned to the ground, the other rides a roller the piston slides, and the platform sits on the top pair
        // the same way. Every stage closes a loop.
        Lift addScissorLift(Host &h, const Vec3 &at, int stages, int crates) {
            const double span = 2.0, rise = 0.6, y0 = 0.3;
            const auto y = [&](int k) { return y0 + rise * k; };
            Lift lift;
            const std::size_t n = std::size_t(stages);
            std::vector<BodyId> a(n), b(n);
            for (int k = 0; k < stages; ++k) {
                a[std::size_t(k)] = bar(h, at + Vec3{0, y(k), 0}, at + Vec3{span, y(k + 1), 0}, 0.06, 40.0);
                b[std::size_t(k)] = bar(h, at + Vec3{span, y(k), 0}, at + Vec3{0, y(k + 1), 0}, 0.06, 40.0);
                h.addJoint(h.jointAt(JointType::Hinge, a[std::size_t(k)], b[std::size_t(k)],
                                     at + Vec3{0.5 * span, 0.5 * (y(k) + y(k + 1)), 0}, kZ));
                lift.parts.insert(lift.parts.end(), {a[std::size_t(k)], b[std::size_t(k)]});
                if (k == 0) continue;
                h.addJoint(
                    h.jointAt(JointType::Hinge, b[std::size_t(k - 1)], a[std::size_t(k)], at + Vec3{0, y(k), 0}, kZ));
                h.addJoint(h.jointAt(JointType::Hinge, a[std::size_t(k - 1)], b[std::size_t(k)],
                                     at + Vec3{span, y(k), 0}, kZ));
            }
            h.addJoint(h.jointAt(JointType::Hinge, hf::BodyId{}, a[0], at + Vec3{0, y0, 0}, kZ));
            const BodyId roller = h.addBox(at + Vec3{span, y0, 0}, Quat{}, Vec3{0.1, 0.1, 0.1}, 20.0);
            hf::JointDesc slide = h.jointAt(JointType::Piston, hf::BodyId{}, roller, at + Vec3{span, y0, 0}, kX);
            slide.maxForce = 2e5;
            slide.speed = 0.3;
            slide.min = -0.6;
            slide.max = 0.2;
            lift.piston = h.addJoint(slide);
            h.addJoint(h.jointAt(JointType::Hinge, roller, b[0], at + Vec3{span, y0, 0}, kZ));

            const double top = y(stages);
            const BodyId deck =
                h.addBox(at + Vec3{0.5 * span, top + 0.12, 0}, Quat{}, Vec3{0.5 * span + 0.3, 0.06, 0.7}, 80.0);
            h.addJoint(h.jointAt(JointType::Hinge, b[std::size_t(stages - 1)], deck, at + Vec3{0, top, 0}, kZ));
            const BodyId runner = h.addBox(at + Vec3{span, top, 0}, Quat{}, Vec3{0.08, 0.05, 0.08}, 5.0);
            h.addJoint(h.jointAt(JointType::Slider, deck, runner, at + Vec3{span, top, 0}, kX));
            h.addJoint(h.jointAt(JointType::Hinge, runner, a[std::size_t(stages - 1)], at + Vec3{span, top, 0}, kZ));
            lift.parts.insert(lift.parts.end(), {roller, deck, runner});
            ignoreAll(h, lift.parts);

            for (int i = 0; i < crates; ++i)
                h.addBox(at + Vec3{0.3 + 0.7 * (i % 3), top + 0.45 + 0.62 * (i / 3), -0.2 + 0.2 * (i % 2)}, Quat{},
                         Vec3{0.3, 0.3, 0.3}, 25.0);
            return lift;
        }

        // Up and down its course, again and again.
        double cycle(int step, double period, double lo, double hi) {
            const double t = step / 60.0;
            return lo + (hi - lo) * (0.5 - 0.5 * std::cos(2.0 * kPi * t / period));
        }

    } // namespace

    Scene scissorLift(int stages) {
        Scene s{"scissor lift", std::make_unique<Host>(), nullptr};
        const Lift lift = addScissorLift(*s.host, Vec3{}, stages, 9);
        s.operate = [lift](Host &h, int i) { h.world.setTarget(lift.piston, cycle(i, 6.0, 0.0, -0.5)); };
        return s;
    }

} // namespace host

namespace host {

    namespace {

        Vec3 perp(const Vec3 &d) { return normalize(Vec3{-d.y, d.x, 0}); }

        // A hydraulic cylinder from a point on one body to a point on another: a barrel and a rod, hinged at either
        // end, and the piston between them.
        JointId addCylinder(Host &h, BodyId onA, const Vec3 &from, BodyId onB, const Vec3 &to, double maxForce,
                            std::vector<BodyId> &parts) {
            const Vec3 u = normalize(to - from), mid = 0.5 * (from + to);
            const double len = length(to - from);
            const BodyId barrel = bar(h, from, mid + 0.15 * len * u, 0.08, 40.0);
            const BodyId rod = bar(h, mid - 0.15 * len * u, to, 0.04, 20.0);
            h.addJoint(h.jointAt(JointType::Hinge, barrel, onA, from, kZ));
            h.addJoint(h.jointAt(JointType::Hinge, rod, onB, to, kZ));
            hf::JointDesc p = h.jointAt(JointType::Piston, barrel, rod, mid, u);
            p.maxForce = maxForce;
            p.speed = 0.5;
            p.min = -0.15 * len;
            p.max = 0.15 * len;
            parts.insert(parts.end(), {barrel, rod});
            return h.addJoint(p);
        }

    } // namespace

    Excavator addExcavator(Host &h, const Vec3 &at) {
        const BodyId base = h.addBox(at + Vec3{0, 0.75, 0}, Quat{}, Vec3{1.0, 0.75, 1.0}, 0.0);
        const Vec3 O1 = at + Vec3{0.5, 1.6, 0}, T1 = at + Vec3{3.5, 3.6, 0}, T2 = at + Vec3{4.7, 1.1, 0};
        const Vec3 stickTop = T1 + 0.6 * normalize(T1 - T2);
        const BodyId boom = bar(h, O1, T1, 0.15, 600.0);
        const BodyId stick = bar(h, stickTop, T2, 0.12, 300.0);
        const BodyId bucket = h.addBox(T2 + Vec3{0.35, -0.3, 0}, Quat{}, Vec3{0.45, 0.06, 0.5}, 150.0);
        std::vector<BodyId> parts{base, boom, stick, bucket};
        h.addJoint(h.jointAt(JointType::Hinge, boom, base, O1, kZ));
        h.addJoint(h.jointAt(JointType::Hinge, stick, boom, T1, kZ));
        h.addJoint(h.jointAt(JointType::Hinge, bucket, stick, T2, kZ));
        const Vec3 boomDir = normalize(T1 - O1), stickDir = normalize(T2 - stickTop);
        Excavator e;
        e.boom = addCylinder(h, base, at + Vec3{1.0, 1.3, 0}, boom, O1 + 1.5 * boomDir, 2e5, parts);
        e.stick = addCylinder(h, boom, O1 + 2.2 * boomDir + 0.35 * perp(boomDir), stick, stickTop, 1e5, parts);
        e.bucket = addCylinder(h, stick, stickTop + 1.2 * stickDir - 0.25 * perp(stickDir), bucket,
                               T2 + Vec3{0.75, 0.5, 0}, 5e4, parts);
        ignoreAll(h, parts);
        return e;
    }

    namespace {

        struct Crane {
            JointId slew, luff;
        };

        // A slewing mast, a luffing jib, and a load on a short chain from the jib's tip.
        Crane addCrane(Host &h, const Vec3 &at, int chainLinks) {
            const BodyId base = h.addBox(at + Vec3{0, 0.5, 0}, Quat{}, Vec3{1.0, 0.5, 1.0}, 0.0);
            const BodyId mast = bar(h, at + Vec3{0, 1.0, 0}, at + Vec3{0, 5.0, 0}, 0.2, 400.0);
            const BodyId jib = bar(h, at + Vec3{0, 5.0, 0}, at + Vec3{6, 5.0, 0}, 0.15, 300.0);
            hf::JointDesc slew = h.jointAt(JointType::Motor, mast, base, at + Vec3{0, 1.0, 0}, kY);
            slew.maxForce = 2e4;
            Crane c;
            c.slew = h.addJoint(slew);
            h.addJoint(h.jointAt(JointType::Hinge, jib, mast, at + Vec3{0, 5.0, 0}, kZ));
            std::vector<BodyId> parts{base, mast, jib};
            c.luff = addCylinder(h, mast, at + Vec3{0, 3.5, 0}, jib, at + Vec3{2.0, 5.0, 0}, 1e5, parts);
            ignoreAll(h, parts);
            BodyId prev = jib;
            Vec3 p = at + Vec3{6, 4.85, 0};
            for (int i = 0; i < chainLinks; ++i) {
                const Vec3 q = p - Vec3{0, 0.25, 0};
                const BodyId link = bar(h, p, q, 0.03, 4.0);
                h.addJoint(h.jointAt(JointType::Ball, prev, link, p, kZ));
                prev = link;
                p = q;
            }
            const BodyId load = h.addBox(p - Vec3{0, 0.4, 0}, Quat{}, Vec3{0.4, 0.4, 0.4}, 150.0);
            h.addJoint(h.jointAt(JointType::Ball, prev, load, p, kZ));
            return c;
        }

        struct Car {
            JointId rearLeft, rearRight, steer;
        };

        // Rear-wheel drive; the front wheels steer on knuckles, through tie rods from a rack the piston slides.
        Car addCar(Host &h, const Vec3 &at) {
            const BodyId chassis = h.addBox(at + Vec3{0, 0.8, 0}, Quat{}, Vec3{1.9, 0.3, 0.9}, 900.0);
            std::vector<BodyId> parts{chassis};
            Car car;
            for (double sz : {1.0, -1.0}) {
                const Vec3 w = at + Vec3{-1.3, 0.45, 0.95 * sz};
                const BodyId wheel = h.addSphere(w, 0.45, 25.0);
                h.shapes.back().rolling = 0.01;
                h.shapes.back().spinning = 0.02;
                hf::JointDesc m = h.jointAt(JointType::Motor, chassis, wheel, w, Vec3{0, 0, -1});
                m.maxForce = 900.0;
                (sz > 0 ? car.rearLeft : car.rearRight) = h.addJoint(m);
                parts.push_back(wheel);
            }
            const double armX = 1.0, rackZ = 0.18;
            const BodyId rack = bar(h, at + Vec3{armX, 0.45, -rackZ}, at + Vec3{armX, 0.45, rackZ}, 0.03, 6.0);
            hf::JointDesc slide = h.jointAt(JointType::Piston, chassis, rack, at + Vec3{armX, 0.45, 0}, kZ);
            slide.maxForce = 3e4;
            slide.speed = 0.4;
            slide.min = -0.11;
            slide.max = 0.11;
            car.steer = h.addJoint(slide);
            parts.push_back(rack);
            for (double sz : {1.0, -1.0}) {
                const Vec3 K = at + Vec3{1.3, 0.45, 0.55 * sz}, W = at + Vec3{1.3, 0.45, 0.95 * sz};
                const Vec3 A = at + Vec3{armX, 0.45, 0.5 * sz}, E = at + Vec3{armX, 0.45, rackZ * sz};
                const BodyId knuckle = bar(h, K, W, 0.04, 12.0);
                const BodyId tie = bar(h, E, A, 0.02, 2.0);
                const BodyId wheel = h.addSphere(W, 0.45, 25.0);
                h.shapes.back().rolling = 0.01;
                h.shapes.back().spinning = 0.02;
                h.addJoint(h.jointAt(JointType::Hinge, knuckle, chassis, K, kY));
                h.addJoint(h.jointAt(JointType::Hinge, wheel, knuckle, W, kZ));
                h.addJoint(h.jointAt(JointType::Hinge, tie, rack, E, kY));
                h.addJoint(h.jointAt(JointType::Hinge, tie, knuckle, A, kY));
                parts.insert(parts.end(), {knuckle, tie, wheel});
            }
            ignoreAll(h, parts);
            return car;
        }

    } // namespace

} // namespace host

namespace host {

    // A pit of loose bodies, an excavator digging a heap, a crane swinging its load into a column of crates, a car
    // circling through cones, a scissor lift raising a stack, a long chain over crates and a fallen ring.
    Scene yard() {
        Scene s{"yard", std::make_unique<Host>(), nullptr};
        Host &h = *s.host;
        Rng rng{3};
        addPile(h, Vec3{0, 0, 0}, 60, 6.0, rng);

        const Excavator dig = addExcavator(h, Vec3{-9, 0, 6});
        for (int i = 0; i < 24; ++i)
            h.addBox(Vec3{-4.6 + 0.65 * (i % 4), 0.3 + 0.6 * (i / 8), 5.2 + 0.6 * ((i / 4) % 2)},
                     hf::angleAxis(rng.in(-0.2, 0.2), kY), Vec3{0.28, 0.28, 0.28}, 30.0);

        // Its load swings over open ground and through a column of crates. Swung into the pit's wall and pile it
        // blows up: the joint-against-contact weakness the solver inherits, not this scene's subject.
        const Crane crane = addCrane(h, Vec3{8, 0, -10}, 12);
        for (int i = 0; i < 6; ++i) h.addBox(Vec3{8, 0.35 + 0.7 * i, -4}, Quat{}, Vec3{0.35, 0.35, 0.35}, 20.0);

        const Car car = addCar(h, Vec3{0, 0, 12});
        for (int i = 0; i < 10; ++i)
            h.addBox(Vec3{-2.0 + 2.2 * i, 0.3, i % 2 ? 9.5 : 14.5}, Quat{}, Vec3{0.2, 0.3, 0.2}, 4.0);

        const Lift lift = addScissorLift(h, Vec3{-9, 0, -7}, 8, 6);
        addChain(h, Vec3{-4, 0, -12}, 120, 1.2);

        s.operate = [dig, crane, car, lift](Host &host, int i) {
            hf::World &w = host.world;
            w.setTarget(dig.boom, cycle(i, 5.0, -0.1, 0.25));
            w.setTarget(dig.stick, cycle(i, 4.0, -0.3, 0.3));
            w.setTarget(dig.bucket, cycle(i, 3.0, -0.2, 0.3));
            w.setTarget(crane.slew, 0.4 * std::sin(i / 60.0 * 0.5)); // eased, as an operator would
            w.setTarget(crane.luff, cycle(i, 7.0, -0.2, 0.2));
            w.setTarget(car.rearLeft, 8.0);
            w.setTarget(car.rearRight, 8.0);
            w.setTarget(car.steer, 0.08);
            w.setTarget(lift.piston, cycle(i, 6.0, 0.0, -0.5));
        };
        return s;
    }

} // namespace host
