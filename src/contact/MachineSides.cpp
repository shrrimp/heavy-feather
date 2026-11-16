#include <algorithm>
#include <cstddef>

#include "math/Dense.hpp"
#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        // This sub-step's response of one machine link, made the first time a contact asks for it.
        int responseFor(World::Impl &I, int mi, int link) {
            MachineData &md = I.machines[std::size_t(mi)];
            if (md.response.size() != md.m.art.links.size()) md.response.assign(md.m.art.links.size(), -1);
            int &at = md.response[std::size_t(link)];
            if (at < 0) {
                LinkResponse R;
                R.machine = mi;
                R.link = link;
                R.off = I.responseArena.size();
                I.responseArena.resize(R.off + 6 * std::size_t(md.m.dofs()));
                md.m.linkResponse(link, I.responseArena.data() + R.off, R.held, R.unheld);
                at = int(I.responses.size());
                I.responses.push_back(R);
            }
            return at;
        }

        double quadratic(const double M[6][6], const Force &a, const Force &b) {
            double s = 0.0;
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c) s += component(a, r) * M[r][c] * component(b, c);
            return s;
        }

        Motion times(const double M[6][6], const Force &f) {
            double out[6];
            for (int r = 0; r < 6; ++r) {
                out[r] = 0.0;
                for (int c = 0; c < 6; ++c) out[r] += M[r][c] * component(f, c);
            }
            return {Vec3{out[0], out[1], out[2]}, Vec3{out[3], out[4], out[5]}};
        }

        int rowsOf(const ContactPrep &p) { return p.turns ? 6 : 3; }

        // The rows' directions as wrenches at x (about O): forces along n, t1, t2, then torques about t1, t2, n.
        void wrenchesAt(const ContactPrep &p, const ContactData &c, const Vec3 &x, double sign, Force *psi) {
            psi[0] = sign * pointForce(x, c.n);
            psi[1] = sign * pointForce(x, p.t1);
            psi[2] = sign * pointForce(x, p.t2);
            psi[3] = Force{sign * p.t1, Vec3{}};
            psi[4] = Force{sign * p.t2, Vec3{}};
            psi[5] = Force{sign * c.n, Vec3{}};
        }

        // The contact rows' mobility through its machine ends: with the loops held, and let go for the floor.
        struct Mobility {
            double held[6] = {}, unheld[6] = {}, held12 = 0.0;
        };

        // An end on a machine link: no solve of its own, the link's response carries it. Its wrenches are the
        // caller's to write.
        MachineSide &linkSide(World::Impl &I, ContactPrep &p, int machine, int link) {
            MachineSide &s = p.sides[p.ns++];
            s.machine = machine;
            s.response = responseFor(I, machine, link);
            s.n = I.machines[std::size_t(machine)].m.dofs();
            return s;
        }

        void addOwn(World::Impl &I, const ContactPrep &p, const MachineSide &s, Mobility &k) {
            const LinkResponse &R = I.responses[std::size_t(s.response)];
            for (int x = 0; x < rowsOf(p); ++x) {
                k.held[x] += quadratic(R.held, s.psi[x], s.psi[x]);
                k.unheld[x] += quadratic(R.unheld, s.psi[x], s.psi[x]);
            }
            k.held12 += quadratic(R.held, s.psi[1], s.psi[2]);
        }

        // Two ends on one machine: a push on either link moves the other, so the pair's mobility has a cross term.
        // Without it, two links a loop locks together would look free to each other.
        void addCross(World::Impl &I, const ContactPrep &p, const MachineSide &sa, const MachineSide &sb, Mobility &k) {
            const LinkResponse &Ra = I.responses[std::size_t(sa.response)];
            const LinkResponse &Rb = I.responses[std::size_t(sb.response)];
            const Articulation &art = I.machines[std::size_t(sa.machine)].m.art;
            double H[6][6]; // a's twist from a unit wrench on b
            for (int col = 0; col < 6; ++col) {
                const double *Z = I.responseArena.data() + Rb.off + std::size_t(col) * std::size_t(sb.n);
                const Motion t = art.twistOf(Ra.link, Z);
                for (int r = 0; r < 6; ++r) H[r][col] = component(t, r);
            }
            for (int x = 0; x < rowsOf(p); ++x) k.held[x] += 2.0 * quadratic(H, sa.psi[x], sb.psi[x]);
            k.held12 += quadratic(H, sa.psi[1], sb.psi[2]) + quadratic(H, sa.psi[2], sb.psi[1]);
        }

        void addMobility(ContactPrep &p, const Mobility &k) {
            const auto floored = [&](int x) { return std::max(k.held[x], kLockedMobility * k.unheld[x]); };
            p.kn += floored(0);
            p.kt1 += floored(1);
            p.kt2 += floored(2);
            p.kt12 += k.held12;
            if (p.turns) {
                p.kr1 += floored(3);
                p.kr2 += floored(4);
                p.ks += floored(5);
            }
        }

        // A belt's end pushes its guide at the axle, so the guide's suspension carries it and its spin doesn't count
        // the belt's speed twice. The carrier takes the push's moment about the axle, and the pull along the belt
        // turns the belt's link against the carrier. Three links of one machine: its rows are written over the
        // machine's whole velocity once, rather than swept a link at a time.
        void addGripRows(World::Impl &I, ContactPrep &p, const ContactData &c) {
            const BeltSlot &belt = I.belts[std::size_t(c.grip.belt)];
            const Group &guide = I.groups[std::size_t(p.b)], &carrier = I.groups[std::size_t(I.carrierGroup(belt))];
            Machine &m = I.machines[std::size_t(guide.machine)].m;
            const Articulation &art = m.art;
            const std::size_t n = std::size_t(m.dofs());
            const Vec3 axle = guide.c + c.grip.axle, lever = p.r2 - c.grip.axle;
            const Vec3 axis = art.links[std::size_t(belt.link)].R * Vec3{1.0, 0.0, 0.0};
            const bool onCarrier = carrier.machine == guide.machine && !carrier.frozen;
            const Vec3 along[3] = {c.n, p.t1, p.t2};
            p.denseMachine = guide.machine;
            p.denseOff = I.responseArena.size();
            I.responseArena.resize(p.denseOff + 6 * n, 0.0);
            double *g = I.responseArena.data() + p.denseOff;
            Mobility k;
            for (int x = 0; x < 3; ++x) {
                double *gx = g + std::size_t(x) * n, *wx = g + std::size_t(3 + x) * n;
                const double pull = dot(along[x], c.grip.way) * belt.radius;
                art.addPointRow(guide.link, axle, along[x], 1.0, gx);
                art.addAngularRow(belt.link, axis, pull, gx);
                if (onCarrier) {
                    const Vec3 turn = cross(lever, along[x]) - pull * axis;
                    const double size = length(turn);
                    if (size > 0.0) art.addAngularRow(carrier.link, turn / size, size, gx);
                }
                k.unheld[x] = m.mobility(gx, wx);
                k.held[x] = dot(gx, wx, n);
            }
            k.held12 = dot(g + n, g + 5 * n, n);
            addMobility(p, k);
        }

    } // namespace

    void addMachineSides(World::Impl &I, ContactPrep &p, const ContactData &c) {
        const Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
        if (machineOf(a) >= 0) {
            MachineSide &s = linkSide(I, p, a.machine, a.link);
            wrenchesAt(p, c, a.c + p.r1 - I.machines[std::size_t(a.machine)].m.art.O, -1.0, s.psi);
        }
        // A belt touching its own machine has no cross terms with its dense rows: there it is an ordinary contact.
        if (c.grip.belt >= 0 && machineOf(a) != b.machine) {
            addGripRows(I, p, c);
        } else if (machineOf(b) >= 0) {
            MachineSide &s = linkSide(I, p, b.machine, b.link);
            wrenchesAt(p, c, b.c + p.r2 - I.machines[std::size_t(b.machine)].m.art.O, 1.0, s.psi);
        }
        // Each machine's mobility is floored as one, its ends' cross terms in it.
        for (int i = 0; i < p.ns; ++i) {
            const int mi = p.sides[i].machine;
            bool seen = false;
            for (int j = 0; j < i; ++j) seen = seen || p.sides[j].machine == mi;
            if (seen) continue;
            Mobility k;
            for (int j = i; j < p.ns; ++j) {
                if (p.sides[j].machine != mi) continue;
                addOwn(I, p, p.sides[j], k);
                for (int l = i; l < j; ++l)
                    if (p.sides[l].machine == mi) addCross(I, p, p.sides[l], p.sides[j], k);
            }
            addMobility(p, k);
        }
    }

} // namespace hf

namespace hf {

    namespace {

        void flushRun(World::Impl &I, int mi, bool pseudo) {
            LinkRun &run = I.runs[pseudo][std::size_t(mi)];
            if (run.response < 0) return;
            const LinkResponse &R = I.responses[std::size_t(run.response)];
            I.machines[std::size_t(mi)].m.applyWrench(I.responseArena.data() + R.off, run.F, pseudo);
            run.response = -1;
        }

        // The run on this link: its twist read once, then kept current by the impulses swept into it.
        LinkRun &runOn(World::Impl &I, int mi, int response, bool pseudo) {
            LinkRun &run = I.runs[pseudo][std::size_t(mi)];
            if (run.response != response) {
                flushRun(I, mi, pseudo);
                const Machine &m = I.machines[std::size_t(mi)].m;
                const LinkResponse &R = I.responses[std::size_t(response)];
                run.response = response;
                run.V = m.art.twistOf(R.link, pseudo ? m.pseudo().data() : m.art.qdot.data());
                run.F = Force{};
            }
            return run;
        }

    } // namespace

    void machineVelocity(World::Impl &I, const ContactPrep &p, bool pseudo, int rows, double *v) {
        for (int k = 0; k < p.ns; ++k) {
            const MachineSide &s = p.sides[k];
            const Motion &V = runOn(I, s.machine, s.response, pseudo).V;
            for (int r = 0; r < rows; ++r) v[r] += dot(V, s.psi[r]);
        }
        if (p.denseMachine < 0) return;
        // A run on the machine has impulses not yet in its velocity: in first, and its twist read again after.
        flushRun(I, p.denseMachine, pseudo);
        const Machine &m = I.machines[std::size_t(p.denseMachine)].m;
        const std::size_t n = std::size_t(m.dofs());
        const double *x = pseudo ? m.pseudo().data() : m.art.qdot.data();
        for (int r = 0; r < std::min(rows, 3); ++r) v[r] += dot(I.responseArena.data() + p.denseOff + std::size_t(r) * n, x, n);
    }

    void machineImpulse(World::Impl &I, const ContactPrep &p, bool pseudo, int rows, const double *l) {
        for (int k = 0; k < p.ns; ++k) {
            const MachineSide &s = p.sides[k];
            Force f{};
            for (int r = 0; r < rows; ++r) f += l[r] * s.psi[r];
            LinkRun &run = runOn(I, s.machine, s.response, pseudo);
            run.F += f;
            run.V += times(I.responses[std::size_t(s.response)].held, f);
        }
        if (p.denseMachine < 0) return;
        flushRun(I, p.denseMachine, pseudo);
        Machine &m = I.machines[std::size_t(p.denseMachine)].m;
        const std::size_t n = std::size_t(m.dofs());
        for (int r = 0; r < std::min(rows, 3); ++r)
            if (l[r] != 0.0) m.applyContact(I.responseArena.data() + p.denseOff + std::size_t(3 + r) * n, l[r], pseudo);
    }

    void flushRuns(World::Impl &I, bool pseudo) {
        for (std::size_t mi = 0; mi < I.runs[pseudo].size(); ++mi) flushRun(I, int(mi), pseudo);
    }

    void flushRuns(World::Impl &I, bool pseudo, const int *machines, int count) {
        for (int k = 0; k < count; ++k) flushRun(I, machines[k], pseudo);
    }

} // namespace hf
