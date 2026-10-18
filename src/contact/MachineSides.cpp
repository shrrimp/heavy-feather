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

        // One end on a machine link: no solve of its own, the link's response carries it.
        void addLinkSide(World::Impl &I, ContactPrep &p, const ContactData &c, const Group &g, const Vec3 &r,
                         double sign, Mobility &k) {
            const Machine &m = I.machines[std::size_t(g.machine)].m;
            MachineSide &s = p.sides[p.ns++];
            s.machine = g.machine;
            s.response = responseFor(I, g.machine, g.link);
            s.n = m.dofs();
            wrenchesAt(p, c, g.c + r - m.art.O, sign, s.psi);
            const LinkResponse &R = I.responses[std::size_t(s.response)];
            for (int x = 0; x < rowsOf(p); ++x) {
                k.held[x] += quadratic(R.held, s.psi[x], s.psi[x]);
                k.unheld[x] += quadratic(R.unheld, s.psi[x], s.psi[x]);
            }
            k.held12 += quadratic(R.held, s.psi[1], s.psi[2]);
        }

        // Both ends on one machine: a push on either link moves the other, so the pair's mobility has a cross term.
        // Without it, two links a loop locks together would look free to each other.
        void addCross(World::Impl &I, const ContactPrep &p, Mobility &k) {
            const MachineSide &sa = p.sides[0], &sb = p.sides[1];
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

    } // namespace

    void addMachineSides(World::Impl &I, ContactPrep &p, const ContactData &c) {
        const Group &a = I.groups[std::size_t(p.a)], &b = I.groups[std::size_t(p.b)];
        if (a.machine >= 0 && a.machine == b.machine) {
            Mobility k;
            addLinkSide(I, p, c, a, p.r1, -1.0, k);
            addLinkSide(I, p, c, b, p.r2, 1.0, k);
            addCross(I, p, k);
            addMobility(p, k);
            return;
        }
        if (a.machine >= 0) {
            Mobility k;
            addLinkSide(I, p, c, a, p.r1, -1.0, k);
            addMobility(p, k);
        }
        if (b.machine >= 0) {
            Mobility k;
            addLinkSide(I, p, c, b, p.r2, 1.0, k);
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
    }

    void flushRuns(World::Impl &I, bool pseudo) {
        for (std::size_t mi = 0; mi < I.runs[pseudo].size(); ++mi) flushRun(I, int(mi), pseudo);
    }

} // namespace hf
