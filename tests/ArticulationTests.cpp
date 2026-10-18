#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "Random.hpp"
#include "dynamics/Articulation.hpp"
#include "math/Dense.hpp"

using namespace hf;

namespace {

    // Link 0 is a floating base or hinged to the world; the rest attach to random earlier links with random
    // joint kinds and frames.
    Articulation randomTree(Random &rnd, int nLinks, bool floating, bool allowSpherical = true) {
        Articulation a;
        a.links.resize(std::size_t(nLinks));
        for (int i = 0; i < nLinks; ++i) {
            Link &L = a.links[std::size_t(i)];
            L.mass = rnd.uni(0.5, 20.0);
            L.inertiaBody = rnd.inertia(L.mass);
            if (i == 0) {
                L.parent = -1;
                L.kind = floating ? JointKind::Free : JointKind::Revolute;
                L.anchorPos = rnd.vec(3.0);
                L.anchorRot = rnd.quat();
            } else {
                L.parent = int(rnd.uni(0, i - 0.001));
                const double k = rnd.uni(0, 1);
                L.kind = (k < 0.5)                      ? JointKind::Revolute
                         : (k < 0.8 || !allowSpherical) ? JointKind::Prismatic
                                                        : JointKind::Spherical;
            }
            L.pp = rnd.vec(1.0);
            L.pc = rnd.vec(1.0);
            L.Qp = rnd.quat();
            L.Qc = rnd.quat();
        }
        a.finalise();
        for (double &x : a.q) x = rnd.uni(-1.5, 1.5);
        for (Quat &r : a.qRot) r = rnd.quat();
        for (double &x : a.qdot) x = rnd.uni(-1.0, 1.0);
        if (floating) {
            a.links[0].c = rnd.vec(5.0);
            a.links[0].R = rnd.quat();
        }
        a.forwardKinematics();
        a.updateKinematics(a.links[0].c);
        return a;
    }

    double potential(const Articulation &a, const Vec3 &g) {
        double V = 0.0;
        for (const Link &L : a.links) V -= L.mass * dot(g, L.c);
        return V;
    }

    // One explicit sub-step of the unconstrained dynamics with joint torques tau.
    void stepFree(Articulation &a, double h, const Vec3 &g, const std::vector<double> &tau) {
        const int n = a.dofs();
        std::vector<double> M, C;
        a.updateKinematics(a.links[0].c);
        a.massMatrix(M);
        a.biasForces(g, {}, C);
        std::vector<double> rhs(a.qdot.size());
        for (std::size_t i = 0; i < rhs.size(); ++i) rhs[i] = (tau.empty() ? 0.0 : tau[i]) - C[i];
        choleskyFactor(M.data(), n, n);
        choleskySolve(M.data(), n, n, rhs.data());
        for (std::size_t i = 0; i < rhs.size(); ++i) a.qdot[i] += h * rhs[i];
        a.integrate(h, a.qdot);
    }

} // namespace

TEST_CASE("maximal <-> reduced round trip, point rows") {
    Random rnd;
    for (int t = 0; t < 30; ++t) {
        Articulation a = randomTree(rnd, 2 + int(rnd.uni(0, 6)), t % 2 == 0);
        const std::vector<double> q0 = a.q;
        const std::vector<Quat> r0 = a.qRot;
        std::vector<Vec3> c, v, w;
        std::vector<Quat> R;
        for (std::size_t i = 0; i < a.links.size(); ++i) {
            c.push_back(a.links[i].c);
            R.push_back(a.links[i].R);
            v.push_back(a.linkVelocity(int(i)));
            w.push_back(a.linkOmega(int(i)));
        }
        Articulation b = a;
        std::fill(b.q.begin(), b.q.end(), 0.0);
        std::fill(b.qdot.begin(), b.qdot.end(), 0.0);
        b.setFromMaximal(c, R, v, w);
        b.updateKinematics(b.links[0].c);

        double dq = 0.0;
        for (std::size_t i = 0; i < q0.size(); ++i) {
            double d = std::fabs(b.q[i] - q0[i]);
            d = std::min(d, std::fabs(d - 2.0 * std::numbers::pi)); // angles come back wrapped
            dq = std::max(dq, d);
        }
        CHECK_MESSAGE(dq < 1e-9, "coordinates come back");

        double dr = 0.0;
        for (std::size_t i = 0; i < r0.size(); ++i) dr = std::max(dr, 1.0 - std::fabs(dot(r0[i], b.qRot[i])));
        CHECK_MESSAGE(dr < 1e-12, "spherical rotations come back");

        double dv = 0.0;
        for (std::size_t i = 0; i < a.links.size(); ++i) {
            dv = std::max(dv, length(b.linkVelocity(int(i)) - v[i]));
            dv = std::max(dv, length(b.linkOmega(int(i)) - w[i]));
            dv = std::max(dv, length(b.links[i].c - c[i]));
        }
        CHECK_MESSAGE(dv < 1e-9, "poses and velocities come back");

        for (int k = 0; k < 5; ++k) {
            const int link = int(rnd.uni(0, double(a.links.size()) - 0.001));
            const Vec3 x = a.links[std::size_t(link)].c + rnd.vec(2.0);
            const Vec3 d = normalize(rnd.vec());
            std::vector<double> row(std::size_t(a.dofs()), 0.0);
            a.addPointRow(link, x, d, 1.0, row.data());
            double rq = 0.0;
            for (std::size_t i = 0; i < row.size(); ++i) rq += row[i] * a.qdot[i];
            const Vec3 vx = pointVelocity(a.twistOf(link, a.qdot.data()), x - a.O);
            CHECK_MESSAGE(std::fabs(rq - dot(d, vx)) < 1e-10, "a point row gives the point's velocity");
        }
    }
}

TEST_CASE("CRBA and RNEA against the Lagrangian") {
    const Vec3 g{0.0, -9.81, 0.0};
    Random rnd;
    for (int t = 0; t < 20; ++t) {
        // A fixed base and no spherical joints: every coordinate is a plain scalar to differentiate.
        Articulation a = randomTree(rnd, 2 + int(rnd.uni(0, 5)), false, false);
        const std::size_t n = std::size_t(a.dofs());
        std::vector<double> M, C;
        a.updateKinematics(a.links[0].c);
        a.massMatrix(M);
        a.biasForces(g, {}, C);

        double keM = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) keM += 0.5 * a.qdot[i] * M[i * n + j] * a.qdot[j];
        CHECK_MESSAGE(std::fabs(keM - a.kineticEnergy()) < 1e-9 * (1.0 + keM), "½ qdotᵀ M qdot = the links' energies");

        const auto at = [&](const std::vector<double> &q) {
            Articulation b = a;
            b.q = q;
            b.forwardKinematics();
            b.updateKinematics(b.links[0].c);
            return b;
        };
        const double eps = 1e-6;

        // d/dt (M qdot) along the motion
        std::vector<double> qp = a.q, qm = a.q;
        for (std::size_t i = 0; i < n; ++i) {
            qp[i] += eps * a.qdot[i];
            qm[i] -= eps * a.qdot[i];
        }
        std::vector<double> Mp, Mm;
        at(qp).massMatrix(Mp);
        at(qm).massMatrix(Mm);
        std::vector<double> Cfd(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) Cfd[i] += (Mp[i * n + j] - Mm[i * n + j]) / (2.0 * eps) * a.qdot[j];

        // - dT/dq + dV/dq
        for (std::size_t k = 0; k < n; ++k) {
            std::vector<double> q1 = a.q, q2 = a.q;
            q1[k] += eps;
            q2[k] -= eps;
            const Articulation b1 = at(q1), b2 = at(q2);
            Cfd[k] -= (b1.kineticEnergy() - b2.kineticEnergy()) / (2.0 * eps);
            Cfd[k] += (potential(b1, g) - potential(b2, g)) / (2.0 * eps);
        }
        double scale = 1.0;
        for (double c : C) scale = std::max(scale, std::fabs(c));
        CHECK_MESSAGE(maxAbsDiff(C, Cfd) < 1e-5 * scale, "bias forces = Lagrangian (Coriolis, centrifugal, gravity)");
    }
}

TEST_CASE("a pendulum keeps its energy") {
    const Vec3 g{0.0, -9.81, 0.0};
    Random rnd;
    Articulation a = randomTree(rnd, 4, false, false);
    for (double &x : a.qdot) x = 0.0;
    a.updateKinematics(a.links[0].c);
    const double e0 = a.kineticEnergy() + potential(a, g);

    double drift = 0.0;
    for (int s = 0; s < 50000; ++s) {
        stepFree(a, 2e-5, g, {});
        if (s % 1000 == 999) {
            a.updateKinematics(a.links[0].c);
            drift = std::max(drift, std::fabs(a.kineticEnergy() + potential(a, g) - e0));
        }
    }
    double scale = 0.0;
    for (const Link &L : a.links) scale += L.mass * 9.81;
    CHECK_MESSAGE(drift < 2e-3 * scale, "a 4-link pendulum keeps its energy over 1 s (symplectic Euler)");
}

TEST_CASE("a floating machine keeps its momentum") {
    const Vec3 g{};
    Random rnd;
    Articulation a = randomTree(rnd, 5, true);

    // About the world origin, not about O.
    const auto momentum = [](const Articulation &b) {
        Force h{};
        for (const Link &L : b.links) {
            const Force hl = L.I * L.V;
            h.f += hl.f;
            h.n += hl.n + cross(b.O, hl.f);
        }
        return h;
    };
    a.updateKinematics(a.links[0].c);
    const Force h0 = momentum(a);

    std::vector<double> tau(std::size_t(a.dofs()), 0.0);
    for (std::size_t i = 6; i < tau.size(); ++i) tau[i] = rnd.uni(-5.0, 5.0); // internal joint torques only
    Articulation half = a;
    for (int s = 0; s < 20000; ++s) stepFree(a, 5e-5, g, tau);
    for (int s = 0; s < 40000; ++s) stepFree(half, 2.5e-5, g, tau);
    half.updateKinematics(half.links[0].c);
    a.updateKinematics(a.links[0].c);
    const Force hh = momentum(half), h1 = momentum(a);
    const double scale = 1.0 + length(h0.f) + length(h0.n);

    // Symplectic Euler keeps momentum to first order in h: halving the step must about halve the drift. A wrong
    // bias term would not converge at all.
    CHECK_MESSAGE(length(hh.f - h0.f) < 0.65 * length(h1.f - h0.f) + 1e-12 * scale,
                  "linear momentum drift is first-order integration error");
    CHECK_MESSAGE(length(h1.f - h0.f) < 5e-2 * scale, "internal torques keep linear momentum, to O(h)");
    CHECK_MESSAGE(length(h1.n - h0.n) < 2e-3 * scale, "internal torques keep angular momentum");
}
