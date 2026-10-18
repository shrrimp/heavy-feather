#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

#include <doctest/doctest.h>

#include "Random.hpp"
#include "solver/Machine.hpp"

using namespace hf;

namespace {

    // A planar four-bar, the shape of every hydraulic loop on the machines. Ground pivots at (0, 0) and
    // (2.4, 0); crank 1 m, coupler 2.6 m, rocker 2.2 m, all hinged about Z: a Grashof crank-rocker, so a motor on
    // the crank can turn it all the way round. The tree is crank and rocker on the world, coupler on the crank;
    // the coupler-rocker hinge is the cut.
    struct FourBar {
        Machine m;
        double crankLen = 1.0, couplerLen = 2.6, rockerLen = 2.2;
        Vec3 groundB{2.4, 0.0, 0.0};
    };

    constexpr double kPi = std::numbers::pi;
    const Vec3 kGravity{0.0, -9.81, 0.0};

    // A hinge about world Z needs a joint frame whose X is Z.
    const Quat kXtoZ = angleAxis(-kPi / 2.0, Vec3{0.0, 1.0, 0.0});

    Link hingeLink(int parent, double mass, double len, double rodInertia) {
        Link L;
        L.parent = parent;
        L.kind = JointKind::Revolute;
        L.Qp = kXtoZ;
        L.Qc = kXtoZ;
        L.pc = Vec3{-0.5 * len, 0, 0};
        L.mass = mass;
        L.inertiaBody = Mat3(0.1, 0, 0, 0, rodInertia, 0, 0, 0, rodInertia);
        return L;
    }

    Vec3 tipOf(const Articulation &a, int link, double len) {
        return a.links[std::size_t(link)].c + a.links[std::size_t(link)].R * Vec3{0.5 * len, 0, 0};
    }

    double closureError(FourBar &fb) {
        return length(tipOf(fb.m.art, 1, fb.couplerLen) - tipOf(fb.m.art, 2, fb.rockerLen));
    }

    FourBar makeFourBar() {
        FourBar fb;
        Articulation &a = fb.m.art;
        a.links = {
            hingeLink(-1, 2.0, fb.crankLen, 2.0 / 12.0 * 1.0),
            hingeLink(0, 3.0, fb.couplerLen, 3.0 / 12.0 * 6.76),
            hingeLink(-1, 1.5, fb.rockerLen, 1.5 / 12.0 * 4.84),
        };
        a.links[1].pp = Vec3{0.5 * fb.crankLen, 0, 0};
        a.links[2].anchorPos = fb.groundB;
        a.finalise();

        // Assemble: the crank at 80 degrees, then the loop solved for the other two angles (elbow up).
        const double th1 = 80.0 * kPi / 180.0;
        const Vec3 P{fb.crankLen * std::cos(th1), fb.crankLen * std::sin(th1), 0.0};
        const Vec3 B = fb.groundB;
        const double d = length(B - P);
        const double alpha = std::atan2(B.y - P.y, B.x - P.x);
        const double cb =
            (fb.couplerLen * fb.couplerLen + d * d - fb.rockerLen * fb.rockerLen) / (2 * fb.couplerLen * d);
        const double th2abs = alpha + std::acos(cb);
        const Vec3 Q = P + fb.couplerLen * Vec3{std::cos(th2abs), std::sin(th2abs), 0.0};
        const double th3 = std::atan2(Q.y - B.y, Q.x - B.x);

        a.q = {th1, th2abs - th1, th3};
        a.forwardKinematics();
        // The hinge frames may measure angles the other way round: then flip all three.
        if (closureError(fb) > 1e-6) {
            a.q = {-th1, -(th2abs - th1), -th3};
            a.forwardKinematics();
        }
        return fb;
    }

    // The cut hinge's rows: 3 linear (tips together) and 2 angular (axes parallel), a 3D hinge's five. A planar
    // loop only needs the in-plane linear pair, so two are redundant by construction.
    void addClosure(FourBar &fb, double h, bool position, bool planarOnly) {
        Articulation &a = fb.m.art;
        const Vec3 pc = tipOf(a, 1, fb.couplerLen), pr = tipOf(a, 2, fb.rockerLen);
        for (int k = 0; k < (planarOnly ? 2 : 3); ++k) {
            Vec3 e{};
            e[k] = 1.0;
            MachineRow r;
            r.target = position ? -0.8 * dot(pc - pr, e) / h : 0.0;
            double *J = position ? fb.m.addPosRow(r) : fb.m.addRow(r);
            a.addPointRow(1, pc, e, 1.0, J);
            a.addPointRow(2, pc, e, -1.0, J);
        }
        if (planarOnly) return;

        const Vec3 u = a.links[1].R * Vec3{0, 0, 1}; // the hinge axis is body Z here
        for (int k = 0; k < 2; ++k) {
            const Vec3 pv = a.links[2].R * (k == 0 ? Vec3{1, 0, 0} : Vec3{0, 1, 0});
            const Vec3 ax = cross(u, pv);
            MachineRow r;
            r.target = position ? -0.8 * dot(u, pv) / h : 0.0;
            double *J = position ? fb.m.addPosRow(r) : fb.m.addRow(r);
            a.addAngularRow(1, ax, 1.0, J);
            a.addAngularRow(2, ax, -1.0, J);
        }
    }

    // Steps the loop under gravity, returning the worst closure gap and the peak kinetic energy.
    void runFourBar(FourBar &fb, double seconds, double h, bool planarOnly, double &worstErr, double &worstKE) {
        std::vector<RowState> st, pst;
        worstErr = 0.0;
        worstKE = 0.0;
        for (int s = 0; s < int(seconds / h); ++s) {
            fb.m.dynamics(h, kGravity, {});
            addClosure(fb, h, false, planarOnly);
            fb.m.prepareRows();
            fb.m.solveRows(st);
            fb.m.buildMobility();
            addClosure(fb, h, true, planarOnly);
            fb.m.solvePositions(pst);
            fb.m.integrate(h);
            worstErr = std::max(worstErr, closureError(fb));
            fb.m.art.updateKinematics(fb.m.art.links[0].c);
            worstKE = std::max(worstKE, fb.m.art.kineticEnergy());
        }
    }

    // A four-bar whose rocker hinge is tilted out of plane by tiltDeg, the authored-misalignment case, driven at
    // 1 rad/s for 2 s. Closed with the cut's three point rows (a ball closure) as one group yielding by
    // groupYield (0: rigid rows). Returns how far the crank turned.
    double drivenTiltedFourBar(double tiltDeg, double groupYield, double &gapOut, bool &finite) {
        FourBar fb = makeFourBar();
        Link &rocker = fb.m.art.links[2];
        rocker.Qp = angleAxis(tiltDeg * kPi / 180.0, Vec3{1.0, 0.0, 0.0}) * rocker.Qp;
        fb.m.art.forwardKinematics();

        const double h = 1.0 / 960.0;
        const double crank0 = fb.m.art.q[0];
        const auto pointRows = [&](bool position) {
            Articulation &a = fb.m.art;
            const Vec3 pc = tipOf(a, 1, fb.couplerLen), pr = tipOf(a, 2, fb.rockerLen);
            for (int k = 0; k < 3; ++k) {
                Vec3 e{};
                e[k] = 1.0;
                MachineRow r;
                r.group = 7;
                r.groupYield = groupYield;
                r.target = position ? -0.8 * dot(pc - pr, e) / h : 0.0;
                double *J = position ? fb.m.addPosRow(r) : fb.m.addRow(r);
                a.addPointRow(1, pc, e, 1.0, J);
                a.addPointRow(2, pc, e, -1.0, J);
            }
        };

        std::vector<RowState> st, pst;
        finite = true;
        gapOut = 0.0;
        for (int s = 0; s < 1920; ++s) {
            fb.m.dynamics(h, kGravity, {});
            pointRows(false);
            MachineRow drive;
            drive.target = 1.0;
            drive.relEps = 1e-7;
            fb.m.addRow(drive)[0] = 1.0;
            fb.m.prepareRows();
            fb.m.solveRows(st);
            pointRows(true);
            fb.m.solvePositions(pst);
            fb.m.integrate(h);
            for (double v : fb.m.art.qdot) finite = finite && std::isfinite(v);
            gapOut = std::max(gapOut, closureError(fb));
        }
        return std::fabs(fb.m.art.q[0] - crank0);
    }

} // namespace

TEST_CASE("a four-bar swings with its loop held") {
    FourBar fb = makeFourBar();
    CHECK_MESSAGE(closureError(fb) < 1e-9, "the four-bar assembles closed");
    double err, ke;
    runFourBar(fb, 3.0, 1.0 / 960.0, true, err, ke);
    CHECK_MESSAGE(err < 1e-4, "3 s under gravity with the loop held (planar closure)");
    CHECK_MESSAGE((ke < 200.0 && std::isfinite(ke)), "and it does not gain energy");
}

TEST_CASE("a 3D hinge's five rows close a planar loop, two of them redundant") {
    FourBar fb = makeFourBar();
    double err, ke;
    runFourBar(fb, 3.0, 1.0 / 960.0, false, err, ke);
    CHECK_MESSAGE(err < 1e-4, "the loop is held");
    CHECK_MESSAGE((ke < 200.0 && std::isfinite(ke)), "and the redundancy does not pump it");
}

TEST_CASE("a capped motor drives the loop") {
    // Asked for 3 rad/s on the crank.
    FourBar fb = makeFourBar();
    const double h = 1.0 / 960.0, cap = 400.0;
    double lastSpeed = 0.0, lastLambda = 0.0, err = 0.0;
    std::vector<RowState> st, pst;
    for (int s = 0; s < 1920; ++s) {
        fb.m.dynamics(h, kGravity, {});
        addClosure(fb, h, false, false);
        MachineRow r;
        r.target = 3.0;
        r.lo = -cap * h;
        r.hi = cap * h;
        fb.m.addRow(r)[0] = 1.0;
        fb.m.prepareRows();
        fb.m.solveRows(st);
        lastLambda = fb.m.lambda().back();
        addClosure(fb, h, true, false);
        fb.m.solvePositions(pst);
        fb.m.integrate(h);
        lastSpeed = fb.m.art.qdot[0];
        err = closureError(fb);
    }
    CHECK_MESSAGE((std::fabs(std::fabs(lastSpeed) - 3.0) < 1e-6 || std::fabs(lastLambda) >= cap * h * (1 - 1e-9)),
                  "the motor holds 3 rad/s, or sits on its cap");
    CHECK_MESSAGE(std::fabs(lastLambda) <= cap * h * (1 + 1e-9), "the motor never exceeds its cap");
    CHECK_MESSAGE(err < 1e-4, "the loop is held while driven");
}

TEST_CASE("a weak motor saturates at its cap") {
    // Too weak to turn the loop against gravity: it must saturate, not cheat.
    FourBar fb = makeFourBar();
    const double h = 1.0 / 960.0, cap = 0.5;
    std::vector<RowState> st, pst;
    double worst = 0.0;
    for (int s = 0; s < 960; ++s) {
        fb.m.dynamics(h, kGravity, {});
        addClosure(fb, h, false, false);
        MachineRow r;
        r.target = 3.0;
        r.lo = -cap * h;
        r.hi = cap * h;
        fb.m.addRow(r)[0] = 1.0;
        fb.m.prepareRows();
        fb.m.solveRows(st);
        worst = std::max(worst, std::fabs(fb.m.lambda().back()));
        addClosure(fb, h, true, false);
        fb.m.solvePositions(pst);
        fb.m.integrate(h);
    }
    CHECK_MESSAGE(worst <= cap * h * (1 + 1e-9), "the motor stays on its cap");
}

TEST_CASE("twin inconsistent actuators") {
    // Two drives on one freedom (the loop ties the crank and rocker angles), commanded slightly inconsistently.
    // Rigid rows would be singular; with the per-row compliance the pair stays finite and capped.
    FourBar fb = makeFourBar();
    const double h = 1.0 / 960.0, cap = 2000.0;
    std::vector<RowState> st, pst;
    double worstLam = 0.0;
    bool finite = true;
    for (int s = 0; s < 960; ++s) {
        fb.m.dynamics(h, kGravity, {});
        addClosure(fb, h, false, false);
        MachineRow a;
        a.target = 1.0;
        a.lo = -cap * h;
        a.hi = cap * h;
        a.relEps = 1e-6;
        fb.m.addRow(a)[0] = 1.0;
        MachineRow b = a;
        b.target = 0.6;
        fb.m.addRow(b)[2] = 1.0; // the rocker, whose speed the loop already fixes
        fb.m.prepareRows();
        fb.m.solveRows(st);
        const std::size_t nr = std::size_t(fb.m.rowCount());
        worstLam = std::max({worstLam, std::fabs(fb.m.lambda()[nr - 1]), std::fabs(fb.m.lambda()[nr - 2])});
        addClosure(fb, h, true, false);
        fb.m.solvePositions(pst);
        fb.m.integrate(h);
        for (double v : fb.m.art.qdot) finite = finite && std::isfinite(v);
    }
    CHECK_MESSAGE(finite, "they stay finite");
    CHECK_MESSAGE(worstLam <= cap * h * (1 + 1e-9), "each stays on or under its cap");
    CHECK_MESSAGE(closureError(fb) < 1e-3, "and the loop stays closed");
}

TEST_CASE("authored misalignment: directional yield") {
    // A loop 1.5° out of plane has no mobility with rigid rows, and an unbounded drive forces it. This harness
    // has no position line search, so the release is the only thing that makes it solvable.
    double gapRigid = 0.0, gap = 0.0;
    bool finiteRigid = true, finite = true;
    drivenTiltedFourBar(1.5, 0.0, gapRigid, finiteRigid);
    const double turned = drivenTiltedFourBar(1.5, 1e-3, gap, finite);
    CHECK_MESSAGE(finite, "the yielding loop stays finite");
    CHECK_MESSAGE(turned > 1.5, "with the weak direction released, the drive turns it");
    // The 2.2 m rocker swings in a plane 1.5° off the loop's, so the cut opens by about 2.2 sin 1.5° ≈ 6 cm.
    CHECK_MESSAGE(gap < 0.1, "the cut opens by the authored error and no more");
    CHECK_MESSAGE((!finiteRigid || gapRigid > 10.0 * gap), "rigid rows do far worse");

    const double planar = drivenTiltedFourBar(0.0, 1e-3, gap, finite);
    CHECK_MESSAGE((planar > 1.9 && gap < 1e-3), "a planar loop with the same yield turns freely and stays closed");
}

TEST_CASE("a dead row changes nothing") {
    // An angular row about an in-plane axis, on a linkage hinged only about Z: the tree can't move along it at
    // all. Its nonsense target must be ignored.
    FourBar a = makeFourBar(), b = makeFourBar();
    a.m.art.qdot = b.m.art.qdot = {0.3, -0.2, 0.1};
    const double h = 1.0 / 960.0;
    std::vector<RowState> sa, sb;
    a.m.dynamics(h, kGravity, {});
    b.m.dynamics(h, kGravity, {});
    addClosure(a, h, false, true);
    addClosure(b, h, false, true);
    MachineRow dead;
    dead.target = 5.0;
    double *J = b.m.addRow(dead);
    b.m.art.addAngularRow(1, Vec3{1.0, 0.0, 0.0}, 1.0, J);
    b.m.art.addAngularRow(2, Vec3{1.0, 0.0, 0.0}, -1.0, J);
    a.m.prepareRows();
    b.m.prepareRows();
    a.m.solveRows(sa);
    b.m.solveRows(sb);
    CHECK_MESSAGE(maxAbsDiff(a.m.art.qdot, b.m.art.qdot) < 1e-9, "the solve is the same with or without it");
}

TEST_CASE("contact mobility W") {
    Random rnd;
    for (int t = 0; t < 10; ++t) {
        FourBar fb = makeFourBar();
        fb.m.art.qdot = {rnd.uni(-1, 1), rnd.uni(-1, 1), rnd.uni(-1, 1)};
        const double h = 1.0 / 960.0;
        fb.m.dynamics(h, kGravity, {});
        addClosure(fb, h, false, true);
        fb.m.prepareRows();
        std::vector<RowState> st;
        fb.m.solveRows(st);
        fb.m.buildMobility();

        const std::size_t n = std::size_t(fb.m.dofs());
        std::vector<double> g(n), w(n), g2(n), w2(n);
        for (double &x : g) x = rnd.uni(-1, 1);
        for (double &x : g2) x = rnd.uni(-1, 1);
        fb.m.mobility(g.data(), w.data());
        fb.m.mobility(g2.data(), w2.data());

        double worst = 0.0;
        for (int r = 0; r < fb.m.rowCount(); ++r) {
            double s = 0.0;
            for (std::size_t i = 0; i < n; ++i) s += fb.m.rowJ(r)[i] * w[i];
            worst = std::max(worst, std::fabs(s));
        }
        CHECK_MESSAGE(worst < 1e-6, "an impulse through W leaves the loop closed (J W g = 0)");

        double a = 0.0, b = 0.0, c = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            a += g2[i] * w[i];
            b += g[i] * w2[i];
            c += g[i] * w[i];
        }
        CHECK_MESSAGE(std::fabs(a - b) < 1e-9 * (1 + std::fabs(a)), "W is symmetric");
        CHECK_MESSAGE(c >= -1e-12, "W is positive semidefinite");
    }
}
