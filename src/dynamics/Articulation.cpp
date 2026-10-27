#include "dynamics/Articulation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "math/Dense.hpp"

namespace hf {

    namespace {

        constexpr Vec3 kX{1.0, 0.0, 0.0}, kY{0.0, 1.0, 0.0}, kZ{0.0, 0.0, 1.0};

        Mat3 worldInertia(const Quat &R, const Mat3 &Ibody) {
            const Mat3 M = toMat3(R);
            return M * Ibody * transpose(M);
        }

        int dofCount(JointKind k) {
            switch (k) {
                case JointKind::Free: return 6;
                case JointKind::Revolute: return 1;
                case JointKind::Prismatic: return 1;
                case JointKind::Spherical: return 3;
            }
            return 0;
        }

    } // namespace

    void Articulation::finalise() {
        ndof_ = 0;
        int nq = 0, nr = 0;
        for (Link &L : links) {
            L.ndof = dofCount(L.kind);
            L.dof0 = ndof_;
            ndof_ += L.ndof;
            if (L.kind == JointKind::Revolute || L.kind == JointKind::Prismatic) L.coord = nq++;
            else if (L.kind == JointKind::Spherical) L.coord = nr++;
            else L.coord = -1;
        }
        q.assign(std::size_t(nq), 0.0);
        qRot.assign(std::size_t(nr), Quat{});
        qdot.assign(std::size_t(ndof_), 0.0);
    }

    Motion Articulation::anchorTwist(const Link &L) const {
        if (L.parent >= 0 || L.kind == JointKind::Free) return Motion{};
        return Motion{L.anchorW, L.anchorV + cross(L.anchorW, O - L.anchorPos)};
    }

    void Articulation::parentPose(int i, Vec3 &c, Quat &R) const {
        const Link &L = links[std::size_t(i)];
        if (L.parent >= 0) {
            c = links[std::size_t(L.parent)].c;
            R = links[std::size_t(L.parent)].R;
        } else {
            c = L.anchorPos;
            R = L.anchorRot;
        }
    }

    Quat Articulation::parentFrame(int i) const {
        Vec3 c;
        Quat R;
        parentPose(i, c, R);
        return R * links[std::size_t(i)].Qp;
    }

    Vec3 Articulation::parentAnchor(int i) const {
        Vec3 c;
        Quat R;
        parentPose(i, c, R);
        return c + R * links[std::size_t(i)].pp;
    }

    void Articulation::setFromMaximal(const std::vector<Vec3> &c, const std::vector<Quat> &R,
                                      const std::vector<Vec3> &v, const std::vector<Vec3> &w) {
        const std::size_t n = links.size();
        for (std::size_t i = 0; i < n; ++i) {
            links[i].c = c[i];
            links[i].R = normalize(R[i]);
        }
        O = n ? c[0] : Vec3{};

        for (std::size_t i = 0; i < n; ++i) {
            Link &L = links[i];
            if (L.kind == JointKind::Free) {
                // O is this link's own centre, so the twist's linear part is its velocity.
                for (int k = 0; k < 3; ++k) {
                    qdot[std::size_t(L.dof0 + k)] = w[i][k];
                    qdot[std::size_t(L.dof0 + 3 + k)] = v[i][k] - cross(w[i], c[i] - O)[k];
                }
                continue;
            }
            Vec3 cP, vP, wP;
            Quat RP;
            if (L.parent >= 0) {
                const std::size_t p = std::size_t(L.parent);
                cP = links[p].c;
                RP = links[p].R;
                vP = v[p];
                wP = w[p];
            } else {
                cP = L.anchorPos;
                RP = L.anchorRot;
                vP = L.anchorV;
                wP = L.anchorW;
            }
            const Quat RPf = RP * L.Qp;
            const Quat Rcf = links[i].R * L.Qc;
            const Vec3 u = RPf * kX;
            switch (L.kind) {
                case JointKind::Revolute: {
                    const Vec3 yC = Rcf * kY, yP = RPf * kY, zP = RPf * kZ;
                    q[std::size_t(L.coord)] = std::atan2(dot(yC, zP), dot(yC, yP));
                    qdot[std::size_t(L.dof0)] = dot(u, w[i] - wP);
                    break;
                }
                case JointKind::Prismatic: {
                    const Vec3 xp = cP + RP * L.pp;
                    const Vec3 xc = links[i].c + links[i].R * L.pc;
                    q[std::size_t(L.coord)] = dot(xc - xp, u);
                    const Vec3 vc = v[i] + cross(w[i], xc - links[i].c);
                    const Vec3 vp = vP + cross(wP, xc - cP);
                    qdot[std::size_t(L.dof0)] = dot(u, vc - vp);
                    break;
                }
                case JointKind::Spherical: {
                    qRot[std::size_t(L.coord)] = normalize(conjugate(RPf) * Rcf);
                    const Vec3 rel = conjugate(RPf) * (w[i] - wP);
                    for (int k = 0; k < 3; ++k) qdot[std::size_t(L.dof0 + k)] = rel[k];
                    break;
                }
                case JointKind::Free: break;
            }
        }
        forwardKinematics();
    }

    void Articulation::forwardKinematics() {
        for (std::size_t i = 0; i < links.size(); ++i) {
            Link &L = links[i];
            if (L.kind == JointKind::Free) continue;
            Vec3 cP;
            Quat RP;
            parentPose(int(i), cP, RP);
            const Quat RPf = RP * L.Qp;
            Vec3 x = cP + RP * L.pp;
            Quat Rcf = RPf;
            switch (L.kind) {
                case JointKind::Revolute: Rcf = RPf * angleAxis(q[std::size_t(L.coord)], kX); break;
                case JointKind::Prismatic: x += (RPf * kX) * q[std::size_t(L.coord)]; break;
                case JointKind::Spherical: Rcf = RPf * qRot[std::size_t(L.coord)]; break;
                case JointKind::Free: break;
            }
            L.R = normalize(Rcf * conjugate(L.Qc));
            L.c = x - L.R * L.pc;
        }
    }

    void Articulation::updateKinematics(const Vec3 &newO) {
        if (!links.empty() && links[0].kind == JointKind::Free) {
            // The base twist's linear part is the velocity of the point at O: moving O moves it.
            const Vec3 w{qdot[0], qdot[1], qdot[2]};
            const Vec3 dv = cross(w, newO - O);
            for (int k = 0; k < 3; ++k) qdot[std::size_t(3 + k)] += dv[k];
        }
        O = newO;

        for (std::size_t i = 0; i < links.size(); ++i) {
            Link &L = links[i];
            switch (L.kind) {
                case JointKind::Free:
                    for (int k = 0; k < 3; ++k) {
                        L.S[k] = Motion{};
                        L.S[k].w[k] = 1.0;
                        L.S[3 + k] = Motion{};
                        L.S[3 + k].v[k] = 1.0;
                    }
                    break;
                case JointKind::Revolute: {
                    const Vec3 u = parentFrame(int(i)) * kX;
                    const Vec3 x = parentAnchor(int(i)) - O;
                    L.S[0] = Motion{u, cross(x, u)};
                    break;
                }
                case JointKind::Prismatic: L.S[0] = Motion{Vec3{}, parentFrame(int(i)) * kX}; break;
                case JointKind::Spherical: {
                    const Quat RPf = parentFrame(int(i));
                    const Vec3 x = parentAnchor(int(i)) - O;
                    for (int k = 0; k < 3; ++k) {
                        Vec3 e{};
                        e[k] = 1.0;
                        e = RPf * e;
                        L.S[k] = Motion{e, cross(x, e)};
                    }
                    break;
                }
            }
            L.I = Inertia::body(L.mass, L.c - O, worldInertia(L.R, L.inertiaBody));
        }
        updateVelocities(qdot.data());

        for (Link &L : links) L.Ic = L.I;
        for (std::size_t i = links.size(); i-- > 0;)
            if (links[i].parent >= 0) links[std::size_t(links[i].parent)].Ic += links[i].Ic;
    }

    void Articulation::updateVelocities(const double *qd) {
        for (Link &L : links) {
            // An anchor moving at a constant twist has no spatial acceleration: its spin's effects come from riding on
            // its twist, in the velocity products below.
            const Motion Vp = L.parent >= 0 ? links[std::size_t(L.parent)].V : anchorTwist(L);
            const Motion Ap = L.parent >= 0 ? links[std::size_t(L.parent)].A : Motion{};
            Motion vj{};
            for (int k = 0; k < L.ndof; ++k) vj += qd[L.dof0 + k] * L.S[k];
            L.V = Vp + vj;
            // S rides on the parent, so dS/dt = Vp × S and the joint adds Vp × (S qdot).
            L.A = Ap + crossM(Vp, vj);
        }
    }

    void Articulation::massMatrix(std::vector<double> &M) const {
        const std::size_t n = std::size_t(ndof_);
        M.assign(n * n, 0.0);
        for (const Link &L : links) {
            for (int a = 0; a < L.ndof; ++a) {
                const Force F = L.Ic * L.S[a];
                const std::size_t ra = std::size_t(L.dof0 + a);
                for (int b = 0; b < L.ndof; ++b) M[ra * n + std::size_t(L.dof0 + b)] = dot(L.S[b], F);
                for (int j = L.parent; j >= 0; j = links[std::size_t(j)].parent) {
                    const Link &P = links[std::size_t(j)];
                    for (int b = 0; b < P.ndof; ++b) {
                        const double val = dot(P.S[b], F);
                        const std::size_t cb = std::size_t(P.dof0 + b);
                        M[ra * n + cb] = val;
                        M[cb * n + ra] = val;
                    }
                }
            }
        }
    }

    bool Articulation::anchorJumps(std::vector<double> &tau) {
        acc_.resize(links.size());
        bool any = false;
        for (std::size_t i = 0; i < links.size(); ++i) {
            Link &L = links[i];
            if (L.parent >= 0) {
                acc_[i] = acc_[std::size_t(L.parent)];
                continue;
            }
            acc_[i] = Motion{};
            if (L.kind == JointKind::Free || (L.anchorJumpV == Vec3{} && L.anchorJumpW == Vec3{})) continue;
            acc_[i] = Motion{L.anchorJumpW, L.anchorJumpV + cross(L.anchorJumpW, O - L.anchorPos)};
            L.anchorJumpV = L.anchorJumpW = Vec3{};
            any = true;
        }
        if (!any) return false;
        tau.assign(std::size_t(ndof_), 0.0);
        for (std::size_t i = 0; i < links.size(); ++i) links[i].F = links[i].I * acc_[i];
        for (std::size_t i = links.size(); i-- > 0;) {
            Link &L = links[i];
            for (int a = 0; a < L.ndof; ++a) tau[std::size_t(L.dof0 + a)] = dot(L.S[a], L.F);
            if (L.parent >= 0) links[std::size_t(L.parent)].F += L.F;
        }
        return true;
    }

    void Articulation::biasForces(const Vec3 &g, const std::vector<Force> &fext, std::vector<double> &C) {
        C.assign(std::size_t(ndof_), 0.0);
        const Motion ag{Vec3{}, g};
        for (std::size_t i = 0; i < links.size(); ++i) {
            Link &L = links[i];
            L.F = L.I * (L.A - L.gravityScale * ag) + crossF(L.V, L.I * L.V) + L.damping * (L.I * L.V);
            if (!fext.empty()) L.F = L.F - fext[i];
        }
        for (std::size_t i = links.size(); i-- > 0;) {
            Link &L = links[i];
            for (int a = 0; a < L.ndof; ++a) C[std::size_t(L.dof0 + a)] = dot(L.S[a], L.F);
            if (L.parent >= 0) links[std::size_t(L.parent)].F += L.F;
        }
    }

    void Articulation::massDiagonal(std::vector<double> &d) const {
        d.assign(std::size_t(ndof_), 0.0);
        for (const Link &L : links)
            for (int a = 0; a < L.ndof; ++a) d[std::size_t(L.dof0 + a)] = dot(L.S[a], L.Ic * L.S[a]);
    }

    // Leaves first: each link's articulated inertia is its own plus what its children pass up once their joints are
    // free to move. A pivot that has vanished pins its joint, as the dense Cholesky does.
    void Articulation::factorMass() {
        const std::size_t n = links.size();
        factor_.resize(n);
        for (std::size_t i = 0; i < n; ++i) factor_[i].IA = ArticulatedInertia::of(links[i].I);
        for (std::size_t i = n; i-- > 0;) {
            const Link &L = links[i];
            Factor &F = factor_[i];
            const int k = L.ndof;
            double D[36];
            for (int a = 0; a < k; ++a) F.U[a] = F.IA * L.S[a];
            for (int a = 0; a < k; ++a)
                for (int b = 0; b < k; ++b) D[a * k + b] = dot(L.S[a], F.U[b]);
            choleskyFactor(D, k, k);
            for (int c = 0; c < k; ++c) {
                double col[6] = {};
                col[c] = 1.0;
                choleskySolve(D, k, k, col);
                for (int r = 0; r < k; ++r) F.Dinv[r][c] = col[r];
            }
            if (L.parent < 0) continue;
            ArticulatedInertia &P = factor_[std::size_t(L.parent)].IA;
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c) {
                    double s = F.IA.a[r][c];
                    for (int a = 0; a < k; ++a)
                        for (int b = 0; b < k; ++b) s -= component(F.U[a], r) * F.Dinv[a][b] * component(F.U[b], c);
                    P.a[r][c] += s;
                }
        }
    }

    // ABA with no velocity and no gravity: the inward pass gathers what each subtree can't absorb, the outward pass
    // turns it into joint accelerations, which are M⁻¹ x.
    void Articulation::solveMass(double *x) const {
        const std::size_t n = links.size();
        bias_.assign(n, Force{});
        acc_.resize(n);
        u_.resize(std::size_t(ndof_));
        for (std::size_t i = n; i-- > 0;) {
            const Link &L = links[i];
            const Factor &F = factor_[i];
            double *u = u_.data() + L.dof0;
            for (int a = 0; a < L.ndof; ++a) u[a] = x[L.dof0 + a] - dot(L.S[a], bias_[i]);
            if (L.parent < 0) continue;
            Force passed = bias_[i];
            for (int a = 0; a < L.ndof; ++a) {
                double du = 0.0;
                for (int b = 0; b < L.ndof; ++b) du += F.Dinv[a][b] * u[b];
                passed += du * F.U[a];
            }
            bias_[std::size_t(L.parent)] += passed;
        }
        for (std::size_t i = 0; i < n; ++i) {
            const Link &L = links[i];
            const Factor &F = factor_[i];
            const Motion ap = L.parent >= 0 ? acc_[std::size_t(L.parent)] : Motion{};
            double r[6];
            for (int a = 0; a < L.ndof; ++a) r[a] = u_[std::size_t(L.dof0 + a)] - dot(ap, F.U[a]);
            Motion ai = ap;
            for (int a = 0; a < L.ndof; ++a) {
                double qdd = 0.0;
                for (int b = 0; b < L.ndof; ++b) qdd += F.Dinv[a][b] * r[b];
                x[L.dof0 + a] = qdd;
                ai += qdd * L.S[a];
            }
            acc_[i] = ai;
        }
    }

    void Articulation::solveAt(int link, const Force *f, int count, double *out) const {
        const std::size_t n = links.size(), nd = std::size_t(ndof_);
        acc_.resize(3 * n);
        u_.resize(3 * nd);
        onPath_.resize(n, 0);
        if (++pathStamp_ == 0) {
            std::fill(onPath_.begin(), onPath_.end(), 0u);
            pathStamp_ = 1;
        }
        Force carried[3] = {};
        for (int i = link; i >= 0; i = links[std::size_t(i)].parent) {
            const Link &L = links[std::size_t(i)];
            const Factor &F = factor_[std::size_t(i)];
            onPath_[std::size_t(i)] = pathStamp_;
            for (int c = 0; c < count; ++c) {
                double *u = u_.data() + std::size_t(c) * nd + std::size_t(L.dof0);
                for (int a = 0; a < L.ndof; ++a) u[a] = dot(L.S[a], f[c]) - dot(L.S[a], carried[c]);
                for (int a = 0; a < L.ndof; ++a) {
                    double du = 0.0;
                    for (int b = 0; b < L.ndof; ++b) du += F.Dinv[a][b] * u[b];
                    carried[c] += du * F.U[a];
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            const Link &L = links[i];
            const Factor &F = factor_[i];
            const bool path = onPath_[i] == pathStamp_;
            for (int c = 0; c < count; ++c) {
                Motion *acc = acc_.data() + std::size_t(c) * n;
                const Motion ap = L.parent >= 0 ? acc[L.parent] : Motion{};
                const double *u = u_.data() + std::size_t(c) * nd + std::size_t(L.dof0);
                double r[6];
                for (int a = 0; a < L.ndof; ++a) r[a] = (path ? u[a] : 0.0) - dot(ap, F.U[a]);
                Motion ai = ap;
                double *x = out + std::size_t(c) * nd + std::size_t(L.dof0);
                for (int a = 0; a < L.ndof; ++a) {
                    double qdd = 0.0;
                    for (int b = 0; b < L.ndof; ++b) qdd += F.Dinv[a][b] * r[b];
                    x[a] = qdd;
                    ai += qdd * L.S[a];
                }
                acc[i] = ai;
            }
        }
    }

    void Articulation::addPointRow(int link, const Vec3 &xWorld, const Vec3 &dir, double scale, double *row) const {
        if (link < 0) return;
        const Vec3 rxd = cross(xWorld - O, dir);
        for (int j = link; j >= 0; j = links[std::size_t(j)].parent) {
            const Link &L = links[std::size_t(j)];
            for (int k = 0; k < L.ndof; ++k) row[L.dof0 + k] += scale * (dot(dir, L.S[k].v) + dot(rxd, L.S[k].w));
        }
    }

    void Articulation::addAngularRow(int link, const Vec3 &axis, double scale, double *row) const {
        if (link < 0) return;
        for (int j = link; j >= 0; j = links[std::size_t(j)].parent) {
            const Link &L = links[std::size_t(j)];
            for (int k = 0; k < L.ndof; ++k) row[L.dof0 + k] += scale * dot(axis, L.S[k].w);
        }
    }

    Motion Articulation::twistOf(int link, const double *qd) const {
        Motion V{};
        int root = link;
        for (int j = link; j >= 0; j = links[std::size_t(j)].parent) {
            const Link &L = links[std::size_t(j)];
            for (int k = 0; k < L.ndof; ++k) V += qd[L.dof0 + k] * L.S[k];
            root = j;
        }
        return V + anchorTwist(links[std::size_t(root)]);
    }

    void Articulation::twists(const double *qd, std::vector<Motion> &V) const {
        V.resize(links.size());
        for (std::size_t i = 0; i < links.size(); ++i) {
            const Link &L = links[i];
            Motion t = L.parent >= 0 ? V[std::size_t(L.parent)] : anchorTwist(L);
            for (int k = 0; k < L.ndof; ++k) t += qd[L.dof0 + k] * L.S[k];
            V[i] = t;
        }
    }

    Vec3 Articulation::linkVelocity(int link) const {
        return pointVelocity(twistOf(link, qdot.data()), links[std::size_t(link)].c - O);
    }

    Vec3 Articulation::linkOmega(int link) const { return twistOf(link, qdot.data()).w; }

    void Articulation::integrate(double h, const std::vector<double> &qd) {
        for (Link &L : links) {
            const std::size_t d = std::size_t(L.dof0);
            switch (L.kind) {
                case JointKind::Free: {
                    const Vec3 w{qd[d], qd[d + 1], qd[d + 2]};
                    const Vec3 vO{qd[d + 3], qd[d + 4], qd[d + 5]};
                    L.c += h * (vO + cross(w, L.c - O));
                    L.R = normalize(rotationVector(h * w) * L.R);
                    break;
                }
                case JointKind::Revolute:
                case JointKind::Prismatic: q[std::size_t(L.coord)] += h * qd[d]; break;
                case JointKind::Spherical: {
                    const Vec3 w{qd[d], qd[d + 1], qd[d + 2]};
                    Quat &r = qRot[std::size_t(L.coord)];
                    r = normalize(rotationVector(h * w) * r);
                    break;
                }
            }
        }
        forwardKinematics();
    }

    double Articulation::kineticEnergy() const {
        double e = 0.0;
        for (const Link &L : links) e += 0.5 * dot(L.V, L.I * L.V);
        return e;
    }

} // namespace hf
