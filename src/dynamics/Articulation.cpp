#include "dynamics/Articulation.hpp"

#include <cmath>
#include <cstddef>

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
            const Motion Vp = L.parent >= 0 ? links[std::size_t(L.parent)].V : Motion{};
            const Motion Ap = L.parent >= 0 ? links[std::size_t(L.parent)].A : Motion{};
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
            Motion vj{};
            for (int k = 0; k < L.ndof; ++k) vj += qdot[std::size_t(L.dof0 + k)] * L.S[k];
            L.V = Vp + vj;
            // S rides on the parent, so dS/dt = Vp × S and the joint adds Vp × (S qdot).
            L.A = Ap + crossM(Vp, vj);
            L.I = Inertia::body(L.mass, L.c - O, worldInertia(L.R, L.inertiaBody));
        }

        for (Link &L : links) L.Ic = L.I;
        for (std::size_t i = links.size(); i-- > 0;)
            if (links[i].parent >= 0) links[std::size_t(links[i].parent)].Ic += links[i].Ic;
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

    void Articulation::biasForces(const Vec3 &g, const std::vector<Force> &fext, std::vector<double> &C) {
        C.assign(std::size_t(ndof_), 0.0);
        const Motion ag{Vec3{}, g};
        for (std::size_t i = 0; i < links.size(); ++i) {
            Link &L = links[i];
            L.F = L.I * (L.A - ag) + crossF(L.V, L.I * L.V);
            if (!fext.empty()) L.F = L.F - fext[i];
        }
        for (std::size_t i = links.size(); i-- > 0;) {
            Link &L = links[i];
            for (int a = 0; a < L.ndof; ++a) C[std::size_t(L.dof0 + a)] = dot(L.S[a], L.F);
            if (L.parent >= 0) links[std::size_t(L.parent)].F += L.F;
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
        for (int j = link; j >= 0; j = links[std::size_t(j)].parent) {
            const Link &L = links[std::size_t(j)];
            for (int k = 0; k < L.ndof; ++k) V += qd[L.dof0 + k] * L.S[k];
        }
        return V;
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
