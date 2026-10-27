#include "solver/Machine.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "math/Dense.hpp"

namespace hf {

    namespace {

        // Eigen-decomposes the symmetric k×k block G (k <= 3) in place by cyclic Jacobi: G ends diagonal and V
        // holds the eigenvectors as columns.
        void symmetricEigen(double G[3][3], double V[3][3], std::size_t k) {
            for (int sweep = 0; sweep < 12; ++sweep)
                for (std::size_t p = 0; p < k; ++p)
                    for (std::size_t q = p + 1; q < k; ++q) {
                        // Converged: smaller than this is rounding. Rotating on would only grind into denormals.
                        if (std::fabs(G[p][q]) <= 1e-15 * (std::fabs(G[p][p]) + std::fabs(G[q][q]))) {
                            G[p][q] = G[q][p] = 0.0;
                            continue;
                        }
                        const double th = (G[q][q] - G[p][p]) / (2.0 * G[p][q]);
                        const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1.0));
                        const double c = 1.0 / std::sqrt(t * t + 1.0), sn = t * c;
                        for (std::size_t r = 0; r < k; ++r) {
                            const double grp = G[r][p], grq = G[r][q];
                            G[r][p] = c * grp - sn * grq;
                            G[r][q] = sn * grp + c * grq;
                        }
                        for (std::size_t r = 0; r < k; ++r) {
                            const double gpr = G[p][r], gqr = G[q][r];
                            G[p][r] = c * gpr - sn * gqr;
                            G[q][r] = sn * gpr + c * gqr;
                        }
                        for (std::size_t r = 0; r < k; ++r) {
                            const double vrp = V[r][p], vrq = V[r][q];
                            V[r][p] = c * vrp - sn * vrq;
                            V[r][q] = sn * vrp + c * vrq;
                        }
                    }
        }

        // Held in the mobility: a row with no bounds of its own, nor any that follow another's.
        bool equality(const MachineRow &r) { return r.follows < 0 && std::isinf(r.lo) && std::isinf(r.hi); }

    } // namespace

    void Machine::dynamics(double h, const Vec3 &gravity, const std::vector<Force> &fext) {
        n_ = art.dofs();
        const std::size_t n = std::size_t(n_);
        art.updateKinematics(art.links.empty() ? Vec3{} : art.links[0].c);
        stats = MachineStats{};
        dense_ = n_ <= denseDofs;
        if (dense_) {
            art.massMatrix(M_);
            diagM_.resize(n);
            for (std::size_t i = 0; i < n; ++i) diagM_[i] = M_[i * n + i];
            L_ = M_;
            stats.pivotsRaised = choleskyFactor(L_.data(), n_, n_);
        } else {
            art.factorMass();
            art.massDiagonal(diagM_);
        }

        // What the machine rides on jumped since its velocities were measured: the joints take the impulse.
        if (art.anchorJumps(tmp_)) {
            solveMass(tmp_.data());
            for (std::size_t i = 0; i < n; ++i) art.qdot[i] -= tmp_[i];
        }

        // The velocity terms at the sub-step's midpoint velocity, found by a few fixed-point passes: taken at the start
        // (explicit), a thin link spinning about its own axis gains energy every step until it explodes.
        const std::vector<double> &qd0 = art.qdot;
        qdotFree_ = qd0;
        qdotMid_.resize(n);
        tmp_.resize(n);
        for (int pass = 0; pass < kMidpointPasses; ++pass) {
            for (std::size_t i = 0; i < n; ++i) qdotMid_[i] = 0.5 * (qd0[i] + qdotFree_[i]);
            art.updateVelocities(qdotMid_.data());
            art.biasForces(gravity, fext, C_);
            for (std::size_t i = 0; i < n; ++i) tmp_[i] = -C_[i];
            solveMass(tmp_.data());
            for (std::size_t i = 0; i < n; ++i) qdotFree_[i] = qd0[i] + h * tmp_[i];
        }
        art.qdot = qdotFree_;
        art.updateVelocities(art.qdot.data());

        pqdot_.assign(n, 0.0);
        clearRows();
        clearPosRows();
    }

    void Machine::solveMass(double *x) const {
        if (dense_) choleskySolve(L_.data(), n_, n_, x);
        else art.solveMass(x);
    }

    // Each row's M⁻¹ Jᵀ.
    void Machine::solveRowsOf(const std::vector<double> &J, std::size_t m, std::vector<double> &Y) const {
        const std::size_t n = std::size_t(n_);
        Y.assign(J.begin(), J.begin() + std::ptrdiff_t(m * n));
        for (std::size_t i = 0; i < m; ++i) solveMass(Y.data() + i * n);
    }

    double Machine::inverseMassAlong(const double *J) const {
        const std::size_t n = std::size_t(n_);
        tmp_.assign(J, J + n);
        solveMass(tmp_.data());
        return dot(J, tmp_.data(), n);
    }

    void Machine::clearRows() {
        rows_.clear();
        J_.clear();
    }

    double *Machine::addRow(const MachineRow &r) {
        rows_.push_back(r);
        J_.resize(rows_.size() * std::size_t(n_), 0.0);
        return J_.data() + (rows_.size() - 1) * std::size_t(n_);
    }

    void Machine::prepareRows() {
        const std::size_t n = std::size_t(n_), m = rows_.size();
        stats.rows = int(m);
        solveRowsOf(J_, m, Y_);
        A_.assign(m * m, 0.0);
        diag_.assign(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = J_.data() + i * n;
            for (std::size_t j = i; j < m; ++j) {
                const double s = dot(Ji, Y_.data() + j * n, n);
                A_[i * m + j] = s;
                A_[j * m + i] = s;
            }
            diag_[i] = A_[i * m + i];
        }
        regularise(A_, diag_, rows_, J_, dead_);
    }

    // Every row gets a compliance, so A is positive definite and a redundant linkage resolves in the minimum-norm
    // sense instead of its rows fighting. No row is dropped and none can explode: an inconsistent loop yields at
    // its cut instead of locking or launching the light links.
    void Machine::regularise(std::vector<double> &A, const std::vector<double> &diag,
                             const std::vector<MachineRow> &rows, const std::vector<double> &J,
                             std::vector<std::uint8_t> &dead) {
        // Relative to the largest diagonal: the floor for a structurally empty row.
        constexpr double kFloor = 1e-9;
        // Mobility below this of a row's Jacobi scale means the tree already satisfies it (a planar loop's cut
        // axis rows, a duplicate pin). Kept, it took huge cancelling multipliers that leaked into light links.
        constexpr double kDead = 1e-6;
        // A Jacobi scale below this of the system's largest is rounding noise, not a row.
        constexpr double kNoise = 1e-12;

        const std::size_t m = rows.size(), n = std::size_t(n_);

        // The Jacobi scale Σ J² / M_kk: the row's effective mass were its dofs uncoupled, its natural size.
        jac_.assign(m, 0.0);
        double jacMax = 0.0;
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = J.data() + i * n;
            for (std::size_t c = 0; c < n; ++c)
                if (Ji[c] != 0.0) jac_[i] += Ji[c] * Ji[c] / std::max(diagM_[c], 1e-300);
            jacMax = std::max(jacMax, jac_[i]);
        }
        dead.assign(m, 0);
        for (std::size_t i = 0; i < m; ++i)
            if (!(diag[i] > kDead * jac_[i]) || !(jac_[i] > kNoise * jacMax)) dead[i] = 1;

        double maxDiag = 0.0;
        for (std::size_t i = 0; i < m; ++i) maxDiag = std::max(maxDiag, diag[i]);
        const double floor = std::max(kFloor * maxDiag, 1e-300);
        for (std::size_t i = 0; i < m; ++i) {
            const double jacobi = rows[i].yield > 0.0 ? jac_[i] : 0.0;
            A[i * m + i] += rows[i].compliance + std::max({rows[i].relEps * diag[i], rows[i].yield * jacobi, floor});
        }

        // Directional yield per group: each eigen-direction of the group's mobility block gets c = y² / g, with
        // y = groupYield × the group's Jacobi scale. A direction the linkage moves freely (g >> y) stays rigid; a
        // weak one, where a slightly out-of-plane loop's closure lives, is released instead of locking the loop.
        for (std::size_t i = 0; i < m; ++i) {
            const int g = rows[i].group;
            if (g < 0 || rows[i].groupYield <= 0.0) continue;
            if (i > 0 && rows[i - 1].group == g) continue; // first row of its group
            std::size_t k = 1;
            while (i + k < m && rows[i + k].group == g && k < 3) ++k;

            // The unregularised block: A off the diagonal, the raw diagonal on it.
            double G[3][3] = {}, V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            for (std::size_t a = 0; a < k; ++a)
                for (std::size_t b = 0; b < k; ++b) G[a][b] = (a == b) ? diag[i + a] : A[(i + a) * m + (i + b)];
            symmetricEigen(G, V, k);

            double gmax = 0.0;
            for (std::size_t a = 0; a < k; ++a) gmax = std::max(gmax, jac_[i + a]);
            if (!(gmax > 0.0)) continue;

            const double y = rows[i].groupYield * gmax;
            for (std::size_t e = 0; e < k; ++e) {
                const double ge = std::max(G[e][e], 1e-300);
                const double ce = std::min(y * y / ge, 1e6 * y); // a dead direction: released, finitely
                for (std::size_t a = 0; a < k; ++a)
                    for (std::size_t b = 0; b < k; ++b) A[(i + a) * m + (i + b)] += ce * V[a][e] * V[b][e];
            }
        }

        // Dead rows leave the system: unit diagonal, no coupling, and a zero right-hand side from the caller.
        for (std::size_t i = 0; i < m; ++i) {
            if (!dead[i]) continue;
            for (std::size_t j = 0; j < m; ++j) {
                A[i * m + j] = 0.0;
                A[j * m + i] = 0.0;
            }
            A[i * m + i] = 1.0;
        }
    }

    void Machine::solveWithBase(const std::vector<double> &base, std::vector<RowState> &state, BoxLcpStats &st) {
        const std::size_t n = std::size_t(n_), m = rows_.size();
        b_.resize(m);
        lo_.resize(m);
        hi_.resize(m);
        follow_.assign(m, Follow{});
        bool follows = false;
        for (std::size_t i = 0; i < m; ++i) {
            const MachineRow &r = rows_[i];
            const double s = dot(J_.data() + i * n, base.data(), n);
            const double compliance = A_[i * m + i] - diag_[i];
            b_[i] = dead_[i] ? 0.0 : s - r.target - compliance * r.anchor;
            lo_[i] = r.lo;
            hi_[i] = r.hi;
            if (r.follows >= 0) {
                follow_[i] = {r.follows, r.ratio};
                follows = true;
            }
        }
        if (!follows) follow_.clear();
        st = solveBoxLcp(A_, b_, lo_, hi_, int(m), lambda_, state, lcp_, follow_);
        art.qdot = base;
        for (std::size_t i = 0; i < m; ++i) {
            const double l = lambda_[i];
            if (l == 0.0) continue;
            const double *Yi = Y_.data() + i * n;
            for (std::size_t c = 0; c < n; ++c) art.qdot[c] += l * Yi[c];
        }
        stats.lcpPasses += st.passes;
        stats.lcpPivots += st.pivots;
        stats.lcpFallbacks += st.fallback ? 1 : 0;
    }

    BoxLcpStats Machine::solveRows(std::vector<RowState> &state) {
        BoxLcpStats st;
        state.resize(rows_.size(), RowState::Free);
        solveWithBase(qdotFree_, state, st);
        return st;
    }

    void Machine::buildMobility() {
        const std::size_t m = rows_.size();
        // Equality rows only. Holding a bounded row here would let the contacts lean on a drive that the next
        // solve finds past its cap: the forklift's braked wheels, released, spun backwards at 90 rad/s.
        freeIdx_.clear();
        for (std::size_t i = 0; i < m; ++i)
            if (equality(rows_[i])) freeIdx_.push_back(int(i));
        const std::size_t k = freeIdx_.size();
        Aff_.resize(k * k);
        for (std::size_t p = 0; p < k; ++p)
            for (std::size_t q = 0; q < k; ++q)
                Aff_[p * k + q] = A_[std::size_t(freeIdx_[p]) * m + std::size_t(freeIdx_[q])];
        if (k) choleskyFactor(Aff_.data(), int(k), int(k));
    }

    // W g = M⁻¹ g − Yᵀ A⁻¹ (J M⁻¹ g) over the held rows.
    double Machine::mobility(const double *g, double *w) const {
        const std::size_t n = std::size_t(n_);
        std::copy(g, g + n, w);
        solveMass(w);
        const double unconstrained = dot(g, w, n);
        holdRows(w);
        return unconstrained;
    }

    // M⁻¹ Jᵀ f for the six wrenches at a link: only the link's path is climbed on the way in.
    void Machine::solveAtLink(int link, const Force *f, double *w) const {
        const std::size_t n = std::size_t(n_);
        if (!dense_) {
            art.solveAt(link, f, 3, w);
            art.solveAt(link, f + 3, 3, w + 3 * n);
            return;
        }
        for (int c = 0; c < 6; ++c) {
            double *wc = w + std::size_t(c) * n;
            std::fill(wc, wc + n, 0.0);
            for (int j = link; j >= 0; j = art.links[std::size_t(j)].parent) {
                const Link &L = art.links[std::size_t(j)];
                for (int a = 0; a < L.ndof; ++a) wc[L.dof0 + a] += dot(L.S[a], f[c]);
            }
            solveMass(wc);
        }
    }

    void Machine::linkResponse(int link, double *Z, double held[6][6], double unheld[6][6]) const {
        const std::size_t n = std::size_t(n_);
        const Force unit[6] = {{Vec3{1, 0, 0}, {}}, {Vec3{0, 1, 0}, {}}, {Vec3{0, 0, 1}, {}},
                               {{}, Vec3{1, 0, 0}}, {{}, Vec3{0, 1, 0}}, {{}, Vec3{0, 0, 1}}};
        // The unheld response first: the link's own twist from it gives the 6×6.
        solveAtLink(link, unit, Z);
        for (int c = 0; c < 6; ++c) {
            const Motion t = art.twistOf(link, Z + std::size_t(c) * n);
            for (int r = 0; r < 6; ++r) unheld[r][c] = component(t, r);
        }
        for (int c = 0; c < 6; ++c) {
            holdRows(Z + std::size_t(c) * n);
            const Motion t = art.twistOf(link, Z + std::size_t(c) * n);
            for (int r = 0; r < 6; ++r) held[r][c] = component(t, r);
        }
    }

    void Machine::holdRows(double *w) const {
        const std::size_t n = std::size_t(n_), k = freeIdx_.size();
        if (k == 0) return;
        free_.resize(k);
        for (std::size_t p = 0; p < k; ++p) {
            free_[p] = dot(J_.data() + std::size_t(freeIdx_[p]) * n, w, n);
        }
        choleskySolve(Aff_.data(), int(k), int(k), free_.data());
        for (std::size_t p = 0; p < k; ++p) {
            const double z = free_[p];
            if (z == 0.0) continue;
            const double *Yp = Y_.data() + std::size_t(freeIdx_[p]) * n;
            for (std::size_t r = 0; r < n; ++r) w[r] -= z * Yp[r];
        }
    }

    void Machine::applyContact(const double *w, double mu, bool pseudo) {
        std::vector<double> &x = pseudo ? pqdot_ : art.qdot;
        for (std::size_t i = 0; i < std::size_t(n_); ++i) x[i] += mu * w[i];
    }

    void Machine::applyWrench(const double *Z, const Force &f, bool pseudo) {
        std::vector<double> &x = pseudo ? pqdot_ : art.qdot;
        const std::size_t n = std::size_t(n_);
        for (int c = 0; c < 6; ++c) {
            const double s = component(f, c);
            if (s == 0.0) continue;
            const double *Zc = Z + std::size_t(c) * n;
            for (std::size_t i = 0; i < n; ++i) x[i] += s * Zc[i];
        }
    }

    void Machine::prepareBounded() {
        const std::size_t n = std::size_t(n_), m = rows_.size();
        bIdx_.clear();
        bOf_.assign(m, -1);
        for (std::size_t i = 0; i < m; ++i)
            if (!equality(rows_[i])) {
                bOf_[i] = int(bIdx_.size());
                bIdx_.push_back(int(i));
            }
        bW_.assign(bIdx_.size() * n, 0.0);
        bK_.assign(bIdx_.size(), 0.0);
        bG_.assign(bIdx_.size(), 0.0);
        bK12_.assign(bIdx_.size(), 0.0);
        for (std::size_t b = 0; b < bIdx_.size(); ++b) {
            const std::size_t i = std::size_t(bIdx_[b]);
            const double *Ji = J_.data() + i * n;
            double *w = bW_.data() + b * n;
            mobility(Ji, w);
            bG_[b] = A_[i * m + i] - diag_[i]; // the compliance the direct solve gave it
            bK_[b] = dot(Ji, w, n) + bG_[b];
        }
        for (std::size_t b = 0; b + 1 < bIdx_.size(); ++b)
            if (rows_[std::size_t(bIdx_[b])].pairedWithNext)
                bK12_[b] = dot(J_.data() + std::size_t(bIdx_[b]) * n, bW_.data() + (b + 1) * n, n);
    }

    // Gauss-Seidel on the bounded rows, against the same compliant system the direct solve answered: its answer
    // is a fixed point here, and only what the contacts changed moves them. A friction pair's 2×2 block is solved
    // exactly and held within its circle: the direct solve could only box it.
    void Machine::sweepBounded() {
        const std::size_t n = std::size_t(n_);
        const auto residual = [&](std::size_t b) {
            const std::size_t i = std::size_t(bIdx_[b]);
            const MachineRow &r = rows_[i];
            return dot(J_.data() + i * n, art.qdot.data(), n) - r.target + bG_[b] * (lambda_[i] - r.anchor);
        };
        const auto apply = [&](std::size_t b, double nl) {
            const std::size_t i = std::size_t(bIdx_[b]);
            const double d = nl - lambda_[i];
            if (d == 0.0) return;
            lambda_[i] = nl;
            const double *w = bW_.data() + b * n;
            for (std::size_t c = 0; c < n; ++c) art.qdot[c] += d * w[c];
        };
        for (std::size_t b = 0; b < bIdx_.size(); ++b) {
            const std::size_t i = std::size_t(bIdx_[b]);
            const MachineRow &r = rows_[i];
            const double cap = r.follows >= 0 ? r.ratio * std::max(0.0, lambda_[std::size_t(r.follows)]) : 0.0;
            if (r.pairedWithNext && b + 1 < bIdx_.size()) {
                const std::size_t j = std::size_t(bIdx_[b + 1]);
                const double r1 = residual(b), r2 = residual(b + 1);
                const double a = bK_[b], c = bK12_[b], d = bK_[b + 1], det = a * d - c * c;
                // A row with no mobility takes nothing: its impulse would be a division by zero.
                double l1 = a > 0.0 ? lambda_[i] - r1 / a : lambda_[i], l2 = d > 0.0 ? lambda_[j] - r2 / d : lambda_[j];
                if (det > 1e-12 * a * d) {
                    l1 = lambda_[i] - (d * r1 - c * r2) / det;
                    l2 = lambda_[j] - (a * r2 - c * r1) / det;
                }
                const double len = std::hypot(l1, l2);
                if (!(len <= cap)) {
                    const double s = std::isfinite(len) && len > 0.0 ? cap / len : 0.0;
                    l1 *= s;
                    l2 *= s;
                }
                apply(b, l1);
                apply(++b, l2);
                continue;
            }
            if (!(bK_[b] > 0.0)) continue;
            const double lo = r.follows >= 0 ? -cap : r.lo, hi = r.follows >= 0 ? cap : r.hi;
            apply(b, std::clamp(lambda_[i] - residual(b) / bK_[b], lo, hi));
        }
    }

    void Machine::pushRow(int row, double target) {
        const int b = row >= 0 && std::size_t(row) < bOf_.size() ? bOf_[std::size_t(row)] : -1;
        if (b < 0) return;
        const std::size_t n = std::size_t(n_), i = std::size_t(row);
        const double v = dot(J_.data() + i * n, art.qdot.data(), n);
        const double k = bK_[std::size_t(b)] - bG_[std::size_t(b)];
        if (!(k > 0.0)) return;
        const double nl = std::clamp(lambda_[i] + (target - v) / k, rows_[i].lo, rows_[i].hi);
        const double d = nl - lambda_[i];
        lambda_[i] = nl;
        const double *w = bW_.data() + std::size_t(b) * n;
        for (std::size_t c = 0; c < n; ++c) art.qdot[c] += d * w[c];
    }

    void Machine::clearPosRows() {
        posRows_.clear();
        posJ_.clear();
    }

    double *Machine::addPosRow(const MachineRow &r) {
        posRows_.push_back(r);
        posJ_.resize(posRows_.size() * std::size_t(n_), 0.0);
        return posJ_.data() + (posRows_.size() - 1) * std::size_t(n_);
    }

    BoxLcpStats Machine::solvePositions(std::vector<RowState> &state) {
        const std::size_t n = std::size_t(n_), m = posRows_.size();
        stats.posRows = int(m);
        BoxLcpStats st;
        if (m == 0) return st;

        solveRowsOf(posJ_, m, posY_);
        posA_.assign(m * m, 0.0);
        posB_.resize(m);
        posLo_.resize(m);
        posHi_.resize(m);
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = posJ_.data() + i * n;
            for (std::size_t j = i; j < m; ++j) {
                const double s = dot(Ji, posY_.data() + j * n, n);
                posA_[i * m + j] = s;
                posA_[j * m + i] = s;
            }
            posB_[i] = dot(Ji, pqdot_.data(), n) - posRows_[i].target;
            posLo_[i] = posRows_[i].lo;
            posHi_[i] = posRows_[i].hi;
        }
        posDiag_.resize(m);
        for (std::size_t i = 0; i < m; ++i) posDiag_[i] = posA_[i * m + i];
        regularise(posA_, posDiag_, posRows_, posJ_, posDead_);
        for (std::size_t i = 0; i < m; ++i)
            if (posDead_[i]) posB_[i] = 0.0;

        state.resize(m, RowState::Free);
        st = solveBoxLcp(posA_, posB_, posLo_, posHi_, int(m), posLambda_, state, lcp_);
        for (std::size_t i = 0; i < m; ++i) {
            const double l = posLambda_[i];
            if (l == 0.0) continue;
            const double *Yi = posY_.data() + i * n;
            for (std::size_t c = 0; c < n; ++c) pqdot_[c] += l * Yi[c];
        }
        stats.lcpPasses += st.passes;
        stats.lcpPivots += st.pivots;
        stats.lcpFallbacks += st.fallback ? 1 : 0;
        return st;
    }

    void Machine::integrate(double h) {
        tmp_.resize(std::size_t(n_));
        for (std::size_t i = 0; i < std::size_t(n_); ++i) tmp_[i] = art.qdot[i] + pqdot_[i];
        art.integrate(h, tmp_);
    }

} // namespace hf
