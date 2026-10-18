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

    } // namespace

    void Machine::dynamics(double h, const Vec3 &gravity, const std::vector<Force> &fext,
                           const std::vector<double> &jointDiag) {
        n_ = art.dofs();
        const std::size_t n = std::size_t(n_);
        art.updateKinematics(art.links.empty() ? Vec3{} : art.links[0].c);
        art.massMatrix(M_);
        for (std::size_t i = 0; i < jointDiag.size() && i < n; ++i) M_[i * n + i] += jointDiag[i];

        L_ = M_;
        stats = MachineStats{};
        stats.pivotsRaised = choleskyFactor(L_.data(), n_, n_);

        Minv_.assign(n * n, 0.0);
        tmp_.assign(n, 0.0);
        for (std::size_t c = 0; c < n; ++c) {
            std::fill(tmp_.begin(), tmp_.end(), 0.0);
            tmp_[c] = 1.0;
            choleskySolve(L_.data(), n_, n_, tmp_.data());
            for (std::size_t r = 0; r < n; ++r) Minv_[r * n + c] = tmp_[r];
        }

        art.biasForces(gravity, fext, C_);
        for (std::size_t i = 0; i < n; ++i) tmp_[i] = -C_[i];
        choleskySolve(L_.data(), n_, n_, tmp_.data());
        qdotFree_.resize(n);
        for (std::size_t i = 0; i < n; ++i) qdotFree_[i] = art.qdot[i] + h * tmp_[i];
        art.qdot = qdotFree_;

        pqdot_.assign(n, 0.0);
        clearRows();
        clearPosRows();
    }

    double Machine::inverseMassAlong(const double *J) const {
        const std::size_t n = std::size_t(n_);
        double s = 0.0;
        for (std::size_t r = 0; r < n; ++r) {
            if (J[r] == 0.0) continue;
            double t = 0.0;
            for (std::size_t c = 0; c < n; ++c) t += Minv_[r * n + c] * J[c];
            s += J[r] * t;
        }
        return s;
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
        Y_.assign(m * n, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = J_.data() + i * n;
            double *Yi = Y_.data() + i * n;
            for (std::size_t c = 0; c < n; ++c) {
                const double jc = Ji[c];
                if (jc == 0.0) continue;
                for (std::size_t r = 0; r < n; ++r) Yi[r] += Minv_[r * n + c] * jc;
            }
        }
        A_.assign(m * m, 0.0);
        diag_.assign(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = J_.data() + i * n;
            for (std::size_t j = i; j < m; ++j) {
                const double *Yj = Y_.data() + j * n;
                double s = 0.0;
                for (std::size_t c = 0; c < n; ++c) s += Ji[c] * Yj[c];
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
                if (Ji[c] != 0.0) jac_[i] += Ji[c] * Ji[c] / std::max(M_[c * n + c], 1e-300);
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
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = J_.data() + i * n;
            double s = 0.0;
            for (std::size_t c = 0; c < n; ++c) s += Ji[c] * base[c];
            b_[i] = dead_[i] ? 0.0 : s - rows_[i].target;
            lo_[i] = rows_[i].lo;
            hi_[i] = rows_[i].hi;
        }
        st = solveBoxLcp(A_, b_, lo_, hi_, int(m), lambda_, state, lcp_);
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
        const std::size_t n = std::size_t(n_), m = rows_.size();
        W_ = Minv_;

        // Equality rows only. Holding a bounded row here would let the contacts lean on a drive that the next
        // solve finds past its cap: the forklift's braked wheels, released, spun backwards at 90 rad/s.
        freeIdx_.clear();
        for (std::size_t i = 0; i < m; ++i)
            if (std::isinf(rows_[i].lo) && std::isinf(rows_[i].hi)) freeIdx_.push_back(int(i));
        const std::size_t k = freeIdx_.size();
        if (k == 0) return;

        Aff_.resize(k * k);
        for (std::size_t p = 0; p < k; ++p)
            for (std::size_t q = 0; q < k; ++q)
                Aff_[p * k + q] = A_[std::size_t(freeIdx_[p]) * m + std::size_t(freeIdx_[q])];
        choleskyFactor(Aff_.data(), int(k), int(k));
        tmp_.resize(k);
        for (std::size_t c = 0; c < n; ++c) {
            for (std::size_t p = 0; p < k; ++p) tmp_[p] = Y_[std::size_t(freeIdx_[p]) * n + c];
            choleskySolve(Aff_.data(), int(k), int(k), tmp_.data());
            for (std::size_t p = 0; p < k; ++p) {
                const double z = tmp_[p];
                if (z == 0.0) continue;
                const double *Yp = Y_.data() + std::size_t(freeIdx_[p]) * n;
                for (std::size_t r = 0; r < n; ++r) W_[r * n + c] -= Yp[r] * z;
            }
        }

        // W is symmetric in exact arithmetic; a contact's effective mass g W g must not depend on rounding.
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t c = r + 1; c < n; ++c) {
                const double a = 0.5 * (W_[r * n + c] + W_[c * n + r]);
                W_[r * n + c] = a;
                W_[c * n + r] = a;
            }
    }

    void Machine::mobility(const double *g, double *w) const {
        const std::size_t n = std::size_t(n_);
        std::fill(w, w + n, 0.0);
        for (std::size_t c = 0; c < n; ++c) {
            const double gc = g[c];
            if (gc == 0.0) continue;
            for (std::size_t r = 0; r < n; ++r) w[r] += W_[r * n + c] * gc;
        }
    }

    void Machine::applyContact(const double *w, double mu) {
        for (std::size_t i = 0; i < std::size_t(n_); ++i) art.qdot[i] += mu * w[i];
    }

    void Machine::prepareBounded() {
        const std::size_t n = std::size_t(n_), m = rows_.size();
        bIdx_.clear();
        for (std::size_t i = 0; i < m; ++i)
            if (!(std::isinf(rows_[i].lo) && std::isinf(rows_[i].hi))) bIdx_.push_back(int(i));
        bW_.assign(bIdx_.size() * n, 0.0);
        bK_.assign(bIdx_.size(), 0.0);
        for (std::size_t b = 0; b < bIdx_.size(); ++b) {
            const std::size_t i = std::size_t(bIdx_[b]);
            const double *Ji = J_.data() + i * n;
            double *w = bW_.data() + b * n;
            mobility(Ji, w);
            double k = 0.0;
            for (std::size_t c = 0; c < n; ++c) k += Ji[c] * w[c];
            k += rows_[i].compliance + std::max(rows_[i].relEps * diag_[i], 1e-300);
            bK_[b] = k;
        }
    }

    void Machine::sweepBounded() {
        const std::size_t n = std::size_t(n_);
        for (std::size_t b = 0; b < bIdx_.size(); ++b) {
            const std::size_t i = std::size_t(bIdx_[b]);
            const MachineRow &r = rows_[i];
            const double *Ji = J_.data() + i * n;
            double v = 0.0;
            for (std::size_t c = 0; c < n; ++c) v += Ji[c] * art.qdot[c];
            const double old = lambda_[i];
            const double nl = std::clamp(old - (v - r.target + r.compliance * old) / bK_[b], r.lo, r.hi);
            const double d = nl - old;
            if (d == 0.0) continue;
            lambda_[i] = nl;
            const double *w = bW_.data() + b * n;
            for (std::size_t c = 0; c < n; ++c) art.qdot[c] += d * w[c];
        }
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

        posY_.assign(m * n, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = posJ_.data() + i * n;
            double *Yi = posY_.data() + i * n;
            for (std::size_t c = 0; c < n; ++c) {
                const double jc = Ji[c];
                if (jc == 0.0) continue;
                for (std::size_t r = 0; r < n; ++r) Yi[r] += Minv_[r * n + c] * jc;
            }
        }
        posA_.assign(m * m, 0.0);
        posB_.resize(m);
        posLo_.resize(m);
        posHi_.resize(m);
        for (std::size_t i = 0; i < m; ++i) {
            const double *Ji = posJ_.data() + i * n;
            for (std::size_t j = i; j < m; ++j) {
                const double *Yj = posY_.data() + j * n;
                double s = 0.0;
                for (std::size_t c = 0; c < n; ++c) s += Ji[c] * Yj[c];
                posA_[i * m + j] = s;
                posA_[j * m + i] = s;
            }
            double s = 0.0;
            for (std::size_t c = 0; c < n; ++c) s += Ji[c] * pqdot_[c];
            posB_[i] = s - posRows_[i].target;
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

    void Machine::applyPseudo(const double *w, double mu) {
        for (std::size_t i = 0; i < std::size_t(n_); ++i) pqdot_[i] += mu * w[i];
    }

    void Machine::integrate(double h) {
        tmp_.resize(std::size_t(n_));
        for (std::size_t i = 0; i < std::size_t(n_); ++i) tmp_[i] = art.qdot[i] + pqdot_[i];
        art.integrate(h, tmp_);
    }

} // namespace hf
