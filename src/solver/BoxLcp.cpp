#include "solver/BoxLcp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "math/Dense.hpp"

namespace hf {

    namespace {

        bool unbounded(double lo, double hi) { return std::isinf(lo) && lo < 0.0 && std::isinf(hi) && hi > 0.0; }

        double boundOf(RowState s, double lo, double hi) { return s == RowState::AtLo ? lo : hi; }

        double rowResidual(const std::vector<double> &A, const std::vector<double> &b, const std::vector<double> &x,
                           int m, int i) {
            return b[std::size_t(i)] + dot(A.data() + std::size_t(i) * std::size_t(m), x.data(), std::size_t(m));
        }

    } // namespace

    BoxLcpStats solveBoxLcp(const std::vector<double> &A, const std::vector<double> &b, const std::vector<double> &lo,
                            const std::vector<double> &hi, int m, std::vector<double> &x, std::vector<RowState> &state,
                            BoxLcpScratch &s, int maxPasses) {
        BoxLcpStats st;
        x.resize(std::size_t(m), 0.0);
        state.resize(std::size_t(m), RowState::Free);
        if (m == 0) return st;

        // The largest |b| is the largest velocity error in play: it scales the tolerance on w.
        double bScale = 1e-12;
        for (int i = 0; i < m; ++i) bScale = std::max(bScale, std::fabs(b[std::size_t(i)]));
        const double wTol = 1e-10 * bScale;

        // A warm start may name a bound that has since become infinite.
        for (int i = 0; i < m; ++i) {
            RowState &r = state[std::size_t(i)];
            if (unbounded(lo[std::size_t(i)], hi[std::size_t(i)])) r = RowState::Free;
            if (r == RowState::AtLo && !std::isfinite(lo[std::size_t(i)])) r = RowState::Free;
            if (r == RowState::AtHi && !std::isfinite(hi[std::size_t(i)])) r = RowState::Free;
        }

        s.w.resize(std::size_t(m));
        int bestCount = std::numeric_limits<int>::max();
        int stall = 0;

        for (int pass = 0; pass < maxPasses; ++pass) {
            ++st.passes;

            s.freeIdx.clear();
            for (int i = 0; i < m; ++i) {
                if (state[std::size_t(i)] == RowState::Free) s.freeIdx.push_back(i);
                else x[std::size_t(i)] = boundOf(state[std::size_t(i)], lo[std::size_t(i)], hi[std::size_t(i)]);
            }
            const int k = static_cast<int>(s.freeIdx.size());
            if (k > 0) {
                s.Aff.resize(std::size_t(k) * std::size_t(k));
                s.rhs.resize(std::size_t(k));
                for (int p = 0; p < k; ++p) {
                    const int i = s.freeIdx[std::size_t(p)];
                    const double *Ai = A.data() + std::size_t(i) * std::size_t(m);
                    double r = -b[std::size_t(i)];
                    for (int j = 0; j < m; ++j)
                        if (state[std::size_t(j)] != RowState::Free) r -= Ai[j] * x[std::size_t(j)];
                    s.rhs[std::size_t(p)] = r;
                    for (int q = 0; q < k; ++q)
                        s.Aff[std::size_t(p) * std::size_t(k) + std::size_t(q)] = Ai[s.freeIdx[std::size_t(q)]];
                }
                choleskyFactor(s.Aff.data(), k, k);
                choleskySolve(s.Aff.data(), k, k, s.rhs.data());
                for (int p = 0; p < k; ++p) x[std::size_t(s.freeIdx[std::size_t(p)])] = s.rhs[std::size_t(p)];
            }

            for (int i = 0; i < m; ++i) s.w[std::size_t(i)] = rowResidual(A, b, x, m, i);
            s.viol.clear();
            for (int i = 0; i < m; ++i) {
                const double l = lo[std::size_t(i)], h = hi[std::size_t(i)], xi = x[std::size_t(i)];
                const double xTol = 1e-12 * (1.0 + std::fabs(xi));
                switch (state[std::size_t(i)]) {
                    case RowState::Free:
                        if (xi < l - xTol || xi > h + xTol) s.viol.push_back(i);
                        break;
                    case RowState::AtLo:
                        if (s.w[std::size_t(i)] < -wTol) s.viol.push_back(i);
                        break;
                    case RowState::AtHi:
                        if (s.w[std::size_t(i)] > wTol) s.viol.push_back(i);
                        break;
                }
            }
            if (s.viol.empty()) return st;

            const int count = static_cast<int>(s.viol.size());
            if (count < bestCount) {
                bestCount = count;
                stall = 0;
            } else {
                ++stall;
            }

            const auto flip = [&](int i) {
                RowState &r = state[std::size_t(i)];
                if (r == RowState::Free) r = x[std::size_t(i)] < lo[std::size_t(i)] ? RowState::AtLo : RowState::AtHi;
                else r = RowState::Free;
                ++st.pivots;
            };
            // Block pivots while the violation count falls, Murty's least index once it stalls.
            if (stall < 3)
                for (int i : s.viol) flip(i);
            else flip(s.viol.front());
        }

        // Out of passes: projected Gauss-Seidel from wherever the pivoting got to.
        st.fallback = true;
        for (int i = 0; i < m; ++i)
            x[std::size_t(i)] = std::clamp(x[std::size_t(i)], lo[std::size_t(i)], hi[std::size_t(i)]);
        for (int sweep = 0; sweep < 200; ++sweep)
            for (int i = 0; i < m; ++i) {
                const double d = A[std::size_t(i) * std::size_t(m) + std::size_t(i)];
                if (!(d > 0.0)) continue;
                const double w = rowResidual(A, b, x, m, i);
                x[std::size_t(i)] = std::clamp(x[std::size_t(i)] - w / d, lo[std::size_t(i)], hi[std::size_t(i)]);
            }

        double worst = 0.0;
        for (int i = 0; i < m; ++i) {
            const double w = rowResidual(A, b, x, m, i);
            const double xi = x[std::size_t(i)];
            RowState &r = state[std::size_t(i)];
            if (xi <= lo[std::size_t(i)]) {
                r = RowState::AtLo;
                worst = std::max(worst, std::max(0.0, -w));
            } else if (xi >= hi[std::size_t(i)]) {
                r = RowState::AtHi;
                worst = std::max(worst, std::max(0.0, w));
            } else {
                r = RowState::Free;
                worst = std::max(worst, std::fabs(w));
            }
        }
        st.residual = worst;
        return st;
    }

} // namespace hf
