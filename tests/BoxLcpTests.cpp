#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <doctest/doctest.h>

#include "Random.hpp"
#include "math/Dense.hpp"
#include "solver/BoxLcp.hpp"

using namespace hf;

namespace {

    // Tries every assignment of Free / AtLo / AtHi and keeps the one that satisfies every condition.
    bool bruteLcp(const std::vector<double> &A, const std::vector<double> &b, const std::vector<double> &lo,
                  const std::vector<double> &hi, int m, std::vector<double> &xOut) {
        const std::size_t M = std::size_t(m);
        int combos = 1;
        for (int i = 0; i < m; ++i) combos *= 3;
        std::vector<RowState> st(M);
        std::vector<double> x(M);
        for (int c = 0; c < combos; ++c) {
            int r = c;
            bool ok = true;
            for (std::size_t i = 0; i < M; ++i) {
                st[i] = RowState(r % 3);
                r /= 3;
                if (st[i] == RowState::AtLo && !std::isfinite(lo[i])) ok = false;
                if (st[i] == RowState::AtHi && !std::isfinite(hi[i])) ok = false;
            }
            if (!ok) continue;

            std::vector<std::size_t> F;
            for (std::size_t i = 0; i < M; ++i) {
                if (st[i] == RowState::Free) F.push_back(i);
                else x[i] = st[i] == RowState::AtLo ? lo[i] : hi[i];
            }
            const std::size_t k = F.size();
            std::vector<double> Aff(k * k), rhs(k);
            for (std::size_t p = 0; p < k; ++p) {
                double s = -b[F[p]];
                for (std::size_t j = 0; j < M; ++j)
                    if (st[j] != RowState::Free) s -= A[F[p] * M + j] * x[j];
                rhs[p] = s;
                for (std::size_t q = 0; q < k; ++q) Aff[p * k + q] = A[F[p] * M + F[q]];
            }
            if (k) {
                choleskyFactor(Aff.data(), int(k), int(k));
                choleskySolve(Aff.data(), int(k), int(k), rhs.data());
            }
            for (std::size_t p = 0; p < k; ++p) x[F[p]] = rhs[p];

            const double tol = 1e-9;
            for (std::size_t i = 0; i < M && ok; ++i) {
                double w = b[i];
                for (std::size_t j = 0; j < M; ++j) w += A[i * M + j] * x[j];
                switch (st[i]) {
                    case RowState::Free: ok = x[i] >= lo[i] - tol && x[i] <= hi[i] + tol; break;
                    case RowState::AtLo: ok = w >= -tol; break;
                    case RowState::AtHi: ok = w <= tol; break;
                }
            }
            if (ok) {
                xOut = x;
                return true;
            }
        }
        return false;
    }

} // namespace

TEST_CASE("boxed LCP, against brute force") {
    const double inf = std::numeric_limits<double>::infinity();
    Random rnd;
    BoxLcpScratch scratch;
    for (int t = 0; t < 300; ++t) {
        const int m = 1 + int(rnd.uni(0, 6.999));
        const std::size_t M = std::size_t(m);
        std::vector<double> B(M * M), A(M * M), b(M), lo(M), hi(M);
        for (double &v : B) v = rnd.uni(-1, 1);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) {
                double s = (i == j) ? 1e-3 : 0.0;
                for (std::size_t k = 0; k < M; ++k) s += B[i * M + k] * B[j * M + k];
                A[i * M + j] = s;
            }
        for (std::size_t i = 0; i < M; ++i) {
            b[i] = rnd.uni(-2, 2);
            const double kind = rnd.uni(0, 1);
            if (kind < 0.25) { // equality
                lo[i] = -inf;
                hi[i] = inf;
            } else if (kind < 0.5) { // contact or limit
                lo[i] = 0.0;
                hi[i] = inf;
            } else { // capped drive
                const double c = rnd.uni(0.01, 1.0);
                lo[i] = -c;
                hi[i] = c;
            }
        }

        std::vector<double> x, ref;
        std::vector<RowState> st;
        const BoxLcpStats s = solveBoxLcp(A, b, lo, hi, m, x, st, scratch);
        CHECK_MESSAGE(!s.fallback, "pivoting finishes without the PGS backstop");
        CHECK_MESSAGE(bruteLcp(A, b, lo, hi, m, ref), "brute force finds the unique solution");
        CHECK_MESSAGE(maxAbsDiff(x, ref) < 1e-7, "pivoting agrees with brute force");

        const BoxLcpStats s2 = solveBoxLcp(A, b, lo, hi, m, x, st, scratch);
        CHECK_MESSAGE((s2.passes == 1 && s2.pivots == 0), "a correct warm start is accepted in one pass");
    }
}
