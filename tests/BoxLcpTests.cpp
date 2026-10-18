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

TEST_CASE("friction rows bounded by the push they ride on") {
    // A box landing on up to four corners, each a push-only normal and two friction rows within mu of it, the rows
    // as soft as a machine's contact rows: four corners are twelve rows on six freedoms.
    const double inf = std::numeric_limits<double>::infinity();
    Random rnd(7);
    BoxLcpScratch scratch;
    int fallbacks = 0;
    for (int t = 0; t < 300; ++t) {
        const double mass = rnd.uni(1.0, 100.0);
        const Mat3 invI = inverse(rnd.inertia(mass));
        const Vec3 v = Vec3{rnd.uni(-2, 2), rnd.uni(-3, -0.1), rnd.uni(-2, 2)}, w = rnd.vec(2.0);
        const int k = 1 + int(rnd.uni(0, 3.999)), m = 3 * k;
        const std::size_t M = std::size_t(m);
        std::vector<Vec3> lin(M), ang(M);
        for (int c = 0; c < k; ++c) {
            const Vec3 r{c & 1 ? 0.5 : -0.5, -0.5, c & 2 ? 0.5 : -0.5};
            const Vec3 n = normalize(Vec3{0, 1, 0} + rnd.vec(0.2));
            const Vec3 t1 = normalize(cross(n, Vec3{0, 0, 1})), t2 = cross(n, t1);
            const Vec3 dirs[3] = {n, t1, t2};
            for (int d = 0; d < 3; ++d) {
                lin[std::size_t(3 * c + d)] = dirs[d];
                ang[std::size_t(3 * c + d)] = cross(r, dirs[d]);
            }
        }
        std::vector<double> A(M * M), b(M), lo(M, 0.0), hi(M, 0.0);
        std::vector<Follow> follow(M);
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < M; ++j)
                A[i * M + j] = dot(lin[i], lin[j]) / mass + dot(ang[i], invI * ang[j]);
            b[i] = dot(lin[i], v) + dot(ang[i], w);
            if (i % 3 == 0) hi[i] = inf;
            else follow[i] = {int(i - i % 3), rnd.uni(0.1, 1.0)};
        }
        for (std::size_t i = 0; i < M; ++i) A[i * M + i] *= i % 3 == 0 ? 1.15 : 1.01;

        std::vector<double> x;
        std::vector<RowState> st;
        fallbacks += solveBoxLcp(A, b, lo, hi, m, x, st, scratch, follow).fallback;

        for (std::size_t i = 0; i < M; ++i) {
            double wi = b[i];
            for (std::size_t j = 0; j < M; ++j) wi += A[i * M + j] * x[j];
            const double tol = 1e-6;
            if (i % 3 == 0) {
                CHECK(x[i] >= -tol);
                CHECK(wi >= -tol);
                CHECK(std::fabs(x[i] * wi) < tol);
                continue;
            }
            const double c = follow[i].ratio * x[i - i % 3];
            CHECK(std::fabs(x[i]) <= c + tol);
            if (x[i] < c - tol && x[i] > -c + tol) CHECK(std::fabs(wi) < tol);
            else if (x[i] >= c - tol && c > tol) CHECK(wi <= tol);
            else if (c > tol) CHECK(wi >= -tol);
        }
    }
    CHECK_MESSAGE(fallbacks <= 3, "the pivoting settles the moving bounds, bar a rare backstop: " << fallbacks);
}
