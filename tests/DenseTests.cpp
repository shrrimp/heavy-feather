#include <cstddef>
#include <vector>

#include <doctest/doctest.h>

#include "Random.hpp"
#include "math/Dense.hpp"

using namespace hf;

TEST_CASE("dense Cholesky") {
    Random rnd;
    for (int t = 0; t < 10; ++t) {
        const int n = 1 + int(rnd.uni(0, 30));
        const std::size_t N = std::size_t(n);
        std::vector<double> B(N * N), A(N * N, 0.0), x(N), b(N, 0.0);
        for (double &v : B) v = rnd.uni(-1, 1);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                double s = (i == j) ? 0.1 : 0.0;
                for (std::size_t k = 0; k < N; ++k) s += B[i * N + k] * B[j * N + k];
                A[i * N + j] = s;
            }
        for (double &v : x) v = rnd.uni(-1, 1);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) b[i] += A[i * N + j] * x[j];

        std::vector<double> L = A;
        CHECK_MESSAGE(choleskyFactor(L.data(), n, n) == 0, "an SPD matrix factors without raised pivots");
        choleskySolve(L.data(), n, n, b.data());
        CHECK_MESSAGE(maxAbsDiff(b, x) < 1e-8, "the solve recovers x");
    }
}
