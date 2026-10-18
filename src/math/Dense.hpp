#pragma once

#include <cstddef>

// Dense symmetric positive (semi)definite solves, row-major with leading dimension ld.

namespace hf {

    // a · b over n entries, in four running sums so each add doesn't wait on the one before.
    inline double dot(const double *a, const double *b, std::size_t n) {
        double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
        std::size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            s0 += a[i] * b[i];
            s1 += a[i + 1] * b[i + 1];
            s2 += a[i + 2] * b[i + 2];
            s3 += a[i + 3] * b[i + 3];
        }
        for (; i < n; ++i) s0 += a[i] * b[i];
        return (s0 + s1) + (s2 + s3);
    }

    // In-place Cholesky of the n×n leading block of A: on return its lower triangle holds L, A = L Lᵀ.
    // A pivot below relFloor × its own diagonal is raised to that floor, as a tiny compliance would; a zero
    // diagonal pins its unknown at zero. Returns how many pivots were raised or pinned.
    int choleskyFactor(double *A, int n, int ld, double relFloor = 1e-12);

    // Solves L Lᵀ x = b in place, with L from choleskyFactor.
    void choleskySolve(const double *L, int n, int ld, double *b);

} // namespace hf
