#pragma once

// Dense symmetric positive (semi)definite solves, row-major with leading dimension ld.

namespace hf {

    // In-place Cholesky of the n×n leading block of A: on return its lower triangle holds L, A = L Lᵀ.
    // A pivot below relFloor × its own diagonal is raised to that floor, as a tiny compliance would; a zero
    // diagonal pins its unknown at zero. Returns how many pivots were raised or pinned.
    int choleskyFactor(double *A, int n, int ld, double relFloor = 1e-12);

    // Solves L Lᵀ x = b in place, with L from choleskyFactor.
    void choleskySolve(const double *L, int n, int ld, double *b);

} // namespace hf
