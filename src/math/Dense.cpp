#include "math/Dense.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace hf {

    int choleskyFactor(double *A, int n, int ld, double relFloor) {
        int raised = 0;
        for (int j = 0; j < n; ++j) {
            double *Aj = A + std::size_t(j) * std::size_t(ld);
            const double orig = Aj[j];
            double d = orig - dot(Aj, Aj, std::size_t(j));
            const double floor = relFloor * std::max(orig, 0.0);
            if (!(orig > 0.0)) {
                Aj[j] = 1.0;
                for (int k = 0; k < j; ++k) Aj[k] = 0.0;
                for (int i = j + 1; i < n; ++i) A[std::size_t(i) * std::size_t(ld) + std::size_t(j)] = 0.0;
                ++raised;
                continue;
            }
            if (!(d > floor)) {
                d = std::max(floor, 1e-300);
                ++raised;
            }
            const double l = std::sqrt(d);
            Aj[j] = l;
            const double inv = 1.0 / l;
            for (int i = j + 1; i < n; ++i) {
                double *Ai = A + std::size_t(i) * std::size_t(ld);
                Ai[j] = (Ai[j] - dot(Ai, Aj, std::size_t(j))) * inv;
            }
        }
        return raised;
    }

    void choleskySolve(const double *L, int n, int ld, double *b) {
        for (int i = 0; i < n; ++i) {
            const double *Li = L + std::size_t(i) * std::size_t(ld);
            b[i] = (b[i] - dot(Li, b, std::size_t(i))) / Li[i];
        }
        // Lᵀ by rows of L: each solved unknown is taken out of the ones above it.
        for (int i = n - 1; i >= 0; --i) {
            const double *Li = L + std::size_t(i) * std::size_t(ld);
            b[i] /= Li[i];
            for (int k = 0; k < i; ++k) b[k] -= Li[k] * b[i];
        }
    }

} // namespace hf
