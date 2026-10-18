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
            double d = orig;
            for (int k = 0; k < j; ++k) d -= Aj[k] * Aj[k];
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
                double s = Ai[j];
                for (int k = 0; k < j; ++k) s -= Ai[k] * Aj[k];
                Ai[j] = s * inv;
            }
        }
        return raised;
    }

    void choleskySolve(const double *L, int n, int ld, double *b) {
        for (int i = 0; i < n; ++i) {
            const double *Li = L + std::size_t(i) * std::size_t(ld);
            double s = b[i];
            for (int k = 0; k < i; ++k) s -= Li[k] * b[k];
            b[i] = s / Li[i];
        }
        for (int i = n - 1; i >= 0; --i) {
            double s = b[i];
            for (int k = i + 1; k < n; ++k) s -= L[std::size_t(k) * std::size_t(ld) + std::size_t(i)] * b[k];
            b[i] = s / L[std::size_t(i) * std::size_t(ld) + std::size_t(i)];
        }
    }

} // namespace hf
