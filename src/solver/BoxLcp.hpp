#pragma once

#include <cstdint>
#include <vector>

namespace hf {

    enum class RowState : std::uint8_t { Free = 0, AtLo = 1, AtHi = 2 };

    struct BoxLcpStats {
        int passes = 0;        // exact solves of the free set
        int pivots = 0;        // rows that changed state
        bool fallback = false; // the pivot budget ran out and the answer is the PGS backstop's
        double residual = 0.0; // worst complementarity violation left, in w units
    };

    struct BoxLcpScratch {
        std::vector<int> freeIdx;
        std::vector<double> Aff, rhs, w;
        std::vector<int> viol;
    };

    // Finds x with lo <= x <= hi and w = A x + b, where w = 0 for x strictly inside its bounds, w >= 0 at lo
    // and w <= 0 at hi. A is m×m row-major, symmetric positive definite, so the answer is unique. Solved by
    // block principal pivoting; x and state come in as the warm start and go out as the answer.
    BoxLcpStats solveBoxLcp(const std::vector<double> &A, const std::vector<double> &b, const std::vector<double> &lo,
                            const std::vector<double> &hi, int m, std::vector<double> &x, std::vector<RowState> &state,
                            BoxLcpScratch &scratch, int maxPasses = 64);

} // namespace hf
