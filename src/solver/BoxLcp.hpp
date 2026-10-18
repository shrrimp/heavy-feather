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
        std::vector<double> Aff, rhs, w, lo, hi;
        std::vector<int> viol;
    };

    // A row whose bounds follow another row's answer, as friction follows its contact's push:
    // |x_i| <= ratio × max(0, x_row). row -1: the row keeps its own bounds.
    struct Follow {
        int row = -1;
        double ratio = 0.0;
    };

    // Finds x with lo <= x <= hi and w = A x + b, where w = 0 for x strictly inside its bounds, w >= 0 at lo
    // and w <= 0 at hi. A is m×m row-major, symmetric positive definite, so the answer is unique. Solved by
    // block principal pivoting; x and state come in as the warm start and go out as the answer. With follow (empty,
    // or one per row) some bounds move with the answer, and the passes go on until they agree with it.
    BoxLcpStats solveBoxLcp(const std::vector<double> &A, const std::vector<double> &b, const std::vector<double> &lo,
                            const std::vector<double> &hi, int m, std::vector<double> &x, std::vector<RowState> &state,
                            BoxLcpScratch &scratch, const std::vector<Follow> &follow = {}, int maxPasses = 64);

} // namespace hf
