#pragma once

#include <cstdint>
#include <limits>
#include <vector>

#include "dynamics/Articulation.hpp"
#include "solver/BoxLcp.hpp"

namespace hf {

    // One row over a machine's generalised velocities: J qdot + compliance λ = target, with lo <= λ <= hi.
    struct MachineRow {
        double target = 0.0; // wanted J qdot (velocity rows) or J pqdot (position rows)
        double lo = -std::numeric_limits<double>::infinity();
        double hi = std::numeric_limits<double>::infinity();
        double compliance = 0.0; // a spring's
        double relEps = 1e-9;    // compliance relative to the row's own diagonal
        double yield = 0.0;      // compliance relative to its Jacobi scale
        // Rows sharing a group (one cut's pin rows, or its axis rows) yield together along their weak directions,
        // by groupYield of the group's Jacobi scale. -1: no group.
        int group = -1;
        double groupYield = 0.0;
        // A friction row: bounded by ±ratio × the multiplier of row follows, whatever the solve finds it to be. Paired
        // with the next row, the two are bounded together, by a circle, wherever the sweeps take them.
        int follows = -1;
        double ratio = 0.0;
        bool pairedWithNext = false;
        // The multiplier the row's compliance measures from. Anchored at its last answer, a soft row gives while it
        // changes and holds rigidly once it doesn't.
        double anchor = 0.0;
    };

    struct MachineStats {
        int rows = 0, posRows = 0;
        int lcpPasses = 0, lcpPivots = 0, lcpFallbacks = 0;
        int pivotsRaised = 0; // mass-matrix pivots the Cholesky had to floor
    };

    // One machine's solve, per sub-step: its articulation's dynamics plus every row the tree doesn't satisfy by
    // construction, solved together and directly. In order:
    //   dynamics, addRow..., prepareRows, solveRows, buildMobility, prepareBounded,
    //   then contacts through mobility() and applyContact() alongside sweepBounded(),
    //   then the position pass, then integrate.
    class Machine {
    public:
        Articulation art;
        MachineStats stats;

        // M factored (dense up to a few dofs, the articulated-body recursion beyond), C (RNEA) at the sub-step's
        // midpoint velocity, and the free velocity qdot + h M⁻¹(-C).
        void dynamics(double h, const Vec3 &gravity, const std::vector<Force> &fext);
        // Machines up to this many dofs are solved through a dense Cholesky of M; larger ones in O(n).
        int denseDofs = 12;

        // Adds a row and returns its J, dofs() zeros to fill in.
        double *addRow(const MachineRow &r);
        int rowCount() const { return int(rows_.size()); }
        const double *rowJ(int i) const { return J_.data() + std::size_t(i) * std::size_t(n_); }

        // J M⁻¹ Jᵀ for a row not added yet: what a spring's frequency needs to become a stiffness.
        double inverseMassAlong(const double *J) const;

        // Y = M⁻¹ Jᵀ and A = J M⁻¹ Jᵀ, regularised. Call once the rows are in.
        void prepareRows();

        // The rows' direct solve from the free velocity. state is the warm start, kept per row by the caller.
        BoxLcpStats solveRows(std::vector<RowState> &state);

        // W = M⁻¹ − Yᵀ A⁻¹ Y over the equality rows: the machine's response to an impulse with them held. Only A's
        // block is factored; W itself is never formed.
        void buildMobility();

        // w = W g. Returns gᵀ M⁻¹ g, the mobility with nothing held.
        double mobility(const double *g, double *w) const;

        // A link's response to the six unit wrenches about O: Z, six rows of dofs() values with the held rows kept,
        // and the link's 6×6 inverse inertia with them held and with nothing held. Any wrench ψ at the link then
        // moves the machine by Σ ψ_c Z_c, its own twist by held ψ, and every contact on one link shares one solve.
        void linkResponse(int link, double *Z, double held[6][6], double unheld[6][6]) const;

        // Applies an impulse of size mu along a generalised direction g, given w = W g: to the velocity, or to the
        // position pass's pseudo-velocity.
        void applyContact(const double *w, double mu, bool pseudo = false);
        // Adds Σ f_c Z_c to the velocity, or to the pseudo-velocity: a wrench at a link, through its response.
        void applyWrench(const double *Z, const Force &f, bool pseudo);

        // Bounded rows (capped drives, limits) are swept through W beside the contacts, so a drive at its cap and
        // the contact loading it agree. prepareBounded once after buildMobility, sweepBounded every iteration.
        void prepareBounded();
        void sweepBounded();
        int boundedRows() const { return int(bIdx_.size()); }
        // One bounded row, pushed on its own toward a velocity: a contact's bounce, after the sub-steps.
        void pushRow(int row, double target);

        // The position pass: its own rows, solved for a pseudo-velocity that integrate spends and then forgets,
        // so correcting drift puts no momentum into the machine.
        void clearPosRows();
        double *addPosRow(const MachineRow &r);
        BoxLcpStats solvePositions(std::vector<RowState> &state);
        const std::vector<double> &pseudo() const { return pqdot_; }

        // The line search on the pseudo-velocity: keepPseudo, then scale to or restore any fraction of it.
        void keepPseudo() { pqdotKept_ = pqdot_; }
        void scalePseudo(double k) {
            if (pqdotKept_.size() != pqdot_.size()) pqdotKept_ = pqdot_;
            for (std::size_t i = 0; i < pqdot_.size(); ++i) pqdot_[i] = k * pqdotKept_[i];
        }
        void restorePseudo(double k = 1.0) {
            if (pqdotKept_.size() == pqdot_.size())
                for (std::size_t i = 0; i < pqdot_.size(); ++i) pqdot_[i] = k * pqdotKept_[i];
        }

        // q along qdot + pqdot, then poses.
        void integrate(double h);

        int dofs() const { return n_; }
        const std::vector<double> &lambda() const { return lambda_; }

    private:
        int n_ = 0;
        bool dense_ = true;
        static constexpr int kMidpointPasses = 3;
        std::vector<double> M_, L_, diagM_, C_, qdotFree_, qdotMid_, pqdot_, pqdotKept_;
        mutable std::vector<double> tmp_, free_;

        std::vector<MachineRow> rows_;
        std::vector<double> J_, Y_, A_, diag_, b_, lo_, hi_, lambda_;
        std::vector<std::uint8_t> dead_;
        std::vector<int> freeIdx_;
        std::vector<double> Aff_;
        std::vector<Follow> follow_;
        std::vector<int> bIdx_, bOf_;      // bounded rows, and each row's place among them (-1: not bounded)
        std::vector<double> bW_, bK_, bG_, bK12_; // their W Jᵀ, effective inverse mass, compliance, pair coupling

        std::vector<MachineRow> posRows_;
        std::vector<double> posJ_, posY_, posA_, posDiag_, posB_, posLo_, posHi_, posLambda_;
        std::vector<std::uint8_t> posDead_;

        BoxLcpScratch lcp_;
        std::vector<double> jac_;

        void solveMass(double *x) const; // x <- M⁻¹ x
        void solveAtLink(int link, const Force *f, double *w) const;
        void holdRows(double *w) const; // w <- w − Yᵀ A⁻¹ J w over the held rows
        void solveRowsOf(const std::vector<double> &J, std::size_t m, std::vector<double> &Y) const;
        void clearRows();
        void solveWithBase(const std::vector<double> &base, std::vector<RowState> &state, BoxLcpStats &st);
        void regularise(std::vector<double> &A, const std::vector<double> &diag, const std::vector<MachineRow> &rows,
                        const std::vector<double> &J, std::vector<std::uint8_t> &dead);
    };

} // namespace hf
