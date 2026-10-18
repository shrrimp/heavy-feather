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

        // M (CRBA), C (RNEA) and the free velocity qdot + h M⁻¹(-C). jointDiag is added to M's diagonal (implicit
        // joint damping or stiffness) and may be empty.
        void dynamics(double h, const Vec3 &gravity, const std::vector<Force> &fext,
                      const std::vector<double> &jointDiag = {});

        // Adds a row and returns its J, dofs() zeros to fill in.
        double *addRow(const MachineRow &r);
        int rowCount() const { return int(rows_.size()); }
        const double *rowJ(int i) const { return J_.data() + std::size_t(i) * std::size_t(n_); }

        // J M⁻¹ Jᵀ for a row not added yet: what a spring's frequency needs to become a stiffness.
        double inverseMassAlong(const double *J) const;
        double unconstrainedMobility(const double *g) const { return inverseMassAlong(g); }

        // Y = M⁻¹ Jᵀ and A = J M⁻¹ Jᵀ, regularised. Call once the rows are in.
        void prepareRows();

        // The rows' direct solve from the free velocity. state is the warm start, kept per row by the caller.
        BoxLcpStats solveRows(std::vector<RowState> &state);

        // W = M⁻¹ − Yᵀ A⁻¹ Y over the equality rows: the machine's response to an impulse with them held.
        void buildMobility();

        // w = W g
        void mobility(const double *g, double *w) const;

        // Applies an impulse of size mu along a generalised direction g, given w = W g.
        void applyContact(const double *w, double mu);

        // Bounded rows (capped drives, limits) are swept through W beside the contacts, so a drive at its cap and
        // the contact loading it agree. prepareBounded once after buildMobility, sweepBounded every iteration.
        void prepareBounded();
        void sweepBounded();
        int boundedRows() const { return int(bIdx_.size()); }

        // The position pass: its own rows, solved for a pseudo-velocity that integrate spends and then forgets,
        // so correcting drift puts no momentum into the machine.
        void clearPosRows();
        double *addPosRow(const MachineRow &r);
        BoxLcpStats solvePositions(std::vector<RowState> &state);
        void applyPseudo(const double *w, double mu);
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
        std::vector<double> M_, L_, Minv_, C_, qdotFree_, pqdot_, pqdotKept_, tmp_;

        std::vector<MachineRow> rows_;
        std::vector<double> J_, Y_, A_, diag_, b_, lo_, hi_, lambda_;
        std::vector<std::uint8_t> dead_;
        std::vector<double> W_;
        std::vector<int> freeIdx_;
        std::vector<double> Aff_;
        std::vector<int> bIdx_;       // bounded rows
        std::vector<double> bW_, bK_; // their W Jᵀ and effective inverse mass

        std::vector<MachineRow> posRows_;
        std::vector<double> posJ_, posY_, posA_, posDiag_, posB_, posLo_, posHi_, posLambda_;
        std::vector<std::uint8_t> posDead_;

        BoxLcpScratch lcp_;
        std::vector<double> jac_;

        void clearRows();
        void solveWithBase(const std::vector<double> &base, std::vector<RowState> &state, BoxLcpStats &st);
        void regularise(std::vector<double> &A, const std::vector<double> &diag, const std::vector<MachineRow> &rows,
                        const std::vector<double> &J, std::vector<std::uint8_t> &dead);
    };

} // namespace hf
