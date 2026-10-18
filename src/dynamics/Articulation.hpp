#pragma once

#include <cstdint>
#include <vector>

#include "hf/math/Mat3.hpp"
#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"
#include "math/Spatial.hpp"

namespace hf {

    enum class JointKind : std::uint8_t {
        Free,      // floating base, the root only: 6 dof
        Revolute,  // about the parent frame's X axis
        Prismatic, // along the parent frame's X axis
        Spherical, // 3 rotations about the anchor
    };

    struct Link {
        int parent = -1; // -1: the world anchor, or nothing for a Free root
        JointKind kind = JointKind::Free;

        // Joint frames: the parent side in the parent body's frame (centre-of-mass relative), the child side in
        // this body's. With no parent, the parent side is in the anchor frame at anchorPos / anchorRot.
        Vec3 pp{}, pc{};
        Quat Qp{}, Qc{};
        Vec3 anchorPos{};
        Quat anchorRot{};

        double mass = 0.0;
        Mat3 inertiaBody{}; // about the centre of mass, body axes

        int dof0 = 0; // first generalised velocity
        int ndof = 0;
        int coord = -1; // into q (revolute, prismatic) or qRot (spherical)

        Vec3 c{}; // world centre of mass
        Quat R{}; // world orientation

        // Per sub-step, about O.
        Motion S[6]; // motion subspace
        Motion V;    // twist
        Motion A;    // bias acceleration, at qddot = 0
        Inertia I;
        Inertia Ic; // composite: this link and everything below it
        Force F;    // RNEA scratch
    };

    // One machine in reduced coordinates: a tree of links, parents before children. Tree joints can't drift
    // apart, because a child's pose is computed from its parent's.
    class Articulation {
    public:
        std::vector<Link> links;
        std::vector<double> q;    // revolute angles, prismatic slides
        std::vector<Quat> qRot;   // spherical joints' relative rotations
        std::vector<double> qdot; // spherical: relative ω in the parent joint frame; free: the root twist about O
        Vec3 O{};                 // the point every spatial vector is about

        int dofs() const { return ndof_; }

        // Lays out dof0 and coord once the links are set. Only links[0] may be Free.
        void finalise();

        // Coordinates and velocities from every link's world pose and velocity. A pose that breaks a tree joint
        // is projected onto it.
        void setFromMaximal(const std::vector<Vec3> &c, const std::vector<Quat> &R, const std::vector<Vec3> &v,
                            const std::vector<Vec3> &w);

        // Link poses from q. A floating root's pose is state, not derived.
        void forwardKinematics();

        // Moves O, then fills every link's S, V, A and inertias. Needs current poses.
        void updateKinematics(const Vec3 &newO);
        // V and A again for other joint velocities, the configuration unchanged. Needs updateKinematics.
        void updateVelocities(const double *qd);

        // The joint-space mass matrix, n×n row-major (CRBA). Needs updateKinematics.
        void massMatrix(std::vector<double> &M) const;

        // C in M qddot + C = tau + Jᵀλ, with gravity g and per-link wrenches about O (RNEA). fext may be empty.
        void biasForces(const Vec3 &g, const std::vector<Force> &fext, std::vector<double> &C);

        // M's diagonal alone, from the composite inertias. Needs updateKinematics.
        void massDiagonal(std::vector<double> &d) const;

        // M factored through articulated-body inertias (ABA), then x <- M⁻¹ x in O(n) for as many right-hand sides
        // as needed. Needs updateKinematics.
        void factorMass();
        void solveMass(double *x) const;
        // M⁻¹ Jᵀ f for up to three wrenches f at one link (about O): out holds count rows of dofs() values. The
        // inward pass only climbs from the link to the root, so this costs one outward pass over the machine.
        void solveAt(int link, const Force *f, int count, double *out) const;

        // Adds scale × (the row mapping qdot to dir · velocity of the point xWorld on link) into row.
        void addPointRow(int link, const Vec3 &xWorld, const Vec3 &dir, double scale, double *row) const;

        // Adds scale × (the row mapping qdot to axis · angular velocity of link) into row.
        void addAngularRow(int link, const Vec3 &axis, double scale, double *row) const;

        Motion twistOf(int link, const double *qd) const;
        // Every link's twist in one pass, parents first.
        void twists(const double *qd, std::vector<Motion> &V) const;
        Vec3 linkVelocity(int link) const; // of its centre of mass
        Vec3 linkOmega(int link) const;

        // Advances q by h along qd (qdot, or qdot plus a position correction), then poses.
        void integrate(double h, const std::vector<double> &qd);

        // ½ qdotᵀ M qdot, from the link twists. Needs updateKinematics.
        double kineticEnergy() const;

    private:
        int ndof_ = 0;

        // One link's share of the factorisation: U = IA S and D⁻¹ = (Sᵀ IA S)⁻¹, k columns for k dofs.
        struct Factor {
            ArticulatedInertia IA;
            Force U[6];
            double Dinv[6][6];
        };
        std::vector<Factor> factor_;
        mutable std::vector<Force> bias_;
        mutable std::vector<Motion> acc_;
        mutable std::vector<double> u_;
        mutable std::vector<std::uint32_t> onPath_;
        mutable std::uint32_t pathStamp_ = 0;

        void parentPose(int link, Vec3 &c, Quat &R) const;
        Vec3 parentAnchor(int link) const;
        Quat parentFrame(int link) const;
    };

} // namespace hf
