#pragma once

#include <memory>

#include "hf/Body.hpp"
#include "hf/Contact.hpp"
#include "hf/Export.hpp"
#include "hf/Joint.hpp"
#include "hf/Types.hpp"

namespace hf {

    struct WorldSettings {
        Vec3 gravity{0.0, -9.81, 0.0};
        int substeps = 16;
        double maxCorrection = 0.05; // m: the furthest a position correction may move a body in one sub-step
    };

    struct WorldStats {
        int bodies = 0, joints = 0, contacts = 0;
        int machines = 0, links = 0, dofs = 0, rows = 0;
        int lcpPivots = 0, lcpFallbacks = 0;
    };

    // Bodies, joints and the host's contacts in; poses out. Jointed bodies are grouped into machines, each solved
    // as one articulation.
    class HF_API World {
    public:
        explicit World(const WorldSettings &settings = {});
        ~World();
        World(World &&) noexcept;
        World &operator=(World &&) noexcept;

        const WorldSettings &settings() const;
        void setSettings(const WorldSettings &settings);

        BodyId addBody(const BodyDesc &desc);
        // Removes the body and every joint on it.
        void removeBody(BodyId body);
        bool contains(BodyId body) const;

        Vec3 position(BodyId body) const; // of the centre of mass
        Quat orientation(BodyId body) const;
        Vec3 velocity(BodyId body) const;
        Vec3 angularVelocity(BodyId body) const;
        void setPose(BodyId body, const Vec3 &position, const Quat &orientation);
        void setVelocity(BodyId body, const Vec3 &velocity, const Vec3 &angularVelocity);

        // A force at a world point, or a torque, for the next step only.
        void applyForce(BodyId body, const Vec3 &force, const Vec3 &point);
        void applyTorque(BodyId body, const Vec3 &torque);

        JointId addJoint(const JointDesc &desc);
        void removeJoint(JointId joint);
        bool contains(JointId joint) const;

        // A motor's angular velocity, or a piston's position.
        void setTarget(JointId joint, double target);
        double target(JointId joint) const;
        double jointPosition(JointId joint) const; // a hinge's angle in rad, a slider's position in m
        double driveForce(JointId joint) const;    // over the last sub-step, N or N·m
        JointInfo jointInfo(JointId joint) const;

        // The contacts for the next step. They are refreshed from their body-local points every sub-step, so the
        // host only detects once per step.
        void clearContacts();
        void addContact(const Contact &contact);

        void step(double dt);

        WorldStats stats() const;

        struct Impl;

    private:
        std::unique_ptr<Impl> impl_;
    };

} // namespace hf
