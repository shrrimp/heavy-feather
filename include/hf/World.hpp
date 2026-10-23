#pragma once

#include <limits>
#include <memory>
#include <span>

#include "hf/Body.hpp"
#include "hf/Contact.hpp"
#include "hf/Export.hpp"
#include "hf/Fault.hpp"
#include "hf/Joint.hpp"
#include "hf/Types.hpp"

namespace hf {

    struct WorldSettings {
        Vec3 gravity{0.0, -9.81, 0.0};
        int substeps = 16;
        // Contact sweeps per sub-step, at least. More converge heavy wheeled machines, whose contacts and joints
        // otherwise settle into a slow jitter, at a cost in large piles.
        int contactIterations = 8;
        double maxCorrection = 0.05; // m: the furthest a position correction may move a body in one sub-step
        bool sleeping = true;        // still bodies stop being stepped until something wakes them
        // m/s: a body whose fastest point moves faster than this at the end of a step is a fault, like a NaN.
        double faultSpeed = std::numeric_limits<double>::infinity();
        // Looks for a NaN after every phase of every sub-step, not only once the step is done, so a fault names the
        // phase it appeared in. For probes and debugging: it costs a pass over the stepped bodies per phase.
        bool checkPhases = false;
    };

    struct WorldStats {
        int bodies = 0, joints = 0, contacts = 0;
        int machines = 0, links = 0, dofs = 0, rows = 0;
        int lcpPivots = 0, lcpFallbacks = 0;
        int stepped = 0, asleep = 0; // bodies
        int islands = 0;             // stepped bodies touching or jointed to each other, apart from the rest
        int faults = 0;              // islands put back this step: see World::faults
        int rejected = 0;            // calls refused for values the world can't use, since the step before this one
    };

    // Bodies, joints and the host's contacts in; poses out. Jointed bodies are grouped into machines, each solved
    // as one articulation.
    //
    // A call given a value the world can't use, a NaN, a negative mass, a normal of no length, refuses it: it returns
    // false (or an invalid handle), changes nothing, and counts in WorldStats::rejected. A call on a handle that is
    // gone returns false too, but isn't counted.
    class HF_API World {
    public:
        explicit World(const WorldSettings &settings = {});
        ~World();
        World(World &&) noexcept;
        World &operator=(World &&) noexcept;

        const WorldSettings &settings() const;
        void setSettings(const WorldSettings &settings);

        // Refused: a dynamic body without a positive mass, or whose inertia isn't symmetric positive definite.
        BodyId addBody(const BodyDesc &desc);
        // Removes the body and every joint on it.
        void removeBody(BodyId body);
        bool contains(BodyId body) const;

        Vec3 position(BodyId body) const; // of the centre of mass
        Quat orientation(BodyId body) const;
        Vec3 velocity(BodyId body) const;
        Vec3 angularVelocity(BodyId body) const;
        bool setPose(BodyId body, const Vec3 &position, const Quat &orientation);
        bool setVelocity(BodyId body, const Vec3 &velocity, const Vec3 &angularVelocity);
        // New mass and inertia (about the centre of mass, body axes), as an edit changes them; the body carries on.
        bool setMassProperties(BodyId body, double mass, const Mat3 &inertia);
        // A static body doesn't move and isn't pushed: one pinned in place.
        void setStatic(BodyId body, bool isStatic);
        // A body still for long enough falls asleep, and isn't stepped until a push, an edit or an awake neighbour
        // wakes it. Its contacts still matter: they tell the world what it rests on, at no cost to the solve.
        bool asleep(BodyId body) const;
        void wake(BodyId body);
        // The body's shape changed: what rested on it, and what touches it in the next step's contacts, wakes.
        void reshaped(BodyId body);

        // A force at a world point, or a torque, for the next step only.
        bool applyForce(BodyId body, const Vec3 &force, const Vec3 &point);
        bool applyTorque(BodyId body, const Vec3 &torque);

        // Refused: a frame not finite, a force, speed, frequency, damping or length negative or NaN, min over max.
        JointId addJoint(const JointDesc &desc);
        void removeJoint(JointId joint);
        bool contains(JointId joint) const;
        // Changes a joint in place. Drives, courses and springs apply from the next step; a new type, body or frame
        // rebuilds its machine and restarts its drive. False if the joint is gone or the new bodies are invalid.
        bool setJoint(JointId joint, const JointDesc &desc);
        JointDesc jointDesc(JointId joint) const;

        // A motor's angular velocity, a piston's position, or a servo's angle.
        bool setTarget(JointId joint, double target);
        // Works a drive like an operator's lever, from -1 to 1: a piston or servo moves at that fraction of its
        // speed, and back at 0 stops where it actually is and holds there; a motor turns at that fraction of its speed.
        bool setLever(JointId joint, double lever);
        double target(JointId joint) const;
        double jointPosition(JointId joint) const; // a hinge's or servo's angle in rad, a slider's position in m
        double driveForce(JointId joint) const;    // over the last sub-step, N or N·m
        JointInfo jointInfo(JointId joint) const;

        // The contacts for the next step. They are refreshed from their body-local points every sub-step, so the
        // host only detects once per step.
        void clearContacts();
        // Refused: a point, depth or normal not finite, a normal far from unit length, a negative friction or lever
        // arm. Restitution is clamped to [0, 1].
        bool addContact(const Contact &contact);

        void step(double dt);

        WorldStats stats() const;
        // The last step's faults, one per island put back.
        std::span<const Fault> faults() const;

        struct Impl;

    private:
        std::unique_ptr<Impl> impl_;
    };

} // namespace hf
