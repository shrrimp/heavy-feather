#pragma once

#include <limits>
#include <memory>
#include <span>

#include "hf/Body.hpp"
#include "hf/Contact.hpp"
#include "hf/Export.hpp"
#include "hf/Fault.hpp"
#include "hf/Joint.hpp"
#include "hf/Rope.hpp"
#include "hf/Types.hpp"

namespace hf {

    struct WorldSettings {
        Vec3 gravity{0.0, -9.81, 0.0};
        int substeps = 16;
        // Each island is solved on its own budget. An island holding a machine (its links, a mechanism contact or a
        // loose spring) sweeps its contacts this many times a sub-step: more converge heavy wheeled machines, whose
        // contacts and joints otherwise settle into a slow jitter.
        int contactIterations = 8;
        // An island of rigid bodies alone: its velocity sweeps, the relaxing sweeps after integrating, and the
        // position sweeps, per sub-step. Sub-steps do the converging there; a pile pays for each sweep in every row.
        int rigidIterations = 1, rigidRelaxIterations = 1, rigidPositionIterations = 2;
        // Sub-steps an island of rigid bodies alone is stepped in, apart from the machines' substeps; 0: as many.
        // Fewer make a pile cheaper (8 steps an awake 600-body pile in two thirds of the time 16 do) but let a tall
        // stack lean: twenty boxes stand at 16, lean at 10 and fall at 8.
        int rigidSubsteps = 0;
        double maxCorrection = 0.05; // m: the furthest a position correction may move a body in one sub-step
        bool sleeping = true;        // still bodies stop being stepped until something wakes them
        // m/s: a body whose fastest point moves faster than this at the end of a step is a fault, like a NaN.
        double faultSpeed = std::numeric_limits<double>::infinity();
        // Looks for a NaN after every phase of every sub-step, not only once the step is done, so a fault names the
        // phase it appeared in. For probes and debugging: it costs a pass over the stepped bodies per phase.
        bool checkPhases = false;
        bool profile = false; // times each phase of a step into WorldStats::profile
    };

    // Where a step's time went, in ms, with WorldSettings::profile on.
    struct StepProfile {
        double total = 0.0;
        double begin = 0.0;     // rebuilds, sleep and islands, drives, the contacts taken in
        double machines = 0.0;  // the machines' own solves: dynamics, joints, drives
        double contacts = 0.0;  // contact rows prepared and swept, with the machines' answers
        double positions = 0.0; // both position passes
        double integrate = 0.0;
        double relax = 0.0;
        double end = 0.0; // restitution, the fault check, warm starts, sleep, what moved
    };

    struct WorldStats {
        int bodies = 0, joints = 0, contacts = 0, ropes = 0;
        int machines = 0, links = 0, dofs = 0, rows = 0;
        int lcpPivots = 0, lcpFallbacks = 0;
        int stepped = 0, asleep = 0; // bodies
        int islands = 0;             // stepped bodies touching or jointed to each other, apart from the rest
        int faults = 0;              // islands put back this step: see World::faults
        int rejected = 0;            // calls refused for values the world can't use, since the step before this one
        // Since the step before this one: whole-world rebuilds of the bodies' groups, and machines built.
        int groupRebuilds = 0, machineRebuilds = 0;
        int rigidRows = 0, machineRows = 0; // swept contact rows in each kind of island, the last sub-step
        int ropeRows = 0, sweptRopes = 0;   // ropes held in a machine's solve, and swept between bodies, last step
        StepProfile profile;
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

        Vec3 position(BodyId body) const;     // of the origin
        Vec3 centerOfMass(BodyId body) const; // world
        Quat orientation(BodyId body) const;
        Vec3 velocity(BodyId body) const;
        Vec3 angularVelocity(BodyId body) const;
        BodyState state(BodyId body) const; // all four at once
        // The bodies whose pose the last step changed: what a host has to write back.
        std::span<const BodyId> moved() const;

        // Puts the body there, and its weld group with it; a machine link alone pulls its machine apart, see
        // moveMachine. A static body carries the machines anchored on it along. A kinematic body moves there over
        // the next step instead. What a pose or a velocity write sets survives the world rebuilding itself.
        bool setPose(BodyId body, const Vec3 &position, const Quat &orientation);
        // A kinematic body moves at this velocity over the next step.
        bool setVelocity(BodyId body, const Vec3 &velocity, const Vec3 &angularVelocity);
        // The body to this pose, and everything it is welded or jointed to with it, by the same rigid move: its
        // machine's joints and drives are left as they were, its velocities turn with it. False for a machine
        // anchored on something that didn't move with it the same way.
        bool moveMachine(BodyId body, const Vec3 &position, const Quat &orientation);
        // One rigid motion for the body's whole machine, about the body's centre of mass: a throw. False for a
        // static body, and for a machine anchored on something.
        bool setMachineVelocity(BodyId body, const Vec3 &velocity, const Vec3 &angularVelocity);
        // A kinematic body is moved by the host alone: each setPose is where it will be once the next step is done,
        // and it gets there over the step at the speed that takes. Nothing pushes it, and it pushes and carries what
        // it touches. Let go, it flies on at that speed. False for a machine link: a machine moves by moveMachine.
        bool setKinematic(BodyId body, bool kinematic);
        bool kinematic(BodyId body) const;
        // New mass and inertia (about the centre of mass, body axes), as an edit changes them; the body carries on.
        bool setMassProperties(BodyId body, double mass, const Mat3 &inertia);
        // And a centre of mass moved within the body, its origin where it was: its joints, its machine's drives and
        // whatever holds onto it carry on, nothing rebuilt for a body alone in its group.
        bool setMassProperties(BodyId body, double mass, const Mat3 &inertia, const Vec3 &centerOfMass);
        // A static body doesn't move and isn't pushed: one pinned in place.
        void setStatic(BodyId body, bool isStatic);
        // See BodyDesc::fixedRotation. Turned on, the body stops turning where it is; off, it carries on from rest.
        void setFixedRotation(BodyId body, bool fixed);
        // Sets the gravity scale for that body, only finite values accepted
        bool setGravityScale(BodyId body, double scale);
        // Sets the damping of that body, damps world movement of the body, negative and NaN rejected
        bool setDamping(BodyId body, double rate);
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
        double jointVelocity(JointId joint) const; // how fast that position changes, in rad/s or m/s
        // Puts a joint at a position, as the operator turning one part of a machine by hand: the part beyond it, away
        // from the machine's root, turns or slides there, its rates and its subtree's stopped, and a servo or piston
        // holds there. For a tree joint of one coordinate (JointInfo::tree): a hinge, motor, servo, angular spring,
        // slider or piston; clamped to the course of the last three. A loop through what moved closes in the
        // following steps.
        bool setJointPosition(JointId joint, double position);
        double driveForce(JointId joint) const; // over the last sub-step, N or N·m
        JointInfo jointInfo(JointId joint) const;

        // A grab of a body at a point in its own axes, from its origin: the mouse's hand on it. It keeps the
        // body awake until released. Refused: a cap or a frequency negative or NaN, a point not finite.
        GrabId grab(BodyId body, const Vec3 &localPoint, const GrabDesc &desc);
        // Where the grabbed point is pulled to; with an orientation, and a torque cap, where the body is turned to.
        bool setGrabTarget(GrabId grab, const Vec3 &point, const Quat *orientation = nullptr);
        void release(GrabId grab);

        // A rope along its route: its ends first and last, and between them whatever it runs round or slides through.
        // Refused: fewer than two points, a body that is gone, a point or axis not finite, an axis of no length at a
        // sheave, a radius, length, stiffness, strength, speed, force or damping negative or NaN, min over max.
        RopeId addRope(const RopeDesc &desc, std::span<const RopePoint> route);
        // How long a rope along this route would be, its bodies where they are now, sheaves included: what to pay out
        // for one that is just taut. 0 for a route the world can't use. Given room, where it would run too, one bend a
        // route point, as ropePath writes them: a host showing a rope between steps, its bodies moved by hand.
        double routeSpan(std::span<const RopePoint> route, std::span<RopeBend> path = {}) const;
        void removeRope(RopeId rope);
        bool contains(RopeId rope) const;
        // Its stiffness, strength and spools from the next step. A new length is paid out at once; a broken rope is
        // mended.
        bool setRope(RopeId rope, const RopeDesc &desc);
        RopeDesc ropeDesc(RopeId rope) const;
        // A new route for it, keeping what is paid out and how hard it pulls: where it catches on an edge, or comes off.
        bool setRopeRoute(RopeId rope, std::span<const RopePoint> route);
        // Pays it out, or takes it in, at once: the host's hand, as an editor stretching it.
        bool setRopeLength(RopeId rope, double length);
        // Works its spools like an operator's lever, from -1 to 1: positive pays the rope out at that fraction of the
        // spools' speed, negative winds it in. Back at 0 the brake holds it where it is.
        bool setSpoolLever(RopeId rope, double lever);
        double ropeLength(RopeId rope) const;  // m: paid out now
        double ropeSpan(RopeId rope) const;    // m: along its route as it lies, sheaves included
        double ropeTension(RopeId rope) const; // N, over the last sub-step
        RopeInfo ropeInfo(RopeId rope) const;
        // The rope through its route as the last step solved it, one bend a route point: the number written.
        int ropePath(RopeId rope, std::span<RopeBend> out) const;
        // The ropes the last step broke. They stay in the world, slack, until the host removes or mends them.
        std::span<const RopeId> brokenRopes() const;

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
