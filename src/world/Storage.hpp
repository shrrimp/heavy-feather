#pragma once

#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

#include "contact/ContactParams.hpp"
#include "hf/Body.hpp"
#include "hf/Joint.hpp"
#include "hf/Rope.hpp"
#include "hf/math/Mat3.hpp"
#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

namespace hf {

    struct BodySlot {
        BodyDesc desc;
        std::uint32_t generation = 0;
        bool alive = false;

        int group = -1;
        Vec3 offset{}; // its centre of mass in the group's frame
        Quat rot{};    // its orientation in the group's frame

        Vec3 force{}, torque{}; // for the next step, about its centre of mass

        // Asleep, it isn't stepped until something wakes it. Carried across steps, and across machine rebuilds.
        bool asleep = false, wakeNext = false;
        int stillSteps = 0;
        bool reshaped = false; // until the next step: whatever its contacts then touch wakes
        int faultStreak = 0;   // consecutive steps its island has been put back

        // Until the end of the next step: the host wrote its pose or its velocity, which a rebuild keeps rather than
        // reading the body back from its group; it was moved by the host, so its contacts and rows start cold.
        bool poseWritten = false, velocityWritten = false, teleported = false;
        bool listed = false; // in World::Impl::written
        bool forced = false; // in World::Impl::forced
        bool kinematic = false;
        bool hasTarget = false; // kinematic: where the host wants it once the next step is done
        Vec3 targetPosition{};
        Quat targetOrientation{};
    };

    struct JointSlot {
        JointDesc desc;
        std::uint32_t generation = 0;
        bool alive = false;
        JointInfo info;

        // Drive state, carried across steps. A position drive's setpoint follows its target at the drive's speed: it
        // moves from targetFrom to targetTo over a step, targetEff along the way.
        double targetFrom = 0.0, targetTo = 0.0, targetRate = 0.0, targetEff = 0.0;
        double lambda = 0.0; // the drive row's last impulse
        // Where the drive was as this step and the last began, measured from the poses.
        double position = 0.0, positionPrev = 0.0;
        double holdPos = 0.0, driveSat = 0.0;
        double lever = 0.0, leverPrev = 0.0;
        bool locked = false, latchPending = true, commanded = true;
    };

    struct GrabSlot {
        GrabDesc desc;
        BodyId body;
        Vec3 local{}; // from the body's centre of mass, in its axes: the host's point less its centre's
        Vec3 target{};
        Quat targetRot{};
        bool turns = false; // pulled toward targetRot too
        std::uint32_t generation = 0;
        bool alive = false;
    };

    // A rope at one point of its route, this step: the group it is on, and where on it, from the group's centre in its
    // frame. Each sub-step: where that is in the world, and the rope there, carried on to the next sub-step. A sheave's
    // tangent points start from where they were, and its angle stays continuous: unwound past nothing, it is off.
    struct RopeBendState {
        int group = 0;
        Vec3 local{}, axis{};
        Vec3 p{}, n{}; // world
        Vec3 in{}, out{};
        Vec3 dir{}; // how the span grows as the point moves: u_in - u_out
        double theta = 0.0;
        bool touching = true;
    };

    struct RopeSlot {
        RopeDesc desc;
        std::vector<RopePoint> route; // every sheave's side settled
        std::vector<RopeBendState> bends;
        std::uint32_t generation = 0;
        bool alive = false;
        bool placed = false; // its sheaves' angles found once, from the route as it lies

        double length = 0.0;     // m paid out
        double lever = 0.0;      // its spools'
        double lengthRate = 0.0; // m/s the spools pay out this step
        double cap = std::numeric_limits<double>::infinity(); // N its row may pull this step
        double force = 0.0;      // N, the last sub-step's tension: the next one's warm start
        double peak = 0.0;       // N, the most this step
        double span = 0.0;       // m along its route
        bool slipping = false;   // at its cap last sub-step: paying out under the load
        bool broken = false;

        // This step: the machine whose solve holds it, -1 swept, -2 not stepped; its row there; its island. And how
        // it was solved the last step it was stepped.
        int machine = -2, row = -1, island = -1, held = -1;
        double length0 = 0.0; // as the step began: what a fault puts back
    };

    // A rope's spools at work: paying out or winding in this step.
    inline bool reeling(const RopeSlot &r) { return r.lengthRate != 0.0; }

    // A piston or a servo: a drive toward a position, with a valve.
    inline bool positionDrive(JointType t) { return t == JointType::Piston || t == JointType::Servo; }
    // m or rad: within this of its target a position drive only nudges itself back; beyond, it travels at its speed.
    constexpr double kServoDeadband = 1e-3;

    // What the solver moves: one body, or a weld group moving as one. Group 0 is the static world.
    struct Group {
        Vec3 c{}; // centre of mass
        Quat q{};
        Vec3 c0{}; // c and q as the step began
        Quat q0{};
        Vec3 turned{}; // the rotation since then, as a rotation vector: refreshed with the contacts' depths
        Vec3 v{}, w{};
        double invMass = 0.0;
        Mat3 inertiaBody{}, invInertiaBody{};
        double radius = 0.0;

        double gravityScale = 1.0;
        double damping = 0.0;

        bool isStatic = true;
        bool fixedRotation = false; // a member's BodyDesc::fixedRotation: invInertiaBody is zero, w stays zero
        // Moved by the host, static to the solver: its velocity is what takes it to its target this step.
        bool kinematic = false;
        Vec3 targetC{};
        Quat targetQ{};

        Vec3 force{}, torque{}; // for this step, about c
        Vec3 pv{}, pw{};        // pseudo velocities: the position pass, one sub-step only
        int machine = -1, link = -1;
        int members = 0;         // its bodies; none, past group 0, and the slot is free for the next body added
        bool frozen = false;     // not stepped this step: asleep, and nothing awake reached it
        int island = -1;         // stepped: which island, the groups touching or jointed to it; -1 otherwise
        double correction = 0.0; // how far the position pass moved it in the last sub-step

        Mat3 invInertiaWorld() const {
            const Mat3 R = toMat3(q);
            return R * invInertiaBody * transpose(R);
        }
    };

    // Moved by the solver this step: neither static nor asleep. And the machine it moves with, if any.
    inline bool moves(const Group &g) { return !g.isStatic && !g.frozen; }
    inline int machineOf(const Group &g) { return g.frozen ? -1 : g.machine; }

    // A row's spring, from a frequency and a damping ratio.
    struct SoftParams {
        double biasRate = 0.0;
        double massCoeff = 1.0;
        double impulseCoeff = 0.0;

        static SoftParams rigid() { return {}; }
        static SoftParams fromHertz(double hertz, double zeta, double h);
    };

    // A spring of frequency hertz and damping ratio zeta on a row of effective inverse mass k, as an implicit soft
    // row: J v + compliance λ = target, with C the row's position error. At 0 Hz it is a rod: as stiff as a contact.
    inline void springRow(double k, double hertz, double zeta, double C, double h, double &compliance, double &target) {
        if (!(hertz > 0.0)) {
            hertz = kContactHertz;
            zeta = kContactDamping;
        }
        const double mEff = k > 1e-12 ? 1.0 / k : 0.0;
        const double w = 2.0 * std::numbers::pi * hertz;
        const double ks = mEff * w * w;
        const double c = 2.0 * mEff * zeta * w;
        const double den = c + h * ks;
        if (!(den > 0.0)) {
            compliance = 0.0;
            target = 0.0;
            return;
        }
        compliance = 1.0 / (h * den);
        target = -(ks / den) * C;
    }

} // namespace hf
