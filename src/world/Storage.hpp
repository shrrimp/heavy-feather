#pragma once

#include <cstdint>
#include <numbers>

#include "hf/Body.hpp"
#include "hf/Joint.hpp"
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
        bool isStatic = true;

        Vec3 force{}, torque{}; // for this step, about c
        Vec3 pv{}, pw{};        // pseudo velocities: the position pass, one sub-step only
        int machine = -1, link = -1;
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
    // row: J v + compliance λ = target, with C the row's position error.
    inline void springRow(double k, double hertz, double zeta, double C, double h, double &compliance, double &target) {
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
