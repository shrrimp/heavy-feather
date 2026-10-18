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
    };

    struct JointSlot {
        JointDesc desc;
        std::uint32_t generation = 0;
        bool alive = false;
        JointInfo info;

        // Drive state, carried across steps. A piston's setpoint follows its target at the drive's speed: it moves
        // from targetFrom to targetTo over a step, targetEff along the way.
        double targetFrom = 0.0, targetTo = 0.0, targetRate = 0.0, targetEff = 0.0;
        double lambda = 0.0; // the drive row's last impulse
        double slidePos = 0.0, slidePrev = 0.0, holdPos = 0.0, driveSat = 0.0;
        bool locked = false, latchPending = true, commanded = true;
    };

    // What the solver moves: one body, or a weld group moving as one. Group 0 is the static world.
    struct Group {
        Vec3 c{}; // centre of mass
        Quat q{};
        Vec3 c0{}; // c and q as the step began
        Quat q0{};
        Vec3 v{}, w{};
        double invMass = 0.0;
        Mat3 invInertiaBody{};
        double radius = 0.0;
        bool isStatic = true;

        Vec3 force{}, torque{}; // for this step, about c
        Vec3 pv{}, pw{};        // pseudo velocities: the position pass, one sub-step only
        int machine = -1, link = -1;

        Mat3 invInertiaWorld() const {
            const Mat3 R = toMat3(q);
            return R * invInertiaBody * transpose(R);
        }
    };

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
