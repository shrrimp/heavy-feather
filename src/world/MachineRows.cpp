#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include "world/MachineJoints.hpp"
#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        constexpr double kInf = std::numeric_limits<double>::infinity();
        constexpr Vec3 kX{1.0, 0.0, 0.0}, kY{0.0, 1.0, 0.0}, kZ{0.0, 0.0, 1.0};

        constexpr double kJointBeta = 0.8;     // fraction of a position error corrected per sub-step
        constexpr double kMaxPosBias = 1000.0; // m/s or rad/s: a pose-only ceiling, it never becomes momentum
        constexpr double kLimitSlop = 0.05;    // m: a course limit's row exists within this of its stop
        constexpr double kGroupYield = 1e-3;   // a cut's directional yield, see Machine::regularise
        // A piston's position target is a servo command, not a pin: its row yields this much relative to its
        // diagonal, so a closure always wins against it and twin cylinders a fraction of a millimetre apart don't
        // fight.
        constexpr double kServoYield = 1e-2;
        constexpr double kServoRate = 0.01;     // m/s: how fast a piston's servo corrects drift near its target
        constexpr double kServoDeadband = 1e-3; // m: beyond this, it closes at its own speed instead
        constexpr double kDriveReached = 1e-4;  // m: a valved piston this close to its target shuts
        constexpr double kDriveSatSlack = 1e-9;
        constexpr double kDriveBackstep = 1e-6; // m: how far a load must push a saturated piston back to shut it

        // Row tags: which joint and which row, stable across steps, so the active set warm-starts.
        enum Slot : std::uint64_t {
            kLin0 = 0,
            kPerp0 = 5,
            kAng0 = 7,
            kDrive = 10,
            kLimLo = 11,
            kLimHi = 12,
            kSpring = 13,
            kPosBase = 32,
        };

        std::uint64_t tagOf(const JointSlot &j, std::uint32_t index, std::uint64_t slot) {
            return (std::uint64_t(j.generation) << 40) | (std::uint64_t(index) << 8) | slot;
        }

    } // namespace

    // Adds a machine's joint rows: velocity rows, position rows, or (for the line search) only their errors.
    struct RowBuilder {
        enum class Mode { Velocity, Position, Measure };

        World::Impl &I;
        MachineData &md;
        Mode mode;
        std::vector<double> *errors = nullptr;
        std::vector<double> scratch;

        bool position() const { return mode != Mode::Velocity; }

        double *add(const MachineRow &row, std::uint64_t tag, double C) {
            switch (mode) {
                case Mode::Velocity: md.velTags.push_back(tag); return md.m.addRow(row);
                case Mode::Position: md.posTags.push_back(tag); return md.m.addPosRow(row);
                case Mode::Measure: break;
            }
            errors->push_back(C);
            scratch.assign(std::size_t(md.m.dofs()), 0.0);
            return scratch.data();
        }

        Vec3 endPoint(const End &e, int link) const {
            if (link >= 0) {
                const Link &L = md.m.art.links[std::size_t(link)];
                return L.c + L.R * e.anchor;
            }
            const Group &g = I.groups[std::size_t(e.group)];
            return g.c + g.q * e.anchor;
        }

        Quat endFrame(const End &e, int link) const {
            return (link >= 0 ? md.m.art.links[std::size_t(link)].R : I.groups[std::size_t(e.group)].q) * e.frame;
        }

        // A slider's position: B's anchor along A's axis. Its rate row: J qdot = d(pos)/dt.
        double sliderPos(const JointRef &r) const {
            return dot(endPoint(r.B, r.linkB) - endPoint(r.A, r.linkA), endFrame(r.A, r.linkA) * kX);
        }

        void sliderRateRow(const JointRef &r, double scale, double *J) const {
            if (r.tree) {
                J[md.m.art.links[std::size_t(r.treeLink)].dof0] += scale * r.sign;
                return;
            }
            const Vec3 pB = endPoint(r.B, r.linkB);
            const Vec3 u = endFrame(r.A, r.linkA) * kX;
            md.m.art.addPointRow(r.linkB, pB, u, scale, J);
            md.m.art.addPointRow(r.linkA, pB, u, -scale, J);
        }

        // A hinge's angle, B relative to A about A's axis, and its rate row.
        double hingeAngle(const JointRef &r) const {
            const Quat qA = endFrame(r.A, r.linkA), qB = endFrame(r.B, r.linkB);
            const Vec3 yA = qA * kY, zA = qA * kZ, yB = qB * kY;
            return std::atan2(dot(yB, zA), dot(yB, yA));
        }

        void hingeRateRow(const JointRef &r, double scale, double *J) const {
            if (r.tree) {
                J[md.m.art.links[std::size_t(r.treeLink)].dof0] += scale * r.sign;
                return;
            }
            const Vec3 u = endFrame(r.A, r.linkA) * kX;
            md.m.art.addAngularRow(r.linkB, u, scale, J);
            md.m.art.addAngularRow(r.linkA, u, -scale, J);
        }

        void closureRows(const JointRef &r);
        void actuatorRows(std::size_t ri);
    };

    // A cut joint's equalities.
    void RowBuilder::closureRows(const JointRef &r) {
        Machine &m = md.m;
        const JointSlot &j = I.joints[r.joint];
        const double h = I.h;
        const Vec3 pA = endPoint(r.A, r.linkA), pB = endPoint(r.B, r.linkB);
        const Quat qA = endFrame(r.A, r.linkA), qB = endFrame(r.B, r.linkB);

        // One group per cut and row family, its pin rows and its axis rows, regularised as blocks.
        const auto add = [&](double C, std::uint64_t slot) -> double * {
            MachineRow row;
            const bool pinRow = slot < kPerp0 + 2;
            row.group = int((std::uint64_t(r.joint) << 1) | (pinRow ? 0u : 1u));
            row.groupYield = kGroupYield;
            if (position()) row.target = std::clamp(-(kJointBeta / h) * C, -kMaxPosBias, kMaxPosBias);
            return this->add(row, tagOf(j, r.joint, (position() ? std::uint64_t(kPosBase) : std::uint64_t(0)) + slot),
                             C);
        };

        const JointType t = j.desc.type;
        if (t == JointType::Ball || t == JointType::Hinge || t == JointType::Motor || t == JointType::AngularSpring) {
            // A cut hinge closes as a ball: its axis rows hold nothing on a planar loop and lock one that is slightly
            // out of plane, where the misaligned cut should show its tilt instead.
            for (int k = 0; k < 3; ++k) {
                Vec3 e{};
                e[k] = 1.0;
                double *J = add(dot(pA - pB, e), kLin0 + std::uint64_t(k));
                m.art.addPointRow(r.linkA, pA, e, 1.0, J);
                m.art.addPointRow(r.linkB, pB, e, -1.0, J);
            }
        }
        if (t == JointType::Slider || t == JointType::Piston) {
            for (int k = 0; k < 2; ++k) {
                const Vec3 tv = qA * (k == 0 ? kY : kZ);
                double *J = add(dot(pA - pB, tv), kPerp0 + std::uint64_t(k));
                // A's material point at B's anchor: sliding adds no lever arm.
                m.art.addPointRow(r.linkA, pB, tv, 1.0, J);
                m.art.addPointRow(r.linkB, pB, tv, -1.0, J);
            }
            Quat qe = qA * conjugate(qB);
            if (qe.w < 0.0) qe = Quat{-qe.w, -qe.x, -qe.y, -qe.z};
            const Vec3 Cang = 2.0 * Vec3{qe.x, qe.y, qe.z};
            for (int k = 0; k < 3; ++k) {
                Vec3 e{};
                e[k] = 1.0;
                double *J = add(Cang[k], kAng0 + std::uint64_t(k));
                m.art.addAngularRow(r.linkA, e, 1.0, J);
                m.art.addAngularRow(r.linkB, e, -1.0, J);
            }
        }
    }

    // Drives, course limits and springs, on a tree dof or a cut joint alike.
    void RowBuilder::actuatorRows(std::size_t ri) {
        Machine &m = md.m;
        const JointRef &r = md.joints[ri];
        JointSlot &j = I.joints[r.joint];
        const JointDesc &d = j.desc;
        const double h = I.h;
        if (mode == Mode::Velocity) md.driveRow[ri] = -1;

        if (d.type == JointType::Motor && mode == Mode::Velocity) {
            // The force cap, and a ceiling on how fast the force may grow: release stays free.
            MachineRow row;
            row.target = d.target;
            row.lo = -d.maxForce * h;
            row.hi = d.maxForce * h;
            if (std::isfinite(d.forceSlew) && d.forceSlew > 0.0) {
                const double ceiling = std::fabs(j.lambda) + d.forceSlew * h * h;
                row.lo = std::max(row.lo, -ceiling);
                row.hi = std::min(row.hi, ceiling);
            }
            row.relEps = 1e-7;
            md.driveRow[ri] = m.rowCount();
            hingeRateRow(r, 1.0, add(row, tagOf(j, r.joint, kDrive), 0.0));
        }

        if (d.type == JointType::Piston || d.type == JointType::Slider) {
            const double pos = sliderPos(r);
            j.slidePos = pos;
            bool driveUnbounded = false;
            if (d.type == JointType::Piston) {
                // The valve: force-capped while the pump drives it, rigid once shut, holding where it stalled if it
                // was mid-travel.
                const bool valved = std::isfinite(d.maxForce);
                const bool hold = valved && j.locked;
                if (hold && j.latchPending && mode == Mode::Velocity) {
                    j.holdPos = pos;
                    j.latchPending = false;
                }
                const double want = !hold ? j.targetEff : j.commanded ? j.holdPos : d.target;
                MachineRow row;
                if (!hold) {
                    row.lo = -d.maxForce * h;
                    row.hi = d.maxForce * h;
                }
                driveUnbounded = !std::isfinite(row.lo) && !std::isfinite(row.hi);
                row.relEps = 1e-7;
                if (position()) {
                    // A servo with a finite flow, not a teleport: inside the deadband it corrects at kServoRate, so
                    // twin cylinders only nudge each other; outside, at its own speed.
                    row.relEps = kServoYield;
                    const double e = pos - want;
                    const double cap = std::fabs(e) > kServoDeadband ? std::max(d.speed, kServoRate) : kServoRate;
                    row.target = std::clamp(-(kJointBeta / h) * e, -cap, cap);
                    sliderRateRow(r, 1.0, add(row, tagOf(j, r.joint, kPosBase + kDrive), e));
                } else {
                    row.target = hold ? 0.0 : j.targetRate;
                    md.driveRow[ri] = m.rowCount();
                    sliderRateRow(r, 1.0, add(row, tagOf(j, r.joint, kDrive), 0.0));
                }
            }
            // Course limits, within kLimitSlop of their stop. Left out while an unbounded drive owns the axis: its
            // target is held to the course, so the two would be the same row.
            if (!driveUnbounded) {
                const double gap[2] = {pos - d.min, d.max - pos};
                for (int i = 0; i < 2; ++i) {
                    if (gap[i] >= kLimitSlop) continue;
                    MachineRow row;
                    row.lo = 0.0;
                    row.hi = kInf;
                    const double dir = i == 0 ? 1.0 : -1.0;
                    const std::uint64_t slot = i == 0 ? kLimLo : kLimHi;
                    if (position()) {
                        if (gap[i] >= 0.0) continue; // nothing to correct until it is past
                        row.target = std::min(-(kJointBeta / h) * gap[i], kMaxPosBias);
                        sliderRateRow(r, dir, add(row, tagOf(j, r.joint, kPosBase + slot), gap[i]));
                    } else {
                        row.target = gap[i] > 0.0 ? -gap[i] / h : 0.0; // arrive on the stop, no faster
                        sliderRateRow(r, dir, add(row, tagOf(j, r.joint, slot), 0.0));
                    }
                }
            }
        }

        if (position()) return;

        if (d.type == JointType::AngularSpring) {
            std::vector<double> J(std::size_t(m.dofs()), 0.0);
            hingeRateRow(r, 1.0, J.data());
            MachineRow row;
            springRow(m.inverseMassAlong(J.data()), d.hertz, d.damping, hingeAngle(r), h, row.compliance, row.target);
            std::copy(J.begin(), J.end(), add(row, tagOf(j, r.joint, kSpring), 0.0));
        }

        if (d.type == JointType::Spring) {
            const Vec3 pA = endPoint(r.A, r.linkA), pB = endPoint(r.B, r.linkB);
            const double len = length(pA - pB);
            if (len > 1e-9) {
                const Vec3 n = (pA - pB) / len;
                std::vector<double> J(std::size_t(m.dofs()), 0.0);
                m.art.addPointRow(r.linkA, pA, n, 1.0, J.data());
                m.art.addPointRow(r.linkB, pB, n, -1.0, J.data());
                MachineRow row;
                springRow(m.inverseMassAlong(J.data()), d.hertz, d.damping, len - d.restLength, h, row.compliance,
                          row.target);
                std::copy(J.begin(), J.end(), add(row, tagOf(j, r.joint, kSpring), 0.0));
            }
        }
    }

    namespace {

        void warmRows(const World::Impl &I, std::vector<RowState> &state, const std::vector<std::uint64_t> &tags) {
            state.resize(tags.size());
            for (std::size_t i = 0; i < tags.size(); ++i) {
                const auto it = I.rowState.find(tags[i]);
                state[i] = it == I.rowState.end() ? RowState::Free : it->second;
            }
        }

        void keepRows(World::Impl &I, const std::vector<RowState> &state, const std::vector<std::uint64_t> &tags) {
            for (std::size_t i = 0; i < tags.size() && i < state.size(); ++i) I.rowState[tags[i]] = state[i];
        }

        void recordDrives(World::Impl &I, MachineData &md) {
            const std::vector<double> &lambda = md.m.lambda();
            for (std::size_t ri = 0; ri < md.joints.size(); ++ri) {
                const int row = md.driveRow[ri];
                if (row < 0 || std::size_t(row) >= lambda.size()) continue;
                JointSlot &j = I.joints[md.joints[ri].joint];
                j.lambda = lambda[std::size_t(row)];
                if (j.desc.type == JointType::Piston)
                    j.driveSat = j.desc.maxForce > 0.0 ? j.lambda / (j.desc.maxForce * I.h) : 0.0;
            }
        }

        void addJointRows(RowBuilder &rb) {
            for (std::size_t ri = 0; ri < rb.md.joints.size(); ++ri) {
                const JointRef &r = rb.md.joints[ri];
                if (!r.tree && !r.duplicate && kinematic(rb.I.joints[r.joint].desc.type)) rb.closureRows(r);
                rb.actuatorRows(ri);
            }
        }

        // The error the line search scores: Σ C² over the position rows, from the current poses.
        double positionResidual(World::Impl &I, MachineData &md) {
            std::vector<double> errors;
            RowBuilder rb{I, md, RowBuilder::Mode::Measure, &errors, {}};
            addJointRows(rb);
            double s = 0.0;
            for (double c : errors) s += c * c;
            return s;
        }

        struct Saved {
            std::vector<double> q;
            std::vector<Quat> qRot;
            Vec3 c0{};
            Quat R0{};
        };

        void save(const Machine &m, Saved &st) {
            st.q = m.art.q;
            st.qRot = m.art.qRot;
            st.c0 = m.art.links[0].c;
            st.R0 = m.art.links[0].R;
        }

        void restore(Machine &m, const Saved &st) {
            m.art.q = st.q;
            m.art.qRot = st.qRot;
            m.art.links[0].c = st.c0;
            m.art.links[0].R = st.R0;
            m.art.forwardKinematics();
        }

    } // namespace

    void World::Impl::solveInternal(MachineData &md) {
        Machine &m = md.m;
        const Vec3 O = m.art.links[0].c;
        for (std::size_t l = 0; l < md.groupOfLink.size(); ++l) {
            const Group &g = groups[std::size_t(md.groupOfLink[l])];
            md.fext[l] = Force{cross(m.art.links[l].c - O, g.force) + g.torque, g.force};
        }
        m.dynamics(h, settings.gravity, md.fext);

        md.velTags.clear();
        md.driveRow.assign(md.joints.size(), -1);
        RowBuilder rb{*this, md, RowBuilder::Mode::Velocity, nullptr, {}};
        addJointRows(rb);
        m.prepareRows();
        warmRows(*this, md.velState, md.velTags);
        m.solveRows(md.velState);
        keepRows(*this, md.velState, md.velTags);
        m.buildMobility();
        m.prepareBounded();
        recordDrives(*this, md);
        writeVelocities(md);
    }

    void World::Impl::solvePositions(MachineData &md) {
        md.m.clearPosRows();
        md.posTags.clear();
        RowBuilder rb{*this, md, RowBuilder::Mode::Position, nullptr, {}};
        addJointRows(rb);
        warmRows(*this, md.posState, md.posTags);
        md.m.solvePositions(md.posState);
        keepRows(*this, md.posState, md.posTags);
    }

    void World::Impl::afterContacts(MachineData &md) {
        recordDrives(*this, md);
        writeVelocities(md);
    }

    void World::Impl::integrate(MachineData &md) {
        Machine &m = md.m;
        // The ceiling on one sub-step's position correction, over the whole machine: scaled, never zeroed, so a
        // correction that is merely too fast arrives late rather than not at all.
        double worst = 0.0;
        for (std::size_t l = 0; l < md.groupOfLink.size(); ++l) {
            const Motion V = m.art.twistOf(int(l), m.pseudo().data());
            const Vec3 vc = pointVelocity(V, m.art.links[l].c - m.art.O);
            const double radius = groups[std::size_t(md.groupOfLink[l])].radius;
            worst = std::max(worst, (length(vc) + length(V.w) * radius) * h);
        }
        m.keepPseudo();
        if (!std::isfinite(worst)) m.scalePseudo(0.0);
        else if (worst > settings.maxCorrection) m.scalePseudo(settings.maxCorrection / worst);
        m.keepPseudo();

        // The safeguard: a correction is kept only if it lowers the position error against the plain dynamics
        // step, otherwise half of it, otherwise none. A set of position rows that can't all be met would otherwise
        // ask for the same correction every sub-step and the machine would creep.
        bool anyPseudo = false;
        for (double x : m.pseudo()) anyPseudo |= (x != 0.0);
        if (anyPseudo) {
            Saved st;
            save(m, st);
            m.integrate(h);
            const double rFull = positionResidual(*this, md);
            restore(m, st);
            m.scalePseudo(0.0);
            m.integrate(h);
            const double rNone = positionResidual(*this, md);
            restore(m, st);
            if (rFull <= rNone) {
                m.restorePseudo();
            } else {
                m.restorePseudo(0.5);
                m.integrate(h);
                const double rHalf = positionResidual(*this, md);
                restore(m, st);
                if (rHalf < rNone) m.restorePseudo(0.5);
                else m.scalePseudo(0.0);
            }
        }
        m.integrate(h);
        writePoses(md);
    }

    // Once per step, so the row set can't change under the solver mid-step. A valved piston (one with a force cap)
    // shuts once it reaches its target, or when it is saturated and the load is pushing it back: a check valve shuts
    // against reverse flow, not against full pressure. The setpoint follows the target at the drive's speed, as
    // commanded motion the velocity drive carries; whatever the target jumped beyond that is reached over later steps.
    void World::Impl::beginDrives(double dt) {
        for (JointSlot &j : joints) {
            if (!j.alive) continue;
            const JointDesc &d = j.desc;
            if (d.type == JointType::Piston && std::isfinite(d.maxForce)) {
                const double moved = j.slidePos - j.slidePrev;
                j.slidePrev = j.slidePos;
                if (std::fabs(d.target - j.slidePos) <= kDriveReached) j.commanded = false;
                const bool wasLocked = j.locked;
                if (!j.commanded) j.locked = true;
                else if (d.target > j.slidePos)
                    j.locked = j.driveSat >= 1.0 - kDriveSatSlack && moved < -kDriveBackstep;
                else j.locked = j.driveSat <= -1.0 + kDriveSatSlack && moved > kDriveBackstep;
                if (j.locked && !wasLocked) j.latchPending = true;
            }
            const double gap = d.target - j.targetFrom;
            const double move = d.speed > 0.0 ? std::clamp(gap, -d.speed * dt, d.speed * dt) : gap;
            j.targetTo = j.targetFrom + move;
            j.targetRate = d.speed > 0.0 ? move / dt : 0.0;
        }
    }

} // namespace hf
