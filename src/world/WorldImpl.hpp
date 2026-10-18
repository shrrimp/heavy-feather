#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "hf/World.hpp"
#include "solver/Machine.hpp"
#include "world/Storage.hpp"

namespace hf {

    // A joint end resolved onto the group it acts on.
    struct End {
        int group = 0;
        bool isStatic = true;
        Vec3 anchor{}; // in the group's frame, from its centre of mass
        Quat frame{};  // in the group's frame
    };

    struct JointRef {
        std::uint32_t joint = 0;
        End A, B;
        int linkA = -1, linkB = -1; // -1: a static end
        bool tree = false;
        bool duplicate = false;
        int treeLink = -1; // the child link of a tree joint
        double sign = 1.0; // the joint's own coordinate = sign × the articulation's
    };

    struct MachineData {
        Machine m;
        std::vector<int> groupOfLink;
        std::vector<JointRef> joints;
        std::vector<Force> fext;
        std::vector<RowState> velState, posState;
        std::vector<std::uint64_t> velTags, posTags;
        std::vector<int> driveRow; // per joint ref: its drive row, or -1
    };

    struct ContactData {
        Contact in;
        int ga = 0, gb = 0; // groups
        Vec3 rA{}, rB{};    // the point from each group's centre, as the step began
        double depth0 = 0.0;
        Vec3 n{};
        double depth = 0.0; // refreshed every sub-step
        double ln = 0.0, lt1 = 0.0, lt2 = 0.0;
        double lnMax = 0.0, vnImpact = 0.0;
    };

    // A contact end on a machine: its rows through the machine's mobility, n doubles each, in the arena as
    // gn, gt1, gt2, then W gn, W gt1, W gt2.
    struct ContactBlock {
        int machine = -1;
        std::size_t off = 0;
        int n = 0;
    };

    struct ContactPrep {
        int contact = 0;
        int a = 0, b = 0; // groups
        Mat3 invIwA{}, invIwB{};
        Vec3 r1{}, r2{}, t1{}, t2{};
        Vec3 r1cn{}, r2cn{}, r1ct1{}, r1ct2{}, r2ct1{}, r2ct2{};
        Vec3 imR1cn{}, imR2cn{}, imR1ct1{}, imR1ct2{}, imR2ct1{}, imR2ct2{};
        double imA = 0.0, imB = 0.0;
        double kn = 0.0, kt1 = 0.0, kt2 = 0.0;
        double kt12 = 0.0; // t1 W t2, for the exact 2×2 friction solve on a machine end
        SoftParams soft;
        double vtarget = 0.0, refLn = 0.0, posBias = 0.0, pseudo = 0.0;
        int nb = 0;
        ContactBlock blocks[2];
    };

    // A spring between bodies that aren't in one machine: a soft row swept with the contacts. Its velocity is
    // vA - vB along n at the anchors; a machine end keeps g then W g in the arena.
    struct SpringPrep {
        int a = 0, b = 0;
        Vec3 n{}, r1{}, r2{}, imR1n{}, imR2n{};
        double imA = 0.0, imB = 0.0;
        double k = 0.0, gamma = 0.0, target = 0.0, lambda = 0.0;
        int nb = 0;
        ContactBlock blocks[2];
    };

    struct ContactKey {
        BodyId a, b;
        std::uint64_t feature = 0;
        bool operator==(const ContactKey &) const = default;
    };

    struct ContactKeyHash {
        std::size_t operator()(const ContactKey &k) const {
            std::uint64_t h = k.feature * 0x9E3779B97F4A7C15ull;
            for (std::uint32_t x : {k.a.index, k.a.generation, k.b.index, k.b.generation})
                h = (h ^ x) * 0x100000001B3ull;
            return std::size_t(h);
        }
    };

    struct HF_NO_EXPORT World::Impl {
        WorldSettings settings;
        std::vector<BodySlot> bodies;
        std::vector<JointSlot> joints;
        std::vector<std::uint32_t> freeBodies, freeJoints;
        std::vector<Group> groups; // 0: the world
        std::vector<MachineData> machines;
        std::vector<JointRef> looseSprings; // springs not inside one machine, solved in the contact sweep
        bool dirty = true;
        double h = 0.0;

        std::vector<ContactData> contacts;
        std::vector<ContactPrep> prep;
        std::vector<SpringPrep> springs;
        std::vector<double> arena;
        std::unordered_map<ContactKey, std::pair<double, Vec3>, ContactKeyHash> warm; // ln, (lt1, lt2, 0)
        std::unordered_map<std::uint64_t, RowState> rowState;
        WorldStats stats;

        BodySlot *body(BodyId id);
        const BodySlot *body(BodyId id) const;
        JointSlot *joint(JointId id);
        const JointSlot *joint(JointId id) const;

        // Groups.cpp
        void rebuildGroups();
        void bodyState(const BodySlot &b, Vec3 &c, Quat &q, Vec3 &v, Vec3 &w) const;

        // Machines.cpp
        void buildMachines();
        void readState(MachineData &md);
        void writeVelocities(MachineData &md);
        void writePoses(MachineData &md);

        // MachineRows.cpp
        void beginDrives(double dt);
        void solveInternal(MachineData &md);
        void afterContacts(MachineData &md);
        void solvePositions(MachineData &md);
        void integrate(MachineData &md);

        // Contacts.cpp
        void beginContacts();
        void prepareContacts();
        void prepareSprings();
        void solveContacts();
        void relaxContacts();
        void solveContactPositions();
        void restitution();
        void storeContacts();

        void step(double dt);
    };

} // namespace hf
