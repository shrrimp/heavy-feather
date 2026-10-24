#pragma once

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>
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

    // What a machine is anchored on: a group that doesn't move with it, its root links' (link >= 0, whose pose the
    // articulation keeps as the anchor's) or, for a floating machine, its closures'. Moved, it carries the machine.
    struct Anchor {
        int group = 0;
        int link = -1;
        Vec3 c{}; // as last followed, when link < 0
        Quat q{};
    };

    struct MachineData {
        Machine m;
        std::vector<int> groupOfLink;
        std::vector<Anchor> anchors;
        bool anchored = false; // any of its joints holds it to something that doesn't move with it
        std::vector<JointRef> joints;
        std::vector<Force> fext;
        std::vector<RowState> velState, posState;
        std::vector<std::uint64_t> velTags, posTags;
        std::vector<int> driveRow; // per joint ref: its drive row, or -1
        // Each row's active set as last solved, sorted by tag: the next solve's warm start. Only the rows of the
        // last solve are kept, so a removed joint leaves nothing behind.
        std::vector<std::pair<std::uint64_t, RowState>> lastVel, lastPos;

        // Per link: its response in this sub-step's cache, or -1.
        std::vector<int> response;
        // The contacts it solves with its joints: between two of its links, or one and the static world.
        std::vector<int> contacts;
        double correction = 0.0; // how far the position pass moved its links in the last sub-step, at most
        bool frozen = false;     // asleep this step, every link of it

        // Kept from step to step, so a warm step allocates nothing.
        std::vector<Motion> twists;
        std::vector<Vec3> c, v, w;
        std::vector<Quat> R;
        std::vector<double> row, errors, savedQ;
        std::vector<Quat> savedQRot;
    };

    struct ContactData {
        Contact in;
        int ga = 0, gb = 0; // groups
        Vec3 rA{}, rB{};    // the point from each group's centre, as the step began
        double depth0 = 0.0;
        Vec3 n{};
        double depth = 0.0; // refreshed every sub-step
        double ln = 0.0, lt1 = 0.0, lt2 = 0.0;
        double lr1 = 0.0, lr2 = 0.0, ls = 0.0; // angular impulses resisting rolling (about t1, t2) and spin (about n)
        double lnMax = 0.0, vnImpact = 0.0;
        // A contact solved in a machine's direct solve: that machine, its first row this sub-step (-1: none), and
        // the friction directions of its rows.
        int machine = -1, row = -1;
        Vec3 t1{}, t2{};
        int island = -1; // its moving end's
        // Sub-steps left of an impact: a gap closed in this one, and the next warm-starts with what stopped it.
        int impact = 0;
        bool fromKept = false; // the world kept it between sleepers: its impulses are its own warm start
        bool absorbed = false; // handed again for a contact the world keeps, which took its fresh look at it
    };

    // A contact end on a machine: its rows through the machine's mobility, n doubles each, in the arena as
    // gn, gt1, gt2, then W gn, W gt1, W gt2, and when the contact resists rolling or spin, the angular rows about t1,
    // t2 and n after them the same way.
    struct ContactBlock {
        int machine = -1;
        std::size_t off = 0;
        int n = 0;
    };

    // Two links of one machine its held rows lock together have ~0 relative mobility through W. Floored at this much
    // of their unconstrained mobility, far below any real mass ratio, a contact or spring between them stays finite.
    constexpr double kLockedMobility = 1e-7;

    // A contact's end on a machine link, through that link's shared response: its wrenches psi for n, t1, t2 and the
    // turn rows about t1, t2, n.
    struct MachineSide {
        int machine = -1;
        int response = -1;
        Force psi[6];
        int n = 0; // the machine's dofs
    };

    // A run of contacts on one link, swept against its twist, written to the machine once the run ends.
    struct LinkRun {
        int response = -1;
        Motion V;
        Force F;
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
        double kt12 = 0.0;  // t1 W t2, for the exact 2×2 friction solve on a machine end
        bool turns = false; // resists rolling or spin
        Vec3 imAr1{}, imAr2{}, imAs{}, imBr1{}, imBr2{}, imBs{}; // each end's I⁻¹ along t1, t2 and n
        double kr1 = 0.0, kr2 = 0.0, ks = 0.0;
        SoftParams soft;
        double vtarget = 0.0, refLn = 0.0, posBias = 0.0, pseudo = 0.0;
        // Which ends the sweeps write: a moving rigid group. A static, sleeping or machine end is only read, so
        // nothing a row computes ever lands in a group that another island reads.
        bool writeA = false, writeB = false;
        int ns = 0;
        MachineSide sides[2];
    };

    // A machine link's response to unit wrenches this sub-step (Machine::linkResponse), shared by its contacts: Z in
    // the response arena at off, and the link's 6×6 inverse inertia with nothing held.
    struct LinkResponse {
        int machine = -1, link = -1;
        std::size_t off = 0;
        double held[6][6] = {}, unheld[6][6] = {};
    };

    // A spring between bodies that aren't in one machine: a soft row swept with the contacts. Its velocity is
    // vA - vB along n at the anchors; a machine end keeps g then W g in the arena.
    struct SpringPrep {
        int a = 0, b = 0;
        Vec3 n{}, r1{}, r2{}, imR1n{}, imR2n{};
        double imA = 0.0, imB = 0.0;
        double k = 0.0, gamma = 0.0, target = 0.0, lambda = 0.0;
        bool writeA = false, writeB = false; // as a contact's
        int nb = 0;
        ContactBlock blocks[2];
    };

    // RigidRows.cpp: what a rigid island's contact rows push. A body that moves this step gets its own slot, and
    // anything else they touch (the ground, a sleeper, a kinematic body) a slot that is only ever read.
    struct SolverBody {
        Vec3 v{}, w{}, pv{}, pw{};
        Mat3 invI{}; // world, as the step began; zero for a body that doesn't move
        double invMass = 0.0;
        int group = 0;
        bool moves = false;
    };

    // A contact in a rigid island, lean: built once a step, its lever arms and masses held for the step, only its
    // depth refreshed each sub-step. Its impulses live here until the step is done.
    struct RigidRow {
        int contact = 0;
        int a = 0, b = 0; // solver bodies
        Vec3 n{}, t1{}, t2{}, r1{}, r2{};
        double kn = 0.0, kt1 = 0.0, kt2 = 0.0;
        double friction = 0.0;
        double ln = 0.0, lt1 = 0.0, lt2 = 0.0, lnMax = 0.0, refLn = 0.0;
        double vtarget = 0.0, posBias = 0.0, pseudo = 0.0, massCoeff = 1.0, impulseCoeff = 0.0;
        int turn = -1;           // into rigidTurns, when it resists rolling or spin
        std::uint8_t active = 0; // swept this sub-step: touching, or a gap it can reach
    };

    struct RigidTurn {
        double kr1 = 0.0, kr2 = 0.0, ks = 0.0;
        double lr1 = 0.0, lr2 = 0.0, ls = 0.0;
        double rolling = 0.0, spinning = 0.0;
    };

    // Sleep.cpp: a link between two dynamic groups, through the host's contacts or a joint, for the flood.
    struct SleepLink {
        int a = 0, b = 0;
        bool touching = false;    // a contact this step
        bool approaching = false; // and closing
        bool conducts = false;    // touching now or a step or two ago: force travels across it
        bool jointed = false;
    };

    // When two bodies last touched, so that a support that stops touching wakes what it held.
    struct Touch {
        std::uint64_t a = 0, b = 0; // each body's index and generation, a < b
        std::uint64_t step = 0;
        bool operator<(const Touch &o) const { return std::tie(a, b) < std::tie(o.a, o.b); }
    };

    struct ContactKey {
        BodyId a, b;
        std::uint64_t feature = 0;
        bool operator==(const ContactKey &) const = default;
        bool operator<(const ContactKey &o) const {
            return std::tie(a.index, a.generation, b.index, b.generation, feature) <
                   std::tie(o.a.index, o.a.generation, o.b.index, o.b.generation, o.feature);
        }
    };

    // Last step's impulses on a contact, found again by its key.
    struct WarmStart {
        ContactKey key;
        double ln = 0.0, lt1 = 0.0, lt2 = 0.0, lr1 = 0.0, lr2 = 0.0, ls = 0.0;
        bool operator<(const WarmStart &o) const { return key < o.key; }
    };

    // MachineSides.cpp: a contact's machine sides. The rows are n, t1, t2 and then the turn rows r1, r2, s; velocities
    // add into v, impulses come in l. pseudo picks the position pass's pseudo-velocity over the velocity. Anything
    // that reads or writes a machine's velocity directly must flush the runs first.
    void addMachineSides(World::Impl &I, ContactPrep &p, const ContactData &c);
    void machineVelocity(World::Impl &I, const ContactPrep &p, bool pseudo, int rows, double *v);
    void machineImpulse(World::Impl &I, const ContactPrep &p, bool pseudo, int rows, const double *l);
    void flushRuns(World::Impl &I, bool pseudo);
    void flushRuns(World::Impl &I, bool pseudo, const int *machines, int count);

    struct HF_NO_EXPORT World::Impl {
        WorldSettings settings;
        std::vector<BodySlot> bodies;
        std::vector<JointSlot> joints;
        std::vector<std::uint32_t> freeBodies, freeJoints;
        std::vector<Group> groups; // 0: the world
        std::vector<MachineData> machines;
        std::vector<JointRef> looseSprings; // springs not inside one machine, solved in the contact sweep
        // Every group and machine rebuilt, as a weld changing needs. Short of that, only the machines of the groups
        // touched: what a joint changing, or a body joining or leaving one, needs.
        bool dirty = true;
        std::vector<int> touched, touchedMachines, freeGroups;
        std::vector<std::uint8_t> machineMark;
        double h = 0.0;

        std::vector<ContactData> contacts;
        std::vector<ContactPrep> prep;
        std::vector<SpringPrep> springs;
        std::vector<double> arena;
        std::vector<LinkResponse> responses;
        std::vector<double> responseArena;
        std::vector<LinkRun> runs[2]; // per machine: on its velocity, and on its pseudo-velocity
        std::vector<WarmStart> warm;  // sorted by key
        // Contacts between bodies that are all asleep or static, kept from step to step with their impulses, so the
        // host hands only pairs with an awake end. Sorted by key.
        std::vector<ContactData> kept, keptScratch, newlyKept;
        // This step's handed contacts by key, and which pairs were handed: open addressing, reused step to step.
        std::vector<int> handedByKey, handedByPair;
        std::vector<std::pair<std::uint64_t, std::uint64_t>> keptPairs; // built when a lost support needs them
        // The rows of machines just rebuilt, so their first solve still warm-starts.
        std::vector<std::pair<std::uint64_t, RowState>> carriedRows;
        std::vector<int> groupDepth;
        std::vector<std::size_t> depthStart;
        std::vector<ContactData> contactScratch;
        WorldStats stats;

        // Faults.cpp: what a faulted island is put back to, and what the host is told.
        struct DriveSave {
            std::uint32_t joint = 0;
            JointSlot slot;
        };
        struct FaultSite {
            int island = -1;
            int group = -1; // the offending group, whose first body the fault names; -1: named already
        };
        std::vector<DriveSave> driveSaves; // the stepped machines' joints, as the step began
        std::vector<Fault> faults;         // this step's, and where each was found
        std::vector<FaultSite> faultSites;
        std::vector<std::uint8_t> islandFaulted;
        std::vector<int> bodyOfGroup, islandBodies;
        int rejected = 0; // calls refused since the last step

        // Moves.cpp: what the host wrote since the last step, and what the last step moved.
        std::vector<std::uint32_t> written;
        std::vector<BodyId> movedBodies;
        int groupRebuilds = 0, machineRebuilds = 0; // since the last step
        std::vector<std::uint8_t> groupMark;

        std::uint64_t stepCount = 0;
        bool reshapedAny = false;
        std::vector<SleepLink> sleepLinks;
        std::vector<Touch> touches, touchScratch;
        std::vector<int> linkStart, linkOf, floodDepth, floodQueue, unitOf, unitSteps, islandOf;
        int islands = 0;
        // This step's moving and kinematic groups, and its stepped machines: what a sub-step loops over.
        std::vector<int> steppedGroups, steppedMachines;
        std::vector<std::uint32_t> forced; // bodies with a force or torque for the next step
        // Contacts.cpp: per island, its rows (in prep), its springs and its machines, and whether it holds a machine.
        std::vector<int> rowStart, springStart, springOrder, machineStart, machineOrder, heads;
        std::vector<std::uint8_t> machineIsland;
        std::vector<int> impactRows; // the rows taking an impact this sub-step, island by island (impactStart)
        std::vector<int> impactStart;
        std::vector<SolverBody> solverBodies;
        std::vector<RigidRow> rigidRows; // island by island (rigidStart)
        std::vector<RigidTurn> rigidTurns;
        std::vector<int> rigidStart, slotOfGroup;
        int rowsOfKind[2] = {0, 0}; // rigid, machine
        std::vector<std::uint8_t> seeded, reached, unitStepped, unitStill, unitAsleep, unitVetoed;

        BodySlot *body(BodyId id);
        const BodySlot *body(BodyId id) const;
        JointSlot *joint(JointId id);
        const JointSlot *joint(JointId id) const;
        double jointPosition(const JointSlot &j) const; // from the bodies' poses

        // Groups.cpp
        void rebuildGroups();
        int addGroup(BodySlot &b); // a body alone: a group of its own, in place
        void freeGroup(int group);
        void setGroupMass(int group, const BodySlot &b);
        void touch(int group); // its machine, if any, and the one it may join, rebuilt before the next step
        void touchJoint(const JointDesc &d);
        void touchJointsOf(BodyId body);   // its joints' machines: it joins or leaves one
        void bodyMassChanged(BodySlot &b); // in place for a body alone in its group, else a rebuild
        template <typename Find>
        void placeWelds(Find &find);
        void bodyState(const BodySlot &b, Vec3 &c, Quat &q, Vec3 &v, Vec3 &w) const;

        // Machines.cpp
        void buildMachines(); // every machine, or only those of the touched groups
        void readState(MachineData &md);
        void writeVelocities(MachineData &md);
        void writePoses(MachineData &md);

        // MachineRows.cpp
        void beginDrives(double dt);
        void beginDrive(JointSlot &j, double dt);
        void solveInternal(MachineData &md);
        void afterContacts(MachineData &md);
        void solvePositions(MachineData &md);
        void integrate(MachineData &md);

        // Contacts.cpp
        void beginContacts();
        void contactRows(MachineData &md);         // its direct contacts' velocity rows, beside its joints'
        void contactPositionRows(MachineData &md); // and their position rows
        void recordContacts(MachineData &md);      // their answers, back into the contacts
        double contactDepth(const ContactData &c) const; // refreshed from the poses, after updateTurned
        void updateTurned();
        void prepareContacts();
        void prepareSprings();
        void solveContacts();
        void relaxContacts();
        void solveContactPositions();
        void restitution();
        void storeContacts();
        void keepResting(); // once the step is done: contacts left between sleepers stay in the world

        // RigidRows.cpp
        void beginRigidRows();             // once a step, after beginContacts
        void prepareRigidRows();           // each sub-step, before the sweeps
        void solveRigidIsland(int island); // its velocity sweeps
        void rigidPositions(int island);
        void relaxRigidIsland(int island);
        void rigidVelocitiesOut(); // the sweeps' answers, back into the groups
        void rigidPositionsOut();
        void rigidRestitution();
        void endRigidRows(); // the impulses, back into the contacts

        // Faults.cpp
        template <typename T>
        T reject(T refused) {
            ++rejected;
            return refused;
        }
        void beginFaults(); // before beginDrives: what a fault puts back
        void afterPhase(FaultPhase phase, int substep);
        void checkFaults(FaultPhase phase); // a fault in each island that has a NaN, or at the end a body too fast
        void rollBack();                    // the faulted islands, as the step found them
        bool faulted(const Group &g) const {
            return g.island >= 0 && std::size_t(g.island) < islandFaulted.size() &&
                   islandFaulted[std::size_t(g.island)];
        }

        // Moves.cpp
        void build(); // the groups and machines, if anything changed them
        void noteWritten(BodySlot &b);
        // Moves the marked groups rigidly, taking the pose (fromC, fromQ) to (toC, toQ): positions turn about fromC,
        // velocities turn too. Their bodies are moved by the host: their contacts and machine rows start cold.
        void moveGroups(const std::vector<std::uint8_t> &mark, const Vec3 &fromC, const Quat &fromQ, const Vec3 &toC,
                        const Quat &toQ);
        void movedByHost(const std::vector<std::uint8_t> &mark); // their bodies written, cold, and awake
        bool followAnchors(MachineData &md); // carried along by an anchor that moved: true if it was
        void carryAnchoredOn(int group);     // every machine anchored on the group, sleeping or not
        void markMachine(int machine, std::vector<std::uint8_t> &mark) const;
        void beginKinematic(double dt); // each kinematic group's velocity, to its target over the step
        void moveKinematic(Group &g, double t, bool last) const;
        void endStepBodies(); // what moved, and the host's writes forgotten

        // Sleep.cpp
        void beginSleep(); // who is stepped: the awake, and what they reach
        void endSleep();   // who falls asleep, and who lost the support they slept on
        void wakePartners(BodyId id); // whatever it touched lately

        void step(double dt);
    };

} // namespace hf
