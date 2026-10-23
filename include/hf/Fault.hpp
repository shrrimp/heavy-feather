#pragma once

#include <cstdint>

#include "hf/Types.hpp"

namespace hf {

    enum class FaultKind : std::uint8_t {
        NonFinite, // a NaN or an infinity
        Speed,     // a body faster than WorldSettings::faultSpeed
    };

    // Where in the step a fault was first seen. Without WorldSettings::checkPhases the world only looks once the step
    // is done, and every fault is seen at EndOfStep.
    enum class FaultPhase : std::uint8_t {
        EndOfStep,
        Internal,         // the machines' own solve: their dynamics, joints and drives
        Contacts,         // the contact sweeps
        AfterContacts,    // their answers back into the machines
        Positions,        // the machines' position pass
        ContactPositions, // the contacts' position pass
        Integrate,
        Relax,
        Restitution,
    };

    // An island that went wrong in a step: a NaN, or a body too fast, appeared in it. Before the host or any other
    // island could see it, the world put the island back as the step found it, at rest, its contacts and drives cold.
    struct Fault {
        std::uint64_t step = 0; // which step, counted from 1
        FaultKind kind = FaultKind::NonFinite;
        FaultPhase phase = FaultPhase::EndOfStep;
        BodyId body;      // the first offending body
        int machine = -1; // its machine; -1: none
        int bodies = 0;   // in the island put back
        int streak = 0;   // consecutive steps this body's island has faulted, this one included
    };

} // namespace hf
