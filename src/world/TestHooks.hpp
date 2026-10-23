#pragma once

#include "hf/Fault.hpp"
#include "hf/World.hpp"

// For the tests only: the library never sets these.
namespace hf::testing {

    // Called after every phase of every sub-step (substep is n for Restitution, after the last), before the world
    // looks for faults there: a test can write a NaN into the state at a known point.
    using PhaseHook = void (*)(World::Impl &world, int substep, FaultPhase phase);
    extern PhaseHook phaseHook;

} // namespace hf::testing
