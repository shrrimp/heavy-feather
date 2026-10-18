#pragma once

#include <functional>
#include <memory>
#include <string>

#include "host/Host.hpp"

// Stress scenes for the determinism, performance and memory tests: a host with its bodies and joints, and what an
// operator does to them at each step.
namespace host {

    struct Scene {
        std::string name;
        std::unique_ptr<Host> host;
        std::function<void(Host &, int step)> operate;

        void step(int i, double dt) {
            if (operate) operate(*host, i);
            host->step(dt);
        }
    };

    Scene pile(int bodies);        // contacts: boxes and balls poured into a walled pit
    Scene chain(int links);        // one long hanging chain, draped over a row of crates
    Scene ring(int links);         // a closed loop of links, one machine with one huge loop
    Scene branches(int depth);     // a branching tree machine, 2^depth - 1 links
    Scene scissorLift(int stages); // a loop per stage, raised and lowered by one piston, a stack on top
    Scene yard();                  // everything at once: machines working in piles and on each other

    // A boom, stick and bucket, each on its own hydraulic cylinder: three loops.
    struct Excavator {
        hf::JointId boom, stick, bucket;
    };
    Excavator addExcavator(Host &h, const hf::Vec3 &at);

} // namespace host
