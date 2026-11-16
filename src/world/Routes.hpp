#pragma once

#include <vector>

#include "hf/Rope.hpp"
#include "world/Storage.hpp"

// A cord laid through its route: a rope from end to end, or a belt round and back to where it started.
namespace hf {

    inline bool sheave(const RopePoint &p) { return p.radius > 0.0; }

    // The cord through its route, each bend's p and n where its point is: each sheave's tangent points and angle,
    // which way the span grows as each point moves, and the span. A sheave starts from its last answer, and its angle
    // stays continuous; one unwound past nothing is off, and the cord runs straight past it. The first time, a
    // sheave's angle is taken in [-π/2, 3π/2): a half turn round it, as a block's, is a half turn.
    //
    // An open route's first and last points are its ends; a closed one has none, and its last point runs on to its
    // first.
    double laySpan(const std::vector<RopePoint> &route, std::vector<RopeBendState> &bends, bool first, bool closed);

} // namespace hf
