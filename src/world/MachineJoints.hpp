#pragma once

#include "hf/Joint.hpp"

namespace hf {

    // A joint an articulation can carry: a tree joint, or a loop's cut held by rows.
    bool kinematic(JointType t);

} // namespace hf
