#pragma once

#include <cstddef>

// Every operator new in the test programs is counted, so a test can ask what a call allocated and what is still
// held.
namespace host {

    struct Allocations {
        std::size_t count = 0; // allocations made so far
        std::size_t live = 0;  // bytes allocated and not freed
    };

    Allocations allocations();

} // namespace host
