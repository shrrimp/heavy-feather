#include <doctest/doctest.h>

#include "hf/Version.hpp"

TEST_CASE("the linked binary is the version its headers describe") {
    const hf::Version v = hf::version();
    CHECK(v.major == HF_VERSION_MAJOR);
    CHECK(v.minor == HF_VERSION_MINOR);
    CHECK(v.patch == HF_VERSION_PATCH);
}
