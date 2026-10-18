#include <hf/hf.hpp>

#include <cstdio>

int main() {
    const hf::Version v = hf::version();
    std::printf("Heavy Feather %d.%d.%d (headers %s)\n", v.major, v.minor, v.patch,
                HF_VERSION_STRING);
    return v.major == HF_VERSION_MAJOR && v.minor == HF_VERSION_MINOR ? 0 : 1;
}
