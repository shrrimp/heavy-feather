#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "hf/math/Mat3.hpp"
#include "hf/math/Quat.hpp"
#include "hf/math/Vec3.hpp"

// Seeded random inputs, one generator per test case so no test depends on another's draws.
struct Random {
    std::mt19937 gen;

    explicit Random(unsigned seed = 12345) : gen(seed) {}

    double uni(double a, double b) { return std::uniform_real_distribution<double>(a, b)(gen); }

    hf::Vec3 vec(double s = 1.0) { return {uni(-s, s), uni(-s, s), uni(-s, s)}; }

    hf::Quat quat() { return hf::normalize(hf::Quat{uni(-1, 1), uni(-1, 1), uni(-1, 1), uni(-1, 1)}); }

    // A box of random proportions, randomly rotated: always a physical inertia.
    hf::Mat3 inertia(double m) {
        const hf::Vec3 e{uni(0.2, 2.0), uni(0.2, 2.0), uni(0.2, 2.0)};
        const hf::Mat3 D(m / 12.0 * (e.y * e.y + e.z * e.z), 0, 0, 0, m / 12.0 * (e.x * e.x + e.z * e.z), 0, 0, 0,
                         m / 12.0 * (e.x * e.x + e.y * e.y));
        const hf::Mat3 R = hf::toMat3(quat());
        return R * D * hf::transpose(R);
    }
};

inline double maxAbsDiff(const std::vector<double> &a, const std::vector<double> &b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::fabs(a[i] - b[i]));
    return m;
}
