#include <cmath>

#include <doctest/doctest.h>

#include "Random.hpp"
#include "math/Spatial.hpp"

using namespace hf;

TEST_CASE("spatial algebra") {
    Random rnd;
    for (int t = 0; t < 20; ++t) {
        const double m = rnd.uni(0.1, 50.0);
        const Vec3 c = rnd.vec(3.0);
        const Mat3 Ic = rnd.inertia(m);
        const Inertia I = Inertia::body(m, c, Ic);
        const Motion V{rnd.vec(2.0), rnd.vec(2.0)};
        const Vec3 vc = pointVelocity(V, c);
        const double keDirect = 0.5 * m * dot(vc, vc) + 0.5 * dot(V.w, Ic * V.w);
        const double keSpatial = 0.5 * dot(V, I * V);

        CHECK_MESSAGE(std::fabs(keDirect - keSpatial) < 1e-9 * (1.0 + keDirect),
                      "kinetic energy about O = about the centre");
        CHECK_MESSAGE(std::fabs(dot(V, crossF(V, I * V))) < 1e-8 * (1.0 + keDirect), "velocity product does no work");
    }
}
