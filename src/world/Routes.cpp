#include "world/Routes.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace hf {

    namespace {

        constexpr double kPi = std::numbers::pi;

        // Where a cord running round a circle of radius r about c, in the plane square to n, touches it: arriving from
        // a, or leaving toward a. Side +1 goes round anticlockwise about n, -1 clockwise. A point inside the circle
        // has no tangent: the rim's nearest point stands in.
        Vec3 tangentPoint(const Vec3 &a, const Vec3 &c, const Vec3 &n, double r, double side, bool arriving) {
            const Vec3 d = (a - c) - dot(a - c, n) * n;
            const double D = length(d);
            if (!(D > r)) return D > 0.0 ? c + (r / D) * d : c;
            const Vec3 e = d / D, f = cross(n, e);
            const double cosPhi = r / D, sinPhi = std::sqrt(std::max(0.0, 1.0 - cosPhi * cosPhi));
            return c + r * (cosPhi * e + (arriving ? side : -side) * sinPhi * f);
        }

        // How far round from where the cord arrives to where it leaves, the way it goes, in (-π, π].
        double arcAngle(const Vec3 &in, const Vec3 &out, const Vec3 &c, const Vec3 &n, double side) {
            const Vec3 a = in - c, b = out - c;
            return side * std::atan2(dot(n, cross(a, b)), dot(a, b));
        }

    } // namespace

    double laySpan(const std::vector<RopePoint> &route, std::vector<RopeBendState> &bends, bool first, bool closed) {
        const int count = int(route.size());
        for (int k = 0; k < count; ++k) {
            RopeBendState &b = bends[std::size_t(k)];
            b.dir = Vec3{};
            if (!sheave(route[std::size_t(k)])) {
                b.in = b.out = b.p;
                b.touching = true;
            } else if (first) {
                b.in = b.out = b.p;
            }
        }
        // The nearest touching point before k and after it; round the loop when it is closed.
        const auto touchingFrom = [&](int k, int step) {
            int reach = count - 1;
            if (!closed) reach = step < 0 ? k : count - 1 - k;
            for (int i = 1; i <= reach; ++i) {
                const int j = (k + step * i + count) % count;
                if (bends[std::size_t(j)].touching) return j;
            }
            return -1;
        };
        const auto prevOut = [&](int k) {
            const int j = touchingFrom(k, -1);
            return j >= 0 ? bends[std::size_t(j)].out : bends[std::size_t(k)].p;
        };
        const auto nextIn = [&](int k) {
            const int j = touchingFrom(k, 1);
            return j >= 0 ? bends[std::size_t(j)].in : bends[std::size_t(k)].p;
        };
        // Neighbouring sheaves depend on each other's tangent points: two passes, each from the last.
        for (int pass = 0; pass < 2; ++pass)
            for (int k = 0; k < count; ++k) {
                const RopePoint &rp = route[std::size_t(k)];
                if (!sheave(rp)) continue;
                RopeBendState &b = bends[std::size_t(k)];
                const double side = double(rp.side);
                // A drum at an end: the cord leaves it, or arrives at it, along one tangent.
                if (!closed && (k == 0 || k == count - 1)) {
                    b.in = b.out = k == 0 ? tangentPoint(nextIn(k), b.p, b.n, rp.radius, side, false)
                                          : tangentPoint(prevOut(k), b.p, b.n, rp.radius, side, true);
                    b.touching = true;
                    continue;
                }
                const Vec3 in = tangentPoint(prevOut(k), b.p, b.n, rp.radius, side, true);
                const Vec3 out = tangentPoint(nextIn(k), b.p, b.n, rp.radius, side, false);
                const double raw = arcAngle(in, out, b.p, b.n, side);
                double theta;
                if (first && pass == 0) theta = raw < -0.5 * kPi ? raw + 2.0 * kPi : raw;
                else theta = raw + 2.0 * kPi * std::round((b.theta - raw) / (2.0 * kPi));
                b.in = in;
                b.out = out;
                b.theta = theta;
                b.touching = theta > 0.0;
            }

        double span = 0.0;
        Vec3 u{};
        const auto segment = [&](int from, int to) {
            const Vec3 seg = bends[std::size_t(to)].in - bends[std::size_t(from)].out;
            const double len = length(seg);
            if (len > 1e-12) u = seg / len; // a segment of no length runs on the way the last one did
            span += len;
            bends[std::size_t(from)].dir -= u;
            bends[std::size_t(to)].dir += u;
        };
        // A closed route starts with the stretch that closes it, from its last touching point round to the first.
        int prev = 0;
        if (closed) {
            prev = -1;
            for (int k = count - 1; k >= 0 && prev < 0; --k)
                if (bends[std::size_t(k)].touching) prev = k;
        }
        for (int k = closed ? 0 : 1; k < count && prev >= 0; ++k) {
            if (!bends[std::size_t(k)].touching) continue;
            segment(prev, k);
            prev = k;
        }
        const int firstArc = closed ? 0 : 1, endArc = closed ? count : count - 1;
        for (int k = firstArc; k < endArc; ++k) {
            const RopeBendState &b = bends[std::size_t(k)];
            if (sheave(route[std::size_t(k)]) && b.touching) span += route[std::size_t(k)].radius * b.theta;
        }
        return span;
    }

} // namespace hf
