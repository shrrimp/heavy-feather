#include "host/Collision.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

using hf::Vec3;

namespace host {

    namespace {

        struct Emitter {
            const Shape &a, &b;
            double margin;
            std::uint64_t base;
            std::vector<hf::Contact> &out;

            void operator()(const Vec3 &point, const Vec3 &normal, double depth, std::uint64_t feature) const {
                if (depth < -margin) return;
                hf::Contact c;
                c.a = a.body;
                c.b = b.body;
                c.point = point;
                c.normal = normal;
                c.depth = depth;
                c.feature = base | feature;
                c.friction = std::sqrt(a.friction * b.friction);
                c.restitution = std::max(a.restitution, b.restitution);
                c.rolling = std::max(a.rolling, b.rolling);
                c.spinning = std::max(a.spinning, b.spinning);
                out.push_back(c);
            }
        };

        Vec3 axis(const Pose &p, int i) { return p.R[i]; }

        void planeBox(const Pose &pl, const Shape &box, const Pose &pb, const Emitter &emit, double sign) {
            const Vec3 n = axis(pl, 1);
            for (int i = 0; i < 8; ++i) {
                const Vec3 local{(i & 1) ? box.half.x : -box.half.x, (i & 2) ? box.half.y : -box.half.y,
                                 (i & 4) ? box.half.z : -box.half.z};
                const Vec3 p = pb.c + pb.R * local;
                emit(p, sign * n, -dot(p - pl.c, n), std::uint64_t(i));
            }
        }

        void planeSphere(const Pose &pl, const Shape &s, const Pose &ps, const Emitter &emit, double sign) {
            const Vec3 n = axis(pl, 1);
            emit(ps.c - n * s.radius, sign * n, s.radius - dot(ps.c - pl.c, n), 0);
        }

        void sphereSphere(const Shape &a, const Pose &pa, const Shape &b, const Pose &pb, const Emitter &emit) {
            const Vec3 d = pb.c - pa.c;
            const double dist = length(d);
            const Vec3 n = dist > 1e-12 ? d / dist : Vec3{0, 1, 0};
            const double depth = a.radius + b.radius - dist;
            emit(pa.c + n * (a.radius - 0.5 * depth), n, depth, 0);
        }

        // The normal goes from the sphere to the box; sign flips it when the box is a.
        void sphereBox(const Shape &s, const Pose &ps, const Shape &box, const Pose &pb, const Emitter &emit,
                       double sign) {
            const Vec3 local = transpose(pb.R) * (ps.c - pb.c);
            const Vec3 clamped{std::clamp(local.x, -box.half.x, box.half.x),
                               std::clamp(local.y, -box.half.y, box.half.y),
                               std::clamp(local.z, -box.half.z, box.half.z)};
            const Vec3 d = clamped - local;
            const double dist = length(d);
            if (dist > 1e-9) {
                emit(pb.c + pb.R * clamped, sign * (pb.R * (d / dist)), s.radius - dist, 0);
                return;
            }
            int best = 0;
            double bestGap = 1e300, side = 1.0;
            for (int i = 0; i < 3; ++i)
                for (double sg : {-1.0, 1.0}) {
                    const double gap = box.half[i] - sg * local[i];
                    if (gap < bestGap) {
                        bestGap = gap;
                        best = i;
                        side = sg;
                    }
                }
            emit(ps.c, -sign * side * axis(pb, best), s.radius + bestGap, 1);
        }

    } // namespace

} // namespace host

namespace host {

    namespace {

        // A face clipped against another box's side planes has at most 8 corners.
        struct Polygon {
            Vec3 p[8];
            int n = 0;
        };

        // Keeps the part of poly where dot(p, axis) <= d.
        void clip(Polygon &poly, const Vec3 &axisDir, double d) {
            Polygon out;
            for (int i = 0; i < poly.n; ++i) {
                const Vec3 &p = poly.p[i], &q = poly.p[(i + 1) % poly.n];
                const double dp = dot(p, axisDir) - d, dq = dot(q, axisDir) - d;
                if (dp <= 0.0 && out.n < 8) out.p[out.n++] = p;
                if ((dp <= 0.0) != (dq <= 0.0) && out.n < 8) out.p[out.n++] = p + (q - p) * (dp / (dp - dq));
            }
            poly = out;
        }

        void closestOnLines(const Vec3 &p, const Vec3 &u, const Vec3 &q, const Vec3 &v, Vec3 &onP, Vec3 &onQ) {
            const Vec3 w = p - q;
            const double a = dot(u, u), b = dot(u, v), c = dot(v, v), d = dot(u, w), e = dot(v, w);
            const double den = a * c - b * b;
            const double s = den > 1e-12 ? (b * e - c * d) / den : 0.0;
            const double t = den > 1e-12 ? (a * e - b * d) / den : e / c;
            onP = p + s * u;
            onQ = q + t * v;
        }

        // Separating axes: three faces each and the nine edge pairs. An edge axis only wins when clearly better,
        // so resting faces keep their face contacts.
        void boxBox(const Shape &a, const Pose &pa, const Shape &b, const Pose &pb, const Emitter &emit) {
            const Vec3 d = pb.c - pa.c;
            const auto overlapOn = [&](const Vec3 &L) {
                double ra = 0.0, rb = 0.0;
                for (int i = 0; i < 3; ++i) {
                    ra += a.half[i] * std::fabs(dot(axis(pa, i), L));
                    rb += b.half[i] * std::fabs(dot(axis(pb, i), L));
                }
                return ra + rb - std::fabs(dot(d, L));
            };

            double bestFace = 1e300, bestEdge = 1e300;
            int face = -1, edgeA = -1, edgeB = -1;
            for (int i = 0; i < 6; ++i) {
                const double o = overlapOn(i < 3 ? axis(pa, i) : axis(pb, i - 3));
                if (o < -emit.margin) return;
                if (o < bestFace) {
                    bestFace = o;
                    face = i;
                }
            }
            Vec3 edgeAxis{};
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    Vec3 L = cross(axis(pa, i), axis(pb, j));
                    const double len = length(L);
                    if (len < 1e-6) continue;
                    L = L / len;
                    const double o = overlapOn(L);
                    if (o < -emit.margin) return;
                    if (o < bestEdge) {
                        bestEdge = o;
                        edgeA = i;
                        edgeB = j;
                        edgeAxis = L;
                    }
                }

            if (edgeA >= 0 && bestEdge < 0.95 * bestFace - 0.005) {
                Vec3 n = edgeAxis;
                if (dot(n, d) < 0.0) n = -n;
                Vec3 ea = pa.c, eb = pb.c;
                for (int k = 0; k < 3; ++k) {
                    if (k != edgeA) ea += (dot(axis(pa, k), n) > 0 ? 1.0 : -1.0) * a.half[k] * axis(pa, k);
                    if (k != edgeB) eb -= (dot(axis(pb, k), n) > 0 ? 1.0 : -1.0) * b.half[k] * axis(pb, k);
                }
                Vec3 onA, onB;
                closestOnLines(ea, axis(pa, edgeA), eb, axis(pb, edgeB), onA, onB);
                emit(0.5 * (onA + onB), n, bestEdge, std::uint64_t(0x300 + edgeA * 4 + edgeB));
                return;
            }

            // The reference box owns the face; the incident box's most opposed face is clipped to it.
            const bool refIsA = face < 3;
            const Shape &rs = refIsA ? a : b, &is = refIsA ? b : a;
            const Pose &rp = refIsA ? pa : pb, &ip = refIsA ? pb : pa;
            const int ri = face % 3;
            Vec3 n = axis(rp, ri);
            if (dot(n, refIsA ? d : -d) < 0.0) n = -n;

            int incAxis = 0;
            double most = -1.0;
            for (int k = 0; k < 3; ++k) {
                const double m = std::fabs(dot(axis(ip, k), n));
                if (m > most) {
                    most = m;
                    incAxis = k;
                }
            }
            const double incSign = dot(axis(ip, incAxis), n) > 0 ? -1.0 : 1.0;
            const int u = (incAxis + 1) % 3, v = (incAxis + 2) % 3;
            const Vec3 fc = ip.c + incSign * is.half[incAxis] * axis(ip, incAxis);
            Polygon poly;
            for (const auto [su, sv] : {std::pair{1.0, 1.0}, {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}})
                poly.p[poly.n++] = fc + su * is.half[u] * axis(ip, u) + sv * is.half[v] * axis(ip, v);
            for (int k = 1; k < 3; ++k) {
                const int t = (ri + k) % 3;
                const Vec3 T = axis(rp, t);
                clip(poly, T, dot(rp.c, T) + rs.half[t]);
                clip(poly, -T, -dot(rp.c, T) + rs.half[t]);
            }
            const Vec3 refPoint = rp.c + rs.half[ri] * n;
            const Vec3 normal = refIsA ? n : -n;
            const std::uint64_t f0 = (refIsA ? 0x100u : 0x200u) + std::uint64_t(face) * 16;
            for (int k = 0; k < poly.n; ++k)
                emit(poly.p[k], normal, -dot(poly.p[k] - refPoint, n), f0 + std::uint64_t(k));
        }

    } // namespace

    double boundingRadius(const Shape &s) {
        switch (s.kind) {
            case ShapeKind::Box: return length(s.half);
            case ShapeKind::Sphere: return s.radius;
            case ShapeKind::Plane: return 1e300;
        }
        return 0.0;
    }

    void collide(const Shape &a, const Pose &pa, const Shape &b, const Pose &pb, double margin,
                 std::uint64_t featureBase, std::vector<hf::Contact> &out) {
        const Emitter emit{a, b, margin, featureBase, out};
        const ShapeKind ka = a.kind, kb = b.kind;
        if (ka == ShapeKind::Plane && kb == ShapeKind::Box) planeBox(pa, b, pb, emit, 1.0);
        else if (ka == ShapeKind::Box && kb == ShapeKind::Plane) planeBox(pb, a, pa, emit, -1.0);
        else if (ka == ShapeKind::Plane && kb == ShapeKind::Sphere) planeSphere(pa, b, pb, emit, 1.0);
        else if (ka == ShapeKind::Sphere && kb == ShapeKind::Plane) planeSphere(pb, a, pa, emit, -1.0);
        else if (ka == ShapeKind::Sphere && kb == ShapeKind::Sphere) sphereSphere(a, pa, b, pb, emit);
        else if (ka == ShapeKind::Sphere && kb == ShapeKind::Box) sphereBox(a, pa, b, pb, emit, 1.0);
        else if (ka == ShapeKind::Box && kb == ShapeKind::Sphere) sphereBox(b, pb, a, pa, emit, -1.0);
        else if (ka == ShapeKind::Box && kb == ShapeKind::Box) boxBox(a, pa, b, pb, emit);
    }

} // namespace host
