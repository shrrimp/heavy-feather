#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#include "world/MachineJoints.hpp"
#include "world/WorldImpl.hpp"

namespace hf {

    namespace {

        constexpr Vec3 kX{1.0, 0.0, 0.0};

        // Which joints stay in the tree: Kruskal takes them in this order, so a loop is cut at the latest. A piston
        // or motor is always a tree joint, its drive then a single bounded dof, and the cut lands on a passive hinge.
        int treePriority(JointType t) {
            switch (t) {
                case JointType::Piston: return 0;
                case JointType::Motor:
                case JointType::Servo: return 1;
                case JointType::AngularSpring: return 2;
                case JointType::Slider: return 3;
                case JointType::Ball: return 4;
                case JointType::Hinge: return 5;
                case JointType::Fixed:
                case JointType::Spring: return 9;
            }
            return 9;
        }

        JointKind kindOf(JointType t) {
            switch (t) {
                case JointType::Slider:
                case JointType::Piston: return JointKind::Prismatic;
                case JointType::Ball: return JointKind::Spherical;
                default: return JointKind::Revolute;
            }
        }

        bool revolute(JointType t) {
            return t == JointType::Hinge || t == JointType::Motor || t == JointType::AngularSpring ||
                   t == JointType::Servo;
        }

        std::size_t findRoot(std::vector<std::size_t> &uf, std::size_t x) {
            while (uf[x] != x) x = uf[x] = uf[uf[x]];
            return x;
        }

    } // namespace

    bool kinematic(JointType t) {
        return t == JointType::Ball || t == JointType::Hinge || t == JointType::Motor ||
               t == JointType::AngularSpring || t == JointType::Servo || t == JointType::Slider ||
               t == JointType::Piston;
    }

    namespace {

        End endOf(const World::Impl &I, BodyId body, const Frame &f) {
            End e;
            if (!body.valid()) {
                e.group = 0;
                e.isStatic = true;
                e.anchor = f.position;
                e.frame = f.rotation;
                return e;
            }
            const BodySlot &b = I.bodies[body.index];
            e.group = b.group;
            e.isStatic = I.groups[std::size_t(b.group)].isStatic;
            e.anchor = b.offset + b.rot * (f.position - b.desc.centerOfMass);
            e.frame = normalize(b.rot * f.rotation);
            return e;
        }

        Vec3 axisOf(const World::Impl &I, const JointRef &r) {
            return normalize((I.groups[std::size_t(r.A.group)].q * r.A.frame) * kX);
        }

    } // namespace

    namespace {

        // Kruskal leaves a loop's cut on whichever passive hinge closes it last. On a loop authored slightly out of
        // plane that locks it, because a tree joint can't tilt. So the cut moves to the loop's outlier hinge, the one
        // least parallel to the others, and the authored error shows up at that one joint.
        template <class NodeFn>
        void selectCuts(const World::Impl &I, std::vector<JointRef> &refs, int K, const NodeFn &node) {
            constexpr double kNearParallel = 10.0 * std::numbers::pi / 180.0;
            constexpr double kWorthCutting = 0.2 * std::numbers::pi / 180.0;
            const auto type = [&](std::size_t e) { return I.joints[refs[e].joint].desc.type; };
            const auto angle = [](const Vec3 &a, const Vec3 &b) {
                return std::acos(std::min(1.0, std::fabs(dot(a, b))));
            };

            for (std::size_t c = 0; c < refs.size(); ++c) {
                if (refs[c].tree || !kinematic(type(c))) continue;

                std::vector<std::vector<std::pair<int, std::size_t>>> adj(std::size_t(K + 1));
                for (std::size_t o = 0; o < refs.size(); ++o) {
                    if (!refs[o].tree) continue;
                    const int a = node(refs[o].A), b = node(refs[o].B);
                    adj[std::size_t(a)].push_back({b, o});
                    adj[std::size_t(b)].push_back({a, o});
                }
                const int from = node(refs[c].A), to = node(refs[c].B);
                std::vector<int> prevNode(std::size_t(K + 1), -1);
                std::vector<std::size_t> prevEdge(std::size_t(K + 1), 0);
                std::vector<int> queue{from};
                prevNode[std::size_t(from)] = from;
                for (std::size_t head = 0; head < queue.size() && prevNode[std::size_t(to)] < 0; ++head)
                    for (auto [nb, o] : adj[std::size_t(queue[head])])
                        if (prevNode[std::size_t(nb)] < 0) {
                            prevNode[std::size_t(nb)] = queue[head];
                            prevEdge[std::size_t(nb)] = o;
                            queue.push_back(nb);
                        }
                if (prevNode[std::size_t(to)] < 0) continue;

                std::vector<std::size_t> loop{c};
                for (int v = to; v != from; v = prevNode[std::size_t(v)]) loop.push_back(prevEdge[std::size_t(v)]);
                std::vector<std::size_t> revs;
                for (std::size_t e : loop)
                    if (revolute(type(e))) revs.push_back(e);
                if (revs.size() < 3) continue;

                const auto score = [&](std::size_t e) {
                    double sum = 0.0;
                    int count = 0;
                    for (std::size_t f : revs) {
                        if (f == e) continue;
                        const double a = angle(axisOf(I, refs[e]), axisOf(I, refs[f]));
                        if (a < kNearParallel) {
                            sum += a;
                            ++count;
                        }
                    }
                    return count ? sum / count : 0.0;
                };
                std::size_t best = c;
                double bestScore = revolute(type(c)) ? score(c) : -1.0;
                for (std::size_t e : revs) {
                    if (e == c || type(e) != JointType::Hinge) continue;
                    const double sc = score(e);
                    if (sc > bestScore + 1e-6) {
                        bestScore = sc;
                        best = e;
                    }
                }
                if (best == c || bestScore < kWorthCutting) continue;
                refs[best].tree = false;
                refs[c].tree = true;
            }
        }

    } // namespace

    namespace {

        void buildMachine(World::Impl &I, MachineData &md, const std::vector<int> &groupsOf,
                          std::vector<JointRef> &refs) {
            // Local nodes: 0..K-1 the machine's groups, K the world.
            std::unordered_map<int, int> local;
            for (std::size_t i = 0; i < groupsOf.size(); ++i) local[groupsOf[i]] = int(i);
            const int K = int(groupsOf.size()), W = K;
            const auto node = [&](const End &e) { return e.isStatic ? W : local.at(e.group); };
            const auto type = [&](const JointRef &r) { return I.joints[r.joint].desc.type; };

            // Kruskal: drives first, then a fixed order, ties broken on the joint so the tree is reproducible.
            std::vector<std::size_t> order(refs.size());
            std::iota(order.begin(), order.end(), std::size_t(0));
            std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
                const int px = treePriority(type(refs[x])), py = treePriority(type(refs[y]));
                if (px != py) return px < py;
                return refs[x].joint < refs[y].joint;
            });
            std::vector<std::size_t> uf(std::size_t(K + 1));
            std::iota(uf.begin(), uf.end(), std::size_t(0));
            for (std::size_t o : order) {
                JointRef &r = refs[o];
                if (!kinematic(type(r))) continue;
                const std::size_t ra = findRoot(uf, std::size_t(node(r.A))), rb = findRoot(uf, std::size_t(node(r.B)));
                if (ra == rb) continue; // closes a loop: a cut
                uf[std::max(ra, rb)] = std::min(ra, rb);
                r.tree = true;
            }
            selectCuts(I, refs, K, node);

            std::vector<std::vector<std::pair<int, std::size_t>>> adj(std::size_t(K + 1));
            bool anchored = false;
            for (std::size_t o = 0; o < refs.size(); ++o) {
                if (!refs[o].tree) continue;
                const int a = node(refs[o].A), b = node(refs[o].B);
                adj[std::size_t(a)].push_back({b, o});
                adj[std::size_t(b)].push_back({a, o});
                anchored |= (a == W || b == W);
            }

            // The root: the world if a tree joint reaches it, else the heaviest group on a floating base.
            int start = W;
            if (!anchored) {
                double best = -1.0;
                for (int i = 0; i < K; ++i) {
                    const double m = 1.0 / std::max(I.groups[std::size_t(groupsOf[std::size_t(i)])].invMass, 1e-300);
                    if (m > best) {
                        best = m;
                        start = i;
                    }
                }
            }

            Articulation &art = md.m.art;
            art.links.clear();
            md.groupOfLink.clear();
            md.anchors.clear();
            std::vector<int> linkOfNode(std::size_t(K + 1), -1);
            const auto addLink = [&](int nd, int parentLink, const JointRef *via, bool parentIsA) {
                const Group &g = I.groups[std::size_t(groupsOf[std::size_t(nd)])];
                Link L;
                L.parent = parentLink;
                L.mass = 1.0 / g.invMass;
                L.inertiaBody = g.fixedRotation ? g.inertiaBody
                                : std::fabs(determinant(g.invInertiaBody)) > 1e-300 ? inverse(g.invInertiaBody)
                                                                                    : Mat3(L.mass * 1e-3);
                if (!via) {
                    L.kind = JointKind::Free;
                } else {
                    L.kind = kindOf(type(*via));
                    const End &P = parentIsA ? via->A : via->B;
                    const End &C = parentIsA ? via->B : via->A;
                    L.pp = P.anchor;
                    L.Qp = P.frame;
                    L.pc = C.anchor;
                    L.Qc = C.frame;
                    if (parentLink < 0) {
                        L.anchorPos = I.groups[std::size_t(P.group)].c;
                        L.anchorRot = I.groups[std::size_t(P.group)].q;
                        md.anchors.push_back({P.group, int(art.links.size())});
                    }
                }
                linkOfNode[std::size_t(nd)] = int(art.links.size());
                art.links.push_back(L);
                md.groupOfLink.push_back(groupsOf[std::size_t(nd)]);
            };

            // Breadth first from the root: parents before children.
            std::vector<std::uint8_t> seen(std::size_t(K + 1), 0);
            std::vector<int> queue;
            const auto visit = [&](int from, int parentLink) {
                for (auto [nb, o] : adj[std::size_t(from)]) {
                    if (seen[std::size_t(nb)] || nb == W) continue;
                    seen[std::size_t(nb)] = 1;
                    JointRef &r = refs[o];
                    const bool parentIsA = node(r.A) == from;
                    addLink(nb, parentLink, &r, parentIsA);
                    r.treeLink = linkOfNode[std::size_t(nb)];
                    r.sign = parentIsA ? 1.0 : -1.0;
                    queue.push_back(nb);
                }
            };
            seen[std::size_t(start)] = 1;
            if (start == W) {
                visit(W, -1);
            } else {
                addLink(start, -1, nullptr, true);
                queue.push_back(start);
            }
            for (std::size_t head = 0; head < queue.size(); ++head)
                visit(queue[head], linkOfNode[std::size_t(queue[head])]);
            art.finalise();

            for (JointRef &r : refs) {
                r.linkA = r.A.isStatic ? -1 : linkOfNode[std::size_t(local.at(r.A.group))];
                r.linkB = r.B.isStatic ? -1 : linkOfNode[std::size_t(local.at(r.B.group))];
            }

            // A floating machine closed onto something static is carried by that instead.
            md.anchored = false;
            for (const JointRef &r : refs) {
                if (!kinematic(type(r)) || !(r.A.isStatic || r.B.isStatic)) continue;
                md.anchored = true;
                const int group = r.A.isStatic ? r.A.group : r.B.group;
                bool known = start == W; // rooted on it already
                for (const Anchor &a : md.anchors) known = known || a.group == group;
                if (!known) {
                    const Group &g = I.groups[std::size_t(group)];
                    md.anchors.push_back({group, -1, g.c, g.q});
                }
            }

            // Two parallel hinges on one pair of links are one axle held at two points. Authored, the second point
            // is never exactly on the first's axis, so as a closure it would lock the loop: the first one on the pair
            // carries the axle and the second emits nothing.
            for (std::size_t x = 0; x < refs.size(); ++x) {
                JointRef &r = refs[x];
                if (r.tree || !revolute(type(r))) continue;
                const Vec3 ax = axisOf(I, r);
                for (std::size_t y = 0; y < refs.size(); ++y) {
                    const JointRef &o = refs[y];
                    if (y == x || o.duplicate || !revolute(type(o))) continue;
                    if (!(o.tree || y < x)) continue;
                    const bool samePair =
                        (o.linkA == r.linkA && o.linkB == r.linkB) || (o.linkA == r.linkB && o.linkB == r.linkA);
                    if (!samePair) continue;
                    if (std::fabs(dot(ax, axisOf(I, o))) > std::cos(5.0 * std::numbers::pi / 180.0)) {
                        r.duplicate = true;
                        break;
                    }
                }
            }
            md.joints = std::move(refs);
        }

    } // namespace

    // Rebuilds the machines of the touched groups, or every machine after a full rebuild of the groups. A machine
    // holding a touched group is taken apart, and its groups and the touched ones are put back together into
    // machines anew; every other machine stays as it was, its rows warm.
    void World::Impl::buildMachines() {
        const std::size_t nG = groups.size();
        const bool all = (touched.empty() && touchedMachines.empty()) || dirty;
        std::vector<std::uint8_t> &redo = groupMark;
        redo.assign(nG, all ? 1 : 0);
        machineMark.assign(machines.size(), all ? 1 : 0);
        if (!all) {
            for (int g : touched) {
                if (std::size_t(g) >= nG || (g > 0 && groups[std::size_t(g)].members == 0)) continue;
                redo[std::size_t(g)] = 1;
                if (groups[std::size_t(g)].machine >= 0) machineMark[std::size_t(groups[std::size_t(g)].machine)] = 1;
            }
            for (int m : touchedMachines)
                if (std::size_t(m) < machines.size()) machineMark[std::size_t(m)] = 1;
            for (std::size_t mi = 0; mi < machines.size(); ++mi)
                if (machineMark[mi])
                    for (int g : machines[mi].groupOfLink)
                        if (g == 0 || groups[std::size_t(g)].members > 0) redo[std::size_t(g)] = 1;
        }
        touched.clear();
        touchedMachines.clear();

        // What the machines taken apart had warm, for the ones built in their place.
        carriedRows.clear();
        for (std::size_t mi = 0; mi < machines.size(); ++mi) {
            if (!machineMark[mi]) continue;
            const MachineData &md = machines[mi];
            carriedRows.insert(carriedRows.end(), md.lastVel.begin(), md.lastVel.end());
            carriedRows.insert(carriedRows.end(), md.lastPos.begin(), md.lastPos.end());
        }
        std::sort(carriedRows.begin(), carriedRows.end(),
                  [](const auto &x, const auto &y) { return x.first < y.first; });
        // The links' responses name machines by index, and the machines that stay keep their indices into them.
        for (const LinkResponse &R : responses) machines[std::size_t(R.machine)].response[std::size_t(R.link)] = -1;
        responses.clear();

        // The machines that stay close up, in order.
        std::vector<int> newIndex(machines.size(), -1);
        std::size_t stay = 0;
        for (std::size_t mi = 0; mi < machines.size(); ++mi) {
            if (machineMark[mi]) continue;
            newIndex[mi] = int(stay);
            if (stay != mi) machines[stay] = std::move(machines[mi]);
            ++stay;
        }
        machines.resize(stay);
        for (std::size_t g = 0; g < nG; ++g) {
            Group &gr = groups[g];
            if (redo[g] || gr.machine < 0) gr.machine = gr.link = -1;
            else gr.machine = newIndex[std::size_t(gr.machine)];
        }
        for (JointSlot &j : joints) {
            if (j.info.machine < 0) continue;
            const int m = newIndex[std::size_t(j.info.machine)];
            if (m < 0) j.info = JointInfo{};
            else j.info.machine = m;
        }

        // Every joint's ends; only those reaching a group being rebuilt make machines, but every spring is sorted
        // anew between the machines and the loose ones.
        looseSprings.clear();
        std::vector<JointRef> refs;
        for (std::uint32_t k = 0; k < joints.size(); ++k) {
            const JointSlot &j = joints[k];
            if (!j.alive || j.desc.type == JointType::Fixed) continue;
            JointRef r;
            r.joint = k;
            r.A = endOf(*this, j.desc.a, j.desc.frameA);
            r.B = endOf(*this, j.desc.b, j.desc.frameB);
            if (r.A.group == r.B.group || (r.A.isStatic && r.B.isStatic)) continue;
            const bool rebuilt =
                (!r.A.isStatic && redo[std::size_t(r.A.group)]) || (!r.B.isStatic && redo[std::size_t(r.B.group)]);
            if (!rebuilt && j.desc.type != JointType::Spring) continue;
            refs.push_back(r);
        }

        // Machines are components over the dynamic-dynamic kinematic joints. A static end anchors one, never merges
        // two: every machine standing on the same ground stays its own.
        std::vector<std::size_t> uf(nG);
        std::iota(uf.begin(), uf.end(), std::size_t(0));
        std::vector<std::uint8_t> inMachine(nG, 0);
        for (const JointRef &r : refs) {
            if (!hf::kinematic(joints[r.joint].desc.type)) continue;
            if (!r.A.isStatic) inMachine[std::size_t(r.A.group)] = 1;
            if (!r.B.isStatic) inMachine[std::size_t(r.B.group)] = 1;
            if (r.A.isStatic || r.B.isStatic) continue;
            const std::size_t a = findRoot(uf, std::size_t(r.A.group)), b = findRoot(uf, std::size_t(r.B.group));
            if (a != b) uf[std::max(a, b)] = std::min(a, b);
        }
        std::vector<int> machineOfRoot(nG, -1);
        std::vector<std::vector<int>> groupsOf;
        for (std::size_t g = 0; g < nG; ++g) {
            if (!inMachine[g] || !redo[g]) continue;
            const std::size_t root = findRoot(uf, g);
            if (machineOfRoot[root] < 0) {
                machineOfRoot[root] = int(groupsOf.size());
                groupsOf.emplace_back();
            }
            groupsOf[std::size_t(machineOfRoot[root])].push_back(int(g));
        }
        // An end's machine once this is done: a new one, one that stayed, or none.
        const auto newMachine = [&](const End &e) {
            if (e.isStatic || !redo[std::size_t(e.group)] || !inMachine[std::size_t(e.group)]) return -1;
            return machineOfRoot[findRoot(uf, std::size_t(e.group))];
        };
        const auto stayed = [&](const End &e) { return e.isStatic ? -1 : groups[std::size_t(e.group)].machine; };

        std::vector<std::vector<JointRef>> refsOf(groupsOf.size());
        for (const JointRef &r : refs) {
            const int na = newMachine(r.A), nb = newMachine(r.B);
            if (hf::kinematic(joints[r.joint].desc.type)) {
                refsOf[std::size_t(na >= 0 ? na : nb)].push_back(r);
                continue;
            }
            // A spring inside one machine (or between it and the static world) is that machine's.
            if ((na >= 0 && (na == nb || r.B.isStatic)) || (nb >= 0 && r.A.isStatic)) {
                refsOf[std::size_t(na >= 0 ? na : nb)].push_back(r);
                continue;
            }
            const int sa = stayed(r.A), sb = stayed(r.B);
            if ((sa >= 0 && (sa == sb || r.B.isStatic)) || (sb >= 0 && r.A.isStatic)) continue; // a machine's own
            looseSprings.push_back(r);
        }

        const std::size_t first = machines.size();
        machines.resize(first + groupsOf.size());
        for (std::size_t k = 0; k < groupsOf.size(); ++k) {
            const std::size_t mi = first + k;
            MachineData &md = machines[mi];
            md = MachineData{};
            buildMachine(*this, md, groupsOf[k], refsOf[k]);
            for (std::size_t l = 0; l < md.groupOfLink.size(); ++l) {
                Group &g = groups[std::size_t(md.groupOfLink[l])];
                g.machine = int(mi);
                g.link = int(l);
            }
            for (const JointRef &r : md.joints) {
                JointInfo &info = joints[r.joint].info;
                info.machine = int(mi);
                info.tree = r.tree;
                info.duplicate = r.duplicate;
            }
            md.fext.assign(md.groupOfLink.size(), Force{});
        }
        machineRebuilds += int(groupsOf.size());
    }

    void World::Impl::writeVelocities(MachineData &md) {
        const Articulation &art = md.m.art;
        art.twists(art.qdot.data(), md.twists);
        for (std::size_t l = 0; l < md.groupOfLink.size(); ++l) {
            Group &g = groups[std::size_t(md.groupOfLink[l])];
            g.v = pointVelocity(md.twists[l], art.links[l].c - art.O);
            g.w = md.twists[l].w;
        }
    }

    void World::Impl::writePoses(MachineData &md) {
        const Articulation &art = md.m.art;
        for (std::size_t l = 0; l < md.groupOfLink.size(); ++l) {
            Group &g = groups[std::size_t(md.groupOfLink[l])];
            g.c = art.links[l].c;
            g.q = art.links[l].R;
        }
        writeVelocities(md);
    }

    // The articulation from the groups' poses and velocities. A tree joint the poses had drifted off is snapped
    // back, and written back so the first contact refresh sees the pose the solver will.
    void World::Impl::readState(MachineData &md) {
        followAnchors(md);
        // A kinematic anchor's pose as it is, its motion as the links' velocities were measured against.
        for (const Anchor &a : md.anchors) {
            if (a.link < 0 || !groups[std::size_t(a.group)].kinematic) continue;
            Link &L = md.m.art.links[std::size_t(a.link)];
            L.anchorPos = groups[std::size_t(a.group)].c;
            L.anchorRot = groups[std::size_t(a.group)].q;
        }
        const std::size_t n = md.groupOfLink.size();
        md.c.resize(n);
        md.R.resize(n);
        md.v.resize(n);
        md.w.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const Group &g = groups[std::size_t(md.groupOfLink[i])];
            md.c[i] = g.c;
            md.R[i] = g.q;
            md.v[i] = g.v;
            md.w[i] = g.w;
        }
        md.m.art.setFromMaximal(md.c, md.R, md.v, md.w);
        writePoses(md);
        // Its motion now: whatever changed since is a jump the joints take at the first dynamics.
        for (const Anchor &a : md.anchors) {
            if (a.link < 0) continue;
            const Group &g = groups[std::size_t(a.group)];
            Link &L = md.m.art.links[std::size_t(a.link)];
            const Vec3 v = g.kinematic ? g.v : Vec3{}, w = g.kinematic ? g.w : Vec3{};
            L.anchorJumpV = v - L.anchorV;
            L.anchorJumpW = w - L.anchorW;
            L.anchorV = v;
            L.anchorW = w;
        }
    }

} // namespace hf
