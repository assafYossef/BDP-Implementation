// Layout.cpp - breadth-first unfolding of the triangles along the dual graph.
#include "bdp/Layout.h"
#include "bdp/Geometry.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <queue>

namespace bdp {

LayoutResult layoutMetric(const Mesh& mesh, const MetricData& d, const Eigen::VectorXd& e) {
    const int nF = int(mesh.number_of_faces());
    std::vector<Eigen::Vector2d> uv(3 * nF, Eigen::Vector2d::Zero());
    std::vector<char> placed(nF, 0);

    // corner index of halfedge h inside its face: corner k = source(h_k)
    auto cornerOf = [&](HIdx h) {
        const FIdx f = mesh.face(h);
        HIdx c = mesh.halfedge(f);
        int k = 0;
        while (c != h) { c = mesh.next(c); ++k; }
        return 3 * int(std::size_t(f)) + k;
    };
    auto len2 = [&](HIdx h) { return e(int(std::size_t(mesh.edge(h)))); };

    // ---- seed: the face of largest area in the target metric, laid out in the
    //      chart  p0 = (0,0), p1 = (|e_01|, 0), p2 with y > 0 (counter-clockwise)
    int seed = 0;
    double bestA = -1;
    for (int f = 0; f < nF; ++f) {
        const double T2 = heronT2(e(d.F2E(f, 0)), e(d.F2E(f, 1)), e(d.F2E(f, 2)));
        if (T2 > bestA) { bestA = T2; seed = f; }
    }
    {
        const double s0 = e(d.F2E(seed, 0)), s1 = e(d.F2E(seed, 1)), s2 = e(d.F2E(seed, 2));
        const double L = std::sqrt(s2);
        const double x = (s1 + s2 - s0) / (2 * L);
        uv[3 * seed + 0] = {0, 0};
        uv[3 * seed + 1] = {L, 0};
        uv[3 * seed + 2] = {x, std::sqrt(std::max(0.0, s1 - x * x))};
        placed[seed] = 1;
    }

    // ---- BFS over the dual graph
    std::queue<int> q;
    q.push(seed);
    while (!q.empty()) {
        const FIdx f(q.front());
        q.pop();
        HIdx h = mesh.halfedge(f);
        for (int k = 0; k < 3; ++k, h = mesh.next(h)) {
            const HIdx ho = mesh.opposite(h);
            if (mesh.is_border(ho)) continue;
            const FIdx g = mesh.face(ho);
            if (placed[std::size_t(g)]) continue;
            // In f, h runs a -> b. In g, ho runs b -> a, and the face reads (b, a, c)
            // counter-clockwise, so c lies to the LEFT of the directed edge b -> a.
            const Eigen::Vector2d pa = uv[cornerOf(h)];
            const Eigen::Vector2d pb = uv[cornerOf(mesh.next(h))];
            const HIdx gb = ho;                       // b -> a   (corner b)
            const HIdx ga = mesh.next(gb);            // a -> c   (corner a)
            const HIdx gc = mesh.next(ga);            // c -> b   (corner c)
            const double Lab2 = len2(gb), Lac2 = len2(ga), Lbc2 = len2(gc);
            // Place c in the frame of b with axis along b -> a (law of cosines):
            //   x = (|bc|^2 + |ab|^2 - |ac|^2) / (2|ab|),   y = sqrt(|bc|^2 - x^2)
            Eigen::Vector2d dir = pa - pb;
            const double L = std::sqrt(Lab2);
            dir.normalize();
            const Eigen::Vector2d nrm(-dir.y(), dir.x());   // +90 degrees (left)
            const double x = (Lbc2 + Lab2 - Lac2) / (2 * L);
            const double y = std::sqrt(std::max(0.0, Lbc2 - x * x));
            uv[cornerOf(gb)] = pb;
            uv[cornerOf(ga)] = pa;
            uv[cornerOf(gc)] = pb + x * dir + y * nrm;
            placed[std::size_t(g)] = 1;
            q.push(int(std::size_t(g)));
        }
    }

    // ---- diagnostics: orientation and edge length agreement
    LayoutResult res;
    res.cornerUV.resize(3 * nF);
    for (int f = 0; f < nF; ++f) {
        const Eigen::Vector2d& a = uv[3 * f];
        const Eigen::Vector2d& b = uv[3 * f + 1];
        const Eigen::Vector2d& c = uv[3 * f + 2];
        const double area2 = (b - a).x() * (c - a).y() - (b - a).y() * (c - a).x();
        if (!(area2 > 0)) ++res.flippedTriangles;
        for (int k = 0; k < 3; ++k) {
            const double s = e(d.F2E(f, k));
            const double l2 = (uv[3 * f + (k + 2) % 3] - uv[3 * f + (k + 1) % 3]).squaredNorm();
            res.maxEdgeLengthError = std::max(res.maxEdgeLengthError, std::abs(l2 - s) / s);
        }
        for (int k = 0; k < 3; ++k)
            res.cornerUV[3 * f + k] = {uv[3 * f + k].x() * d.scale, uv[3 * f + k].y() * d.scale};
    }
    return res;
}

}  // namespace bdp
