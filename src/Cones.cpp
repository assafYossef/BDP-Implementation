// Cones.cpp - curvature prescription: cone file, automatic cones, Gauss-Bonnet.
#include "bdp/Cones.h"
#include "bdp/Geometry.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace bdp {

std::vector<Cone> loadConeFile(const std::string& path, const LoadedMesh& lm, bool oneBased) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open cone file " + path);
    // OBJ index -> compacted mesh index
    std::vector<int> meshIndex;
    if (!lm.objIndexOfVertex.empty()) {
        const int maxObj = *std::max_element(lm.objIndexOfVertex.begin(), lm.objIndexOfVertex.end());
        meshIndex.assign(maxObj + 1, -1);
        for (std::size_t i = 0; i < lm.objIndexOfVertex.size(); ++i) meshIndex[lm.objIndexOfVertex[i]] = int(i);
    }
    std::vector<Cone> cones;
    std::string line;
    while (std::getline(in, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.resize(hash);
        std::istringstream ls(line);
        int vid; double deg;
        if (!(ls >> vid >> deg)) continue;
        if (oneBased) --vid;
        int v = vid;
        if (!meshIndex.empty()) v = (vid >= 0 && vid < int(meshIndex.size())) ? meshIndex[vid] : -1;
        if (v < 0 || v >= int(lm.mesh.number_of_vertices()))
            throw std::runtime_error("cone vertex " + std::to_string(vid) + " is not a mesh vertex");
        cones.push_back({v, deg * kPi / 180.0});
    }
    return cones;
}

std::vector<Cone> automaticCones(const MetricData& d, const Topology& topo) {
    std::vector<Cone> cones;
    const int chi = topo.eulerChar;
    if (!topo.closed() || chi == 0) return cones;
    const int nCones = 4 * std::abs(chi);
    const double k0 = (chi > 0 ? 1.0 : -1.0) * kPi / 2.0;   // +-90 degrees each
    const int nV = int(d.V.rows());

    // Candidates: the 25% of vertices whose original curvature has the right sign
    // and largest magnitude (negative cones also need valence >= 5 to open up an
    // angle of 450 degrees without extreme distortion).
    std::vector<int> cand(nV);
    std::iota(cand.begin(), cand.end(), 0);
    const double sgn = chi > 0 ? 1.0 : -1.0;
    std::sort(cand.begin(), cand.end(), [&](int a, int b) { return sgn * d.Korig(a) > sgn * d.Korig(b); });
    std::vector<int> pool;
    for (int v : cand) {
        if (chi < 0 && d.valence(v) < 5) continue;
        pool.push_back(v);
        if (int(pool.size()) >= std::max(nCones * 4, nV / 4)) break;
    }
    if (int(pool.size()) < nCones) throw std::runtime_error("not enough candidate vertices for automatic cones");

    // Farthest point sampling inside the pool, seeded with the most curved vertex.
    std::vector<double> dist(pool.size(), std::numeric_limits<double>::infinity());
    int next = 0;
    for (int c = 0; c < nCones; ++c) {
        const int v = pool[next];
        cones.push_back({v, k0});
        double best = -1;
        for (std::size_t i = 0; i < pool.size(); ++i) {
            dist[i] = std::min(dist[i], (d.V.row(pool[i]) - d.V.row(v)).norm());
            if (dist[i] > best) { best = dist[i]; next = int(i); }
        }
    }
    std::cout << "[cones] automatic placement: " << nCones << " cones of " << k0 * 180.0 / kPi
              << " deg (chi = " << chi << ")\n";
    return cones;
}

void setConesOnMesh(Mesh& mesh, const std::vector<Cone>& cones) {
    clearCones(mesh);
    auto [isCone, coneAngle] = coneMaps(mesh);
    for (const Cone& c : cones) {
        const VIdx v(static_cast<Mesh::size_type>(c.vertex));
        const double full = CGAL::is_border(v, mesh) ? kPi : 2.0 * kPi;
        isCone[v] = 1;
        coneAngle[v] = full - c.curvature;      // K0 = full - (sum of angles)
    }
}

void assignCurvatureConstraints(MetricData& d, const Topology& topo) {
    const int nV = int(d.V.rows());
    if (d.isCone.size() != nV || d.coneAngle.size() != nV)
        throw std::runtime_error("cone vectors not initialised (readConeProperties)");
    std::vector<int> vc;
    std::vector<double> k0;
    int nCones = 0;
    for (int v = 0; v < nV; ++v) {
        // Section 5: K(v) = 2 pi - alpha_v (interior), pi - alpha_v (boundary).
        // A cone prescribes alpha_v = coneAngle(v); other interior vertices are flat.
        const double full = d.isBnd(v) ? kPi : 2.0 * kPi;
        if (d.isCone(v)) {
            vc.push_back(v);
            k0.push_back(full - d.coneAngle(v));
            ++nCones;
        } else if (!d.isBnd(v)) {
            vc.push_back(v);
            k0.push_back(0.0);
        }
    }
    d.Vc = Eigen::Map<Eigen::VectorXi>(vc.data(), Eigen::Index(vc.size()));
    d.K0 = Eigen::Map<Eigen::VectorXd>(k0.data(), Eigen::Index(k0.size()));

    // Discrete Gauss-Bonnet: if every vertex is constrained the prescription must
    // sum to 2 pi chi, otherwise no metric can satisfy it.
    if (int(vc.size()) == nV) {
        const double total = d.K0.sum(), target = 2.0 * kPi * topo.eulerChar;
        if (std::abs(total - target) > 1e-6)
            throw std::runtime_error("cone curvatures sum to " + std::to_string(total * 180 / kPi) +
                                     " deg but Gauss-Bonnet requires 360*chi = " +
                                     std::to_string(360 * topo.eulerChar) +
                                     " deg (check the cone angles and the cone index base)");
    }
    std::cout << "[cones] " << nCones << " cone(s), " << vc.size()
              << " constrained vertices, " << (nV - int(vc.size())) << " free boundary vertices\n";
}

}  // namespace bdp
