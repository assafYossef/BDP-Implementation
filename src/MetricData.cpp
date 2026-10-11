// MetricData.cpp - CGAL topology/geometry -> Eigen arrays of the metric problem.
#include "bdp/MetricData.h"
#include "bdp/Geometry.h"

#include <Eigen/Dense>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace bdp {

MetricData buildMetricData(const LoadedMesh& lm) {
    const Mesh& mesh = lm.mesh;
    const int nV = int(mesh.number_of_vertices());
    const int nE = int(mesh.number_of_edges());
    const int nF = int(mesh.number_of_faces());
    MetricData d;

    // ---------------- connectivity ---------------------------------------
    d.V.resize(nV, 3);
    for (VIdx v : mesh.vertices()) {
        const auto& p = mesh.point(v);
        d.V.row(std::size_t(v)) << p.x(), p.y(), p.z();
    }
    d.EV.resize(nE, 2);
    for (EIdx e : mesh.edges()) {
        const HIdx h = mesh.halfedge(e);
        d.EV(std::size_t(e), 0) = int(std::size_t(mesh.source(h)));
        d.EV(std::size_t(e), 1) = int(std::size_t(mesh.target(h)));
    }
    // Face f with halfedges h0 -> h1 -> h2 (h0 = halfedge(f)):
    //   corner k = source(h_k),   edge opposite corner k = edge(h_{k+1})
    d.F.resize(nF, 3);
    d.F2E.resize(nF, 3);
    for (FIdx f : mesh.faces()) {
        HIdx h[3];
        h[0] = mesh.halfedge(f);
        h[1] = mesh.next(h[0]);
        h[2] = mesh.next(h[1]);
        for (int k = 0; k < 3; ++k) {
            d.F(std::size_t(f), k)   = int(std::size_t(mesh.source(h[k])));
            d.F2E(std::size_t(f), k) = int(std::size_t(mesh.edge(h[(k + 1) % 3])));
        }
    }
    d.isBnd = Eigen::VectorXi::Zero(nV);
    d.valence = Eigen::VectorXi::Zero(nV);
    for (VIdx v : mesh.vertices()) {
        d.isBnd(std::size_t(v)) = CGAL::is_border(v, mesh) ? 1 : 0;
        d.valence(std::size_t(v)) = int(mesh.degree(v));
    }

    // ---------------- original metric (normalised scale) -------------------
    // Squared lengths e0 of the input. The mesh is scaled so that the mean edge
    // length is 1: this conditions the numerics; distortion is measured relative
    // to this scaled copy, which only changes the global scale of the UVs (undone
    // by `scale` when the layout is exported).
    Eigen::VectorXd len(nE);
    for (int e = 0; e < nE; ++e) len(e) = (d.V.row(d.EV(e, 0)) - d.V.row(d.EV(e, 1))).norm();
    d.scale = len.mean();
    if (!(d.scale > 0)) throw std::runtime_error("degenerate mesh (zero edge lengths)");
    d.e0 = (len / d.scale).array().square();

    // ---------------- per-face linear maps (Section 4.3) --------------------
    d.cA.resize(nF, 3); d.cV.resize(nF, 3); d.cW.resize(nF, 3);
    d.wF.resize(nF);
    for (int f = 0; f < nF; ++f) {
        const double s0 = d.e0(d.F2E(f, 0)), s1 = d.e0(d.F2E(f, 1)), s2 = d.e0(d.F2E(f, 2));
        // Isometric chart g_i of the original triangle (Fig. 3):
        //   p0 = (0,0), p1 = (sqrt(s2), 0), p2 = (x, y), y > 0
        //   |p2 - p0|^2 = s1,  |p2 - p1|^2 = s0  ->  x = (s1 + s2 - s0) / (2 sqrt(s2))
        const double L2 = std::sqrt(s2);
        const double x  = (s1 + s2 - s0) / (2.0 * L2);
        const double y2 = s1 - x * x;
        if (!(y2 > 0)) throw std::runtime_error("degenerate triangle " + std::to_string(f) + " in the input mesh");
        const double y = std::sqrt(y2);
        const Eigen::Vector2d p0(0, 0), p1(L2, 0), p2(x, y);
        // Edge vectors opposite each corner, in the chart.
        const Eigen::Vector2d dv[3] = {p2 - p1, p0 - p2, p1 - p0};
        // s_k = d_k^T M d_k = m11 dx^2 + 2 m12 dx dy + m22 dy^2   ->  Q m = s
        Eigen::Matrix3d Q;
        for (int k = 0; k < 3; ++k) Q.row(k) << dv[k].x() * dv[k].x(), 2.0 * dv[k].x() * dv[k].y(), dv[k].y() * dv[k].y();
        const Eigen::Matrix3d Qi = Q.inverse();      // rows: m11, m12, m22 as functions of s
        d.cA.row(f) = 0.5 * (Qi.row(0) + Qi.row(2)); // A        = (m11 + m22) / 2
        d.cV.row(f) = 0.5 * (Qi.row(0) - Qi.row(2)); // 2Re(eta) = (m11 - m22) / 2
        d.cW.row(f) = Qi.row(1);                     // 2Im(eta) = m12
        d.wF(f) = 0.5 * L2 * y;                      // area of the original triangle
    }
    d.wF /= d.wF.sum();                              // w_i = Area_i / Area(S)  (Section 8.1)

    d.Korig = vertexCurvature(d.F, d.F2E, d.e0, d.isBnd);
    readConeProperties(mesh, d);

    // ---------------- reference metric for V5 (optional) --------------------
    if (lm.hasCornerUV) {
        d.Eref.resize(nF, 3);
        bool ok = true;
        for (int f = 0; f < nF && ok; ++f) {
            // UV of each corner of the mesh face, looked up in the soup triangle
            Eigen::Vector2d uv[3];
            for (int k = 0; k < 3; ++k) {
                int j = 0;
                while (j < 3 && lm.faceSoupVertices[f][j] != d.F(f, k)) ++j;
                if (j == 3) { ok = false; break; }
                uv[k] = Eigen::Vector2d(lm.faceCornerUV[f][j][0], lm.faceCornerUV[f][j][1]);
            }
            if (!ok) break;
            for (int k = 0; k < 3; ++k)   // edge opposite corner k joins corners k+1, k+2
                d.Eref(f, k) = (uv[(k + 2) % 3] - uv[(k + 1) % 3]).squaredNorm();
        }
        if (ok) {
            // The reference map is compared with the metric in normalised units;
            // rescale its lengths so that its total area matches the input area.
            double areaRef = 0, area3D = 0;
            for (int f = 0; f < nF; ++f) {
                areaRef += 0.25 * std::sqrt(std::max(0.0, heronT2(d.Eref(f, 0), d.Eref(f, 1), d.Eref(f, 2))));
                area3D  += 0.25 * std::sqrt(std::max(0.0, heronT2(d.e0(d.F2E(f, 0)), d.e0(d.F2E(f, 1)), d.e0(d.F2E(f, 2)))));
            }
            if (areaRef > 0) d.Eref *= area3D / areaRef;
        } else {
            d.Eref.resize(0, 0);
        }
    }
    return d;
}

void readConeProperties(const Mesh& mesh, MetricData& d) {
    const int nV = int(mesh.number_of_vertices());
    d.isCone = Eigen::VectorXi::Zero(nV);
    d.coneAngle = Eigen::VectorXd::Zero(nV);
    const auto isCone = mesh.property_map<VIdx, unsigned char>(kIsConeProperty);
    const auto angle  = mesh.property_map<VIdx, double>(kConeAngleProperty);
    if (!isCone || !angle) return;
    for (VIdx v : mesh.vertices()) {
        if (!(*isCone)[v]) continue;
        d.isCone(std::size_t(v)) = 1;
        d.coneAngle(std::size_t(v)) = (*angle)[v];
    }
}

namespace {
template <typename M>
void writeMatrix(const std::filesystem::path& p, const M& m, int offset = 0) {
    std::ofstream out(p);
    if (!out) throw std::runtime_error("cannot write " + p.string());
    out.precision(17);
    for (Eigen::Index i = 0; i < m.rows(); ++i) {
        for (Eigen::Index j = 0; j < m.cols(); ++j) out << (j ? " " : "") << m(i, j) + offset;
        out << '\n';
    }
}
}  // namespace

void dumpMetricData(const MetricData& d, const std::string& dir, const std::string& variant) {
    namespace fs = std::filesystem;
    fs::create_directories(dir);
    const fs::path p(dir);
    writeMatrix(p / "F.txt", d.F, 1);
    writeMatrix(p / "F2E.txt", d.F2E, 1);
    writeMatrix(p / "EV.txt", d.EV, 1);
    writeMatrix(p / "e0.txt", d.e0);
    writeMatrix(p / "cA.txt", d.cA);
    writeMatrix(p / "cV.txt", d.cV);
    writeMatrix(p / "cW.txt", d.cW);
    writeMatrix(p / "wF.txt", d.wF);
    writeMatrix(p / "isBnd.txt", d.isBnd);
    writeMatrix(p / "Vc.txt", d.Vc, 1);
    writeMatrix(p / "K0.txt", d.K0);
    writeMatrix(p / "isCone.txt", d.isCone);
    writeMatrix(p / "coneAngle.txt", d.coneAngle);
    if (d.Eref.size() > 0) writeMatrix(p / "Eref.txt", d.Eref);
    std::ofstream(p / "variant.txt") << variant << '\n';
    std::cout << "[dump] problem written to " << dir << "  (run bdp_run_from_dump('" << dir << "') in MATLAB)\n";
}

}  // namespace bdp
