// MetricData.h - translation of the CGAL mesh into the dense Eigen arrays that
// define the optimisation problem in the space of metrics (Sections 4-6).
//
// Variables: the squared edge lengths e = (d_1, ..., d_|E|) of the mesh (the
// discrete metric). Everything the MATLAB solver needs is a function of the
// connectivity plus a few per-face LINEAR maps:
//
//   For triangle i with squared edge lengths s = (s_0, s_1, s_2) (s_k opposite
//   corner k), the pull-back metric tensor M_i of the affine map
//   "original triangle -> triangle with lengths s" is linear in s:
//        s_k = d_k^T M_i d_k        (d_k = edge vector in a planar chart of
//                                     the ORIGINAL triangle, Section 4.3)
//   Writing M = [m11 m12; m12 m22], the three equations Q m = s with
//   Q_k = [dx^2, 2 dx dy, dy^2] are inverted once per face:  m = Q^{-1} s.
//   From m we form the chart-invariant coordinates of Section 3.2.1:
//        A      = (m11 + m22)/2   = cA . s      (= ||J||_F^2 / 2)
//        2Re(eta) = (m11 - m22)/2 = cV . s
//        2Im(eta) = m12           = cW . s
//   so that sigma_a^2 = A + 2|eta|, sigma_b^2 = A - 2|eta| (Eqs. (1)-(2)).
//   These rows cA, cV, cW (|F| x 3) reproduce the closed forms of Section 4.3
//   and Appendix A (up to the normalisation A = 1 for the identity map).
#pragma once

#include "bdp/Mesh.h"

#include <Eigen/Core>
#include <string>

namespace bdp {

struct MetricData {
    // ---- topology (0-based here; converted to 1-based when sent to MATLAB)
    Eigen::MatrixXi F;      // |F| x 3   corner k of face f = source of k-th halfedge
    Eigen::MatrixXi F2E;    // |F| x 3   edge opposite corner k
    Eigen::MatrixXi EV;     // |E| x 2   edge end points
    Eigen::VectorXi isBnd;  // |V|       1 = boundary vertex
    Eigen::VectorXi valence;
    // ---- geometry
    Eigen::MatrixXd V;      // |V| x 3   vertex positions (original units)
    double scale = 1.0;     // positions were divided by `scale` (mean edge length := 1)
    Eigen::VectorXd e0;     // |E|       squared edge lengths of the (normalised) input
    Eigen::MatrixXd cA, cV, cW;   // |F| x 3 linear maps s -> A, 2Re(eta), 2Im(eta)
    Eigen::VectorXd wF;     // |F|       normalised original areas (sum = 1)
    Eigen::VectorXd Korig;  // |V|       angle defect of the input metric
    // ---- cone singularities, translated from the CGAL vertex property maps
    //      "v:bdp_is_cone" / "v:bdp_cone_angle" (see Mesh.h)
    Eigen::VectorXi isCone;     // |V|  1 = cone vertex
    Eigen::VectorXd coneAngle;  // |V|  prescribed cone angle (radians), 0 where not a cone
    // ---- curvature prescription (Section 5), derived from the cones:
    //      K0 = 2 pi - coneAngle (interior) or pi - coneAngle (boundary) at cones,
    //      K0 = 0 at the other interior vertices, free boundary elsewhere
    Eigen::VectorXi Vc;     // constrained vertices
    Eigen::VectorXd K0;     // prescribed curvatures (radians)
    // ---- optional reference metric for variant V5: per face, squared lengths
    //      (opposite corner k) of a reference map (e.g. existing UVs)
    Eigen::MatrixXd Eref;   // |F| x 3, empty if not available
};

// Builds F, F2E, EV, e0, cA, cV, cW, wF, isBnd, valence, Korig, and reads the
// cone property maps of the mesh (readConeProperties). Vc/K0 are filled by
// assignCurvatureConstraints (Cones.h).
MetricData buildMetricData(const LoadedMesh& lm);

// CGAL -> Eigen: copies the vertex property maps "v:bdp_is_cone" and
// "v:bdp_cone_angle" into d.isCone / d.coneAngle (zeros if the maps are absent).
void readConeProperties(const Mesh& mesh, MetricData& d);

// Writes every array as a whitespace separated text file <dir>/<name>.txt with
// 1-based indices, readable by matlab/bdp_run_from_dump.m.
void dumpMetricData(const MetricData& d, const std::string& dir, const std::string& variant);

}  // namespace bdp
