// Layout.h - realisation of a (flat) discrete metric in the plane.
//
// Section 7.2 / [Springborn et al. 2008, Sec. 3.3]: once the metric is flat
// away from the cones, the triangles are laid out one by one while traversing
// a spanning tree of the dual graph. The non-tree edges form the cut graph; it
// passes through every cone, so the map is locally injective (flatness at every
// vertex that is interior to the cut disk) and the transition across the seams
// are rigid motions.
#pragma once

#include "bdp/Mesh.h"
#include "bdp/MetricData.h"

#include <array>
#include <vector>

namespace bdp {

struct LayoutResult {
    std::vector<std::array<double, 2>> cornerUV;  // 3|F|, corner k of face f at 3f+k
    int flippedTriangles = 0;                     // faces with non-positive signed area
    double maxEdgeLengthError = 0.0;              // max | |uv_i - uv_j|^2 - e_ij | / e_ij
};

// Lays out the metric e (normalised units) and rescales the UVs to the units of
// the input mesh (multiplying by d.scale).
LayoutResult layoutMetric(const Mesh& mesh, const MetricData& d, const Eigen::VectorXd& e);

}  // namespace bdp
