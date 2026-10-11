// Cones.h - prescription of the target Gaussian curvature (Section 5).
//
// The method flattens the metric everywhere except at a prescribed set of cone
// vertices whose curvatures must satisfy discrete Gauss-Bonnet:
//        sum_v K(v) = 2 pi chi = 2 pi (2 - 2g - m).
// Cone placement itself is not part of the paper (it uses the cones of the
// benchmark of Myles et al. 2014); we support a cone file and a simple
// automatic heuristic for closed surfaces.
#pragma once

#include "bdp/Mesh.h"
#include "bdp/MetricData.h"

#include <string>
#include <vector>

namespace bdp {

// Cone as produced by the cone file / the automatic placement. Cones read from
// `c` records of the OBJ file go directly into the mesh property maps (Mesh.h).
struct Cone {
    int vertex;        // mesh vertex index
    double curvature;  // prescribed angle defect K0 (radians)
};

// Cone file: one cone per line "<vertex index> <curvature in degrees>", e.g.
//   12  90      (a cone of angle 270 degrees at vertex 12)
//   57 -90      (a cone of angle 450 degrees)
// '#' starts a comment. Indices are 0-based OBJ vertex indices (use
// `--cones-one-based` for 1-based files).
std::vector<Cone> loadConeFile(const std::string& path, const LoadedMesh& lm, bool oneBased);

// Heuristic for closed surfaces when no cone file is given:
//   chi > 0 : 4*chi cones of +90 deg  at well separated vertices of high positive
//             Gaussian curvature (8 cones for a sphere),
//   chi < 0 : 4|chi| cones of -90 deg at well separated vertices of very negative
//             curvature and valence >= 5,
//   chi = 0 : no cones (a torus admits a flat metric).
std::vector<Cone> automaticCones(const MetricData& d, const Topology& topo);

// Writes cones into the CGAL vertex property maps of the mesh (replacing any
// cones already stored there). The curvature K0 is converted into the cone angle
//   angle = 2 pi - K0 (interior vertex),   angle = pi - K0 (boundary vertex).
void setConesOnMesh(Mesh& mesh, const std::vector<Cone>& cones);

// Fills d.Vc / d.K0 from the cone vectors d.isCone / d.coneAngle (which come
// from the CGAL property maps): every interior vertex is constrained (flat,
// K0 = 0, unless it is a cone); boundary vertices are free unless they carry a
// cone. Checks Gauss-Bonnet when every vertex is constrained. Throws on error.
void assignCurvatureConstraints(MetricData& d, const Topology& topo);

}  // namespace bdp
