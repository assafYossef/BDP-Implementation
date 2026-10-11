// Mesh.h - CGAL mesh type, OBJ input/output and topological summary.
//
// The 3D surface S of the paper (Section 4.1) is stored in a CGAL::Surface_mesh,
// a halfedge data structure that answers all topological queries we need:
// corner/edge incidences, boundary detection, dual-graph traversal for the layout.
//
// Cone singularities (Section 5) are attached to the mesh itself as two dynamic
// vertex property maps, so that every source of cones (the `c` records of the
// OBJ file, a separate cone file, the automatic placement) ends up in one place:
//
//     "v:bdp_is_cone"     unsigned char   1 if the vertex is a cone
//     "v:bdp_cone_angle"  double          prescribed cone angle (radians), i.e. the
//                                         target sum of the corner angles at the vertex
//
// MetricData.cpp converts these maps into Eigen vectors for the MATLAB solver.
#pragma once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace bdp {

using Kernel = CGAL::Simple_cartesian<double>;
using Point3 = Kernel::Point_3;
using Mesh   = CGAL::Surface_mesh<Point3>;
using VIdx   = Mesh::Vertex_index;
using HIdx   = Mesh::Halfedge_index;
using EIdx   = Mesh::Edge_index;
using FIdx   = Mesh::Face_index;

// ---- cone property maps on the CGAL mesh ------------------------------------------
inline const char* const kIsConeProperty    = "v:bdp_is_cone";
inline const char* const kConeAngleProperty = "v:bdp_cone_angle";
using IsConeMap    = Mesh::Property_map<VIdx, unsigned char>;
using ConeAngleMap = Mesh::Property_map<VIdx, double>;

// Returns the cone property maps, creating them (no cones, angle 0) if needed.
std::pair<IsConeMap, ConeAngleMap> coneMaps(Mesh& mesh);
// Removes every cone from the maps (they stay attached to the mesh).
void clearCones(Mesh& mesh);
// Number of vertices flagged as cones (0 if the maps do not exist).
int countCones(const Mesh& mesh);

// ---- OBJ reading -------------------------------------------------------------------
// Index base of the vertex index in `c <vertex> <angle/pi>` records.
//   Auto: 0-based if an index 0 occurs, 1-based if the largest possible 1-based index
//         (= number of `v` records) occurs, otherwise 1-based like OBJ face indices.
enum class ConeIndexBase { Auto, Zero, One };

struct ObjReadOptions {
    ConeIndexBase coneIndexBase = ConeIndexBase::Auto;
    bool verbose = true;
};

// Non-geometric content of the file, kept for logging / diagnostics.
struct ObjMetadata {
    std::vector<std::pair<std::string, std::string>> headers;  // "#scale_uv", "#num_v", ...
    int vLines = 0, vtLines = 0, fLines = 0;
    int coneLines = 0;          // `c` records
    int seamLines = 0;          // `sm` records (ignored)
    int ignoredLines = 0;       // vn, g, o, s, usemtl, mtllib, unknown tags
    int declaredNumV = -1;      // value of #num_v if present
    std::string coneIndexBaseUsed;   // "0-based" / "1-based" (when cones were read)
};

// Result of reading an OBJ file.
struct LoadedMesh {
    Mesh mesh;                  // carries the cone property maps (see above)
    ObjMetadata meta;
    // Optional texture coordinates of an existing (reference) parametrization,
    // one (u,v) per face corner, ordered like the polygon soup triangle
    // (i.e. like the face's vertices in the file). Used by variant V5.
    bool hasCornerUV = false;
    std::vector<std::array<std::array<double, 2>, 3>> faceCornerUV;  // per face
    std::vector<std::array<int, 3>> faceSoupVertices;                // per face, mesh vertex ids
    // compacted vertex index -> index of the vertex in the OBJ file (0-based)
    std::vector<int> objIndexOfVertex;
};

// Topological summary of S: Euler characteristic chi = V - E + F = 2 - 2g - m.
struct Topology {
    int nV = 0, nE = 0, nF = 0;
    int boundaryLoops = 0;   // m
    int eulerChar = 0;       // chi
    int genus = 0;           // g
    bool closed() const { return boundaryLoops == 0; }
};

// Reads an OBJ file and builds the CGAL halfedge mesh. Throws std::runtime_error.
//
//   v  x y z             vertex
//   vt u v               texture coordinate (reference map for variant V5)
//   f  a b c ...         facet; "a", "a/t", "a//n", "a/t/n", negative = relative;
//                        polygons are fan-triangulated
//   c  vid angle         cone singularity: vertex index and cone ANGLE in multiples
//                        of pi (e.g. "c 17 1.5" = 270 degree cone). Stored in the
//                        cone property maps of the mesh.
//   sm ...               seam definition: ignored (the layout computes its own cut)
//   #scale_uv, #translate_uv, #cone_grid_res, #num_v   metadata headers: logged
//   anything else        ignored (vn, g, o, s, usemtl, comments, ...)
//
// Vertices not referenced by any facet are dropped (their cones too, with a warning).
LoadedMesh loadObj(const std::string& path, const ObjReadOptions& options = ObjReadOptions());

Topology computeTopology(const Mesh& mesh);

// Writes the mesh with one texture coordinate per face corner (3|F| x 2, face
// major, corner k of face f = source of the k-th halfedge of f). Cones of the
// mesh are written back as `c` records (1-based vertex, angle / pi).
void writeObjWithUV(const std::string& path, const Mesh& mesh,
                    const std::vector<std::array<double, 2>>& cornerUV);

}  // namespace bdp
