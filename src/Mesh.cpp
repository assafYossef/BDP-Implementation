// Mesh.cpp - OBJ parsing, CGAL mesh construction and topology queries.
#include "bdp/Mesh.h"

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace bdp {

namespace {

constexpr double kPiMesh = 3.14159265358979323846;

[[noreturn]] void parseError(const std::string& path, int lineNo, const std::string& what) {
    throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": " + what);
}

// Parses one "f" token "v", "v/vt", "v//vn" or "v/vt/vn" (1-based, negative = relative).
// Returns false on malformed input.
bool parseFaceToken(const std::string& tok, int nV, int nVT, int& vi, int& ti) {
    vi = -1; ti = -1;
    auto toIndex = [](const std::string& s, int n, int& out) {
        if (s.empty()) return false;
        std::size_t pos = 0;
        long long i = 0;
        try { i = std::stoll(s, &pos); } catch (...) { return false; }
        if (pos != s.size() || i == 0) return false;
        out = int(i > 0 ? i - 1 : n + i);        // OBJ: negative indices are relative to the end
        return true;
    };
    const auto s1 = tok.find('/');
    if (!toIndex(tok.substr(0, s1), nV, vi)) return false;
    if (s1 != std::string::npos) {
        const auto s2 = tok.find('/', s1 + 1);
        const std::string t = tok.substr(s1 + 1, s2 == std::string::npos ? std::string::npos : s2 - s1 - 1);
        if (!t.empty() && !toIndex(t, nVT, ti)) return false;
    }
    return true;
}

// Metadata headers that are logged (all other comments are skipped silently).
bool isKnownHeader(const std::string& key) {
    return key == "scale_uv" || key == "translate_uv" || key == "cone_grid_res" || key == "num_v";
}

struct RawCone {
    long long index;     // as written in the file
    double angleOverPi;  // cone angle / pi
    int line;
};

}  // namespace

// ------------------------------------------------------------------------------------
//  Cone property maps
// ------------------------------------------------------------------------------------
std::pair<IsConeMap, ConeAngleMap> coneMaps(Mesh& mesh) {
    IsConeMap isCone = mesh.add_property_map<VIdx, unsigned char>(kIsConeProperty, 0).first;
    ConeAngleMap angle = mesh.add_property_map<VIdx, double>(kConeAngleProperty, 0.0).first;
    return {isCone, angle};
}

void clearCones(Mesh& mesh) {
    auto [isCone, angle] = coneMaps(mesh);
    for (VIdx v : mesh.vertices()) { isCone[v] = 0; angle[v] = 0.0; }
}

int countCones(const Mesh& mesh) {
    const auto isCone = mesh.property_map<VIdx, unsigned char>(kIsConeProperty);
    if (!isCone) return 0;
    int n = 0;
    for (VIdx v : mesh.vertices()) n += (*isCone)[v] ? 1 : 0;
    return n;
}

// ------------------------------------------------------------------------------------
//  OBJ reader
// ------------------------------------------------------------------------------------
LoadedMesh loadObj(const std::string& path, const ObjReadOptions& options) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);

    LoadedMesh out;
    ObjMetadata& meta = out.meta;
    std::vector<Point3> points;
    std::vector<std::array<double, 2>> uvs;
    std::vector<std::vector<std::size_t>> polys;      // triangles of the soup
    std::vector<std::array<int, 3>> polyUV;           // vt per triangle corner (-1 = none)
    std::vector<RawCone> rawCones;

    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        // tolerate Windows line endings and a UTF-8 byte order mark
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (lineNo == 1 && line.size() >= 3 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);

        std::istringstream ls(line);
        std::string tag;
        if (!(ls >> tag)) continue;                    // empty line

        if (tag[0] == '#') {
            // ---- comments and metadata headers ("#scale_uv 1.0 1.0", "# num_v 1234", ...)
            std::string key = tag.substr(1);
            if (key.empty()) ls >> key;
            std::string value;
            std::getline(ls, value);
            const auto b = value.find_first_not_of(" \t");
            value = b == std::string::npos ? std::string() : value.substr(b);
            if (isKnownHeader(key)) {
                meta.headers.emplace_back(key, value);
                if (key == "num_v") {
                    try { meta.declaredNumV = std::stoi(value); } catch (...) {}
                }
            }
        } else if (tag == "v") {
            double x, y, z;
            if (!(ls >> x >> y >> z)) parseError(path, lineNo, "malformed vertex record");
            points.emplace_back(x, y, z);
            ++meta.vLines;
        } else if (tag == "vt") {
            double u = 0, v = 0;
            ls >> u >> v;
            uvs.push_back({u, v});
            ++meta.vtLines;
        } else if (tag == "f") {
            std::vector<int> fv, ft;
            std::string tok;
            while (ls >> tok) {
                int vi, ti;
                if (!parseFaceToken(tok, int(points.size()), int(uvs.size()), vi, ti))
                    parseError(path, lineNo, "malformed facet index '" + tok + "'");
                fv.push_back(vi);
                ft.push_back(ti);
            }
            if (fv.size() < 3) parseError(path, lineNo, "facet with fewer than 3 vertices");
            for (std::size_t k = 1; k + 1 < fv.size(); ++k) {   // fan triangulation
                polys.push_back({std::size_t(fv[0]), std::size_t(fv[k]), std::size_t(fv[k + 1])});
                polyUV.push_back({ft[0], ft[k], ft[k + 1]});
            }
            ++meta.fLines;
        } else if (tag == "c") {
            // ---- cone singularity: "c <vertex index> <cone angle in multiples of pi>"
            RawCone c;
            if (!(ls >> c.index >> c.angleOverPi))
                parseError(path, lineNo, "malformed cone record (expected: c <vertex> <angle/pi>)");
            if (!(c.angleOverPi > 0)) parseError(path, lineNo, "cone angle must be positive");
            c.line = lineNo;
            rawCones.push_back(c);
            ++meta.coneLines;
        } else if (tag == "sm") {
            ++meta.seamLines;                          // seams: our layout computes its own cut
        } else {
            ++meta.ignoredLines;                       // vn, g, o, s, usemtl, mtllib, ...
        }
    }
    if (polys.empty()) throw std::runtime_error("no faces found in " + path);

    if (options.verbose) {
        std::cout << "[obj] " << meta.vLines << " v, " << meta.fLines << " f, " << meta.vtLines << " vt, "
                  << meta.coneLines << " c (cones), " << meta.seamLines << " sm (seams, ignored)";
        if (meta.ignoredLines) std::cout << ", " << meta.ignoredLines << " other records ignored";
        std::cout << "\n";
        for (const auto& [k, v] : meta.headers) std::cout << "[obj] header #" << k << " = " << v << "\n";
    }
    if (meta.declaredNumV >= 0 && meta.declaredNumV != meta.vLines)
        std::cerr << "[obj] warning: #num_v = " << meta.declaredNumV << " but the file has "
                  << meta.vLines << " vertices\n";

    // ---- drop vertices not referenced by any face (they would be isolated
    //      "interior" vertices with curvature 2*pi and break the curvature constraints)
    //      Kept vertices retain the relative order of the file, so that mesh vertex i
    //      is OBJ vertex i whenever every vertex is used.
    std::vector<char> used(points.size(), 0);
    for (const auto& p : polys)
        for (const auto v : p) {
            if (v >= points.size()) throw std::runtime_error("face index out of range in " + path);
            used[v] = 1;
        }
    std::vector<int> newIndex(points.size(), -1);
    std::vector<Point3> usedPoints;
    for (std::size_t v = 0; v < points.size(); ++v)
        if (used[v]) {
            newIndex[v] = int(usedPoints.size());
            usedPoints.push_back(points[v]);
            out.objIndexOfVertex.push_back(int(v));
        }
    if (usedPoints.size() != points.size())
        std::cerr << "[obj] warning: " << points.size() - usedPoints.size()
                  << " vertices are used by no facet and were dropped\n";
    for (auto& p : polys)
        for (auto& v : p) v = std::size_t(newIndex[v]);

    // ---- polygon soup -> halfedge mesh. If the soup is not a valid oriented
    //      2-manifold it is repaired/oriented by CGAL (the corner UVs are then
    //      no longer aligned with the faces and are discarded).
    bool uvValid = !uvs.empty();
    for (const auto& t : polyUV)
        for (int k = 0; k < 3; ++k) uvValid = uvValid && t[k] >= 0 && t[k] < int(uvs.size());

    bool repaired = false;
    if (!PMP::is_polygon_soup_a_polygon_mesh(polys)) {
        std::cerr << "[mesh] polygon soup is not a valid oriented manifold mesh: repairing/orienting\n";
        PMP::repair_polygon_soup(usedPoints, polys);
        PMP::orient_polygon_soup(usedPoints, polys);
        uvValid = false;
        repaired = true;
        out.objIndexOfVertex.clear();   // indices were renumbered
    }
    PMP::polygon_soup_to_polygon_mesh(usedPoints, polys, out.mesh);
    if (!CGAL::is_triangle_mesh(out.mesh)) throw std::runtime_error("mesh is not a triangle mesh");

    // ---- connectivity: the method parametrizes one connected surface
    auto fcc = out.mesh.add_property_map<FIdx, std::size_t>("f:cc", 0).first;
    const std::size_t ncc = PMP::connected_components(out.mesh, fcc);
    out.mesh.remove_property_map(fcc);
    if (ncc != 1)
        throw std::runtime_error("mesh has " + std::to_string(ncc) +
                                 " connected components; please provide a single connected surface");

    // ---- cones: OBJ vertex index -> CGAL vertex, stored in the vertex property maps
    auto [isCone, coneAngle] = coneMaps(out.mesh);
    if (!rawCones.empty()) {
        // index base of the `c` records
        bool zeroBased;
        long long minIdx = rawCones[0].index, maxIdx = rawCones[0].index;
        for (const RawCone& c : rawCones) { minIdx = std::min(minIdx, c.index); maxIdx = std::max(maxIdx, c.index); }
        switch (options.coneIndexBase) {
            case ConeIndexBase::Zero: zeroBased = true; break;
            case ConeIndexBase::One:  zeroBased = false; break;
            default:
                if (minIdx == 0) zeroBased = true;                                  // 0 only exists 0-based
                else if (maxIdx == (long long)points.size()) zeroBased = false;     // n only exists 1-based
                else {
                    zeroBased = false;                                              // OBJ convention
                    std::cout << "[obj] cone indices are ambiguous; assuming 1-based like OBJ facets "
                                 "(override with --cone-index-base 0)\n";
                }
        }
        meta.coneIndexBaseUsed = zeroBased ? "0-based" : "1-based";

        // after a repair the vertices were renumbered: find them by position
        std::map<std::array<double, 3>, VIdx> byPosition;
        if (repaired)
            for (VIdx v : out.mesh.vertices()) {
                const Point3& p = out.mesh.point(v);
                byPosition[{p.x(), p.y(), p.z()}] = v;
            }

        int stored = 0;
        for (const RawCone& c : rawCones) {
            const long long obj = zeroBased ? c.index : c.index - 1;
            if (obj < 0 || obj >= (long long)points.size())
                parseError(path, c.line, "cone vertex index " + std::to_string(c.index) + " out of range (" +
                                         meta.coneIndexBaseUsed + ", " + std::to_string(points.size()) + " vertices)");
            VIdx v;
            if (repaired) {
                const Point3& p = points[std::size_t(obj)];
                const auto it = byPosition.find({p.x(), p.y(), p.z()});
                if (it == byPosition.end()) {
                    std::cerr << "[obj] warning: cone at line " << c.line << " is on a vertex removed by the repair; skipped\n";
                    continue;
                }
                v = it->second;
            } else {
                if (newIndex[std::size_t(obj)] < 0) {
                    std::cerr << "[obj] warning: cone at line " << c.line << " is on a vertex used by no facet; skipped\n";
                    continue;
                }
                v = VIdx(Mesh::size_type(newIndex[std::size_t(obj)]));
            }
            if (isCone[v])
                std::cerr << "[obj] warning: vertex " << c.index << " has several cone records; using line " << c.line << "\n";
            isCone[v] = 1;
            coneAngle[v] = c.angleOverPi * kPiMesh;    // multiples of pi -> radians
            ++stored;
        }
        if (options.verbose)
            std::cout << "[obj] " << stored << " cone(s) attached to the CGAL mesh (" << meta.coneIndexBaseUsed
                      << " indices, property maps " << kIsConeProperty << " / " << kConeAngleProperty << ")\n";
    }

    // Faces of a mesh built from a valid soup keep the soup order -> keep the UVs.
    out.faceSoupVertices.resize(polys.size());
    for (std::size_t f = 0; f < polys.size(); ++f)
        for (int k = 0; k < 3; ++k) out.faceSoupVertices[f][k] = int(polys[f][k]);
    if (uvValid && out.mesh.number_of_faces() == polyUV.size()) {
        out.hasCornerUV = true;
        out.faceCornerUV.resize(polyUV.size());
        for (std::size_t f = 0; f < polyUV.size(); ++f)
            for (int k = 0; k < 3; ++k) out.faceCornerUV[f][k] = uvs[polyUV[f][k]];
    }
    return out;
}

Topology computeTopology(const Mesh& mesh) {
    Topology t;
    t.nV = int(mesh.number_of_vertices());
    t.nE = int(mesh.number_of_edges());
    t.nF = int(mesh.number_of_faces());
    // count boundary loops by walking border halfedges
    std::unordered_set<std::size_t> seen;
    for (HIdx h : mesh.halfedges()) {
        if (!mesh.is_border(h) || seen.count(std::size_t(h))) continue;
        ++t.boundaryLoops;
        HIdx c = h;
        do {
            seen.insert(std::size_t(c));
            c = mesh.next(c);
        } while (c != h);
    }
    t.eulerChar = t.nV - t.nE + t.nF;
    t.genus = (2 - t.eulerChar - t.boundaryLoops) / 2;
    return t;
}

void writeObjWithUV(const std::string& path, const Mesh& mesh,
                    const std::vector<std::array<double, 2>>& cornerUV) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out.precision(17);
    out << "# Bounded distortion parametrization (Chien, Levi, Weber 2016)\n";
    for (VIdx v : mesh.vertices()) {
        const auto& p = mesh.point(v);
        out << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
    }
    for (const auto& uv : cornerUV) out << "vt " << uv[0] << ' ' << uv[1] << '\n';
    // cone singularities: "c <vertex (1-based)> <cone angle / pi>"
    const auto isCone = mesh.property_map<VIdx, unsigned char>(kIsConeProperty);
    const auto angle  = mesh.property_map<VIdx, double>(kConeAngleProperty);
    if (isCone && angle)
        for (VIdx v : mesh.vertices())
            if ((*isCone)[v])
                out << "c " << std::size_t(v) + 1 << ' ' << (*angle)[v] / 3.14159265358979323846 << '\n';
    for (FIdx f : mesh.faces()) {
        out << 'f';
        HIdx h = mesh.halfedge(f);
        for (int k = 0; k < 3; ++k, h = mesh.next(h))
            out << ' ' << std::size_t(mesh.source(h)) + 1 << '/' << 3 * std::size_t(f) + k + 1;
        out << '\n';
    }
}

}  // namespace bdp
