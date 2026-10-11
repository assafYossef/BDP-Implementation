// main.cpp - Bounded Distortion Parametrization in the Space of Metrics
//            (Chien, Levi, Weber - SIGGRAPH Asia 2016).
//
// Pipeline
//   1. read an .obj file and build the CGAL halfedge mesh (topology),
//   2. prescribe the target curvature (flat everywhere except at cones),
//   3. translate topology + geometry into Eigen arrays: squared edge lengths e0,
//      face->edge incidences and the per-face linear maps s -> (A, Re eta, Im eta),
//   4. push them into MATLAB through the Engine API; MATLAB runs Sequential
//      Convex Programming with a custom Newton (log-barrier) solver and returns
//      the optimised squared edge lengths e,
//   5. lay the flat cone metric out in the plane and visualise everything in
//      Polyscope (input mesh, distortion heat maps, flattened metric).
#include "bdp/Cones.h"
#include "bdp/Geometry.h"
#include "bdp/Layout.h"
#include "bdp/MatlabBridge.h"
#include "bdp/Mesh.h"
#include "bdp/MetricData.h"
#ifdef BDP_HAVE_POLYSCOPE
#include "bdp/Viewer.h"
#endif

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#ifndef BDP_MATLAB_DIR
#define BDP_MATLAB_DIR "matlab"
#endif

using namespace bdp;

namespace {

struct Args {
    std::string mesh;
    std::string cones;
    bool conesOneBased = false;
    ObjReadOptions obj;
    std::string out;
    std::string dumpDir;
    std::string loadMetric;
    std::string matlabDir = BDP_MATLAB_DIR;
    bool matlabShared = false;
    bool gui = true;
    bool headless = false;
    std::string screenshot;
    SolverOptions solver;
};

void usage() {
    std::cout <<
R"(usage: bdp <mesh.obj> [options]

  --variant V1R|V1|V2|V3|V4|V5  energy (paper Sec. 8.1), default V1R
                               V1R max tau + 0.01 * average (best worst case; = reference code)
                               V1 max isometric tau      V2 area-weighted beta
                               V3 Green-Lagrange          V4 max |eta| (conformal)
                               V5 projection of the mesh's own UVs (vt) to bounded distortion
  --avg-weight <w>           V1R: weight of the average term (default 0.01)
  --solver auto|mosek|cvx|coneprog|newton   convex subproblem solver. auto (default): MOSEK
                               called directly (also found inside the CVX folder), else CVX with
                               its default solver, else MATLAB coneprog, else built-in Newton
  --Ca <x>  --Cb <x>  --CK <x>  distortion bounds sigma_a^2<=Ca, sigma_b^2>=Cb, K^2<=CK
                               (default: variant's values; 'inf' / 0 disables a bound)
  --cones <file>             cone file: "<vertex> <curvature in degrees>" per line (0-based);
                             overrides the `c` records of the OBJ file
  --cones-one-based          vertex indices in the cone file are 1-based
  --cone-index-base auto|0|1 index base of the OBJ `c <vertex> <angle/pi>` records (default auto)
  --max-iters <n>            SCP iterations (default 500)
  --no-cetm                  disable the CETM acceleration (Sec. 8.2)
  --cetm-threshold <deg>     L_inf curvature error that triggers CETM (default 1 deg)
  --constraint auto|hard|soft  linearised curvature as equality (hard), as the paper's
                               L_inf penalty (soft), or hard with soft fallback (auto)
  --out <file.obj>           output mesh with UVs (default <input>_bdp.obj)
  --dump <dir>               write the problem matrices as text (for bdp_run_from_dump.m)
  --load-metric <file>       skip the optimisation; read squared edge lengths (one per line)
  --matlab-dir <dir>         folder with bdp_solve.m (default: the project's matlab/)
  --matlab-shared            attach to a shared MATLAB session (matlab.engine.shareEngine)
  --no-gui                   do not open the Polyscope viewer
  --screenshot <file.png>    render the viewer once into an image instead of opening a window
  --headless                 use Polyscope's EGL backend (no display needed)
  --quiet                    less solver output
)";
}

double parseBound(const std::string& s) {
    if (s == "inf" || s == "Inf" || s == "INF") return std::numeric_limits<double>::infinity();
    return std::stod(s);
}

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string k = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + k);
            return argv[++i];
        };
        if (k == "-h" || k == "--help") { usage(); std::exit(0); }
        else if (k == "--variant") a.solver.variant = next();
        else if (k == "--avg-weight") a.solver.avgWeight = std::stod(next());
        else if (k == "--solver") a.solver.solver = next();
        else if (k == "--Ca") a.solver.Ca = parseBound(next());
        else if (k == "--Cb") a.solver.Cb = parseBound(next());
        else if (k == "--CK") a.solver.CK = parseBound(next());
        else if (k == "--cones") a.cones = next();
        else if (k == "--cones-one-based") a.conesOneBased = true;
        else if (k == "--cone-index-base") {
            const std::string b = next();
            if (b == "0") a.obj.coneIndexBase = ConeIndexBase::Zero;
            else if (b == "1") a.obj.coneIndexBase = ConeIndexBase::One;
            else if (b == "auto") a.obj.coneIndexBase = ConeIndexBase::Auto;
            else throw std::runtime_error("--cone-index-base expects auto, 0 or 1");
        }
        else if (k == "--max-iters") a.solver.maxIters = std::stoi(next());
        else if (k == "--no-cetm") a.solver.useCETM = false;
        else if (k == "--cetm-threshold") a.solver.cetmThresholdDeg = std::stod(next());
        else if (k == "--constraint") a.solver.curvatureConstraint = next();
        else if (k == "--out") a.out = next();
        else if (k == "--dump") a.dumpDir = next();
        else if (k == "--load-metric") a.loadMetric = next();
        else if (k == "--matlab-dir") a.matlabDir = next();
        else if (k == "--matlab-shared") a.matlabShared = true;
        else if (k == "--no-gui") a.gui = false;
        else if (k == "--screenshot") a.screenshot = next();
        else if (k == "--headless") a.headless = true;
        else if (k == "--quiet") a.solver.verbose = false;
        else if (!k.empty() && k[0] == '-') throw std::runtime_error("unknown option " + k);
        else a.mesh = k;
    }
    if (a.mesh.empty()) { usage(); throw std::runtime_error("no input mesh given"); }
    if (a.out.empty()) {
        std::filesystem::path p(a.mesh);
        a.out = (p.parent_path() / (p.stem().string() + "_bdp.obj")).string();
    }
    return a;
}

Eigen::VectorXd readVector(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<double> v;
    double x;
    while (in >> x) v.push_back(x);
    return Eigen::Map<Eigen::VectorXd>(v.data(), Eigen::Index(v.size()));
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Args args = parseArgs(argc, argv);

        // ---------------- 1) mesh & topology (CGAL) ----------------------------
        LoadedMesh lm = loadObj(args.mesh, args.obj);   // also attaches `c` cones as property maps
        const Topology topo = computeTopology(lm.mesh);
        std::cout << "[mesh] " << args.mesh << ": |V|=" << topo.nV << " |E|=" << topo.nE
                  << " |F|=" << topo.nF << "  chi=" << topo.eulerChar << " genus=" << topo.genus
                  << " boundary loops=" << topo.boundaryLoops << "\n";

        // ---------------- 2) CGAL -> Eigen: the metric problem -------------------
        MetricData data = buildMetricData(lm);
        std::cout << "[data] mean edge length " << data.scale << " (normalised to 1)"
                  << (data.Eref.size() ? ", reference UVs found (usable by V5)" : "") << "\n";
        if (args.solver.variant == "V5" && data.Eref.size() == 0)
            std::cout << "[data] V5 without usable vt coordinates: the reference map is the identity\n";

        // ---------------- 3) curvature prescription (Section 5) -----------------
        // Every cone source ends up in the CGAL vertex property maps of the mesh;
        // priority: --cones file  >  `c` records of the OBJ  >  automatic placement.
        if (!args.cones.empty()) {
            if (countCones(lm.mesh) > 0)
                std::cout << "[cones] --cones given: replacing the " << countCones(lm.mesh)
                          << " cone(s) of the OBJ file\n";
            setConesOnMesh(lm.mesh, loadConeFile(args.cones, lm, args.conesOneBased));
        } else if (countCones(lm.mesh) > 0) {
            std::cout << "[cones] using the " << countCones(lm.mesh) << " cone(s) of the OBJ file\n";
        } else {
            setConesOnMesh(lm.mesh, automaticCones(data, topo));
        }
        readConeProperties(lm.mesh, data);            // CGAL property maps -> Eigen vectors
        assignCurvatureConstraints(data, topo);       // Eigen cone vectors -> Vc, K0

        if (!args.dumpDir.empty()) dumpMetricData(data, args.dumpDir, args.solver.variant);

        // ---------------- 4) optimisation in MATLAB -----------------------------
        Eigen::VectorXd e;
        Eigen::VectorXd history;
        if (!args.loadMetric.empty()) {
            e = readVector(args.loadMetric);
            if (e.size() != data.e0.size()) throw std::runtime_error("metric file has wrong length");
            std::cout << "[solve] metric loaded from " << args.loadMetric << "\n";
        } else {
#ifdef BDP_HAVE_MATLAB
            MatlabBridge matlab(args.matlabDir, args.matlabShared);
            const SolverResult r = matlab.solve(data, args.solver);
            e = r.e;
            history = r.curvErrHistory;
            std::cout << "[solve] " << r.message << " in " << r.seconds << " s (" << r.solverUsed << ")"
                      << (r.usedCETM ? " (CETM finish)" : "") << "\n";
#else
            if (!args.dumpDir.empty()) {
                std::cout << "[solve] built without MATLAB: run bdp_run_from_dump('" << args.dumpDir
                          << "') in MATLAB, then rerun with --load-metric " << args.dumpDir << "/e_opt.txt\n";
                return 0;
            }
            throw std::runtime_error("built without the MATLAB Engine: use --dump / --load-metric");
#endif
        }

        // ---------------- 5) evaluation, layout, export ---------------------------
        const Eigen::VectorXd K_in  = vertexCurvature(data.F, data.F2E, data.e0, data.isBnd);
        const Eigen::VectorXd K_out = vertexCurvature(data.F, data.F2E, e, data.isBnd);
        Eigen::VectorXd errIn = Eigen::VectorXd::Zero(topo.nV), errOut = Eigen::VectorXd::Zero(topo.nV);
        for (int i = 0; i < data.Vc.size(); ++i) {
            errIn(data.Vc(i))  = std::abs(K_in(data.Vc(i)) - data.K0(i));
            errOut(data.Vc(i)) = std::abs(K_out(data.Vc(i)) - data.K0(i));
        }
        const Distortion D = computeDistortion(data, e);
        const LayoutResult L = layoutMetric(lm.mesh, data, e);
        std::cout << std::setprecision(6)
                  << "[result] max |K-K0| = " << errOut.maxCoeff() * 180 / kPi << " deg"
                  << "  max sqrt(tau) = " << D.sqrtTau.maxCoeff()
                  << "  max k = " << D.k.maxCoeff()
                  << "  sigma in [" << D.sigmaB.minCoeff() << ", " << D.sigmaA.maxCoeff() << "]"
                  << "  flipped = " << L.flippedTriangles << "\n";
        writeObjWithUV(args.out, lm.mesh, L.cornerUV);
        std::cout << "[result] parametrized mesh written to " << args.out << "\n";

        // ---------------- quality of the solution ----------------------------------
        const QualityMetrics Q = computeQuality(data, e, D, L.flippedTriangles);
        const std::string report = formatQuality(Q);
        std::cout << "\n" << report << "\n";
        {
            std::filesystem::path qp(args.out);
            qp.replace_filename(qp.stem().string() + "_quality.txt");
            std::ofstream qf(qp);
            qf << "mesh    : " << args.mesh << "\nvariant : " << args.solver.variant << "\n\n" << report;
            std::cout << "[result] quality report written to " << qp.string() << "\n";
        }

        // ---------------- 6) visualisation ------------------------------------------
        if (args.gui) {
#ifdef BDP_HAVE_POLYSCOPE
            ViewerInput vin;
            vin.data = &data; vin.eOpt = &e; vin.distortion = &D; vin.layout = &L;
            vin.curvatureErrorInitial = errIn; vin.curvatureErrorFinal = errOut;
            vin.curvErrHistory = history;
            vin.quality = Q;
            vin.title = "BDP " + args.solver.variant + " - " + std::filesystem::path(args.mesh).filename().string();
            vin.screenshot = args.screenshot;
            vin.headless = args.headless;
            showViewer(vin);
#else
            std::cout << "[viewer] built without Polyscope\n";
#endif
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 1;
    }
}
