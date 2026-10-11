// Viewer.cpp - Polyscope visualisation.
#include "bdp/Viewer.h"

#include "imgui.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <algorithm>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <vector>

namespace bdp {

void showViewer(const ViewerInput& in) {
    const MetricData& d = *in.data;
    const Distortion& D = *in.distortion;
    const LayoutResult& L = *in.layout;
    const int nF = int(d.F.rows());

    polyscope::options::programName = "Bounded Distortion Parametrization";
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::init(in.headless ? "openGL3_egl" : "");

    // ---------------- 1) the original 3D mesh -------------------------------------
    auto* surf = polyscope::registerSurfaceMesh("input mesh", d.V, d.F);
    // per-corner UVs of the optimised flat cone metric -> checkerboard texture
    Eigen::MatrixXd cornerUV(3 * nF, 2);
    for (int c = 0; c < 3 * nF; ++c) cornerUV.row(c) << L.cornerUV[c][0], L.cornerUV[c][1];
    auto* pq = surf->addParameterizationQuantity("parametrization (flat cone metric)", cornerUV,
                                                 polyscope::ParamCoordsType::WORLD);
    pq->setStyle(polyscope::ParamVizStyle::CHECKER);
    pq->setCheckerSize(d.scale);   // about one mean edge length per checker
    pq->setEnabled(true);

    // Heat maps as in the paper (Sec. 9.1): dark blue = no distortion,
    // red = sqrt(tau) >= 4  /  k >= 0.4.
    surf->addFaceScalarQuantity("isometric distortion sqrt(tau)", D.sqrtTau)
        ->setColorMap("jet")->setMapRange({1.0, 4.0});
    surf->addFaceScalarQuantity("conformal distortion k", D.k)
        ->setColorMap("jet")->setMapRange({0.0, 0.4});
    surf->addFaceScalarQuantity("s1 = sigma_a (max stretch)", D.sigmaA)->setColorMap("jet");
    surf->addFaceScalarQuantity("s2 = sigma_b (min stretch)", D.sigmaB)->setColorMap("jet");
    surf->addFaceScalarQuantity("symmetric Dirichlet (4 = isometric)", D.symDir)
        ->setColorMap("jet")->setMapRange({4.0, 8.0});
    surf->addFaceScalarQuantity("log area ratio", D.logArea, polyscope::DataType::SYMMETRIC)
        ->setColorMap("coolwarm");
    surf->addVertexScalarQuantity("|K - K0| input metric (rad)", in.curvatureErrorInitial)->setColorMap("reds");
    surf->addVertexScalarQuantity("|K - K0| optimised metric (rad)", in.curvatureErrorFinal)->setColorMap("reds");

    // ---------------- 2) cones --------------------------------------------------------
    std::vector<int> coneIdx;                          // vertices flagged in the cone property map
    for (int v = 0; v < d.isCone.size(); ++v)
        if (d.isCone(v)) coneIdx.push_back(v);
    if (!coneIdx.empty()) {
        Eigen::MatrixXd P(coneIdx.size(), 3);
        Eigen::VectorXd angleOverPi(coneIdx.size()), Kdeg(coneIdx.size());
        for (std::size_t i = 0; i < coneIdx.size(); ++i) {
            const int v = coneIdx[i];
            P.row(i) = d.V.row(v);
            angleOverPi(i) = d.coneAngle(v) / kPi;
            Kdeg(i) = ((d.isBnd(v) ? kPi : 2.0 * kPi) - d.coneAngle(v)) * 180.0 / kPi;
        }
        auto* pc = polyscope::registerPointCloud("cones", P);
        pc->setPointRadius(0.6 * d.scale, false);
        pc->addScalarQuantity("cone curvature (deg)", Kdeg, polyscope::DataType::SYMMETRIC)
            ->setColorMap("coolwarm")->setEnabled(true);
        pc->addScalarQuantity("cone angle / pi", angleOverPi);
    }

    // ---------------- 3) the flattened metric: planar layout of the cut mesh ---------
    Eigen::MatrixXd UV2(3 * nF, 2);
    Eigen::MatrixXi Fuv(nF, 3);
    for (int f = 0; f < nF; ++f)
        for (int k = 0; k < 3; ++k) {
            UV2.row(3 * f + k) << L.cornerUV[3 * f + k][0], L.cornerUV[3 * f + k][1];
            Fuv(f, k) = 3 * f + k;
        }
    // shift the layout next to the 3D model
    const double xmax3 = d.V.col(0).maxCoeff();
    const double umin = UV2.col(0).minCoeff();
    UV2.col(0).array() += xmax3 - umin + 0.25 * (d.V.colwise().maxCoeff() - d.V.colwise().minCoeff()).norm();
    UV2.col(1).array() -= UV2.col(1).mean() - d.V.col(1).mean();
    auto* flat = polyscope::registerSurfaceMesh2D("flattened metric (UV layout)", UV2, Fuv);
    flat->setEdgeWidth(0.5);
    flat->addFaceScalarQuantity("isometric distortion sqrt(tau)", D.sqrtTau)
        ->setColorMap("jet")->setMapRange({1.0, 4.0})->setEnabled(true);
    flat->addFaceScalarQuantity("conformal distortion k", D.k)->setColorMap("jet")->setMapRange({0.0, 0.4});

    // ---------------- 4) summary panel ------------------------------------------------
    const QualityMetrics Q = in.quality;
    std::vector<float> hist(in.curvErrHistory.size());
    for (int i = 0; i < int(hist.size()); ++i) hist[i] = float(std::log10(std::max(1e-16, in.curvErrHistory(i))));
    polyscope::state::userCallback = [=]() {
        ImGui::TextUnformatted(in.title.c_str());
        ImGui::Separator();
        ImGui::TextUnformatted("solution quality          value     perfect");
        ImGui::Text("tau    worst stretch    %10.4f    1", Q.tau);
        ImGui::Text("tauSc  after best scale %10.4f    1   (scale x%.4g)", Q.tauSc, Q.bestScale);
        ImGui::Text("Kmax   worst shape      %10.4f    1", Q.Kmax);
        ImGui::Text("symDir area-weighted    %10.4f    4   (excess %.3g)", Q.symDir, Q.symDirExcess);
        ImGui::Text("max|K-K*| (deg)         %10.3e    0", Q.maxCurvErrDeg);
        ImGui::Text("valid                   %10d    1", Q.valid ? 1 : 0);
        ImGui::Separator();
        ImGui::Text("degenerate triangles : %d", Q.degenerateTriangles);
        ImGui::Text("flipped triangles    : %d", L.flippedTriangles);
        ImGui::Text("layout length error  : %.2e", L.maxEdgeLengthError);
        if (!hist.empty())
            ImGui::PlotLines("log10 |K-K0|_inf", hist.data(), int(hist.size()), 0, nullptr,
                             FLT_MAX, FLT_MAX, ImVec2(0, 80));
    };
    if (!in.screenshot.empty()) {
        polyscope::screenshot(in.screenshot, false);
        std::cout << "[viewer] screenshot written to " << in.screenshot << "\n";
        return;
    }
    polyscope::show();
}

}  // namespace bdp
