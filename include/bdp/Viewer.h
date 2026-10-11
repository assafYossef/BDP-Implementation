// Viewer.h - Polyscope visualisation of the input mesh, the distortion heat maps
// and the flattened metric (the planar layout of the cut surface).
#pragma once

#include "bdp/Geometry.h"
#include "bdp/Layout.h"
#include "bdp/MetricData.h"

#include <string>

namespace bdp {

struct ViewerInput {
    const MetricData* data = nullptr;
    const Eigen::VectorXd* eOpt = nullptr;     // optimised squared lengths
    const Distortion* distortion = nullptr;
    const LayoutResult* layout = nullptr;
    Eigen::VectorXd curvatureErrorInitial;     // |K - K0| per vertex, input metric
    Eigen::VectorXd curvatureErrorFinal;       // |K - K0| per vertex, optimised metric
    Eigen::VectorXd curvErrHistory;            // SCP convergence (for the UI panel)
    QualityMetrics quality;                    // tau, tauSc, Kmax, symDir, max|K-K*|, valid
    std::string title;
    std::string screenshot;    // if set: render one frame to this .png and return (no window loop)
    bool headless = false;     // use the EGL backend (no display needed, e.g. for screenshots)
};

void showViewer(const ViewerInput& in);

}  // namespace bdp
