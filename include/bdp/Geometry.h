// Geometry.h - evaluation of a discrete metric (squared edge lengths) in C++:
// corner angles, angle-defect curvature and per-triangle distortion. These are
// cheap O(|F|) evaluations used for reporting and visualisation; all heavy
// optimisation runs in MATLAB.
#pragma once

#include <Eigen/Core>

#include <string>

namespace bdp {

struct MetricData;

constexpr double kPi = 3.14159265358979323846;

// T^2 = 2ab + 2bc + 2ca - a^2 - b^2 - c^2 = (4 Area)^2 for squared lengths a, b, c.
inline double heronT2(double a, double b, double c) {
    return 2.0 * (a * b + b * c + c * a) - a * a - b * b - c * c;
}

// Corner angles (|F| x 3) by the law of cosines in squared lengths (Eq. (4)):
//   theta_k = atan2(T, s_i + s_j - s_k).
Eigen::MatrixXd cornerAngles(const Eigen::MatrixXi& F2E, const Eigen::VectorXd& e);

// Discrete Gaussian curvature (Section 5):  K(v) = 2pi - sum theta (interior),
//                                            K(v) =  pi - sum theta (boundary).
Eigen::VectorXd vertexCurvature(const Eigen::MatrixXi& F, const Eigen::MatrixXi& F2E,
                                const Eigen::VectorXd& e, const Eigen::VectorXi& isBnd);

// Per-face distortion of metric e relative to the original metric (Section 3.2.1).
struct Distortion {
    Eigen::VectorXd sigmaA, sigmaB;   // singular values of the Jacobian
    Eigen::VectorXd tau;              // max(sigma_a^2, 1/sigma_b^2)
    Eigen::VectorXd sqrtTau;          // plotted in the paper's figures
    Eigen::VectorXd K;                // conformal distortion sigma_a / sigma_b
    Eigen::VectorXd k;                // small dilatation (K-1)/(K+1) in [0,1)
    Eigen::VectorXd logArea;          // log(sigma_a sigma_b) = log area ratio
    Eigen::VectorXd symDir;           // s1^2 + s2^2 + 1/s1^2 + 1/s2^2 (symmetric Dirichlet, 4 = isometric)
};
Distortion computeDistortion(const MetricData& d, const Eigen::VectorXd& e);

// Quality summary of a solution. Everything derives from the two stretches of
// each triangle, s1 = sigma_a >= s2 = sigma_b (singular values of the Jacobian;
// 1 = unchanged, 1.3 = stretched by 30%, 0.8 = squashed to 80%).
//
//   tau         max_i max(s1, 1/s2)                 worst stretch              perfect 1
//               (= sqrt of the paper's tau = max(sigma_a^2, 1/sigma_b^2))
//   tauSc       min_c max_i max(c s1, 1/(c s2))     worst stretch after the    perfect 1
//               = sqrt(max s1 / min s2),            best global rescaling
//               attained at c* = 1/sqrt(max s1 min s2)   (free fix: tauSc <= tau)
//   Kmax        max_i s1/s2                         worst shape change         perfect 1
//   symDir      sum_i w_i (s1^2+s2^2+1/s1^2+1/s2^2) area-weighted average      perfect 4
//               (w_i = original area fraction);     symDirExcess = symDir - 4  perfect 0
//   maxCurvErr  max_{v in Vc} |K(v) - K*(v)| (deg)  flat where it should be    perfect 0
//   valid       every triangle is a real triangle   (strict triangle inequality in the
//               metric, s2 > 0) and no triangle is flipped in the planar layout   perfect 1
struct QualityMetrics {
    double tau = 0, tauSc = 0, bestScale = 1, Kmax = 0;
    double symDir = 0, symDirExcess = 0;
    double maxCurvErrDeg = 0;
    bool valid = false;
    int degenerateTriangles = 0;   // triangle inequality violated / zero area in the metric
    int flippedTriangles = 0;      // non-positive orientation in the layout (-1 = not checked)
};
QualityMetrics computeQuality(const MetricData& d, const Eigen::VectorXd& e,
                              const Distortion& D, int flippedTriangles = -1);

// Human-readable table (also written to <output>_quality.txt).
std::string formatQuality(const QualityMetrics& q);

}  // namespace bdp
