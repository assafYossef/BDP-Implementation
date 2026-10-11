// Geometry.cpp - angles, curvature and distortion of a discrete metric.
#include "bdp/Geometry.h"
#include "bdp/MetricData.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>

namespace bdp {

Eigen::MatrixXd cornerAngles(const Eigen::MatrixXi& F2E, const Eigen::VectorXd& e) {
    const Eigen::Index nF = F2E.rows();
    Eigen::MatrixXd th(nF, 3);
    for (Eigen::Index f = 0; f < nF; ++f) {
        const double s[3] = {e(F2E(f, 0)), e(F2E(f, 1)), e(F2E(f, 2))};
        const double T = std::sqrt(std::max(0.0, heronT2(s[0], s[1], s[2])));   // 4 * Area
        for (int k = 0; k < 3; ++k)
            th(f, k) = std::atan2(T, s[(k + 1) % 3] + s[(k + 2) % 3] - s[k]);
    }
    return th;
}

Eigen::VectorXd vertexCurvature(const Eigen::MatrixXi& F, const Eigen::MatrixXi& F2E,
                                const Eigen::VectorXd& e, const Eigen::VectorXi& isBnd) {
    const Eigen::MatrixXd th = cornerAngles(F2E, e);
    Eigen::VectorXd K(isBnd.size());
    for (Eigen::Index v = 0; v < K.size(); ++v) K(v) = isBnd(v) ? kPi : 2.0 * kPi;
    for (Eigen::Index f = 0; f < F.rows(); ++f)
        for (int k = 0; k < 3; ++k) K(F(f, k)) -= th(f, k);
    return K;
}

Distortion computeDistortion(const MetricData& d, const Eigen::VectorXd& e) {
    const Eigen::Index nF = d.F.rows();
    Distortion D;
    D.sigmaA.resize(nF); D.sigmaB.resize(nF); D.tau.resize(nF); D.sqrtTau.resize(nF);
    D.K.resize(nF); D.k.resize(nF); D.logArea.resize(nF); D.symDir.resize(nF);
    for (Eigen::Index f = 0; f < nF; ++f) {
        double A = 0, v = 0, w = 0;
        for (int k = 0; k < 3; ++k) {
            const double s = e(d.F2E(f, k));
            A += d.cA(f, k) * s;   // A        (linear in s, Section 4.3)
            v += d.cV(f, k) * s;   // 2 Re(eta)
            w += d.cW(f, k) * s;   // 2 Im(eta)
        }
        const double r   = std::hypot(v, w);              // 2|eta|
        const double sa2 = A + r;                         // sigma_a^2   (Eq. 1)
        const double sb2 = A - r;                         // sigma_b^2   (Eq. 2)
        const double inf = std::numeric_limits<double>::infinity();
        D.sigmaA(f) = std::sqrt(std::max(sa2, 0.0));
        D.sigmaB(f) = std::sqrt(std::max(sb2, 0.0));
        D.tau(f) = sb2 > 0 ? std::max(sa2, 1.0 / sb2) : inf;
        D.sqrtTau(f) = std::sqrt(D.tau(f));
        D.K(f) = sb2 > 0 ? std::sqrt(sa2 / sb2) : inf;    // K^2 = sigma_a^2/sigma_b^2 (Eq. 3)
        D.k(f) = sb2 > 0 ? (D.K(f) - 1.0) / (D.K(f) + 1.0) : 1.0;
        D.logArea(f) = sb2 > 0 ? 0.5 * std::log(sa2 * sb2) : -inf;
        // symmetric Dirichlet  s1^2 + s2^2 + 1/s1^2 + 1/s2^2 = tr(M) + tr(M^-1) = 2A (1 + 1/det M)
        D.symDir(f) = sb2 > 0 ? sa2 + sb2 + 1.0 / sa2 + 1.0 / sb2 : inf;
    }
    return D;
}

QualityMetrics computeQuality(const MetricData& d, const Eigen::VectorXd& e,
                              const Distortion& D, int flippedTriangles) {
    QualityMetrics q;
    const Eigen::Index nF = d.F.rows();
    // ---- validity: every triangle must be a real triangle in the new metric
    for (Eigen::Index f = 0; f < nF; ++f) {
        const double T2 = heronT2(e(d.F2E(f, 0)), e(d.F2E(f, 1)), e(d.F2E(f, 2)));
        if (!(T2 > 0) || !(D.sigmaB(f) > 0) || !std::isfinite(D.sigmaA(f))) ++q.degenerateTriangles;
    }
    q.flippedTriangles = flippedTriangles;
    q.valid = q.degenerateTriangles == 0 && flippedTriangles <= 0;

    // ---- stretch measures (s1 = sigma_a, s2 = sigma_b)
    const double s1max = D.sigmaA.maxCoeff();
    const double s2min = D.sigmaB.minCoeff();
    q.tau = 0;
    q.Kmax = 0;
    for (Eigen::Index f = 0; f < nF; ++f) {
        const double s1 = D.sigmaA(f), s2 = D.sigmaB(f);
        q.tau  = std::max(q.tau, s2 > 0 ? std::max(s1, 1.0 / s2) : std::numeric_limits<double>::infinity());
        q.Kmax = std::max(q.Kmax, s2 > 0 ? s1 / s2 : std::numeric_limits<double>::infinity());
    }
    // Scaling the whole map by c multiplies every s1, s2 by c. max(c s1max, 1/(c s2min))
    // is smallest when both terms are equal: c* = 1/sqrt(s1max s2min), tauSc = sqrt(s1max/s2min).
    if (s2min > 0) {
        q.bestScale = 1.0 / std::sqrt(s1max * s2min);
        q.tauSc = std::sqrt(s1max / s2min);
    } else {
        q.bestScale = 1.0;
        q.tauSc = std::numeric_limits<double>::infinity();
    }
    // ---- area-weighted symmetric Dirichlet (weights = original area fractions, sum 1)
    q.symDir = d.wF.dot(D.symDir);
    q.symDirExcess = q.symDir - 4.0;

    // ---- curvature error on the constrained vertices
    const Eigen::VectorXd K = vertexCurvature(d.F, d.F2E, e, d.isBnd);
    double err = 0;
    for (Eigen::Index i = 0; i < d.Vc.size(); ++i) err = std::max(err, std::abs(K(d.Vc(i)) - d.K0(i)));
    q.maxCurvErrDeg = err * 180.0 / kPi;
    return q;
}

std::string formatQuality(const QualityMetrics& q) {
    std::ostringstream o;
    o << std::setprecision(6);
    o << "solution quality (per triangle: s1 >= s2 = stretches, 1 = unchanged)\n";
    o << "  tau        = " << std::setw(12) << q.tau   << "   worst stretch max(s1, 1/s2)                 (perfect 1)\n";
    o << "  tauSc      = " << std::setw(12) << q.tauSc << "   same after optimal global scale x" << q.bestScale << "  (perfect 1)\n";
    o << "  Kmax       = " << std::setw(12) << q.Kmax  << "   worst shape change max s1/s2                (perfect 1)\n";
    o << "  symDir     = " << std::setw(12) << q.symDir << "   area-weighted mean s1^2+s2^2+1/s1^2+1/s2^2 (perfect 4), excess "
      << q.symDirExcess << "\n";
    o << "  max|K-K*|  = " << std::setw(12) << q.maxCurvErrDeg << "   worst curvature error in degrees          (perfect 0)\n";
    o << "  valid      = " << std::setw(12) << (q.valid ? 1 : 0) << "   " << q.degenerateTriangles
      << " degenerate triangles, ";
    if (q.flippedTriangles >= 0) o << q.flippedTriangles << " flipped in the layout";
    else o << "layout not checked";
    o << "   (perfect 1)\n";
    return o.str();
}

}  // namespace bdp
