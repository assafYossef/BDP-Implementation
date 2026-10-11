// MatlabBridge.h - inter-process communication with MATLAB through the
// MATLAB Engine API for C++ (MatlabEngine.hpp / MatlabDataArray.hpp).
//
// Division of labour (see README):
//   C++    : builds the topology and every matrix of the problem (MetricData),
//            pushes them into the MATLAB workspace, receives the result.
//   MATLAB : runs the SCP loop with the custom Newton barrier solver
//            (matlab/bdp_solve.m) - sparse factorisations and linear systems.
#pragma once

#include "bdp/MetricData.h"

#include <Eigen/Core>
#include <memory>
#include <string>

namespace bdp {

struct SolverOptions {
    std::string variant = "V1R";      // V1R (default), V1 ... V5 (Section 8.1)
    double avgWeight = 0.01;          // V1R: weight of the average term
    std::string solver = "auto";      // auto | mosek | cvx | coneprog | newton
    double Ca = -1, Cb = -1, CK = -1; // < 0: variant default; use "inf" to disable
    int maxIters = 500;
    bool useCETM = true;
    double cetmThresholdDeg = 1.0;
    std::string curvatureConstraint = "auto";   // auto | hard | soft
    bool verbose = true;
};

struct SolverResult {
    Eigen::VectorXd e;               // optimised squared edge lengths (normalised units)
    Eigen::VectorXd curvErrHistory;  // ||K - K0||_inf per SCP iteration
    Eigen::VectorXd energyHistory;
    bool converged = false;
    bool usedCETM = false;
    double seconds = 0;
    std::string message;
    std::string solverUsed;          // mosek / coneprog / newton
};

#ifdef BDP_HAVE_MATLAB
class MatlabBridge {
public:
    // Starts a new MATLAB process, or attaches to a shared session
    // (started in MATLAB with `matlab.engine.shareEngine`) if `shared` is true.
    MatlabBridge(const std::string& matlabCodeDir, bool shared);
    ~MatlabBridge();

    // ---- workspace transfer (Eigen is column major, like MATLAB) ----------
    void put(const std::string& name, const Eigen::MatrixXd& M);
    void put(const std::string& name, const Eigen::MatrixXi& M, int offset = 0);
    void putScalar(const std::string& name, double v);
    void putString(const std::string& name, const std::string& s);
    Eigen::MatrixXd getMatrix(const std::string& name);
    double getScalar(const std::string& name);
    std::string getString(const std::string& name);

    // Evaluates a MATLAB statement; MATLAB's console output is forwarded to stdout.
    void eval(const std::string& statement);

    // Full pipeline: push MetricData, call bdp_solve, fetch e and the history.
    SolverResult solve(const MetricData& d, const SolverOptions& o);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
#endif

}  // namespace bdp
