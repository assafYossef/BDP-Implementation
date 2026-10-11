// MatlabBridge.cpp - MATLAB Engine API for C++ integration.
#include "bdp/MatlabBridge.h"

#ifdef BDP_HAVE_MATLAB

#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <streambuf>
#include <vector>

namespace bdp {

namespace me = matlab::engine;
namespace md = matlab::data;

namespace {

// Forwards MATLAB's UTF-16 console output (fprintf in bdp_solve) to a std::ostream
// so the progress of the SCP iterations is visible from the C++ process.
class ForwardingBuffer : public std::basic_streambuf<char16_t> {
public:
    explicit ForwardingBuffer(std::ostream& os) : os_(os) {}
protected:
    int_type overflow(int_type ch) override {
        if (!traits_type::eq_int_type(ch, traits_type::eof())) put(char16_t(ch));
        return traits_type::not_eof(ch);
    }
    std::streamsize xsputn(const char16_t* s, std::streamsize n) override {
        for (std::streamsize i = 0; i < n; ++i) put(s[i]);
        return n;
    }
private:
    void put(char16_t c) {
        os_.put(c < 128 ? char(c) : '?');
        if (c == u'\n') os_.flush();
    }
    std::ostream& os_;
};

std::u16string u16(const std::string& s) { return me::convertUTF8StringToUTF16String(s); }

}  // namespace

struct MatlabBridge::Impl {
    std::unique_ptr<me::MATLABEngine> engine;
    md::ArrayFactory factory;
    std::shared_ptr<ForwardingBuffer> out = std::make_shared<ForwardingBuffer>(std::cout);
    std::shared_ptr<ForwardingBuffer> err = std::make_shared<ForwardingBuffer>(std::cerr);
};

MatlabBridge::MatlabBridge(const std::string& matlabCodeDir, bool shared) : impl_(new Impl) {
    if (shared) {
        const auto names = me::findMATLAB();
        if (names.empty()) throw std::runtime_error("no shared MATLAB session found (run matlab.engine.shareEngine)");
        std::cout << "[matlab] connecting to shared session\n";
        impl_->engine = me::connectMATLAB(names[0]);
    } else {
        std::cout << "[matlab] starting MATLAB engine (this can take a few seconds)...\n";
        impl_->engine = me::startMATLAB({u"-nodesktop", u"-nosplash"});
    }
    // Make the solver functions (matlab/*.m) visible to MATLAB.
    eval("addpath('" + matlabCodeDir + "');");
}

MatlabBridge::~MatlabBridge() = default;

void MatlabBridge::eval(const std::string& statement) {
    try {
        impl_->engine->eval(u16(statement), impl_->out, impl_->err);
    } catch (const me::MATLABExecutionException& ex) {
        throw std::runtime_error(std::string("MATLAB error: ") + ex.what());
    } catch (const me::MATLABSyntaxException& ex) {
        throw std::runtime_error(std::string("MATLAB syntax error: ") + ex.what());
    }
}

void MatlabBridge::put(const std::string& name, const Eigen::MatrixXd& M) {
    // Eigen and MATLAB share the column-major layout -> straight copy.
    md::TypedArray<double> A = impl_->factory.createArray(
        md::ArrayDimensions{std::size_t(M.rows()), std::size_t(M.cols())}, M.data(), M.data() + M.size());
    impl_->engine->setVariable(u16(name), std::move(A));
}

void MatlabBridge::put(const std::string& name, const Eigen::MatrixXi& M, int offset) {
    // Indices are shipped as doubles (MATLAB's sparse()/indexing default type).
    Eigen::MatrixXd D = M.cast<double>().array() + double(offset);
    put(name, D);
}

void MatlabBridge::putScalar(const std::string& name, double v) {
    impl_->engine->setVariable(u16(name), impl_->factory.createScalar(v));
}

void MatlabBridge::putString(const std::string& name, const std::string& s) {
    impl_->engine->setVariable(u16(name), impl_->factory.createCharArray(s));
}

Eigen::MatrixXd MatlabBridge::getMatrix(const std::string& name) {
    md::Array a = impl_->engine->getVariable(u16(name));
    if (a.getType() != md::ArrayType::DOUBLE) throw std::runtime_error("MATLAB variable " + name + " is not double");
    md::TypedArray<double> t(std::move(a));
    const auto dims = t.getDimensions();
    const std::size_t r = dims.size() > 0 ? dims[0] : 0, c = dims.size() > 1 ? dims[1] : 1;
    Eigen::MatrixXd M(r, c);
    std::size_t i = 0;
    for (const double v : t) M.data()[i++] = v;   // column-major iteration
    return M;
}

double MatlabBridge::getScalar(const std::string& name) {
    const Eigen::MatrixXd M = getMatrix(name);
    if (M.size() != 1) throw std::runtime_error("MATLAB variable " + name + " is not a scalar");
    return M(0, 0);
}

std::string MatlabBridge::getString(const std::string& name) {
    md::CharArray c(impl_->engine->getVariable(u16(name)));
    return c.toAscii();
}

SolverResult MatlabBridge::solve(const MetricData& d, const SolverOptions& o) {
    // ---- 1) push every matrix built by the C++ side into the MATLAB workspace.
    //         Topology arrays become 1-based (MATLAB indexing).
    put("bdp_F",   d.F,   1);
    put("bdp_F2E", d.F2E, 1);
    put("bdp_EV",  d.EV,  1);
    put("bdp_e0",  Eigen::MatrixXd(d.e0));
    put("bdp_cA",  d.cA);
    put("bdp_cV",  d.cV);
    put("bdp_cW",  d.cW);
    put("bdp_wF",  Eigen::MatrixXd(d.wF));
    put("bdp_isBnd", Eigen::MatrixXi(d.isBnd));
    put("bdp_Vc",  Eigen::MatrixXi(d.Vc), 1);
    put("bdp_K0",  Eigen::MatrixXd(d.K0));
    // cone singularities: Eigen copies of the CGAL vertex property maps
    // "v:bdp_is_cone" / "v:bdp_cone_angle" (|V| x 1 each, angle in radians)
    put("bdp_isCone",    Eigen::MatrixXi(d.isCone));
    put("bdp_coneAngle", Eigen::MatrixXd(d.coneAngle));
    put("bdp_Eref", d.Eref.size() ? d.Eref : Eigen::MatrixXd(0, 3));

    // ---- 2) options
    putString("bdp_variant", o.variant);
    putString("bdp_ccmode", o.curvatureConstraint);
    putString("bdp_solver", o.solver);
    auto bound = [&](const char* name, double v) {
        if (v < 0) eval(std::string("bdp_opts.") + name + " = [];");          // variant default
        else { putScalar("bdp_tmp", v); eval(std::string("bdp_opts.") + name + " = bdp_tmp;"); }
    };
    eval("bdp_opts = bdp_default_options();"
         "bdp_opts.variant = bdp_variant;"
         "bdp_opts.curvatureConstraint = bdp_ccmode;"
         "bdp_opts.solver = bdp_solver;");
    putScalar("bdp_tmp", o.avgWeight);         eval("bdp_opts.avgWeight = bdp_tmp;");
    bound("Ca", o.Ca);
    bound("Cb", o.Cb);
    bound("CK", o.CK);
    putScalar("bdp_tmp", o.maxIters);          eval("bdp_opts.maxIters = bdp_tmp;");
    putScalar("bdp_tmp", o.useCETM ? 1 : 0);   eval("bdp_opts.useCETM = logical(bdp_tmp);");
    putScalar("bdp_tmp", o.cetmThresholdDeg);  eval("bdp_opts.cetmThreshold = bdp_tmp*pi/180;");
    putScalar("bdp_tmp", o.verbose ? 1 : 0);   eval("bdp_opts.verbose = logical(bdp_tmp);");

    // ---- 3) run the SCP / Newton solver inside MATLAB
    eval("bdp_data = struct('F',bdp_F,'F2E',bdp_F2E,'EV',bdp_EV,'e0',bdp_e0,"
         "'cA',bdp_cA,'cV',bdp_cV,'cW',bdp_cW,'wF',bdp_wF,'isBnd',bdp_isBnd,"
         "'Vc',bdp_Vc,'K0',bdp_K0,'Eref',bdp_Eref,"
         "'isCone',bdp_isCone,'coneAngle',bdp_coneAngle);");
    eval("[bdp_e, bdp_info] = bdp_solve(bdp_data, bdp_opts);");

    // ---- 4) fetch the optimised squared edge lengths and diagnostics
    eval("bdp_hist_err = bdp_info.curvErr; bdp_hist_E = bdp_info.energy;"
         "bdp_conv = double(bdp_info.converged); bdp_cetm = double(bdp_info.usedCETM);"
         "bdp_time = bdp_info.time; bdp_msg = bdp_info.message;"
         "bdp_quality = bdp_info.quality; bdp_solver_used = bdp_info.solver;");   // tau, tauSc, Kmax, symDir, ... (see bdp_quality.m)
    SolverResult r;
    r.e = getMatrix("bdp_e");
    r.curvErrHistory = getMatrix("bdp_hist_err");
    r.energyHistory = getMatrix("bdp_hist_E");
    r.converged = getScalar("bdp_conv") != 0;
    r.usedCETM = getScalar("bdp_cetm") != 0;
    r.seconds = getScalar("bdp_time");
    r.message = getString("bdp_msg");
    r.solverUsed = getString("bdp_solver_used");
    if (r.e.size() != d.e0.size()) throw std::runtime_error("MATLAB returned a metric of wrong size");
    return r;
}

}  // namespace bdp

#endif  // BDP_HAVE_MATLAB
