function [e, info] = bdp_solve(data, opts)
%BDP_SOLVE Bounded distortion flattening of a discrete metric (squared edge lengths).
%
%   [e, info] = bdp_solve(data, opts)
%
%   Implements Chien, Levi & Weber, "Bounded Distortion Parametrization in the
%   Space of Metrics" (SIGGRAPH Asia 2016), with the CVX/Mosek SOCP solver of the
%   paper replaced by a custom Newton-based log-barrier interior point solver
%   (bdp_barrier_newton) built on sparse LDL' factorisations of the Newton/KKT
%   systems.
%
%   INPUT  data  (built by the C++ process and pushed through the MATLAB Engine):
%     F    |F|x3  triangle vertex indices (1-based)
%     F2E  |F|x3  F2E(f,k) = edge opposite corner k of face f (1-based)
%     EV   |E|x2  edge end points (1-based)
%     e0   |E|x1  squared edge lengths of the input mesh (original metric)
%     cA, cV, cW  |F|x3  linear maps from the 3 squared lengths of a face to
%                 A, 2Re(eta), 2Im(eta) (see bdp_distortion / MetricData.cpp)
%     wF   |F|x1  normalised original areas  (sum = 1)
%     isBnd|V|x1  1 for boundary vertices
%     Vc   |Vc|x1 vertices with curvature constraints (1-based)
%     K0   |Vc|x1 prescribed curvatures (angle defects, radians)
%     Eref |F|x3  (optional, variant V5) reference squared lengths per face
%     isCone    |V|x1 (optional) 1 at cone vertices     } Eigen copies of the CGAL
%     coneAngle |V|x1 (optional) cone angle (radians)   } vertex property maps
%                 (from the `c <vertex> <angle/pi>` records of the OBJ file, a cone
%                 file or the automatic placement); used to cross-check K0.
%   OUTPUT
%     e    |E|x1  optimised squared edge lengths (flat cone metric if converged)
%     info        convergence history, per-face distortion of the result and
%                 info.quality (see bdp_quality: tau, tauSc, Kmax, symDir,
%                 maxCurvErrDeg, valid)
%
%   Optimisation problem (Section 6.1, Eqs. (5)-(9)):
%       min_e  E(e)
%       s.t.   sigma_a^2(e) <= Ca, sigma_b^2(e) >= Cb, K^2(e) <= CK   (every face)
%              K_j(e) = K0_j                                          (j in Vc)
%   All constraints but the curvature are convex second-order cones in e.
%
%   Sequential Convex Programming (Section 7):
%       l = 0,1,2,...
%         J_l = dK/de(e_l)                       (Appendix C)
%         b   = K0 - K(e_l) + J_l e_l
%         e*  = argmin E(e) + lambda_l || J_l e - b ||_inf   s.t. bounds   ('soft', Sec. 7.2)
%               lambda_l = 100 E(e_l) / ||K0 - K(e_l)||_inf
%           or  argmin E(e)  s.t. J_l e = b and the bounds             ('hard', Eq. (12))
%           ('auto': soft far from flat, hard close to flat, soft if hard is infeasible)
%         e_{l+1} = e_l + alpha (e* - e_l), alpha halved until ||K - K0||_inf decreases
%   Because the bounded distortion set is convex and both e_l and e* are inside
%   it, EVERY iterate e_l satisfies the distortion bounds.
%   Once ||K-K0||_inf < 1 degree, CETM (Newton on conformal log-scale factors) is
%   used to finish the flattening (Section 8.2).

if nargin < 2, opts = bdp_default_options(); end
opts = merge_defaults(opts, bdp_default_options());
tStart = tic;

data = prepare_data(data);
[opts, bnd] = resolve_bounds(opts, data);
variant = upper(opts.variant);
data.avgWeight = opts.avgWeight;

e = data.e0(:) * bnd.startScale;
if opts.verbose
    fprintf('[bdp] |V|=%d |E|=%d |F|=%d |Vc|=%d cones=%d variant=%s  Ca=%g Cb=%g CK=%g\n', ...
        data.nV, data.nE, data.nF, numel(data.Vc), data.nCones, variant, bnd.Ca, bnd.Cb, bnd.CK);
end

info.curvErr = []; info.energy = []; info.alpha = []; info.newtonSteps = [];
info.usedCETM = false; info.converged = false; info.message = '';
cetmThreshold = opts.cetmThreshold;
lastHardOk = false;
[solverName, solverFcn, solverReason] = select_solver(opts.solver);
info.solver = solverName;
if opts.verbose
    fprintf('[bdp] subproblem solver: %s  (%s)\n', solverName, solverReason);
end

for l = 0:opts.maxIters
    [K, ~, T] = bdp_curvature(e, data);
    r   = K(data.Vc) - data.K0;
    err = max([abs(r); 0]);
    El  = bdp_energy(e, data, variant);
    info.curvErr(end+1,1) = err;
    info.energy(end+1,1)  = El;
    if opts.verbose
        fprintf('[bdp] iter %3d  |K-K0|_inf = %10.4e (%8.4f deg)  E = %12.6e\n', ...
            l, err, err*180/pi, El);
    end
    if err <= opts.tolFlat
        info.converged = true; info.message = 'flat metric reached'; break;
    end

    % ---------------- CETM acceleration (Section 8.2) ---------------------
    if opts.useCETM && err < cetmThreshold
        [ec, okC, errC] = bdp_cetm_newton(e, data, opts);
        accept = okC;
        if okC && opts.cetmRespectBounds
            accept = satisfies_bounds(ec, data, bnd, opts.boundTol);
        end
        if opts.verbose
            fprintf('[bdp]   CETM: success=%d  |K-K0|_inf=%9.3e  bounds ok=%d\n', okC, errC, accept);
        end
        if accept
            e = ec; info.usedCETM = true; info.converged = true;
            info.curvErr(end+1,1) = errC;
            info.energy(end+1,1)  = bdp_energy(e, data, variant);
            info.message = 'flat metric reached (CETM finish)';
            break;
        end
        cetmThreshold = err / 4;        % try again once the error dropped 4x
    end
    if l == opts.maxIters
        info.message = 'maximum number of iterations reached'; break;
    end

    % ---------------- Linearised curvature constraint (Section 7.1) --------
    %   K(e) ~ K(e_l) + J_l (e - e_l) = K0   <=>   J_l e = b,  b = K0 - K(e_l) + J_l e_l
    %   (J_l e_l = 0 by scale invariance, Eq. (13); we keep the term for round-off.)
    JK = bdp_curvature_jacobian(e, data, T);
    J  = JK(data.Vc, :);
    b  = data.K0 - K(data.Vc) + J * e;
    lambda = opts.lambdaFactor * max(El, opts.energyFloor) / err;     % Section 7.2

    % ---------------- Convex subproblem -------------------------------------
    %   'hard': J e = b as an equality. Exact, and the curvature error of e* is
    %           then purely second order (this is what the reference script does).
    %   'soft': the paper's  E + lambda ||J e - b||_inf  (always feasible).
    %   'auto': hard, falling back to soft when J e = b misses the bounded
    %           distortion set (this is exactly the situation of Section 7.2).
    %   Solver: a conic interior point solver (MOSEK / coneprog, see
    %   bdp_socp_solve) when available, otherwise our Newton barrier method.
    mode = lower(opts.curvatureConstraint);
    usedMode = ''; stepInfo = '';
    eStar = [];
    while ~strcmp(solverName, 'newton')
        modes = {'hard', 'soft'};
        if strcmp(mode, 'soft'), modes = {'soft'}; elseif strcmp(mode, 'hard'), modes = {'hard'}; end
        crashed = false;
        for im = 1:numel(modes)
            try
                cp = bdp_socp_build(e, data, bnd, variant, J, b, lambda, modes{im});
                [xs, sst] = bdp_socp_solve(cp, solverFcn, opts);
            catch solverErr
                crashed = true; break;
            end
            usedMode = modes{im};
            stepInfo = sprintf('%s %.2fs (%s it) %s', solverName, sst.time, num2str(sst.iterations), sst.status);
            info.newtonSteps(end+1,1) = sst.iterations;
            if strcmp(sst.status, 'optimal')
                eStar = xs(cp.idxE); break;
            end
        end
        if ~crashed, break; end
        % e.g. MOSEK installed but not licensed: mosek -> cvx -> coneprog -> newton
        nextSolver = next_solver(solverName);
        warning('bdp:conic', 'solver "%s" failed (%s); switching to "%s"', solverName, solverErr.message, nextSolver);
        solverName = nextSolver; solverFcn = nextSolver; info.solver = nextSolver;
    end
    if isempty(eStar) && ~strcmp(solverName, 'newton')
        info.message = 'conic solver failed on the subproblem'; break;
    end
    if strcmp(solverName, 'newton')
        st.eqInfeasible = true;
        tryHard = strcmp(mode, 'hard') || ...
                  (strcmp(mode, 'auto') && (err < opts.hardBelow || lastHardOk));
        if tryHard
            [prob, x0] = build_subproblem(e, data, bnd, variant, J, b, lambda, 'hard');
            [x, st] = bdp_barrier_newton(prob, x0, opts.newton);
            usedMode = 'hard';
            lastHardOk = ~st.eqInfeasible;
        end
        if st.eqInfeasible && ~strcmp(mode, 'hard')
            [prob, x0] = build_subproblem(e, data, bnd, variant, J, b, lambda, 'soft');
            [x, st] = bdp_barrier_newton(prob, x0, opts.newton);
            usedMode = 'soft';
        end
        eStar = x(1:data.nE);
        info.newtonSteps(end+1,1) = st.newtonSteps;
        stepInfo = sprintf('newton steps=%3d  barrier gap=%9.2e', st.newtonSteps, st.gap);
    end

    % ---------------- Backtracking on the curvature error (Section 7.2) -----
    alpha = 1; accepted = false;
    for h = 0:opts.maxHalvings
        eTry = e + alpha * (eStar - e);
        [Kt, ~, ~, valid] = bdp_curvature(eTry, data);
        if valid && max(abs(Kt(data.Vc) - data.K0)) < err
            accepted = true; break;
        end
        alpha = alpha / 2;
    end
    info.alpha(end+1,1) = alpha * accepted;
    if opts.verbose
        fprintf('[bdp]   %s  %s  step=%g\n', usedMode, stepInfo, alpha*accepted);
    end
    if ~accepted
        info.message = 'line search failed to decrease the curvature error (stalled)';
        break;
    end
    e = eTry;
end

% ---------------- V4: optimal global scale ------------------------------------
% max|eta| is 1-homogeneous in e while angles (curvatures) and K are scale
% invariant, so for a fixed flat shape the best scale puts min sigma_b^2 = Cb.
% (Needed in particular after a CETM finish, which freezes the current scale.)
if strcmp(variant, 'V4')
    D = bdp_distortion(e, data);
    c = bnd.Cb / min(D.sb2);
    if c < 1 && satisfies_bounds(c * e, data, bnd, opts.boundTol)
        e = c * e;
    end
end

% ---------------- Final report ---------------------------------------------
D = bdp_distortion(e, data);
K = bdp_curvature(e, data);
info.K = K;
info.finalEnergy = bdp_energy(e, data, variant);
info.tau = D.tau; info.k = D.k; info.sigmaA = D.sigmaA; info.sigmaB = D.sigmaB;
info.Kconf = D.K;
info.boundsSatisfied = satisfies_bounds(e, data, bnd, opts.boundTol);
info.bounds = bnd;
info.iterations = numel(info.alpha);
info.quality = bdp_quality(e, data);   % tau, tauSc, Kmax, symDir, maxCurvErrDeg, valid
info.time = toc(tStart);
if opts.verbose
    q = info.quality;
    fprintf('[bdp] done: %s  (%d SCP iterations, %.2f s)\n', info.message, info.iterations, info.time);
    fprintf('[bdp] quality                     value   perfect\n');
    fprintf('[bdp]   tau    worst stretch %12.6g   1\n', q.tau);
    fprintf('[bdp]   tauSc  after scale   %12.6g   1   (best global scale x%.4g)\n', q.tauSc, q.bestScale);
    fprintf('[bdp]   Kmax   worst shape   %12.6g   1\n', q.Kmax);
    fprintf('[bdp]   symDir area-avg      %12.6g   4   (excess %.4g)\n', q.symDir, q.symDirExcess);
    fprintf('[bdp]   max|K-K*| (deg)      %12.3e   0\n', q.maxCurvErrDeg);
    fprintf('[bdp]   valid                %12d   1   (%d degenerate triangles)\n', q.valid, q.degenerateTriangles);
    fprintf('[bdp]   distortion bounds satisfied: %d\n', info.boundsSatisfied);
end
end

% =============================================================================
function [name, fcn, reason] = select_solver(s)
% 'auto' -> MOSEK directly (also found inside the CVX folder) -> CVX (with its
% default solver, e.g. MOSEK) -> coneprog -> our Newton barrier solver.
if isa(s, 'function_handle'), name = func2str(s); fcn = s; reason = 'user-supplied function'; return; end
[name, reason] = bdp_find_solver(s);
fcn = name;
end

function n = next_solver(name)
% fallback chain when a solver crashes (e.g. no licence)
isOctave = exist('OCTAVE_VERSION', 'builtin') > 0;
chain = {'mosek', 'cvx', 'coneprog', 'newton'};
avail = [exist('mosekopt', 'file') > 0, exist('cvx_begin', 'file') > 0, ~isOctave && exist('coneprog', 'file') > 0, true];
k = find(strcmp(chain, name));
if isempty(k), k = 0; end
n = 'newton';
for i = k+1:numel(chain)
    if avail(i), n = chain{i}; return; end
end
end

function data = prepare_data(data)
fld = {'F','F2E','EV','e0','cA','cV','cW','wF','isBnd','Vc','K0'};
for i = 1:numel(fld)
    if ~isfield(data, fld{i}), error('bdp:data', 'data.%s is missing', fld{i}); end
    data.(fld{i}) = double(data.(fld{i}));
end
data.e0 = data.e0(:); data.wF = data.wF(:); data.isBnd = data.isBnd(:);
data.Vc = data.Vc(:); data.K0 = data.K0(:);
data.nF = size(data.F, 1);
data.nE = numel(data.e0);
data.nV = max([max(data.F(:)), numel(data.isBnd)]);
% Linear operators (|F| x |E|): A = GA e, 2Re(eta) = GV e, 2Im(eta) = GW e
rows = repmat((1:data.nF)', 1, 3);
data.GA = sparse(rows(:), data.F2E(:), data.cA(:), data.nF, data.nE);
data.GV = sparse(rows(:), data.F2E(:), data.cV(:), data.nF, data.nE);
data.GW = sparse(rows(:), data.F2E(:), data.cW(:), data.nF, data.nE);
% Cones: K0 must equal the angle defect implied by the prescribed cone angles,
%   K0(v) = 2*pi - coneAngle(v) (interior)  or  pi - coneAngle(v) (boundary)  (Section 5)
data.nCones = 0;
if isfield(data, 'isCone') && ~isempty(data.isCone)
    isCone = logical(data.isCone(:));
    data.nCones = nnz(isCone);
    if isfield(data, 'coneAngle') && numel(data.coneAngle) == numel(isCone)
        fullAngle = 2*pi - pi*data.isBnd(:);
        Kc = fullAngle - double(data.coneAngle(:));
        Kfull = nan(numel(isCone), 1);
        Kfull(data.Vc) = data.K0;
        bad = isCone & ~(abs(Kfull - Kc) < 1e-9);
        if any(bad)
            error('bdp:cones', '%d cone vertices have K0 inconsistent with their cone angle', nnz(bad));
        end
    end
end
% Reference metric quantities for V5 (A0, 2Re eta0, 2Im eta0 per face, Eq. (10))
if isfield(data, 'Eref') && ~isempty(data.Eref)
    Er = double(data.Eref);
    data.A0 = sum(data.cA .* Er, 2);
    data.v0 = sum(data.cV .* Er, 2);
    data.w0 = sum(data.cW .* Er, 2);
else           % identity reference map: V5 reduces to (half) the Green-Lagrange energy
    data.A0 = ones(data.nF,1); data.v0 = zeros(data.nF,1); data.w0 = zeros(data.nF,1);
end
end

function [opts, bnd] = resolve_bounds(opts, data) %#ok<INUSD>
% Variant defaults (Section 8.1 of the paper)
switch upper(opts.variant)
    case {'V1','V1R'}, d = [Inf, 0, Inf];             % bounds implied by min max tau
    case {'V2','V3'}, d = [20^2, 1/20^2, 20^2];
    case 'V4', d = [20^2, 1, Inf];
        % Paper: only sigma_b^2 >= 1 (PD + no collapse to the zero metric). max|eta|
        % is invariant along e -> c e for conformal metrics, so without an upper
        % bound the barrier problem is unbounded below (its central path runs to
        % infinite scale). A loose sigma_a^2 <= 20^2 makes it well posed and is
        % inactive at the optimum (which prefers sigma_b = 1).
    case 'V5', d = [3.5^2, 1/3.5^2, Inf];
    otherwise, error('bdp:variant', 'Unknown variant %s', opts.variant);
end
bnd.Ca = d(1); bnd.Cb = d(2); bnd.CK = d(3);
if ~isempty(opts.Ca), bnd.Ca = opts.Ca; end
if strcmpi(opts.variant, 'V4') && ~isfinite(bnd.Ca)
    warning('bdp:V4', 'V4 needs a finite Ca for a bounded barrier problem; using Ca = 1e4.');
    bnd.Ca = 1e4;
end
if ~isempty(opts.Cb), bnd.Cb = opts.Cb; end
if ~isempty(opts.CK), bnd.CK = opts.CK; end
% The bounded distortion set must have a non-empty interior (Section 7.2: Ca > Cb)
if ~(bnd.Ca > bnd.Cb) || bnd.CK <= 1 || bnd.Cb < 0
    error('bdp:bounds', 'Distortion bounds must satisfy Ca > Cb >= 0 and CK > 1.');
end
% Starting metric: the original metric (identity map, sigma_a = sigma_b = 1),
% uniformly rescaled by c if needed so it lies STRICTLY inside the bounds
% (e.g. V4 uses Cb = 1, which puts the identity on the boundary). A uniform
% scale leaves all angles, hence all curvatures, unchanged.
c = 1;
if bnd.Cb >= 1                         % move just above the lower bound
    c = 1.05 * bnd.Cb;
    if isfinite(bnd.Ca), c = min(c, sqrt(bnd.Cb * bnd.Ca)); end
elseif bnd.Ca <= 1                     % move just below the upper bound
    c = bnd.Ca / 1.05;
    if bnd.Cb > 0, c = max(c, sqrt(bnd.Cb * bnd.Ca)); end
end
bnd.startScale = c;
end

function ok = satisfies_bounds(e, data, bnd, tol)
[~, ~, valid] = bdp_angles(e, data.F2E);
D = bdp_distortion(e, data);
ok = valid && all(D.sb2 > 0) ...
     && all(D.sa2 <= bnd.Ca * (1 + tol)) ...
     && all(D.sb2 >= bnd.Cb * (1 - tol)) ...
     && all(D.sa2 ./ D.sb2 <= bnd.CK * (1 + tol));
end

function [prob, x0] = build_subproblem(e, data, bnd, variant, J, b, lambda, mode)
% Builds the SOCP of Section 7 as a barrier problem for bdp_barrier_newton.
%
%   Variables  x = [ e (|E|) ; aux (variant specific) ; t (1, soft mode only) ]
%   mode 'hard':  J e = b  is passed as the equality  Aeq x = beq.
%   mode 'soft':  t is the epigraph variable of the L_inf curvature penalty:
%       lambda ||J e - b||_inf  ->  lambda t,   t >= |J_j e - b_j|  for all j
%   which is a 1-D second-order cone  t >= sqrt(r_j^2),  r = J e - b.
nE = data.nE; nF = data.nF; nC = size(J,1);
D  = bdp_distortion(e, data);
switch variant
    case {'V1','V4'}, nAux = 1;            % global epigraph variable s
    case 'V2',        nAux = 2*nF;         % per-face epigraph variables p_i, q_i
    case 'V1R',       nAux = 1 + 2*nF;     % s (max) + p_i, q_i (average)
    otherwise,        nAux = 0;
end
soft = strcmp(mode, 'soft');
n   = nE + nAux + soft;
it  = n;                                  % index of t (soft mode)
Z   = @(m) sparse(m, n - nE);
GA  = [data.GA, Z(nF)];
GV  = [data.GV, Z(nF)];
GW  = [data.GW, Z(nF)];
sel = @(rows, cols, m) sparse(rows, cols, 1, m, n);

fams = {};
lin  = {};
dense = [];          % variables coupled to (almost) every row -> bordered KKT solve
c    = zeros(n, 1);
quad = {};
x0   = [e; zeros(nAux, 1); zeros(soft, 1)];

% ---- distortion bounds (6)-(8):  u >= |2 eta| = sqrt(v^2 + w^2) ----------------
%   sigma_a^2 <= Ca  <=>  Ca - A        >= 2|eta|
%   sigma_b^2 >= Cb  <=>  A - Cb        >= 2|eta|
%   K^2       <= CK  <=>  A (CK-1)/(CK+1) >= 2|eta|
if isfinite(bnd.Ca), fams{end+1} = fam(-GA, bnd.Ca*ones(nF,1), GV, GW); end
if bnd.Cb > 0,       fams{end+1} = fam( GA, -bnd.Cb*ones(nF,1), GV, GW); end
if isfinite(bnd.CK), fams{end+1} = fam(((bnd.CK-1)/(bnd.CK+1)) * GA, zeros(nF,1), GV, GW); end
hasLower = bnd.Cb > 0 || isfinite(bnd.CK);

switch variant
    case 'V1'
        % E = max_i tau_i = min s  s.t.  sigma_a^2 <= s  and  sigma_b^2 >= 1/s
        %   s - A     >= 2|eta|   (linear SOC)
        %   A - 1/s   >= 2|eta|   (SOC with concave u)
        is = nE + 1;
        S  = sel((1:nF)', is*ones(nF,1), nF);
        fams{end+1} = fam(S - GA, zeros(nF,1), GV, GW);
        fams{end+1} = fam_inv(GA, zeros(nF,1), GV, GW, is*ones(nF,1));
        c(is) = 1;
        x0(is) = 1.05 * max(D.tau) + 1e-9;
        dense(end+1) = is;
    case 'V1R'
        % E = max_i tau_i + avgWeight * 0.5 * sum_i w_i beta_i
        %   (reference script, method 50: max tau + 0.01 * average isometric distortion)
        is = nE + 1;
        S  = sel((1:nF)', is*ones(nF,1), nF);
        fams{end+1} = fam(S - GA, zeros(nF,1), GV, GW);                  % sigma_a^2 <= s
        fams{end+1} = fam_inv(GA, zeros(nF,1), GV, GW, is*ones(nF,1));   % sigma_b^2 >= 1/s
        ip = nE + 1 + (1:nF)'; iq = nE + 1 + nF + (1:nF)';
        fams{end+1} = fam(sel((1:nF)', ip, nF) - GA, zeros(nF,1), GV, GW);  % p_i >= sigma_a^2
        fams{end+1} = fam_inv(GA, zeros(nF,1), GV, GW, iq);                 % q_i >= 1/sigma_b^2
        c(is) = 1;
        c(ip) = 0.5 * data.avgWeight * data.wF;
        c(iq) = 0.5 * data.avgWeight * data.wF;
        x0(is) = 1.05 * max(D.tau) + 1e-9;
        x0(ip) = 1.05 * D.sa2 + 1e-9;
        x0(iq) = 1.05 ./ D.sb2 + 1e-9;
        dense(end+1) = is;
    case 'V2'
        % E = sum w_i (p_i + q_i),  p_i >= sigma_a^2,  q_i >= 1/sigma_b^2
        ip = nE + (1:nF)'; iq = nE + nF + (1:nF)';
        fams{end+1} = fam(sel((1:nF)', ip, nF) - GA, zeros(nF,1), GV, GW);
        fams{end+1} = fam_inv(GA, zeros(nF,1), GV, GW, iq);
        c(ip) = data.wF; c(iq) = data.wF;
        x0(ip) = 1.05 * D.sa2 + 1e-9;
        x0(iq) = 1.05 ./ D.sb2 + 1e-9;
    case 'V3'
        % E_GL = sum w_i [2 (A-1)^2 + 2 (2Re eta)^2 + 2 (2Im eta)^2]   (convex quadratic)
        quad{end+1} = qterm(GA, ones(nF,1), 2*data.wF);
        quad{end+1} = qterm(GV, zeros(nF,1), 2*data.wF);
        quad{end+1} = qterm(GW, zeros(nF,1), 2*data.wF);
    case 'V4'
        % E = max |eta_i| = min s  s.t.  2s >= 2|eta_i|
        is = nE + 1;
        fams{end+1} = fam(2*sel((1:nF)', is*ones(nF,1), nF), zeros(nF,1), GV, GW);
        c(is) = 1;
        dense(end+1) = is;
        % start well inside the cone (not at its apex: at the identity eta = 0)
        x0(is) = 1.05 * max(D.etaAbs) + 0.05 * median(D.A);
    case 'V5'
        % E_GLD = sum w_i [(A-A0)^2 + 4|eta-eta0|^2]   (Eq. (10))
        quad{end+1} = qterm(GA, data.A0, data.wF);
        quad{end+1} = qterm(GV, data.v0, data.wF);
        quad{end+1} = qterm(GW, data.w0, data.wF);
end
if ~hasLower && ~any(strcmp(variant, {'V1','V1R','V2'}))
    % keep the metric positive definite (=> triangle inequalities, Section 4.3)
    fams{end+1} = fam(GA, -1e-8*ones(nF,1), GV, GW);
end

% ---- linearised curvature constraint (Section 7.1 / 7.2) --------------------
if ~soft
    Aeq = [J, sparse(nC, n - nE)];
    beq = b;
    if numel(data.Vc) == data.nV && nC > 0
        % closed surface with every curvature prescribed: Gauss-Bonnet makes the
        % rows of J sum to zero, so one equation is redundant -> drop it
        Aeq = Aeq(1:end-1, :); beq = beq(1:end-1);
    end
    prob.Aeq = Aeq; prob.beq = beq;
elseif nC > 0
    % t - (J e - b) >= 0   and   t + (J e - b) >= 0   (two linear barriers)
    T1 = sel((1:nC)', it*ones(nC,1), nC);
    JJ = [J, sparse(nC, n-nE)];
    lin{1} = struct('G', T1 - JJ, 'g0',  b);
    lin{2} = struct('G', T1 + JJ, 'g0', -b);
    c(it) = lambda;
    x0(it) = 1.1 * max(abs(J*e - b)) + 1e-12;
    dense(end+1) = it;
    prob.gapRef = max(bdp_energy(e, data, variant), 1e-3);   % accuracy w.r.t. E, not lambda*t
else
    x0(it) = 1;                       % no curvature constraints: pin the unused t
    quad{end+1} = qterm(sel(1, it, 1), 1, 1);
end

prob.n = n; prob.c = c; prob.quad = quad; prob.fams = fams; prob.lin = lin; prob.dense = dense;
end

function F = fam(Gu, u0, Gv, Gw, v0)
if nargin < 5, v0 = zeros(size(Gv,1),1); end
F.Gu = Gu; F.u0 = u0; F.Gv = Gv; F.v0 = v0; F.Gw = Gw; F.w0 = zeros(size(Gw,1),1);
F.inv = false; F.invIdx = []; F.m = size(Gu,1);
end

function F = fam_inv(Gu, u0, Gv, Gw, invIdx)
F = fam(Gu, u0, Gv, Gw);
F.inv = true; F.invIdx = invIdx;
end

function Q = qterm(G, y, w)
Q.G = G; Q.y = y; Q.w = w;
end

function o = merge_defaults(o, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(o, f{i})
        o.(f{i}) = d.(f{i});
    elseif isstruct(d.(f{i})) && isstruct(o.(f{i}))
        o.(f{i}) = merge_defaults(o.(f{i}), d.(f{i}));
    end
end
end
