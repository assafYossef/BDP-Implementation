function opts = bdp_default_options()
%BDP_DEFAULT_OPTIONS Default parameters of the Bounded-Distortion Parametrization solver.
%
%   Reference: E. Chien, Z. Levi, O. Weber, "Bounded Distortion Parametrization
%   in the Space of Metrics", ACM TOG 35(6) (SIGGRAPH Asia 2016).
%
%   Variants (Section 8.1 of the paper):
%     'V1R' E = max_i tau_i + avgWeight * 0.5 * sum_i w_i beta_i      (no bounds; DEFAULT)
%           = the objective of the reference implementation (method 50): best worst case
%     'V1'  E = max_i tau_i,  tau = max(sigma_a^2, 1/sigma_b^2)       (no bounds by default)
%     'V2'  E = sum_i w_i beta_i, beta = sigma_a^2 + 1/sigma_b^2        (Ca=20^2, Cb=1/20^2, CK=20^2)
%     'V3'  E = sum_i w_i ||M_i - I||_F^2  (Green-Lagrange)              (Ca=20^2, Cb=1/20^2, CK=20^2)
%     'V4'  E = max_i |eta_i|              (conformal)                   (Cb = 1, loose Ca = 20^2)
%     'V5'  E = sum_i w_i [(A_i-A0_i)^2 + 4|eta_i-eta0_i|^2] (projection) (Ca=3.5^2, Cb=1/3.5^2)
%
%   Distortion bounds (inequalities (6)-(8) of the paper):
%     Ca : sigma_a^2 <= Ca     Cb : sigma_b^2 >= Cb     CK : K^2 = sigma_a^2/sigma_b^2 <= CK
%   Leave a bound empty ([]) to use the variant default; set it to Inf (Ca, CK)
%   or 0 (Cb) to disable it.

opts.variant            = 'V1R';
opts.avgWeight          = 0.01;       % V1R: weight of the average term (reference: 0.01)
opts.Ca                 = [];
opts.Cb                 = [];
opts.CK                 = [];

% ---- Convex subproblem solver ------------------------------------------------
%   'auto'     : MOSEK called directly if usable (also looked up inside the CVX
%                folder, where CVX keeps its own copy off the MATLAB path),
%                else CVX with its default solver (e.g. MOSEK, as the reference code),
%                else MATLAB's coneprog (Optimization Toolbox R2020b+),
%                else our Newton log-barrier solver (bdp_barrier_newton).
%   'mosek' | 'cvx' | 'coneprog' | 'newton' | function handle @(cp) -> [x, status]
%   opts.cvxSolver : optional solver name passed to cvx_solver (default: CVX's default)
opts.solver             = 'auto';

% ---- Sequential Convex Programming (Section 7) -----------------------------
opts.maxIters           = 500;        % paper: at most 500 SCP iterations
opts.tolFlat            = 1e-9;       % L_inf curvature error (radians) declared "flat"
opts.curvatureConstraint = 'auto';   % 'hard' (J e = b), 'soft' (paper Sec. 7.2) or 'auto'
opts.hardBelow          = 5*pi/180;   % 'auto': try the hard form once |K-K0|_inf < 5 deg
opts.lambdaFactor       = 100;        % lambda_l = 100 * E(e_l) / ||K0 - K(e_l)||_inf (Sec. 7.2)
opts.energyFloor        = 1e-2;       % guards lambda_l when E(e_l) = 0 (e.g. V3/V4 at the identity)
opts.maxHalvings        = 30;         % backtracking halvings of the SCP line search (Sec. 7.2)

% ---- CETM acceleration (Section 8.2) ---------------------------------------
opts.useCETM            = true;
opts.cetmThreshold      = pi/180;     % switch to CETM once the L_inf error < 1 degree
opts.cetmMaxIters       = 100;
opts.cetmRespectBounds  = true;       % reject a CETM result that violates the distortion bounds
opts.boundTol           = 1e-6;       % relative tolerance for that check

% ---- Newton / log-barrier interior point solver for each convex subproblem --
opts.newton.maxIters    = 60;         % Newton steps per centering
opts.newton.maxOuter    = 40;         % barrier parameter updates
opts.newton.tol         = 1e-9;       % stop centering when (Newton decrement)^2/2 < tol
opts.newton.muFactor    = 20;         % barrier parameter growth
opts.newton.gapTol      = 1e-6;       % stop when m/mu < gapTol * (|f0| + 1e-6)
opts.newton.verbose     = false;

opts.verbose            = true;
end
