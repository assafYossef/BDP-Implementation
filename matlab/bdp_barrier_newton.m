function [x, st] = bdp_barrier_newton(prob, x, o)
%BDP_BARRIER_NEWTON Custom Newton-based interior point solver for one SCP subproblem.
%
%   This replaces the CVX/Mosek SOCP solve of the paper (Section 7, problem (11)
%   and its soft variant of Section 7.2). Every convex subproblem we generate has
%   the generic form
%
%     minimize   f0(x) = c'x + sum_q sum_r wq_r ((Gq x - yq)_r)^2
%     subject to u_r(x) >= sqrt(v_r(x)^2 + w_r(x)^2)    (second-order cones, prob.fams)
%                a_r(x) = (G x + g0)_r > 0              (linear inequalities, prob.lin)
%                Aeq x = beq                            (linearised curvature, optional)
%
%   where v, w, a are affine in x and u is affine in x, optionally minus 1/x_aux
%   (a concave term, used for sigma_b^2 >= 1/s in variants V1/V2).
%
%   LOG-BARRIER METHOD (Boyd & Vandenberghe, Ch. 11):
%       phi(x) = - sum log(u^2 - v^2 - w^2) - sum log(a)
%       for mu = mu0, mu0*muFactor, ...
%           x <- argmin  F_mu(x) = mu*f0(x) + phi(x)   s.t. Aeq x = beq   (Newton)
%           stop when m/mu < gapTol*(|f0|+1e-6),   m = barrier degree
%   The SOC barrier is convex also for the concave u = A - 1/x_aux, because
%   -log(u^2-|z|^2) = -log(u - |z|^2/u) - log(u) and |z|^2/u is convex and
%   decreasing in u.
%
%   EQUALITY CONSTRAINTS are handled by Newton's method with infeasible start
%   (Boyd & Vandenberghe 10.3): starting from x = e_l (which does NOT satisfy
%   the linearised curvature Aeq x = beq) the KKT step drives the primal
%   residual to zero; a backtracking line search on the norm of the primal-dual
%   residual keeps x strictly inside the cones. If the residual cannot be
%   eliminated (the linearised constraint does not meet the bounded distortion
%   set) st.eqInfeasible is set and the caller falls back to the soft L_inf
%   penalty of Section 7.2.
%
%   NEWTON SYSTEM. With y = (u, v, w) = G x + h for every cone row, the Hessian
%   of F_mu is  H = Hd + G' P G,  P block diagonal (3x3 per cone row):
%       P_r = grad^2[-log(u^2-v^2-w^2)] = (2/D^2)(2 Jy (Jy)' - D J),
%       D = u^2 - v^2 - w^2,   J = diag(1,-1,-1).
%   For rows far from the cone boundary (D >= 1e-3 u^2) G_r' P_r G_r is added
%   to Hd directly (normal equations). For rows close to the boundary or to the
%   apex (e.g. the epigraph cones p_i >= sigma_a^2 at the optimum, where the
%   condition number of P_r explodes) the AUGMENTED form is used instead:
%
%       [ Hd    G'      Aeq' ] [dx]   [ -g ]
%       [ G   -P^{-1}    0   ] [ z] = [  0 ]          z = P G dx
%       [ Aeq   0        0   ] [nu]   [ -r ]
%
%   with the closed-form inverse  P_r^{-1} = y y' - (D/2) J  (no division,
%   exact even at the apex), and P^{-1} = diag(a.^2) for linear inequalities.
%   The sparse symmetric indefinite system is factorised by MA57 (MATLAB's
%   sparse ldl; UMFPACK LU under Octave) and polished by iterative refinement.

if nargin < 3, o = struct(); end
o = fill_defaults(o);
hasEq = isfield(prob, 'Aeq') && ~isempty(prob.Aeq);
if ~hasEq, prob.Aeq = sparse(0, prob.n); prob.beq = zeros(0, 1); end
nEq = size(prob.Aeq, 1);
if isfield(prob, 'dense'), dense = prob.dense(:); else, dense = []; end

mdeg = 0;
for k = 1:numel(prob.fams), mdeg = mdeg + 2*prob.fams{k}.m; end
for k = 1:numel(prob.lin),  mdeg = mdeg + size(prob.lin{k}.G, 1); end

H0 = objective_hessian(prob);
if ~barrier_feasible(prob, x)
    error('bdp:infeasibleStart', 'Starting point of the barrier method is not strictly feasible.');
end
eqTol  = 1e-11 * max(1, full(max(abs(prob.Aeq) * abs(x))));
eqFeas = ~hasEq || max(abs(prob.Aeq * x - prob.beq)) <= eqTol;
nu     = zeros(nEq, 1);

% ---- initial barrier parameter: least-squares fit of mu*g0 + gphi = 0 --------------
%   mu0 = -(g0' Hphi^{-1} gphi) / (g0' Hphi^{-1} g0)   (B&V 11.3.1, projected on Aeq dx = 0)
[f0, g0] = objective(prob, x);
[~, gP, K] = newton_system(prob, x, 0, H0, o.nearTol);
z0 = solve_kkt(K, g0, prob.n, zeros(nEq,1), dense);
mu = -(z0' * gP) / max(z0' * g0, realmin);
if ~isfinite(mu) || mu <= 0
    mu = mdeg / max(abs(f0), 1e-8);
end
mu = max(mu, 1e-8);

st.newtonSteps = 0; st.outer = 0; st.eqInfeasible = false;
for outer = 1:o.maxOuter
    st.outer = outer;
    % ------------------------- centering step -------------------------------
    for it = 1:o.maxIters
        [f0, g0] = objective(prob, x);
        [phi, gP, K] = newton_system(prob, x, mu, H0, o.nearTol);
        g = mu * g0 + gP;
        st.newtonSteps = st.newtonSteps + 1;
        if ~eqFeas
            % ---------- infeasible-start Newton step (B&V Alg. 10.2) ----------
            rp = prob.Aeq * x - prob.beq;
            [dxn, nuPlus] = solve_kkt(K, -g, prob.n, -rp, dense);
            dnu = nuPlus - nu;
            rn0 = norm([g + prob.Aeq' * nu; rp]);
            alpha = 1;
            while ~barrier_feasible(prob, x + alpha*dxn) && alpha > 1e-20
                alpha = alpha / 2;
            end
            while alpha > 1e-20
                xn = x + alpha * dxn;
                [~, g0n] = objective(prob, xn);
                gn = mu * g0n + barrier_gradient(prob, xn);
                rn = norm([gn + prob.Aeq' * (nu + alpha*dnu); prob.Aeq * xn - prob.beq]);
                if rn <= (1 - 0.01*alpha) * rn0, break; end
                alpha = alpha / 2;
            end
            if o.verbose > 1
                fprintf('         it %2d  [infeasible start]  |r_p|=%9.3e  alpha=%g\n', ...
                        it, max(abs(rp)), alpha);
            end
            if alpha <= 1e-20 || (it >= o.maxInfeasibleSteps)
                st.eqInfeasible = true; st.mu = mu; st.f0 = f0; st.gap = mdeg/mu;
                return;
            end
            x  = x + alpha * dxn;
            nu = nu + alpha * dnu;
            if alpha == 1 || max(abs(prob.Aeq * x - prob.beq)) <= eqTol
                eqFeas = true;
            end
            continue;
        end
        % ---------- feasible Newton step:  H dx + Aeq' nu = -g,  Aeq dx = 0 ---------
        [dx, nu] = solve_kkt(K, -g, prob.n, zeros(nEq,1), dense);
        dec2 = -(g' * dx);                         % squared Newton decrement
        if dec2 / 2 <= o.tol || ~isfinite(dec2)
            break;
        end
        alpha = 1;
        while ~barrier_feasible(prob, x + alpha*dx) && alpha > 1e-20
            alpha = alpha / 2;
        end
        if dec2 > 0.1     % outside the quadratic convergence region -> Armijo
            F = mu * f0 + phi;
            while alpha >= 1e-20
                xn = x + alpha*dx;
                if barrier_feasible(prob, xn)
                    Fn = mu * objective(prob, xn) + barrier_value(prob, xn);
                    if Fn <= F - 0.01 * alpha * dec2, break; end
                end
                alpha = alpha / 2;
            end
        end
        if o.verbose > 1
            fprintf('         it %2d  dec2=%9.3e  alpha=%g\n', it, dec2, alpha);
        end
        if alpha < 1e-20, break; end
        x = x + alpha * dx;
    end
    if ~eqFeas               % ran out of Newton steps before reaching Aeq x = beq
        st.eqInfeasible = true; st.mu = mu; st.f0 = f0; st.gap = mdeg/mu;
        return;
    end
    f0 = objective(prob, x);
    gap = mdeg / mu;
    if o.verbose
        fprintf('      barrier outer %2d  mu=%9.3e  f0=%12.6e  gap=%9.2e  newton=%d\n', ...
                outer, mu, f0, gap, st.newtonSteps);
    end
    % relative duality-gap test; the caller may supply the scale of the true
    % energy (soft mode: f0 = E + lambda*t is dominated by the penalty term)
    if isfield(prob, 'gapRef'), ref = prob.gapRef; else, ref = abs(f0); end
    if gap < o.gapTol * (ref + 1e-6)
        break;
    end
    mu = mu * o.muFactor;
end
st.mu = mu; st.f0 = f0; st.gap = mdeg/mu;
end

% =============================================================================
function o = fill_defaults(o)
d.maxIters = 60; d.maxOuter = 40; d.tol = 1e-9; d.muFactor = 20;
d.gapTol = 1e-6; d.verbose = false; d.nearTol = 1e-3; d.maxInfeasibleSteps = 30;
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(o, f{i}), o.(f{i}) = d.(f{i}); end
end
end

function [f, g] = objective(prob, x)
% f0(x) = c'x + sum_q sum_r wq_r (Gq x - yq)_r^2
f = prob.c' * x;
if nargout > 1, g = prob.c; end
for q = 1:numel(prob.quad)
    Q = prob.quad{q};
    r = Q.G * x - Q.y;
    f = f + sum(Q.w .* r.^2);
    if nargout > 1, g = g + 2 * (Q.G' * (Q.w .* r)); end
end
end

function H = objective_hessian(prob)
% constant Hessian of the quadratic part: 2 sum_q Gq' diag(wq) Gq
H = sparse(prob.n, prob.n);
for q = 1:numel(prob.quad)
    Q = prob.quad{q};
    H = H + 2 * (Q.G' * dg(Q.w) * Q.G);
end
end

function [u, v, w, Gu, d2u] = family_values(F, x)
u = F.Gu * x + F.u0;
Gu = F.Gu;
d2u = [];
if F.inv
    xa = x(F.invIdx);
    u  = u - 1 ./ xa;                                   % u = A - 1/x_aux
    if nargout > 3
        Gu  = Gu + sparse(1:F.m, F.invIdx, 1 ./ xa.^2, F.m, size(F.Gu,2));
        d2u = -2 ./ xa.^3;
    end
end
v = F.Gv * x + F.v0;
w = F.Gw * x + F.w0;
end

function D = cone_det(u, v, w)
% D = u^2 - v^2 - w^2 computed as (u - rho)(u + rho): accurate near the boundary
rho = sqrt(v.^2 + w.^2);
D = (u - rho) .* (u + rho);
end

function ok = barrier_feasible(prob, x)
ok = all(isfinite(x));
for k = 1:numel(prob.lin)
    if ~ok, return; end
    ok = all(prob.lin{k}.G * x + prob.lin{k}.g0 > 0);
end
for k = 1:numel(prob.fams)
    if ~ok, return; end
    F = prob.fams{k};
    if F.inv && any(x(F.invIdx) <= 0), ok = false; return; end
    [u, v, w] = family_values(F, x);
    ok = all(u > 0) && all(cone_det(u, v, w) > 0);
end
end

function phi = barrier_value(prob, x)
phi = 0;
for k = 1:numel(prob.lin)
    a = prob.lin{k}.G * x + prob.lin{k}.g0;
    if any(a <= 0), phi = Inf; return; end
    phi = phi - sum(log(a));
end
for k = 1:numel(prob.fams)
    [u, v, w] = family_values(prob.fams{k}, x);
    D = cone_det(u, v, w);
    if any(D <= 0) || any(u <= 0), phi = Inf; return; end
    phi = phi - sum(log(D));
end
end

function g = barrier_gradient(prob, x)
g = zeros(prob.n, 1);
for k = 1:numel(prob.lin)
    L = prob.lin{k};
    g = g - L.G' * (1 ./ (L.G * x + L.g0));
end
for k = 1:numel(prob.fams)
    F = prob.fams{k};
    [u, v, w, Gu] = family_values(F, x);
    D = cone_det(u, v, w);
    g = g - (Gu' * (2*u ./ D) - F.Gv' * (2*v ./ D) - F.Gw' * (2*w ./ D));
end
end

function [phi, g, K] = newton_system(prob, x, mu, H0, nearTol)
% Barrier value and gradient, and the augmented Newton matrix K of
% F_mu = mu*f0 + phi (see the header). K's last block rows are Aeq.
n  = prob.n;
phi = 0;
g  = zeros(n, 1);
Hd = mu * H0;
Gblk = {}; Pblk = {};
% ---- linear inequalities  a = G x + g0 > 0:   augmented, P^{-1} = diag(a.^2) --------
for k = 1:numel(prob.lin)
    L = prob.lin{k};
    a = L.G * x + L.g0;
    phi = phi - sum(log(a));
    g   = g - L.G' * (1 ./ a);
    Gblk{end+1} = L.G;                         %#ok<AGROW>
    Pblk{end+1} = dg(a.^2);                    %#ok<AGROW>
end
% ---- second-order cones  y = (u, v, w) ------------------------------------------------
for k = 1:numel(prob.fams)
    F = prob.fams{k};
    [u, v, w, Gu, d2u] = family_values(F, x);
    D = cone_det(u, v, w);
    phi = phi - sum(log(D));
    % gradient:  -(2u grad u - 2v grad v - 2w grad w) / D
    g = g - (Gu' * (2*u ./ D) - F.Gv' * (2*v ./ D) - F.Gw' * (2*w ./ D));
    % curvature of the nonlinear u = A - 1/x_aux:  -(2 u u''/D) on the aux diagonal
    if F.inv
        Hd = Hd + sparse(F.invIdx, F.invIdx, -2 * u .* d2u ./ D, n, n);
    end
    near = D < nearTol * u.^2;                 % close to the cone boundary / apex
    far  = ~near;
    if any(far)        % ---- normal equations: G' P G with P = (2/D^2)(2Jy(Jy)' - DJ)
        s  = 2 ./ D(far).^2;
        uf = u(far); vf = v(far); wf = w(far); Df = D(far);
        Gu_ = Gu(far,:); Gv_ = F.Gv(far,:); Gw_ = F.Gw(far,:);
        Hd = Hd + Gu_' * dg(s .* (uf.^2 + vf.^2 + wf.^2)) * Gu_ ...
                + Gv_' * dg(s .* (2*vf.^2 + Df)) * Gv_ ...
                + Gw_' * dg(s .* (2*wf.^2 + Df)) * Gw_;
        Cuv = Gu_' * dg(-2 * s .* uf .* vf) * Gv_;
        Cuw = Gu_' * dg(-2 * s .* uf .* wf) * Gw_;
        Cvw = Gv_' * dg( 2 * s .* vf .* wf) * Gw_;
        Hd = Hd + Cuv + Cuv' + Cuw + Cuw' + Cvw + Cvw';
    end
    if any(near)       % ---- augmented rows: P^{-1} = y y' - (D/2) J, ordering [u; v; w]
        un = u(near); vn = v(near); wn = w(near); Dn = D(near);
        m = numel(un); r = (1:m)';
        I = [r; r+m; r+2*m; r; r; r+m; r+m; r+2*m; r+2*m];
        J = [r; r+m; r+2*m; r+m; r+2*m; r; r+2*m; r; r+m];
        V = [un.^2 - Dn/2; vn.^2 + Dn/2; wn.^2 + Dn/2; un.*vn; un.*wn; un.*vn; vn.*wn; un.*wn; vn.*wn];
        Gblk{end+1} = [Gu(near,:); F.Gv(near,:); F.Gw(near,:)];   %#ok<AGROW>
        Pblk{end+1} = sparse(I, J, V, 3*m, 3*m);                   %#ok<AGROW>
    end
end
if isempty(Gblk)
    G = sparse(0, n); Pinv = sparse(0, 0);
else
    G = vertcat(Gblk{:}); Pinv = blkdiag(Pblk{:});
end
na = size(G, 1); nEq = size(prob.Aeq, 1);
K = [Hd,        G',               prob.Aeq';
     G,         -Pinv,            sparse(na, nEq);
     prob.Aeq,  sparse(nEq, na),  sparse(nEq, nEq)];
end

function [z, nu] = solve_kkt(K, r, n, req, dense)
% Solve  K [z; y; nu] = [r; 0; req]   and return z and the multipliers nu.
%
% Factorisation: MATLAB - sparse symmetric indefinite LDL' (MA57, MC64 scaling,
% AMD ordering), Octave - UMFPACK LU (no sparse LDL' there). Iterative
% refinement with the same factors follows: K's entries span many orders of
% magnitude by design (P^{-1} ~ 1e-12 for cones at their apex), so refinement
% is what guarantees an accurate Newton direction. Condition-number warnings
% are meaningless here and silenced.
%
% DENSE VARIABLES (the global epigraph variable s of V1/V4 and the L_inf
% penalty variable t of the soft mode) couple with every face / every
% curvature row and would destroy the sparsity of the factors. They are
% eliminated by a bordered (Schur complement) solve:
%     [Kss Ksd; Kds Kdd] [xs; xd] = [bs; bd]
%     (Kdd - Kds Kss^{-1} Ksd) xd = bd - Kds Kss^{-1} bs,   xs = Kss^{-1}(bs - Ksd xd)
% which costs one sparse factorisation of Kss and (1 + #dense) solves.
if nargin < 5, dense = []; end
N = size(K, 1); nEq = numel(req);
rhs = [r; zeros(N - n - nEq, 1); req];
ws = warning('off', 'all');
if isempty(dense)
    solve = factorize(K);
else
    sp  = true(N, 1); sp(dense) = false; sp = find(sp);
    solveS = factorize(K(sp, sp));
    B  = K(sp, dense);
    Z  = solveS(full(B));
    Sc = full(K(dense, dense)) - full(B' * Z);
    solve = @(b) bordered(b, sp, dense, solveS, B, Z, Sc);
end
sol = solve(rhs);
nb = norm(rhs);
for k = 1:5
    res = rhs - K * sol;
    if norm(res) <= 1e-14 * nb, break; end
    sol = sol + solve(res);
end
if any(~isfinite(sol))
    delta = 1e-14 * max(1, full(max(abs(diag(K)))));
    sol = (K + blkdiag(delta * speye(n), -delta * speye(N - n))) \ rhs;
end
warning(ws);
z  = sol(1:n);
nu = sol(N-nEq+1:N);
end

function x = bordered(b, sp, dense, solveS, B, Z, Sc)
y  = solveS(b(sp));
xd = Sc \ (b(dense) - B' * y);
x  = zeros(size(b));
x(sp) = y - Z * xd;
x(dense) = xd;
end

function solve = factorize(K)
persistent isOctave
if isempty(isOctave), isOctave = exist('OCTAVE_VERSION', 'builtin') > 0; end
solve = [];
if ~isOctave
    try
        [L, Dm, P, S] = ldl((K + K') / 2);     % P' S K S P = L D L'  (MA57)
        solve = @(b) S * (P * (L' \ (Dm \ (L \ (P' * (S * b))))));
    catch
        solve = [];
    end
end
if isempty(solve)
    [L, U, P, Q, R] = lu(K);                   % P (R \ K) Q = L U   (UMFPACK)
    solve = @(b) Q * (U \ (L \ (P * (R \ b))));
end
end

function S = dg(d)
S = spdiags(d(:), 0, numel(d), numel(d));
end
