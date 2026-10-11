function [x, st] = bdp_socp_solve(cp, solver, solverOpts)
%BDP_SOCP_SOLVE Solve a conic subproblem built by bdp_socp_build.
%
%   [x, st] = bdp_socp_solve(cp, solver)
%
%   solver : 'mosek'     MOSEK's interior point optimizer through mosekopt (the
%                        solver used by the reference implementation via CVX).
%                        MOSEK >= 10 uses affine conic constraints (prob.f/g/accs);
%                        MOSEK 8/9 uses cones on auxiliary variables.
%            'cvx'       CVX with its currently selected solver (cvx_solver, e.g.
%                        MOSEK) - the route used by the reference script. Slower
%                        than calling MOSEK directly (CVX re-models every problem)
%                        but works with CVX's bundled, CVX-only MOSEK licence.
%            'coneprog'  MATLAB Optimization Toolbox (R2020b+) primal-dual
%                        interior point second-order cone solver.
%            function handle  @(cp) -> [x, status] for any other conic solver.
%
%   st.status : 'optimal' | 'infeasible' | 'failed'
%   st.time, st.iterations, st.message
%
%   Both MOSEK and coneprog are primal-dual interior point methods with
%   Nesterov-Todd scaling and Mehrotra predictor-corrector steps: they need
%   ~20-50 factorisations per subproblem, compared with the hundreds of Newton
%   steps of the pure log-barrier path following in bdp_barrier_newton.

if nargin < 3, solverOpts = struct(); end
t0 = tic;
st.iterations = NaN; st.message = '';
if isa(solver, 'function_handle')
    [x, status] = solver(cp);
    st.status = status;
else
    switch lower(solver)
        case 'mosek',    [x, st] = solve_mosek(cp, solverOpts, st);
        case 'coneprog', [x, st] = solve_coneprog(cp, solverOpts, st);
        case 'cvx',      [x, st] = solve_cvx(cp, solverOpts, st);
        otherwise, error('bdp:solver', 'Unknown conic solver "%s"', solver);
    end
end
if strcmp(st.status, 'optimal') && (isempty(x) || any(~isfinite(x)))
    st.status = 'failed';
end
st.time = toc(t0);
end

% =============================================================================
function [x, st] = solve_mosek(cp, so, st)
persistent sc
if isempty(sc)
    [~, res] = mosekopt('symbcon echo(0)');
    sc = res.symbcon;
end
n = cp.n;
prob.c   = cp.c;
prob.a   = [cp.Aeq; cp.Ain];
if isempty(prob.a), prob.a = sparse(0, n); end
prob.blc = [cp.beq; -inf(size(cp.bin))];
prob.buc = [cp.beq; cp.bin];
prob.blx = cp.lb;
prob.bux = cp.ub;
K = numel(cp.dims);
if isfield(sc, 'MSK_DOMAIN_QUADRATIC_CONE')
    % MOSEK >= 10: affine conic constraints  F x + g in Q^{d1} x Q^{d2} x ...
    prob.f = cp.F;
    prob.g = cp.g;
    prob.accs = reshape([repmat(double(sc.MSK_DOMAIN_QUADRATIC_CONE), 1, K); cp.dims(:)'], 1, []);
else
    % MOSEK 8/9: cones on variables -> auxiliary y = F x + g
    m = size(cp.F, 1);
    prob.c   = [cp.c; zeros(m, 1)];
    prob.a   = [prob.a, sparse(size(prob.a, 1), m); cp.F, -speye(m)];
    prob.blc = [prob.blc; -cp.g];
    prob.buc = [prob.buc; -cp.g];
    prob.blx = [cp.lb; -inf(m, 1)];
    prob.bux = [cp.ub;  inf(m, 1)];
    prob.cones.type   = repmat(double(sc.MSK_CT_QUAD), 1, K);
    prob.cones.sub    = n + (1:m);
    prob.cones.subptr = cumsum([1; cp.dims(1:end-1)])';
end
param = struct();
if isfield(so, 'mosekParam'), param = so.mosekParam; end
[rcode, res] = mosekopt('minimize echo(0) info', prob, param);
x = []; st.status = 'failed';
if isfield(res, 'sol') && isfield(res.sol, 'itr')
    sol = res.sol.itr;
    x = sol.xx(1:n);
    solsta = sol.solsta; prosta = sol.prosta;
    st.message = sprintf('%s / %s (rcode %d)', solsta, prosta, rcode);
    if any(strcmp(solsta, {'OPTIMAL', 'NEAR_OPTIMAL'}))
        st.status = 'optimal';
    elseif ~isempty(strfind(solsta, 'INFEASIBLE')) || ~isempty(strfind(prosta, 'INFEASIBLE'))
        st.status = 'infeasible';
    elseif strcmp(solsta, 'UNKNOWN') && all(isfinite(x))
        st.status = 'optimal';                       % stalled close to optimal: accept
    end
end
if isfield(res, 'info') && isfield(res.info, 'MSK_IINF_INTPNT_ITER')
    st.iterations = res.info.MSK_IINF_INTPNT_ITER;
end
end

% =============================================================================
function [x, st] = solve_coneprog(cp, so, st)
K = numel(cp.dims);
ptr = cumsum([1; cp.dims(:)]);
Ft = cp.F';                                          % column slicing is fast in CSC
socs = cell(1, K);
denseD = false;
for k = 1:K
    head = ptr(k);
    tail = ptr(k)+1:ptr(k+1)-1;
    % ||F_tail x + g_tail|| <= F_head x + g_head   <=>   ||A x - b|| <= d' x - gamma
    d = Ft(:, head);
    if ~denseD
        try
            socs{k} = secondordercone(Ft(:, tail)', -cp.g(tail), d, -cp.g(head));
            continue;
        catch
            denseD = true;                           % this MATLAB release wants a full d
        end
    end
    socs{k} = secondordercone(Ft(:, tail)', -cp.g(tail), full(d), -cp.g(head));
end
socs = [socs{:}];
opts = optimoptions('coneprog', 'Display', 'off');
if isfield(so, 'coneprogOptions'), opts = so.coneprogOptions; end
lb = cp.lb; ub = cp.ub;
if all(isinf(lb)), lb = []; end
if all(isinf(ub)), ub = []; end
[x, ~, exitflag, output] = coneprog(cp.c, socs, cp.Ain, cp.bin, cp.Aeq, cp.beq, lb, ub, opts);
st.iterations = output.iterations;
st.message = sprintf('exitflag %d: %s', exitflag, output.message);
if exitflag == 1
    st.status = 'optimal';
elseif exitflag == -2
    st.status = 'infeasible';
elseif exitflag == 0 && ~isempty(x) && all(isfinite(x))
    st.status = 'optimal';                           % iteration limit, usable point
else
    st.status = 'failed';
end
end

% =============================================================================
function [x, st] = solve_cvx(cp, so, st)
% Same conic data, stated as a CVX model:  ||tail_k|| <= head_k  for every cone.
ptr = cumsum([1; cp.dims(:)]);
s3 = find(cp.dims == 3);  h3 = ptr(s3);              % the many 3-dimensional cones
big = find(cp.dims ~= 3);                            % e.g. the V3/V5 energy cone
Fh = cp.F(h3, :);   gh = cp.g(h3);
F1 = cp.F(h3+1, :); g1 = cp.g(h3+1);
F2 = cp.F(h3+2, :); g2 = cp.g(h3+2);
n = cp.n;
hasLb = any(isfinite(cp.lb)); hasUb = any(isfinite(cp.ub));
cvx_begin quiet
    if isfield(so, 'cvxSolver') && ~isempty(so.cvxSolver)
        cvx_solver(so.cvxSolver);
    end
    variable x(n)
    minimize( cp.c' * x )
    subject to
        if ~isempty(cp.Aeq), cp.Aeq * x == cp.beq; end
        if ~isempty(cp.Ain), cp.Ain * x <= cp.bin; end
        if ~isempty(s3), norms([F1 * x + g1, F2 * x + g2], 2, 2) <= Fh * x + gh; end
        for k = big(:)'
            rows = ptr(k)+1:ptr(k+1)-1;
            norm(cp.F(rows, :) * x + cp.g(rows)) <= cp.F(ptr(k), :) * x + cp.g(ptr(k));
        end
        if hasLb, x(isfinite(cp.lb)) >= cp.lb(isfinite(cp.lb)); end
        if hasUb, x(isfinite(cp.ub)) <= cp.ub(isfinite(cp.ub)); end
cvx_end
st.message = cvx_status;
if exist('cvx_slvitr', 'var'), st.iterations = cvx_slvitr; end
if any(strcmp(cvx_status, {'Solved', 'Inaccurate/Solved'}))
    st.status = 'optimal';
elseif ~isempty(strfind(cvx_status, 'Infeasible'))
    st.status = 'infeasible';
else
    st.status = 'failed';
end
x = full(x);
end
