function [e, ok, err] = bdp_cetm_newton(e0, data, opts)
%BDP_CETM_NEWTON Finish the flattening with a conformal (CETM) Newton solve.
%
%   Section 8.2 of the paper: once the SCP iterates are within 1 degree of the
%   prescribed curvature, the metric is made exactly flat by Conformal
%   Equivalence of Triangle Meshes [Springborn, Schroeder & Pinkall 2008].
%
%   Conformal change of the metric by vertex log-scale factors u:
%       l_ij <- exp((u_i + u_j)/2) l_ij      <=>      e_ij <- exp(u_i + u_j) e_ij
%   (squared lengths, so the change is LINEAR in log space).
%
%   We solve K(u) = K0 on the constrained vertices with Newton's method. CETM
%   shows K(u) - K0 is the gradient of a convex energy, whose Hessian is the
%   cotangent Laplacian (/2). Using the chain rule with the Jacobian of
%   Appendix C, the Hessian is obtained exactly:
%       dK/du = JK(e) * diag(e) * B,      B(edge, v) = 1 if v is an end point,
%   since d e_ij / d u_i = e_ij. Unconstrained vertices (free boundary) keep u=0.
%   If all vertices are constrained, the energy is invariant to a global shift
%   of u (Gauss-Bonnet), so one vertex is pinned (u = 0) to remove the kernel.
%
%   Line search: backtracking on ||K - K0||_2 while keeping all triangle
%   inequalities valid (CETM fails exactly when they break, Section 9.2).

if nargin < 3, opts = bdp_default_options(); end
e0 = e0(:);
Vc = data.Vc(:); K0 = data.K0(:);
nV = data.nV;
free = Vc;
K0f  = K0;
if numel(Vc) == nV                       % closed surface, all curvatures prescribed
    free = Vc(2:end); K0f = K0(2:end);
end
B = sparse([(1:data.nE)'; (1:data.nE)'], [data.EV(:,1); data.EV(:,2)], 1, data.nE, nV);

u = zeros(nV, 1);
e = e0;
[K, ~, T, valid] = bdp_curvature(e, data);
ok = false;
if ~valid, err = Inf; return; end
res = K(free) - K0f;
for it = 1:opts.cetmMaxIters
    nr = norm(res);
    if max(abs(res)) < 1e-12, break; end
    JK = bdp_curvature_jacobian(e, data, T);
    H  = JK(free, :) * spdiags(e, 0, data.nE, data.nE) * B(:, free);   % dK/du
    H  = (H + H') / 2;
    du = zeros(nV, 1);
    du(free) = -(H \ res);
    alpha = 1; improved = false;
    while alpha > 1e-10
        un = u + alpha * du;
        en = e0 .* exp(B * un);                % e_ij = exp(u_i + u_j) e0_ij
        [Kn, ~, Tn, vn] = bdp_curvature(en, data);
        if vn
            rn = Kn(free) - K0f;
            if norm(rn) < nr, improved = true; break; end
        end
        alpha = alpha / 2;
    end
    if ~improved, break; end
    u = un; e = en; T = Tn; res = rn;
end
K = bdp_curvature(e, data);
err = max(abs(K(Vc) - K0));
ok = err < max(opts.tolFlat, 1e-10);
end
