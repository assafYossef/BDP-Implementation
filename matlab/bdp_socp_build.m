function cp = bdp_socp_build(e, data, bnd, variant, J, b, lambda, mode, dropRedundant)
%BDP_SOCP_BUILD The SCP subproblem as a standard second-order cone program.
%
%   cp = bdp_socp_build(e, data, bnd, variant, J, b, lambda, mode, dropRedundant)
%
%   Builds the convex subproblem of Section 7 (problem (11), or its soft form of
%   Section 7.2) in the solver-independent conic form
%
%       minimize    c' x
%       subject to  Aeq x  = beq
%                   Ain x <= bin
%                   lb <= x <= ub
%                   F_k x + g_k  in  Q^{d_k}   (second-order cones, k = 1..numel(dims))
%
%   where Q^d = { y in R^d : y(1) >= norm(y(2:d)) }. The rows of F/g are the
%   cones stacked one after the other (dims(k) rows each, head row first): this
%   is exactly the "affine conic constraint" layout of MOSEK (prob.f, prob.g,
%   prob.accs) and is converted to secondordercone objects for coneprog.
%
%   Variables  x = [ e (|E|) ; z (|F|) ; aux... ]
%     z_i >= 2|eta_i| = norm(v_i, w_i)              cone [z_i; v_i; w_i]
%   With z, every distortion quantity becomes affine or a rotated cone:
%     sigma_a^2 <= A + z,   sigma_b^2 >= A - z     (exact at the optimum: every
%                                                   use of z is "worse when bigger")
%     1/sigma_b^2 <= q   <=>  (A - z) q >= 1        rotated cone, written as the
%                                                   quadratic cone [A-z+q; 2; A-z-q]
%   Distortion bounds (6)-(8) are then LINEAR in (e, z):
%     A + z <= Ca,   A - z >= Cb,   z <= (CK-1)/(CK+1) A
%
%   Objectives (Section 8.1):
%     V1   min t,                       A + z <= t,  q <= t
%     V1R  min t + avgWeight/2 sum w (A + z + q)        (reference script, method 50)
%     V2   min sum w (A + z + q)
%     V3   min sum w [2(A-1)^2 + 2 v^2 + 2 w^2]   -> epigraph s with one rotated cone
%     V4   min s,                       z <= 2 s
%     V5   min sum w [(A-A0)^2 + (v-v0)^2 + (w-w0)^2]   -> epigraph s
%   Curvature (Section 7.1/7.2):
%     mode 'hard':  J e = b
%     mode 'soft':  -r <= J e - b <= r,   objective + lambda * r
%
%   cp.idxE gives the position of e inside x.

if nargin < 9, dropRedundant = true; end
nE = data.nE; nF = data.nF; nC = size(J, 1);
GA = data.GA; GV = data.GV; GW = data.GW;      % |F| x |E|
w = data.wF(:);
soft = strcmp(mode, 'soft');

% ---------------- variable layout -------------------------------------------------
iz = nE + (1:nF)';            n = nE + nF;
useQ = any(strcmp(variant, {'V1','V1R','V2'}));
useT = any(strcmp(variant, {'V1','V1R','V4'}));
useS = any(strcmp(variant, {'V3','V5'}));
iq = []; it = []; is = []; ir = [];
if useQ, iq = n + (1:nF)'; n = n + nF; end
if useT, it = n + 1;       n = n + 1;  end
if useS, is = n + 1;       n = n + 1;  end
if soft && nC > 0, ir = n + 1; n = n + 1; end

Z  = @(m) sparse(m, n - nE);
Ae = [GA, Z(nF)];                                  % A_i   as rows over x
Ve = [GV, Z(nF)];                                  % 2Re eta
We = [GW, Z(nF)];                                  % 2Im eta
Sz = sparse(1:nF, iz, 1, nF, n);                   % z_i
c  = zeros(n, 1);
Ain = {}; bin = {};
Fc = {};  gc = {};  dims = {};

% ---------------- cone  z_i >= norm(v_i, w_i) ------------------------------------
Fc{end+1} = interleave({Sz, Ve, We});  gc{end+1} = zeros(3*nF, 1);  dims{end+1} = 3*ones(nF,1);

% ---------------- rotated cones  (A - z) q >= 1 -----------------------------------
if useQ
    Sq = sparse(1:nF, iq, 1, nF, n);
    head = Ae - Sz + Sq;   tail2 = Ae - Sz - Sq;
    Fc{end+1} = interleave({head, sparse(nF, n), tail2});
    gc{end+1} = reshape([zeros(1,nF); 2*ones(1,nF); zeros(1,nF)], [], 1);
    dims{end+1} = 3*ones(nF,1);
end

% ---------------- objective -------------------------------------------------------
switch variant
    case {'V1','V1R'}
        St = sparse(1:nF, it, 1, nF, n);
        Ain{end+1} = Ae + Sz - St;  bin{end+1} = zeros(nF,1);      % sigma_a^2 <= t
        Ain{end+1} = Sq - St;       bin{end+1} = zeros(nF,1);      % 1/sigma_b^2 <= t
        c(it) = 1;
        if strcmp(variant, 'V1R')
            aw = 0.5 * data.avgWeight;
            c = c + aw * ((w' * Ae)' + (w' * Sz)' + (w' * Sq)');
        end
    case 'V2'
        c = (w' * Ae)' + (w' * Sz)' + (w' * Sq)';
    case 'V4'
        St = sparse(1:nF, it, 2, nF, n);
        Ain{end+1} = Sz - St;  bin{end+1} = zeros(nF,1);           % 2|eta| <= 2 s
        c(it) = 1;
    case {'V3','V5'}
        if strcmp(variant, 'V3')
            R = [spdiags(sqrt(2*w),0,nF,nF) * Ae; spdiags(sqrt(2*w),0,nF,nF) * Ve; spdiags(sqrt(2*w),0,nF,nF) * We];
            r0 = [sqrt(2*w) .* ones(nF,1); zeros(2*nF,1)];          % residual = R x - r0
        else
            R = [spdiags(sqrt(w),0,nF,nF) * Ae; spdiags(sqrt(w),0,nF,nF) * Ve; spdiags(sqrt(w),0,nF,nF) * We];
            r0 = [sqrt(w) .* data.A0; sqrt(w) .* data.v0; sqrt(w) .* data.w0];
        end
        % ||R x - r0||^2 <= s   <=>   [1 + s; 1 - s; 2 (R x - r0)] in Q
        Ss = sparse(1, is, 1, 1, n);
        Fc{end+1} = [Ss; -Ss; 2*R];  gc{end+1} = [1; 1; -2*r0];  dims{end+1} = size(R,1) + 2;
        c(is) = 1;
    otherwise
        error('bdp:variant', 'Unknown variant %s', variant);
end

% ---------------- distortion bounds (linear in e, z) ------------------------------
if isfinite(bnd.Ca), Ain{end+1} = Ae + Sz;  bin{end+1} = bnd.Ca * ones(nF,1); end
if bnd.Cb > 0,       Ain{end+1} = Sz - Ae;  bin{end+1} = -bnd.Cb * ones(nF,1); end
if isfinite(bnd.CK), Ain{end+1} = Sz - ((bnd.CK-1)/(bnd.CK+1)) * Ae;  bin{end+1} = zeros(nF,1); end
if ~useQ && ~(bnd.Cb > 0) && ~isfinite(bnd.CK)
    % keep the metric positive definite (=> triangle inequalities, Section 4.3)
    Ain{end+1} = Sz - Ae;  bin{end+1} = -1e-8 * ones(nF,1);
end

% ---------------- linearised curvature --------------------------------------------
Aeq = sparse(0, n); beq = zeros(0, 1);
if nC > 0
    JJ = [J, sparse(nC, n - nE)];
    if ~soft
        Aeq = JJ; beq = b;
        if dropRedundant && numel(data.Vc) == data.nV
            Aeq = Aeq(1:end-1, :); beq = beq(1:end-1);   % Gauss-Bonnet: one row is redundant
        end
    else
        Sr = sparse(1:nC, ir, 1, nC, n);
        Ain{end+1} = JJ - Sr;   bin{end+1} = b;          %  J e - b <= r
        Ain{end+1} = -JJ - Sr;  bin{end+1} = -b;         % -(J e - b) <= r
        c(ir) = lambda;
    end
end

cp.n = n; cp.c = c;
cp.Aeq = Aeq; cp.beq = beq;
if isempty(Ain), cp.Ain = sparse(0, n); cp.bin = zeros(0,1);
else, cp.Ain = vertcat(Ain{:}); cp.bin = vertcat(bin{:}); end
cp.F = vertcat(Fc{:}); cp.g = vertcat(gc{:}); cp.dims = vertcat(dims{:});
cp.lb = -inf(n, 1); cp.ub = inf(n, 1);
cp.idxE = (1:nE)';
cp.e0 = e;                                         % current iterate (for warm starts)
cp.mode = mode;
end

function M = interleave(blocks)
% Stack k blocks (each m x n) so that row j of every block is adjacent:
% [B1(1,:); B2(1,:); ...; B1(2,:); B2(2,:); ...]  -> cones of dimension k.
k = numel(blocks); m = size(blocks{1}, 1);
M = vertcat(blocks{:});
perm = reshape(reshape(1:k*m, m, k)', [], 1);
M = M(perm, :);
end
