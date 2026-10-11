function JK = bdp_curvature_jacobian(e, data, T)
%BDP_CURVATURE_JACOBIAN Sparse Jacobian dK/de of ALL vertex curvatures w.r.t. the
%squared edge lengths (Appendix C of the paper).
%
%   JK is |V| x |E|; the rows of the constrained vertices are JK(data.Vc,:).
%
%   Derivation (per triangle, corner k with opposite squared length s_k and the
%   two other corners i, j, with (k,i,j) a cyclic permutation):
%
%     theta_k = atan2(T, s_i + s_j - s_k),  T = 4*Area
%
%     d theta_k / d s_k = 1 / T                              (= 1/(4 Area))
%     d theta_k / d s_i = -cot(theta_j) / (2 s_i)
%     d theta_k / d s_j = -cot(theta_i) / (2 s_j)
%
%   with cot(theta_j) = (s_k + s_i - s_j) / T. Since K(v) = const - sum theta,
%   dK/ds = -sum over the corners at v of d theta/ds. This reproduces the two
%   expressions of Appendix C:
%     dK/da_i = -1/(4 Area_i),
%     dK/db_i = (a_i + b_i - b_{i+1})/(2 b_i T_i) + (a_{i-1} + b_i - b_{i-1})/(2 b_i T_{i-1}).
%
%   Scale invariance K(c e) = K(e) implies JK * e = 0 (used in Eq. (13)).

if nargin < 3
    [~, T] = bdp_angles(e, data.F2E);
end
s = e(data.F2E);
nF = size(s,1);
cyc = [1 2 3; 2 3 1; 3 1 2];      % rows: (k, i, j)

I = zeros(nF, 9); J = zeros(nF, 9); V = zeros(nF, 9);
col = 0;
for r = 1:3
    k = cyc(r,1); i = cyc(r,2); j = cyc(r,3);
    sk = s(:,k); si = s(:,i); sj = s(:,j);
    cot_j = (sk + si - sj) ./ T;          % angle at corner j (opposite s_j)
    cot_i = (sk + sj - si) ./ T;          % angle at corner i (opposite s_i)
    dth_dsk = 1 ./ T;
    dth_dsi = -cot_j ./ (2*si);
    dth_dsj = -cot_i ./ (2*sj);
    vk = data.F(:,k);
    col = col + 1; I(:,col) = vk; J(:,col) = data.F2E(:,k); V(:,col) = -dth_dsk;
    col = col + 1; I(:,col) = vk; J(:,col) = data.F2E(:,i); V(:,col) = -dth_dsi;
    col = col + 1; I(:,col) = vk; J(:,col) = data.F2E(:,j); V(:,col) = -dth_dsj;
end
JK = sparse(I(:), J(:), V(:), data.nV, data.nE);
end
