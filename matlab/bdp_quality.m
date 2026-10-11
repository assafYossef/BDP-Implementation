function Q = bdp_quality(e, data)
%BDP_QUALITY Quality summary of a metric e (squared edge lengths).
%
%   Q = bdp_quality(e, data)      (data as passed to bdp_solve)
%
%   Everything derives from the two stretches of each triangle, the singular
%   values s1 = sigma_a >= s2 = sigma_b of the Jacobian (1 = unchanged,
%   1.3 = stretched by 30%, 0.8 = squashed to 80%):
%
%     Q.tau          max_i max(s1, 1/s2)              worst stretch            perfect 1
%                    (= sqrt of the paper's tau = max(sigma_a^2, 1/sigma_b^2))
%     Q.tauSc        min_c max_i max(c s1, 1/(c s2))  worst stretch after the  perfect 1
%                    = sqrt(max s1 / min s2)          best global rescaling
%     Q.bestScale    the optimal c = 1/sqrt(max s1 * min s2)
%     Q.Kmax         max_i s1/s2                      worst shape change       perfect 1
%     Q.symDir       sum_i w_i (s1^2+s2^2+1/s1^2+1/s2^2), w_i = original area
%                    fraction (sum 1)                 average distortion       perfect 4
%     Q.symDirExcess symDir - 4                                                perfect 0
%     Q.maxCurvErrDeg max over constrained vertices of |K - K*| in degrees     perfect 0
%     Q.valid        every triangle is a real triangle (strict triangle
%                    inequality, s2 > 0)                                       perfect 1
%     Q.degenerateTriangles  number of faces violating that
%
%   Scaling the whole map by c multiplies s1 and s2 by c; max(c*max s1, 1/(c*min s2))
%   is smallest when both terms are equal, which gives the closed form for tauSc.
%   (Flipped triangles of the planar layout are checked in C++, see Layout.cpp.)

if ~isfield(data, 'GA')                     % allow calls with raw (C++) data
    data.nF = size(data.F, 1); data.nE = numel(data.e0);
    data.nV = max([max(data.F(:)), numel(data.isBnd)]);
end
[~, T, okTri] = bdp_angles(e, data.F2E);
D  = bdp_distortion(e, data);
s1 = D.sigmaA; s2 = D.sigmaB;

bad = ~(T > 0) | ~(s2 > 0) | ~isfinite(s1);
Q.degenerateTriangles = nnz(bad);
Q.valid = okTri && Q.degenerateTriangles == 0;

Q.tau  = max(max(s1, 1 ./ s2));
Q.Kmax = max(s1 ./ s2);
if min(s2) > 0
    Q.bestScale = 1 / sqrt(max(s1) * min(s2));
    Q.tauSc = sqrt(max(s1) / min(s2));
else
    Q.bestScale = 1; Q.tauSc = Inf;
end
Q.symDir = sum(data.wF(:) .* (s1.^2 + s2.^2 + 1 ./ s1.^2 + 1 ./ s2.^2));
Q.symDirExcess = Q.symDir - 4;

K = bdp_curvature(e, data);
if isempty(data.Vc)
    Q.maxCurvErrDeg = 0;
else
    Q.maxCurvErrDeg = max(abs(K(data.Vc(:)) - data.K0(:))) * 180 / pi;
end
end
