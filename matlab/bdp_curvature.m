function [K, theta, T, valid] = bdp_curvature(e, data)
%BDP_CURVATURE Discrete Gaussian curvature (angle defect) at every vertex.
%
%   Section 5 of the paper:
%       K(v) = 2*pi - alpha_v   (interior vertex)
%       K(v) =   pi - alpha_v   (boundary vertex)
%   where alpha_v is the sum of the corner angles incident to v, computed
%   from the squared edge lengths by the law of cosines (bdp_angles).

[theta, T, valid] = bdp_angles(e, data.F2E);
alpha = accumarray(data.F(:), theta(:), [data.nV, 1]);   % angle sums
K = (2*pi - pi*data.isBnd(:)) - alpha;
end
