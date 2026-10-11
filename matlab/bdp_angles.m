function [theta, T, valid] = bdp_angles(e, F2E)
%BDP_ANGLES Corner angles of every triangle from the squared edge lengths.
%
%   e     : |E| x 1 squared edge lengths (the discrete metric, our variables).
%   F2E   : |F| x 3, F2E(f,k) = index of the edge OPPOSITE to corner k of face f.
%
%   theta : |F| x 3 corner angles, theta(f,k) is the angle at corner k.
%   T     : |F| x 1, T = sqrt(2ab + 2bc + 2ca - a^2 - b^2 - c^2) = 4 * Area
%           (Heron's formula written in squared lengths, Section 4.3 of the paper).
%   valid : true iff every triangle satisfies the (strict) triangle inequality.
%
%   Law of cosines in squared lengths (Eq. (4) of the paper):
%       cos(theta_k) = (s_j + s_i - s_k) / (2 sqrt(s_i s_j)),
%   and since sin(theta_k) = 2 Area / sqrt(s_i s_j) = T / (2 sqrt(s_i s_j)),
%       theta_k = atan2(T, s_i + s_j - s_k),
%   which is numerically far more robust than acos for needle triangles.

s  = e(F2E);                       % s(f,k) = squared length of edge opposite corner k
a  = s(:,1); b = s(:,2); c = s(:,3);
T2 = 2*(a.*b + b.*c + c.*a) - a.^2 - b.^2 - c.^2;
valid = all(T2 > 0) && all(s(:) > 0);
T  = sqrt(max(T2, 0));

theta = [atan2(T, b + c - a), ...
         atan2(T, c + a - b), ...
         atan2(T, a + b - c)];
end
