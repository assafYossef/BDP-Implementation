function D = bdp_distortion(e, data)
%BDP_DISTORTION Per-triangle distortion of the metric e w.r.t. the original metric.
%
%   For every triangle the pull-back metric tensor M_i = J_i' J_i (J_i the
%   Jacobian of the affine map from the original triangle to the triangle with
%   squared lengths e) is a LINEAR function of the three squared edge lengths
%   (Section 4.3). C++ precomputes, per face, the row vectors cA, cV, cW with
%
%       A_i      = cA_i . s_i           ( = (M11+M22)/2 = ||J||_F^2 / 2 )
%       2Re(eta) = cV_i . s_i           ( = (M11-M22)/2 )
%       2Im(eta) = cW_i . s_i           ( = M12 )
%
%   so that  M_i = A_i I + [2Re eta, 2Im eta; 2Im eta, -2Re eta]  (Section 3.2.1).
%   The eigenvalues of M are the squared singular values of J (Eqs. (1)-(3)):
%       sigma_a^2 = A + 2|eta|,   sigma_b^2 = A - 2|eta|,   K^2 = sigma_a^2 / sigma_b^2.

s = e(data.F2E);
D.A   = sum(data.cA .* s, 2);
D.v   = sum(data.cV .* s, 2);          % 2 Re(eta)
D.w   = sum(data.cW .* s, 2);          % 2 Im(eta)
D.r   = sqrt(D.v.^2 + D.w.^2);         % 2 |eta|
D.etaAbs = D.r / 2;
D.sa2 = D.A + D.r;                     % sigma_a^2   (Eq. 1)
D.sb2 = D.A - D.r;                     % sigma_b^2   (Eq. 2)
D.sigmaA = sqrt(max(D.sa2, 0));
D.sigmaB = sqrt(max(D.sb2, 0));
D.tau  = max(D.sa2, 1 ./ D.sb2);       % tau  = max(sigma_a^2, 1/sigma_b^2)
D.beta = D.sa2 + 1 ./ D.sb2;           % beta = sigma_a^2 + 1/sigma_b^2
D.K    = sqrt(D.sa2 ./ D.sb2);         % conformal distortion K = sigma_a / sigma_b
D.k    = (D.K - 1) ./ (D.K + 1);       % small dilatation k in [0,1)
D.GL   = 2*(D.A - 1).^2 + 2*D.r.^2;    % ||M - I||_F^2 = 2((A-1)^2 + 4|eta|^2)
end
