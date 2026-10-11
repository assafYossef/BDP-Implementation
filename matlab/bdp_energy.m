function E = bdp_energy(e, data, variant)
%BDP_ENERGY Value of the (convex) objective of the selected variant at metric e.
%   Used for the dynamic penalty weight lambda_l of Section 7.2 and for logging.

D = bdp_distortion(e, data);
w = data.wF;
switch upper(variant)
    case 'V1'                                   % max_i tau_i
        E = max(D.tau);
    case 'V1R'                                  % max tau + avgWeight * 0.5 * sum w beta
        aw = 0.01; if isfield(data, 'avgWeight'), aw = data.avgWeight; end
        E = max(D.tau) + aw * 0.5 * sum(w .* D.beta);
    case 'V2'                                   % sum_i w_i beta_i
        E = sum(w .* D.beta);
    case 'V3'                                   % sum_i w_i ||M_i - I||_F^2
        E = sum(w .* D.GL);
    case 'V4'                                   % max_i |eta_i|
        E = max(D.etaAbs);
    case 'V5'                                   % sum_i w_i [(A-A0)^2 + 4|eta-eta0|^2]
        E = sum(w .* ((D.A - data.A0).^2 + (D.v - data.v0).^2 + (D.w - data.w0).^2));
    otherwise
        error('bdp:variant', 'Unknown variant %s', variant);
end
end
