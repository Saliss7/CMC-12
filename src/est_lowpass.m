% est_lowpass.m — Baseline estimator: numerical differentiation + low-pass filter.
% Owner: Pessoa B
%
% Common interface:
%   [xhat, P] = est_lowpass(xhat_prev, P_prev, z_k, u_k, p)
%
% Inputs
%   xhat_prev (4×1) — previous state estimate
%   P_prev    (4×4) — previous covariance (not used here; passed for API compat)
%   z_k       (2×1) — measurement [cart_pos; pole_angle]
%   u_k       (1×1) — control input (not used; passed for API compat)
%   p         — params struct (uses p.dt and p.alpha_lp)
%
% Outputs
%   xhat (4×1) — updated state estimate
%   P    (4×4) — covariance (identity × large value — this estimator has no model)

function [xhat, P] = est_lowpass(xhat_prev, P_prev, z_k, u_k, p) %#ok<INUSD>

% Dropped measurement (S4, NaN sentinel — see sensor_model.m): this filter has
% no dynamic model to predict with, so the only sensible fallback is to hold
% the last estimate (freeze) until a real measurement arrives again.
if any(~isfinite(z_k))
    xhat = xhat_prev;
    P    = P_prev;
    return
end

% Low-pass filter coefficient (0 < alpha < 1); closer to 1 = smoother, more lag
alpha = getfield_default(p, 'alpha_lp', 0.8);

% Positions: filter measurement directly
x_filt     = alpha * xhat_prev(1) + (1 - alpha) * z_k(1);
theta_filt = alpha * xhat_prev(3) + (1 - alpha) * z_k(2);

% Velocities: numerical differentiation of filtered positions
x_dot_est     = (x_filt     - xhat_prev(1)) / p.dt;
theta_dot_est = (theta_filt - xhat_prev(3)) / p.dt;

xhat = [x_filt; x_dot_est; theta_filt; theta_dot_est];

% Covariance placeholder (no probabilistic model)
P = P_prev;

end

% --- helper ------------------------------------------------------------------
function v = getfield_default(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end
