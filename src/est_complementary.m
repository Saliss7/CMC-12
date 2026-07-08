% est_complementary.m — Complementary filter (high-pass model + low-pass meas.).
% Owner: Pessoa B
%
% Common interface:
%   [xhat, P] = est_complementary(xhat_prev, P_prev, z_k, u_k, p)

function [xhat, P] = est_complementary(xhat_prev, P_prev, z_k, u_k, p) %#ok<INUSD>

% Complementary gain (0 = trust model only, 1 = trust measurement only)
alpha = getfield_default(p, 'alpha_comp', 0.02);

% Model prediction: forward-Euler with zero acceleration assumption
x_pred     = xhat_prev(1) + p.dt * xhat_prev(2);
theta_pred = xhat_prev(3) + p.dt * xhat_prev(4);

% Fuse prediction (high-frequency) with measurement (low-frequency)
x_fused     = (1 - alpha) * x_pred     + alpha * z_k(1);
theta_fused = (1 - alpha) * theta_pred + alpha * z_k(2);

% Velocity from fusion difference
x_dot_est     = (x_fused     - xhat_prev(1)) / p.dt;
theta_dot_est = (theta_fused - xhat_prev(3)) / p.dt;

xhat = [x_fused; x_dot_est; theta_fused; theta_dot_est];
P    = P_prev;

end

% --- helper ------------------------------------------------------------------
function v = getfield_default(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end
