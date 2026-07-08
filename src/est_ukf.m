% est_ukf.m — Unscented Kalman Filter (sigma-point / unscented transform).
% Owner: Pessoa C
%
% Common interface:
%   [xhat, P] = est_ukf(xhat_prev, P_prev, z_k, u_k, p)

function [xhat, P] = est_ukf(xhat_prev, P_prev, z_k, u_k, p)

n  = 4;                   % state dimension
Q  = p.Q;
R  = p.R;
C  = [1 0 0 0; 0 0 1 0];  % observation matrix (linear)

% --- UKF tuning parameters ---
alpha = getfield_default(p, 'ukf_alpha', 1e-3);
beta  = getfield_default(p, 'ukf_beta',  2);
kappa = getfield_default(p, 'ukf_kappa', 0);
lambda = alpha^2 * (n + kappa) - n;

% Weights
Wm = [lambda/(n+lambda), repmat(1/(2*(n+lambda)), 1, 2*n)];
Wc = Wm;
Wc(1) = Wc(1) + (1 - alpha^2 + beta);

% --- Generate sigma points ---
try
    S = chol((n + lambda) * P_prev, 'lower');
catch
    % Regularise if not positive definite
    P_prev = P_prev + 1e-6 * eye(n);
    S = chol((n + lambda) * P_prev, 'lower');
end

sigmas = [xhat_prev, xhat_prev(:,ones(1,n)) + S, ...
                     xhat_prev(:,ones(1,n)) - S];   % 4 × (2n+1)

% --- Propagate sigma points through nonlinear dynamics ---
sigmas_pred = zeros(n, 2*n+1);
for i = 1:2*n+1
    sigmas_pred(:,i) = rk4_step(@(t,x) plant_cartpole(t, x, u_k, p), ...
                                 0, sigmas(:,i), p.dt);
end

% --- Predicted mean and covariance ---
x_pred = sigmas_pred * Wm';
P_pred = Q;
for i = 1:2*n+1
    d = sigmas_pred(:,i) - x_pred;
    P_pred = P_pred + Wc(i) * (d * d');
end

% --- Measurement update (linear C, so standard KF suffices) ---
z_pred = C * x_pred;
S_inn  = C * P_pred * C' + R;
K      = P_pred * C' / S_inn;
xhat   = x_pred + K * (z_k - z_pred);
P      = (eye(n) - K * C) * P_pred;

end


% --- helpers -----------------------------------------------------------------
function v = getfield_default(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
