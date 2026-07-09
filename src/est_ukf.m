% est_ukf.m — Unscented Kalman Filter (sigma-point / unscented transform).
% Owner: Pessoa C
%
% Common interface:
%   [xhat, P]      = est_ukf(xhat_prev, P_prev, z_k, u_k, p)
%   [xhat, P, dbg] = est_ukf(...)   % optional debug/consistency output
%
% Inputs
%   xhat_prev (4x1) — previous state estimate [x; x_dot; theta; theta_dot]
%   P_prev    (4x4) — previous covariance
%   z_k       (2x1) — measurement [cart_pos; pole_angle]; NaN entries → dropped
%   u_k       (1x1) — control input applied over the step
%   p               — params struct (uses p.dt, p.Q/p.R or p.Q_ukf/p.R_ukf,
%                     and optional p.ukf_alpha/p.ukf_beta/p.ukf_kappa)
%
% Outputs
%   xhat (4x1) — updated state estimate
%   P    (4x4) — updated covariance (symmetric)
%   dbg        — struct with innovation info for NIS analysis:
%                  .nu (2x1) innovation, .S (2x2) innovation cov, .z_pred (2x1)
%                (empty when the measurement is dropped / update skipped)

function [xhat, P, dbg] = est_ukf(xhat_prev, P_prev, z_k, u_k, p)

n  = 4;                                   % state dimension
Q  = getfield_default(p, 'Q_ukf', p.Q);   % per-filter override, else shared Q
R  = getfield_default(p, 'R_ukf', p.R);
C  = [1 0 0 0; 0 0 1 0];                  % observation matrix (linear)

dbg = [];   % populated below only when an update is performed

% --- UKF tuning parameters (Van der Merwe scaled unscented transform) ---
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
    Ssp = chol((n + lambda) * P_prev, 'lower');
catch
    % Regularise if not positive definite
    P_prev = P_prev + 1e-6 * eye(n);
    Ssp = chol((n + lambda) * P_prev, 'lower');
end

sigmas = [xhat_prev, xhat_prev(:,ones(1,n)) + Ssp, ...
                     xhat_prev(:,ones(1,n)) - Ssp];   % 4 x (2n+1)

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
    dsig   = sigmas_pred(:,i) - x_pred;
    P_pred = P_pred + Wc(i) * (dsig * dsig');
end
P_pred = (P_pred + P_pred') / 2;      % keep symmetric

% --- Measurement update ---
% Skip the correction when the measurement is missing (dropout / low rate):
% coast on the prediction. See sensor NaN convention (S4).
if any(~isfinite(z_k))
    xhat = x_pred;
    P    = P_pred;
    return
end

% The measurement model h(x) = [x1; x3] is LINEAR, so propagating the sigma
% points through h and re-computing the innovation covariance would give
% exactly C*P_pred*C' + R. Using C directly is therefore exact, not an
% approximation — and cheaper.
z_pred = C * x_pred;
S_inn  = C * P_pred * C' + R;
K      = P_pred * C' / S_inn;
nu     = z_k - z_pred;                % innovation
xhat   = x_pred + K * nu;

P    = (eye(n) - K * C) * P_pred;
P    = (P + P') / 2;                  % symmetrise safety net

if nargout > 2
    dbg.nu     = nu;
    dbg.S      = S_inn;
    dbg.z_pred = z_pred;
end

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
