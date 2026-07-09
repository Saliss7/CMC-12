% est_ekf.m — Extended Kalman Filter with online Jacobian computation.
% Owner: Pessoa C
%
% Common interface:
%   [xhat, P]      = est_ekf(xhat_prev, P_prev, z_k, u_k, p)
%   [xhat, P, dbg] = est_ekf(...)   % optional debug/consistency output
%
% Inputs
%   xhat_prev (4x1) — previous state estimate [x; x_dot; theta; theta_dot]
%   P_prev    (4x4) — previous covariance
%   z_k       (2x1) — measurement [cart_pos; pole_angle]; NaN entries → dropped
%   u_k       (1x1) — control input applied over the step
%   p               — params struct (uses p.dt, p.Q/p.R or p.Q_ekf/p.R_ekf)
%
% Outputs
%   xhat (4x1) — updated state estimate
%   P    (4x4) — updated covariance (symmetric)
%   dbg        — struct with innovation info for NIS analysis:
%                  .nu (2x1) innovation, .S (2x2) innovation cov, .z_pred (2x1)
%                (empty when the measurement is dropped / update skipped)
%
% Notes
%   * Predict step integrates the full nonlinear dynamics (RK4); the covariance
%     is propagated with the analytical Jacobian evaluated online at xhat_prev.
%   * The measurement Jacobian C is exact (measurements are linear in the
%     state), so no linearisation is needed there.

function [xhat, P, dbg] = est_ekf(xhat_prev, P_prev, z_k, u_k, p)

n = 4;
Q = getfield_default(p, 'Q_ekf', p.Q);   % per-filter override, else shared Q
R = getfield_default(p, 'R_ekf', p.R);
C = [1 0 0 0;
     0 0 1 0];   % measurement Jacobian (linear, no approximation needed)

dbg = [];   % populated below only when an update is performed

% --- Predict: RK4 integration of nonlinear dynamics ---
x_pred = rk4_step(@(t,x) plant_cartpole(t, x, u_k, p), 0, xhat_prev, p.dt);

% --- Linearise: continuous Jacobian A = df/dx at xhat_prev ---
F  = jacobian_f(xhat_prev, u_k, p);   % 4x4 continuous Jacobian
% Discretise via matrix exponential of the frozen linearisation
Fd = expm(F * p.dt);

P_pred = Fd * P_prev * Fd' + Q;
P_pred = (P_pred + P_pred') / 2;      % keep symmetric

% --- Update ---
% Skip the correction when the measurement is missing (dropout / low rate):
% coast on the prediction. See sensor NaN convention (S4).
if any(~isfinite(z_k))
    xhat = x_pred;
    P    = P_pred;
    return
end

S    = C * P_pred * C' + R;
K    = P_pred * C' / S;
nu   = z_k - C * x_pred;              % innovation
xhat = x_pred + K * nu;

% Joseph form: numerically stable and stays PSD even with mistuned R / large P0
ImKC = eye(n) - K * C;
P    = ImKC * P_pred * ImKC' + K * R * K';
P    = (P + P') / 2;                  % symmetrise safety net

if nargout > 2
    dbg.nu     = nu;
    dbg.S      = S;
    dbg.z_pred = C * x_pred;
end

end


% --- helpers -----------------------------------------------------------------
function v = getfield_default(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

% --- RK4 (duplicated here so est_ekf.m is self-contained) -------------------
function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
