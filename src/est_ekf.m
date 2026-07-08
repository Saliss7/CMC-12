% est_ekf.m — Extended Kalman Filter with online Jacobian computation.
% Owner: Pessoa C
%
% Common interface:
%   [xhat, P] = est_ekf(xhat_prev, P_prev, z_k, u_k, p)

function [xhat, P] = est_ekf(xhat_prev, P_prev, z_k, u_k, p)

Q = p.Q;
R = p.R;
C = [1 0 0 0;
     0 0 1 0];   % measurement Jacobian (linear, no approximation needed)

% --- Predict: RK4 integration of nonlinear dynamics ---
x_pred = rk4_step(@(t,x) plant_cartpole(t, x, u_k, p), 0, xhat_prev, p.dt);

% --- Linearise: compute A = df/dx at xhat_prev ---
F = jacobian_f(xhat_prev, u_k, p);   % 4×4 continuous Jacobian
% Discretise via 1st-order approximation (use expm for better accuracy)
Fd = eye(4) + F * p.dt;

P_pred = Fd * P_prev * Fd' + Q;

% --- Update ---
S    = C * P_pred * C' + R;
K    = P_pred * C' / S;
xhat = x_pred + K * (z_k - C * x_pred);
P    = (eye(4) - K * C) * P_pred;

end


% ---------------------------------------------------------------------------
% jacobian_f — Analytical Jacobian of plant_cartpole w.r.t. state x
%
%   F = jacobian_f(x, u, p)    →  4×4 matrix
%
% TODO (Pessoa C): derive all partial derivatives symbolically.
% ---------------------------------------------------------------------------
function F = jacobian_f(x, u, p) %#ok<INUSD>

% Placeholder: identity — replace with analytical Jacobian
F = zeros(4);

end


% --- RK4 (duplicated here so est_ekf.m is self-contained) -------------------
function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
