% est_ukf.m — Unscented Kalman Filter (scaled sigma points). 
% Owner: Pessoa C
%   [xhat, P, dbg] = est_ukf(xhat_prev, P_prev, z_k, u_k, p)
% z_k = [cart_pos; pole_angle], NaN = dropped sample. Optional dbg (.nu/.S/.z_pred) for NIS.

function [xhat, P, dbg] = est_ukf(xhat_prev, P_prev, z_k, u_k, p)

n  = 4;
Q  = getfield_default(p, 'Q_ukf', p.Q);
R  = getfield_default(p, 'R_ukf', p.R);
C  = [1 0 0 0; 0 0 1 0];
dbg = [];

alpha = getfield_default(p, 'ukf_alpha', 1e-3);
beta  = getfield_default(p, 'ukf_beta',  2);
kappa = getfield_default(p, 'ukf_kappa', 0);
lambda = alpha^2 * (n + kappa) - n;

Wm = [lambda/(n+lambda), repmat(1/(2*(n+lambda)), 1, 2*n)];
Wc = Wm;
Wc(1) = Wc(1) + (1 - alpha^2 + beta);

try
    Ssp = chol((n + lambda) * P_prev, 'lower');
catch
    P_prev = P_prev + 1e-6 * eye(n);   % nudge back to positive definite
    Ssp = chol((n + lambda) * P_prev, 'lower');
end

sigmas = [xhat_prev, xhat_prev(:,ones(1,n)) + Ssp, ...
                     xhat_prev(:,ones(1,n)) - Ssp];

sigmas_pred = zeros(n, 2*n+1);
for i = 1:2*n+1
    sigmas_pred(:,i) = rk4_step(@(t,x) plant_cartpole(t, x, u_k, p), ...
                                 0, sigmas(:,i), p.dt);
end

x_pred = sigmas_pred * Wm';
P_pred = Q;
for i = 1:2*n+1
    dsig   = sigmas_pred(:,i) - x_pred;
    P_pred = P_pred + Wc(i) * (dsig * dsig');
end
P_pred = (P_pred + P_pred') / 2;

if any(~isfinite(z_k))   % dropped measurement: coast on the prediction
    xhat = x_pred;
    P    = P_pred;
    return
end

% h(x) = [x1; x3] is linear, so the update is the plain KF update with C.
z_pred = C * x_pred;
S_inn  = C * P_pred * C' + R;
K      = P_pred * C' / S_inn;
nu     = z_k - z_pred;
xhat   = x_pred + K * nu;

P    = (eye(n) - K * C) * P_pred;
P    = (P + P') / 2;

if nargout > 2
    dbg.nu     = nu;
    dbg.S      = S_inn;
    dbg.z_pred = z_pred;
end

end


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
