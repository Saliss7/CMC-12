% est_ekf.m — Extended Kalman Filter. 
% Owner: Pessoa C
%   [xhat, P, dbg] = est_ekf(xhat_prev, P_prev, z_k, u_k, p)
% z_k = [cart_pos; pole_angle], NaN = dropped sample. Optional dbg (.nu/.S/.z_pred) for NIS.

function [xhat, P, dbg] = est_ekf(xhat_prev, P_prev, z_k, u_k, p)

n = 4;
Q = getfield_default(p, 'Q_ekf', p.Q);
R = getfield_default(p, 'R_ekf', p.R);
C = [1 0 0 0;
     0 0 1 0];
dbg = [];

% Predict: nonlinear mean (RK4), covariance via Jacobian frozen at xhat_prev.
x_pred = rk4_step(@(t,x) plant_cartpole(t, x, u_k, p), 0, xhat_prev, p.dt);
Fd = expm(jacobian_f(xhat_prev, u_k, p) * p.dt);
P_pred = Fd * P_prev * Fd' + Q;
P_pred = (P_pred + P_pred') / 2;

if any(~isfinite(z_k))   % dropped measurement: coast on the prediction
    xhat = x_pred;
    P    = P_pred;
    return
end

S    = C * P_pred * C' + R;
K    = P_pred * C' / S;
nu   = z_k - C * x_pred;
xhat = x_pred + K * nu;

ImKC = eye(n) - K * C;   % Joseph form: stays symmetric PSD under mistuned R / large P0
P    = ImKC * P_pred * ImKC' + K * R * K';
P    = (P + P') / 2;

if nargout > 2
    dbg.nu     = nu;
    dbg.S      = S;
    dbg.z_pred = C * x_pred;
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
