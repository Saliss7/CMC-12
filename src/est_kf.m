% est_kf.m — Linear Kalman Filter (linearised around upright equilibrium).
% Owner: Pessoa B
%
% Common interface:
%   [xhat, P] = est_kf(xhat_prev, P_prev, z_k, u_k, p)
%
% Note: p must contain the discrete-time system matrices p.Ad, p.Bd, p.Cd,
%       computed via linearise_upright.m (called once inside params.m, so
%       any p returned by params() already has them).

function [xhat, P] = est_kf(xhat_prev, P_prev, z_k, u_k, p)

n = 4;

if ~isfield(p, 'Ad')
    error('params missing Ad/Bd/Cd — run linearise_upright() first.');
end

Ad = p.Ad;  Bd = p.Bd;  Cd = p.Cd;
Q  = p.Q;   R  = p.R;

% --- Predict ---
x_pred = Ad * xhat_prev + Bd * u_k;
P_pred = Ad * P_prev * Ad' + Q;
P_pred = (P_pred + P_pred') / 2;      % keep symmetric

% Dropped measurement (S4, NaN sentinel — see sensor_model.m): coast on the
% prediction, same convention as est_ekf.m/est_ukf.m.
if any(~isfinite(z_k))
    xhat = x_pred;
    P    = P_pred;
    return
end

% --- Update ---
S    = Cd * P_pred * Cd' + R;
K    = P_pred * Cd' / S;
xhat = x_pred + K * (z_k - Cd * x_pred);

% Joseph form: numerically stable and stays PSD even with mistuned R / large P0
% (same convention as est_ekf.m).
ImKC = eye(n) - K * Cd;
P    = ImKC * P_pred * ImKC' + K * R * K';
P    = (P + P') / 2;                  % symmetrise safety net

end

% linearise_upright() used to be a local (private) function here, but params.m
% needs to call it too, and MATLAB local functions aren't visible outside the
% file they're defined in — moved to its own file: linearise_upright.m.
