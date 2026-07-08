% est_kf.m — Linear Kalman Filter (linearised around upright equilibrium).
% Owner: Pessoa B
%
% Common interface:
%   [xhat, P] = est_kf(xhat_prev, P_prev, z_k, u_k, p)
%
% Note: p must contain the discrete-time system matrices p.Ad, p.Bd, p.Cd
%       computed offline via linearise_upright() below (or stored in params).

function [xhat, P] = est_kf(xhat_prev, P_prev, z_k, u_k, p)

if ~isfield(p, 'Ad')
    error('params missing Ad/Bd/Cd — run linearise_upright() first.');
end

Ad = p.Ad;  Bd = p.Bd;  Cd = p.Cd;
Q  = p.Q;   R  = p.R;

% --- Predict ---
x_pred = Ad * xhat_prev + Bd * u_k;
P_pred = Ad * P_prev * Ad' + Q;

% --- Update ---
S   = Cd * P_pred * Cd' + R;
K   = P_pred * Cd' / S;
xhat = x_pred + K * (z_k - Cd * x_pred);
P    = (eye(4) - K * Cd) * P_pred;

end


% ---------------------------------------------------------------------------
% linearise_upright  — call once to populate p.Ad, p.Bd, p.Cd
%
%   p = linearise_upright(p)
%
% TODO (Pessoa B): derive A_c, B_c analytically from EOM, then discretise.
% ---------------------------------------------------------------------------
function p = linearise_upright(p)

M = p.M;  m = p.m;  L = p.L;  g = p.g;  b = p.b;
I = m * L^2 / 3;     % moment of inertia (uniform rod)
D = (M + m) * (I + m*L^2) - (m*L)^2;   % denominator

% Continuous-time A and B at theta=0 (upright)
% TODO (Pessoa B): fill in A_c, B_c from the linearised EOM
A_c = zeros(4);   % placeholder
B_c = zeros(4,1); % placeholder

C_c = [1 0 0 0;
       0 0 1 0];   % observe x and theta

% Discretise (zero-order hold)
sys   = ss(A_c, B_c, C_c, 0);
sysd  = c2d(sys, p.dt, 'zoh');
p.Ad  = sysd.A;
p.Bd  = sysd.B;
p.Cd  = sysd.C;

end
