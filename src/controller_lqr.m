% controller_lqr.m — LQR stabilising controller (linearised around upright).
% Owner: Pessoa A
%
% Usage:
%   u = controller_lqr(xhat, p)
%
% Inputs
%   xhat — estimated state (4×1), from any estimator
%   p    — params struct from params.m (must contain field K_lqr after design)
%
% Output
%   u — control force on cart (N), saturated to p.u_max

function u = controller_lqr(xhat, p)

% TODO (Pessoa A): design K_lqr offline via lqr() and store in params.
% Suggested cost matrices (tune as needed):
%   Q_lqr = diag([10, 1, 100, 1]);
%   R_lqr = 0.1;
%   [K, ~, ~] = lqr(A_lin, B_lin, Q_lqr, R_lqr);

if ~isfield(p, 'K_lqr')
    error('params.K_lqr not set — run lqr design first.');
end

u_max = 20;   % N, adjust to match hardware/sim limits
u = -p.K_lqr * xhat;
u = max(-u_max, min(u_max, u));

end
