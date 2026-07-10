% controller_lqr.m — LQR stabilising controller (linearised around upright).
% Owner: Matheus Felipe
%
% Applies the state-feedback control law u = -K*xhat that regulates the
% cart-pole to the upright equilibrium (X = 0). The optimal gain K is designed
% offline (linearisation around theta = 0 + lqr) and stored in params as
% p.K_lqr — see params.m. The linearisation, controllability check and choice
% of Q/R are documented in the report (main.tex).
%
% Usage:
%   u = controller_lqr(xhat, p)
%
% Inputs
%   xhat — estimated state (4x1) [cart_pos; cart_vel; pole_angle; pole_ang_vel]
%   p    — params struct from params.m (must contain field K_lqr)
%
% Output
%   u — control force on cart (N), saturated to +/- p.u_max

function u = controller_lqr(xhat, p)

% Wrap the pole angle to [-pi, pi] so the regulation error is the true angular
% deviation from upright (avoids a huge spurious -K*theta after full rotations)
xw    = xhat;
xw(3) = mod(xhat(3) + pi, 2*pi) - pi;

u = -p.K_lqr * xw;
u_max = p.u_max;
u = max(-u_max, min(u_max, u));

end
