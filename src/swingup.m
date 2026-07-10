% swingup.m — Energy-based swing-up controller.
% Owner: Matheus Felipe
%
% Switches to LQR (controller_lqr) when the pole is near upright.
%
% Usage:
%   [u, mode] = swingup(xhat, p)
%
% Inputs
%   xhat — estimated state (4×1)
%   p    — params struct from params.m
%
% Outputs
%   u    — control force (N)
%   mode — 'swingup' | 'lqr'

function [u, mode] = swingup(xhat, p)

theta     = xhat(3);
theta_dot = xhat(4);

near_upright = (abs(theta) < deg2rad(15)) && (abs(theta_dot) < 2.0);

if near_upright
    u    = controller_lqr(xhat, p);
    mode = 'lqr';
else
    I     = p.I + p.m*p.L^2;
    E     = 0.5*I*theta_dot^2 + p.m*p.g*p.L*cos(theta);
    Ed    = p.m*p.g*p.L;
    Etil  = E - Ed;

    u    = -p.k_swing * theta_dot * Etil * cos(theta);
    u    = max(-p.u_max, min(p.u_max, u));
    mode = 'swingup';
end

end
