% swingup.m — Energy-based swing-up controller.
% Owner: Pessoa A
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

% Switch to LQR when pole is within ~15 deg of upright and slow
near_upright = (abs(theta) < deg2rad(15)) && (abs(theta_dot) < 2.0);

if near_upright
    u    = controller_lqr(xhat, p);
    mode = 'lqr';
else
    % TODO (Pessoa A): implement energy-based swing-up
    % Classic Åström–Furuta method:
    %   E_ref = 2 m g L   (energy of upright equilibrium)
    %   E     = 0.5 (I + m L^2) theta_dot^2 - m g L cos(theta)
    %   u     = k_e * (E - E_ref) * sign(theta_dot * cos(theta))
    u    = 0;   % placeholder
    mode = 'swingup';
end

end
