% plant_cartpole.m — Nonlinear cart-pole dynamics (continuous-time ODE).
% Owner: Pessoa A
%
% Usage:
%   xdot = plant_cartpole(t, x, u, p)
%
% Inputs
%   t  — current time (s), unused but required by ode45 signature
%   x  — state [cart_pos; cart_vel; pole_angle; pole_ang_vel]  (4×1)
%   u  — horizontal force on cart (N)
%   p  — params struct from params.m
%
% Output
%   xdot — time derivative of x  (4×1)

function xdot = plant_cartpole(~, x, u, p)

% Unpack state
% x(1) = cart position    x
% x(2) = cart velocity    x_dot
% x(3) = pole angle       theta  (0 = upright)
% x(4) = pole ang. vel.   theta_dot

% TODO (Pessoa A): implement nonlinear EOM
% Reference equations of motion (Lagrangian formulation):
%
%   (M + m) x_ddot  +  b x_dot  -  m L theta_ddot cos(theta)
%       + m L theta_dot^2 sin(theta) = u
%
%   (I + m L^2) theta_ddot  -  m g L sin(theta)
%       = m L x_ddot cos(theta)
%
% where I = m L^2 / 3  for a uniform rod.

xdot = zeros(4,1);  % placeholder — replace with actual dynamics

end
