% plant_cartpole.m — Nonlinear cart-pole dynamics (continuous-time ODE).
% Owner: Matheus Felipede
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
%
% Expected fields in p:
%   p.M  — cart mass (kg)
%   p.m  — pole mass (kg)
%   p.L  — distance from pivot to pole center of mass (m)
%   p.b  — viscous friction coefficient on the cart (N/(m/s))
%   p.g  — gravitational acceleration (m/s^2)
%   p.I  — pole inertia about its CM (set in params.m; uniform rod = m*L^2/3)
function xdot = plant_cartpole(~, x, u, p)
% Unpack state
% x(1) = cart position    x
% x(2) = cart velocity    x_dot
% x(3) = pole angle       theta  (0 = upright)
% x(4) = pole ang. vel.   theta_dot

x_dot     = x(2);
theta     = x(3);
theta_dot = x(4);

M = p.M;
m = p.m;
L = p.L;
b = p.b;
g = p.g;
I = p.I;

a = M + m;
h = m * L * cos(theta);
J = I + m * L^2;

r1 = u - b*x_dot - m*L*theta_dot^2*sin(theta);
r2 = m*g*L*sin(theta);

Delta = a*J - h^2;

x_ddot     = ( J*r1 + h*r2 ) / Delta;
theta_ddot = ( h*r1 + a*r2 ) / Delta;

xdot = [ x_dot;
         x_ddot;
         theta_dot;
         theta_ddot ];
end