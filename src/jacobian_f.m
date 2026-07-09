% jacobian_f.m — Analytical Jacobian of plant_cartpole w.r.t. state x.
% Owner: Pessoa C (used by est_ekf.m; standalone file so it can be unit-tested)
%
%   F = jacobian_f(x, u, p)    →  4x4 matrix (continuous-time df/dx)
%
% Mirrors the equations of motion documented in plant_cartpole.m's header
% (KEEP IN SYNC). Cart-pole, Lagrangian form, uniform rod I = m L^2 / 3,
% theta = 0 upright. Solving the 2x2 system for [x_ddot; theta_ddot]:
%
%   a = M+m,  d = I+mL^2,  k = mL
%   det   = a*d - k^2 cos^2(th)
%   rhs1  = u - b*xd - k*thd^2*sin(th)
%   rhs2  = g*k*sin(th)
%   x_ddot     = ( d*rhs1 + k*cos(th)*rhs2 ) / det
%   theta_ddot = ( a*rhs2 + k*cos(th)*rhs1 ) / det
%
% State x = [pos; xd; th; thd]. Dynamics are independent of pos, so column 1
% of F is zero (translation invariance).

function F = jacobian_f(x, u, p)

M = p.M;  m = p.m;  L = p.L;  g = p.g;  b = p.b;

xd  = x(2);
th  = x(3);
thd = x(4);

a = M + m;
d = m*L^2/3 + m*L^2;      % I + m L^2, with I = m L^2 / 3
k = m * L;

c = cos(th);
s = sin(th);

det  = a*d - k^2 * c^2;
rhs1 = u - b*xd - k*thd^2*s;
rhs2 = g*k*s;

N1 = d*rhs1 + k*c*rhs2;   % numerator of x_ddot
N2 = a*rhs2 + k*c*rhs1;   % numerator of theta_ddot

% Partials of the building blocks
ddet_dth   = 2*k^2 * c * s;         % d(det)/d(th)
drhs1_dxd  = -b;
drhs1_dth  = -k*thd^2 * c;
drhs1_dthd = -2*k*thd*s;
drhs2_dth  = g*k*c;

% d(N1)/d(.)
dN1_dxd  = d*drhs1_dxd;
dN1_dthd = d*drhs1_dthd;
dN1_dth  = d*drhs1_dth + k*(-s*rhs2 + c*drhs2_dth);

% d(N2)/d(.)
dN2_dxd  = k*c*drhs1_dxd;
dN2_dthd = k*c*drhs1_dthd;
dN2_dth  = a*drhs2_dth + k*(-s*rhs1 + c*drhs1_dth);

% Assemble (quotient rule for the theta column)
F = zeros(4);
F(1,:) = [0, 1, 0, 0];
F(3,:) = [0, 0, 0, 1];

F(2,2) = dN1_dxd / det;
F(2,4) = dN1_dthd / det;
F(2,3) = (dN1_dth*det - N1*ddet_dth) / det^2;

F(4,2) = dN2_dxd / det;
F(4,4) = dN2_dthd / det;
F(4,3) = (dN2_dth*det - N2*ddet_dth) / det^2;

end
