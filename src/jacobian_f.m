% jacobian_f.m — analytical df/dx (4x4) of the cart-pole. 
% Owner: Pessoa C
%   F = jacobian_f(x, u, p),  x = [pos; xd; th; thd]
% Derivatives of plant_cartpole.m's EOM (KEEP IN SYNC). Independent of pos, so column 1 is 0.

function F = jacobian_f(x, u, p)

M = p.M;  m = p.m;  L = p.L;  g = p.g;  b = p.b;

xd  = x(2);
th  = x(3);
thd = x(4);

a = M + m;
d = m*L^2/3 + m*L^2;      % I + m L^2
k = m * L;

c = cos(th);
s = sin(th);

det  = a*d - k^2 * c^2;
rhs1 = u - b*xd - k*thd^2*s;
rhs2 = g*k*s;

N1 = d*rhs1 + k*c*rhs2;
N2 = a*rhs2 + k*c*rhs1;

ddet_dth   = 2*k^2 * c * s;
drhs1_dxd  = -b;
drhs1_dth  = -k*thd^2 * c;
drhs1_dthd = -2*k*thd*s;
drhs2_dth  = g*k*c;

dN1_dxd  = d*drhs1_dxd;
dN1_dthd = d*drhs1_dthd;
dN1_dth  = d*drhs1_dth + k*(-s*rhs2 + c*drhs2_dth);

dN2_dxd  = k*c*drhs1_dxd;
dN2_dthd = k*c*drhs1_dthd;
dN2_dth  = a*drhs2_dth + k*(-s*rhs1 + c*drhs1_dth);

F = zeros(4);
F(1,:) = [0, 1, 0, 0];
F(3,:) = [0, 0, 0, 1];

F(2,2) = dN1_dxd / det;
F(2,4) = dN1_dthd / det;
F(2,3) = (dN1_dth*det - N1*ddet_dth) / det^2;   % quotient rule (det depends on th)

F(4,2) = dN2_dxd / det;
F(4,4) = dN2_dthd / det;
F(4,3) = (dN2_dth*det - N2*ddet_dth) / det^2;

end
