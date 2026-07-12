% linearise_upright.m — Continuous linearisation about the upright equilibrium,
% discretised (ZOH) into p.Ad, p.Bd, p.Cd for est_kf.m.
% Owner: Pessoa B
%
% Usage:
%   p = linearise_upright(p)
%
% Standalone file (like jacobian_f.m is for Pessoa C) so params.m can call it
% and so it can be unit-tested on its own (tests/test_classical_filters.m).
%
% Continuous-time A_c, B_c are the analytical linearisation of the cart-pole
% EOM about the upright equilibrium (theta=0, u=0), derived in the report
% (Seção "Linearização em torno da vertical"). They are re-derived here (not
% imported from jacobian_f.m, which is Pessoa C's file for the EKF) to keep
% this file self-contained; jacobian_f(zeros(4,1),0,p) must match A_c below —
% this is exactly what tests/test_classical_filters.m checks.

function p = linearise_upright(p)

M = p.M;  m = p.m;  L = p.L;  g = p.g;  b = p.b;
I = m * L^2 / 3;     % moment of inertia (uniform rod)
D = (M + m) * (I + m*L^2) - (m*L)^2;   % denominator (Delta_0 in the report)

% Continuous-time A and B at theta=0 (upright), from the small-angle
% linearisation of the EOM (cos(theta)~1, sin(theta)~theta, thetadot^2*sin~0):
A_c = [0                    1        0                  0;
       0   -(I + m*L^2)*b/D     (m^2*L^2*g)/D            0;
       0                    0        0                  1;
       0        -(m*L*b)/D     (M+m)*m*g*L/D             0];

B_c = [0; (I + m*L^2)/D; 0; (m*L)/D];

C_c = [1 0 0 0;
       0 0 1 0];   % observe x and theta

% Discretise (zero-order hold)
sys   = ss(A_c, B_c, C_c, 0);
sysd  = c2d(sys, p.dt, 'zoh');
p.Ad  = sysd.A;
p.Bd  = sysd.B;
p.Cd  = sysd.C;

end
