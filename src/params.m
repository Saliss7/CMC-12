% params.m — Physical and simulation parameters shared by all modules.
% ALL files must use this struct; do not hardcode constants elsewhere.
%
% Cart-pole convention
%   State:       x = [cart_pos; cart_vel; pole_angle; pole_ang_vel]
%   theta = 0   → pole pointing UPWARD (stabilisation target)
%   theta = pi  → pole pointing downward (natural rest)
%   Measurement: z = [cart_pos; pole_angle]
%   Input:       u  (horizontal force on cart, N)

function p = params()

% --- Physical ---
p.M   = 1.0;    % cart mass (kg)
p.m   = 0.3;    % pole mass (kg)
p.L   = 0.5;    % half-pole length to CoM (m)
p.g   = 9.81;   % gravity (m/s^2)
p.b   = 0.1;    % cart viscous friction (N·s/m)
p.I   = p.m*p.L^2/3;   % pole inertia about its CM (uniform rod) — single source of truth

% --- Simulation ---
p.dt  = 0.01;   % sample period (s)
p.Tf  = 10.0;   % default simulation horizon (s)

% --- Sensor noise (std dev) ---
p.sigma_x     = 0.005;   % cart position encoder (m)
p.sigma_theta = 0.005;   % pole angle encoder (rad)

% --- Estimator common priors ---
p.x0hat = zeros(4,1);           % initial state estimate
p.P0    = diag([0.1, 1, 0.1, 1]); % initial covariance

% --- Process / measurement noise covariances ---
p.Q = diag([1e-4, 1e-3, 1e-4, 1e-3]);   % process noise
p.R = diag([p.sigma_x^2, p.sigma_theta^2]); % measurement noise

% --- LQR design (linearização em torno da vertical: theta = 0, u = 0) ---
% Reaproveita o Jacobiano analítico (jacobian_f.m) para a matriz de estados A.
A = jacobian_f([0;0;0;0], 0, p);

% Matriz de entrada B = df/du (só as linhas de aceleração dependem de u)
d     = p.I + p.m*p.L^2;             % I + m L^2
k     = p.m*p.L;                     % m L
a     = p.M + p.m;                   % M + m
Delta = a*d - k^2;                   % determinante em theta = 0
B = [0; d/Delta; 0; k/Delta];

% Verificação de controlabilidade (par (A,B) deve ter posto completo)
if rank(ctrb(A, B)) < 4
    error('params:notControllable', ...
          'par (A,B) não controlável — verifique os parâmetros físicos.');
end

% Pesos LQR pela regra de Bryson: Q_ii = 1/desvio_max^2, R = 1/u_max^2
%   x_max = 0.5 m | xdot_max = 2 m/s | theta_max = 10 deg | thetadot_max = 2 rad/s
Q_lqr = diag([1/0.5^2, 1/2^2, 1/deg2rad(10)^2, 1/2^2]);
R_lqr = 1/20^2;

p.K_lqr = lqr(A, B, Q_lqr, R_lqr);   % ganho ótimo de realimentação (1x4)
p.u_max = 20;                        % limite de força do atuador (N)

% --- Swing-up energético ---
p.k_swing = 5;   % ganho de bombeamento de energia (sintonizar: ~1 a 20)

end
