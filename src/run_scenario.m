% run_scenario.m — Main simulation harness.
% Owner: Pessoa A
%
% Runs one scenario with one estimator and returns a log struct.
%
% Usage:
%   log = run_scenario(estimator_fn, scenario, p)
%
% Inputs
%   estimator_fn — function handle with the common interface:
%                    [xhat, P] = estimator_fn(xhat_prev, P_prev, z_k, u_k, p)
%   scenario     — struct with scenario overrides (see scenarios/S*.m)
%   p            — params struct from params.m
%
% Output
%   log — struct with fields:
%     .t       (N×1)   time vector
%     .x_true  (N×4)   true state
%     .xhat    (N×4)   estimated state
%     .P       (4×4×N) covariance history
%     .u       (N×1)   control input
%     .z       (N×2)   measurements
%     .mode    (N×1)   cell array: 'swingup' | 'lqr'

function log = run_scenario(estimator_fn, scenario, p)

Nsteps = round(p.Tf / p.dt);

% Pre-allocate log arrays
log.t      = zeros(Nsteps, 1);
log.x_true = zeros(Nsteps, 4);
log.xhat   = zeros(Nsteps, 4);
log.P      = zeros(4, 4, Nsteps);
log.u      = zeros(Nsteps, 1);
log.z      = zeros(Nsteps, 2);
log.mode   = cell(Nsteps, 1);

% Initial conditions
x_true = scenario.x0_true;     % true initial state (4×1)
xhat   = p.x0hat;
P      = p.P0;

for k = 1:Nsteps

    t_k = (k-1) * p.dt;

    % 1 — Measure
    z_k = sensor_model(x_true, p, scenario);

    % 2 — Estimate
    [xhat, P] = estimator_fn(xhat, P, z_k, u_k_prev, p);   %#ok u_k_prev set below on k>1

    % 3 — Control
    [u_k, mode_k] = swingup(xhat, p);

    % 4 — Propagate true plant (RK4)
    x_true = rk4_step(@(t,x) plant_cartpole(t, x, u_k, p), t_k, x_true, p.dt);

    % 5 — Log
    log.t(k)        = t_k;
    log.x_true(k,:) = x_true';
    log.xhat(k,:)   = xhat';
    log.P(:,:,k)    = P;
    log.u(k)        = u_k;
    log.z(k,:)      = z_k';
    log.mode{k}     = mode_k;

    u_k_prev = u_k;
end

end

% --- RK4 integrator step -----------------------------------------------------
function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
