% run_scenario.m — Main simulation harness.
% Owner: Matheus Felipe
%
% Runs one scenario with one estimator in closed loop and returns a log struct.
%
% Usage:
%   log = run_scenario(estimator_fn, scenario, p)
%
% Inputs
%   estimator_fn — function handle with the common interface:
%                    [xhat, P]        = estimator_fn(xhat_prev, P_prev, z_k, u_k, p)
%                  or, optionally (EKF/UKF), a 3rd debug output:
%                    [xhat, P, dbg]   = estimator_fn(...)
%                  where dbg.nu is the innovation and dbg.S its covariance.
%   scenario     — struct with scenario overrides (see scenarios/S*.m)
%   p            — params struct from params.m
%
% Output
%   log — struct with fields:
%     .t        (N×1)   time vector
%     .x_true   (N×4)   true state (time-aligned with xhat/z)
%     .xhat     (N×4)   estimated state
%     .P        (4×4×N) covariance history
%     .u        (N×1)   control input applied over [t_k, t_{k+1}]
%     .z        (N×2)   measurements (NaN row when the measurement was dropped)
%     .dropped  (N×1)   logical: true when the measurement was dropped (S4)
%     .innov    (N×2)   innovation nu (NaN if estimator does not provide it)
%     .S        (2×2×N) innovation covariance (NaN if not provided)
%     .mode     (N×1)   cell array: 'swingup' | 'lqr'
%     .cpu_time scalar  accumulated estimator wall-clock time (s)

function log = run_scenario(estimator_fn, scenario, p)

% --- Apply scenario overrides onto an effective params struct ----------------
% Scenarios may override the simulation horizon, the estimator priors and the
% noise covariances SEEN BY THE ESTIMATOR (see scenarios/S*.m). Physical
% parameters and the designed LQR gain are never overridden.
pe = p;
if isfield(scenario, 'Tf'),    pe.Tf    = scenario.Tf;    end
if isfield(scenario, 'x0hat'), pe.x0hat = scenario.x0hat; end
if isfield(scenario, 'P0'),    pe.P0    = scenario.P0;    end
if isfield(scenario, 'Q'),     pe.Q     = scenario.Q;     end
if isfield(scenario, 'R'),     pe.R     = scenario.R;     end

Nsteps = round(pe.Tf / pe.dt);

% --- Pre-allocate log arrays -------------------------------------------------
log.t      = zeros(Nsteps, 1);
log.x_true = zeros(Nsteps, 4);
log.xhat   = zeros(Nsteps, 4);
log.P      = zeros(4, 4, Nsteps);
log.u      = zeros(Nsteps, 1);
log.z       = zeros(Nsteps, 2);
log.dropped = false(Nsteps, 1);
log.innov   = nan(Nsteps, 2);
log.S       = nan(2, 2, Nsteps);
log.mode    = cell(Nsteps, 1);

% Does the estimator also return a debug struct (3rd output)? Non-breaking:
% 2-output estimators (KF/lowpass/complementary) still work. EKF/UKF expose
% dbg.nu (innovation) and dbg.S, enabling the NIS test in compute_metrics.
provides_dbg = abs(nargout(estimator_fn)) >= 3;

% --- Initial conditions ------------------------------------------------------
x_true   = scenario.x0_true;   % true initial state (4×1)
xhat     = pe.x0hat;           % estimator prior mean
P        = pe.P0;              % estimator prior covariance
u_k_prev = 0;                  % no control applied before the first step

cpu_time = 0;                  % accumulated estimator wall-clock time

for k = 1:Nsteps

    t_k = (k-1) * pe.dt;

    % 1 — Measure the true state at t_k (z_k = [] when the reading is dropped)
    z_k     = sensor_model(x_true, pe, scenario);
    dropped = isempty(z_k);

    % 2 — Estimate: predict with the previous control u_{k-1}, update with z_k.
    %     On a dropout (z_k == []) the estimator must do a prediction-only step.
    t0 = tic;
    if provides_dbg
        [xhat, P, dbg] = estimator_fn(xhat, P, z_k, u_k_prev, pe);
    else
        [xhat, P]      = estimator_fn(xhat, P, z_k, u_k_prev, pe);
        dbg = [];
    end
    cpu_time = cpu_time + toc(t0);

    % 3 — Control from the current estimate (swing-up → LQR)
    [u_k, mode_k] = swingup(xhat, pe);

    % 4 — Log at t_k, BEFORE propagating, so x_true, xhat and z are time-aligned
    log.t(k)        = t_k;
    log.x_true(k,:) = x_true';
    log.xhat(k,:)   = xhat';
    log.P(:,:,k)    = P;
    log.u(k)        = u_k;
    log.dropped(k)  = dropped;
    if dropped, log.z(k,:) = [NaN, NaN]; else, log.z(k,:) = z_k'; end
    log.mode{k}     = mode_k;
    if isstruct(dbg)
        if isfield(dbg, 'nu'), log.innov(k,:) = dbg.nu(:)'; end
        if isfield(dbg, 'S'),  log.S(:,:,k)   = dbg.S;      end
    end

    % 5 — Propagate the true plant from t_k to t_{k+1} (RK4) with u_k
    x_true   = rk4_step(@(t,x) plant_cartpole(t, x, u_k, pe), t_k, x_true, pe.dt);
    u_k_prev = u_k;
end

log.cpu_time = cpu_time;

end

% --- RK4 integrator step -----------------------------------------------------
function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
