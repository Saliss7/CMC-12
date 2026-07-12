% compute_metrics.m — Compute comparison metrics from a simulation log.
% Owner: Matheus Felipe
%
% Usage:
%   metrics = compute_metrics(log, p)
%
% Inputs
%   log — struct returned by run_scenario
%   p   — params struct
%
% Output
%   metrics — struct with fields:
%     .rmse        (1×4)  RMSE per state [x, x_dot, theta, theta_dot]
%     .nees        (N×1)  Normalised Estimation Error Squared  (chi2, dof = 4)
%     .nis         (N×1)  Normalised Innovation Squared        (chi2, dof = 2)
%     .nees_mean   scalar mean NEES  (≈ 4 if the filter is consistent)
%     .nis_mean    scalar mean NIS   (≈ 2 if the filter is consistent)
%     .t_converge  scalar time (s) at which the estimation error enters and
%                  stays within tolerance for a rolling window; NaN if never
%     .cpu_time    scalar estimator wall-clock time (s), from run_scenario

function metrics = compute_metrics(log, p)

err    = log.x_true - log.xhat;      % estimation error (N×4)
Nsteps = size(log.x_true, 1);

% --- RMSE per state ---------------------------------------------------------
metrics.rmse = sqrt(mean(err.^2, 1));

% --- NEES: e' P^{-1} e  (chi-square with 4 dof under consistency) -----------
nees = zeros(Nsteps, 1);
for k = 1:Nsteps
    e_k     = err(k,:)';
    P_k     = log.P(:,:,k);
    nees(k) = e_k' * (P_k \ e_k);
end
metrics.nees      = nees;
metrics.nees_mean = mean(nees);

% --- NIS: nu' S^{-1} nu  (chi-square with 2 dof under consistency) ----------
% Uses the innovation (nu) and its covariance (S) logged by run_scenario for
% filters that expose them (EKF/UKF). Steps without an update (dropout) or from
% filters that do not provide the innovation stay NaN.
nis = nan(Nsteps, 1);
if isfield(log, 'innov') && isfield(log, 'S')
    for k = 1:Nsteps
        nu  = log.innov(k,:)';
        S_k = log.S(:,:,k);
        if all(isfinite(nu)) && all(isfinite(S_k(:)))
            nis(k) = nu' * (S_k \ nu);
        end
    end
end
metrics.nis      = nis;
metrics.nis_mean = mean(nis(isfinite(nis)));   % NaN if no valid steps

% --- Convergence time (rolling-window) --------------------------------------
% First instant after which every state error stays within tolerance for a
% window of p.conv_window seconds. Per-state tolerances and window come from
% params if present, otherwise sensible defaults are used.
if isfield(p, 'conv_tol')
    tol = p.conv_tol(:)';
else
    tol = [0.05, 0.10, deg2rad(2), 0.10];   % [x, x_dot, theta, theta_dot]
end
if isfield(p, 'conv_window')
    W = max(1, round(p.conv_window / p.dt));
else
    W = round(1.0 / p.dt);                  % 1 s window by default
end

within             = all(abs(err) <= tol, 2);   % N×1: all states within tol
metrics.t_converge = NaN;
for k = 1:(Nsteps - W + 1)
    if all(within(k:k+W-1))
        metrics.t_converge = log.t(k);
        break;
    end
end

% --- CPU time (measured inside run_scenario) --------------------------------
if isfield(log, 'cpu_time')
    metrics.cpu_time = log.cpu_time;
else
    metrics.cpu_time = NaN;
end

end
