% compute_metrics.m — Compute comparison metrics from a simulation log.
% Owner: Pessoa A
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
%     .nees        (N×1)  Normalised Estimation Error Squared
%     .nis         (N×1)  Normalised Innovation Squared
%     .t_converge  scalar time to first convergence (s), NaN if never
%     .cpu_time    scalar wall-clock time of estimator calls (s)

function metrics = compute_metrics(log, p)

err = log.x_true - log.xhat;

% RMSE per state
metrics.rmse = sqrt(mean(err.^2, 1));

% NEES: (e_k)' * P_k^{-1} * (e_k)  —  chi2(4) under consistency
Nsteps = size(log.x_true, 1);
nees   = zeros(Nsteps, 1);
for k = 1:Nsteps
    e_k      = err(k,:)';
    P_k      = log.P(:,:,k);
    nees(k)  = e_k' * (P_k \ e_k);
end
metrics.nees = nees;

% NIS: TODO (Pessoa A / B / C) — requires innovation sequence from estimator
% For now, fill with NaN as placeholder
metrics.nis = nan(Nsteps, 1);

% Convergence: first time all state RMSEs drop below threshold (rolling 1 s)
% TODO: implement rolling-window convergence detection
metrics.t_converge = NaN;

% CPU time is measured externally by run_scenario or Monte Carlo wrapper
metrics.cpu_time = NaN;

end
