% test_classical_filters.m — Unit/smoke tests for est_lowpass.m,
% est_complementary.m, est_kf.m and linearise_upright.m.
% Owner: Pessoa B
%
% Run from the repo root:
%   matlab -batch "addpath('src','tests'); test_classical_filters"
%
% Checks
%   1. linearise_upright.m (p.Ad, p.Bd) is cross-checked against an
%      independent ZOH discretisation built from jacobian_f.m — Pessoa C's
%      analytical Jacobian, developed independently for the EKF. Agreement
%      confirms the two linearisations of the same plant (done by two
%      different people) are consistent, and that c2d/ZOH was applied
%      correctly.
%   2. The discretised A_d keeps the open-loop instability of the upright
%      equilibrium (one eigenvalue with |eig|>1), consistent with the
%      Routh-Hurwitz analysis in the report.
%   3. Smoke test: all three classical estimators run on a synthetic
%      trajectory, stay finite, and — the three filters had NO handling for
%      dropped (NaN) measurements before this pass — now correctly coast
%      through them without producing NaN/Inf (see the fixes in
%      est_lowpass.m, est_complementary.m and est_kf.m).
%
% This file ALSO regenerates the two figures used in the report (Seção
% "Estimadores de Estado Clássicos"), saved to results/. It uses a small
% closed-loop demo (state feedback with the LQR gain already published in
% the report) purely to produce a representative near-upright trajectory —
% this is NOT Pessoa A's controller_lqr.m/swingup.m (still TODO at the time
% of writing) and should not be confused with it; see the comment above
% demo_closed_loop_run() below.

function test_classical_filters()

fprintf('=== test_classical_filters ===\n');
rng(0);            % reproducible measurement noise
p = params();

npass = 0; nfail = 0;

% ---------------------------------------------------------------------------
% Test 1 — linearise_upright.m cross-checked against jacobian_f.m (Pessoa C)
% ---------------------------------------------------------------------------
M = p.M; m = p.m; L = p.L; g = p.g; b = p.b; %#ok<NASGU>
I = m*L^2/3;
D = (M+m)*(I+m*L^2) - (m*L)^2;

A_c_ref = jacobian_f(zeros(4,1), 0, p);      % Pessoa C's independent Jacobian
B_c_ref = [0; (I+m*L^2)/D; 0; (m*L)/D];      % analytical B (report, Seção LQR)

Maug   = [A_c_ref, B_c_ref; zeros(1,5)];
Md_ref = expm(Maug * p.dt);
Ad_ref = Md_ref(1:4,1:4);
Bd_ref = Md_ref(1:4,5);

err_Ad = max(abs(p.Ad(:) - Ad_ref(:)));
err_Bd = max(abs(p.Bd(:) - Bd_ref(:)));
tol = 1e-8;
if err_Ad < tol && err_Bd < tol
    fprintf('[PASS] linearise_upright vs jacobian_f cross-check: err_Ad=%.2e, err_Bd=%.2e\n', err_Ad, err_Bd);
    npass = npass + 1;
else
    fprintf('[FAIL] linearise_upright vs jacobian_f cross-check: err_Ad=%.2e, err_Bd=%.2e (tol %.0e)\n', err_Ad, err_Bd, tol);
    nfail = nfail + 1;
end

% ---------------------------------------------------------------------------
% Test 2 — discrete model keeps the upright instability (|eig(Ad)|>1)
% ---------------------------------------------------------------------------
eig_Ad = eig(p.Ad);
if any(abs(eig_Ad) > 1 + 1e-9)
    fprintf('[PASS] A_d preserves the unstable mode (max |eig|=%.4f)\n', max(abs(eig_Ad)));
    npass = npass + 1;
else
    fprintf('[FAIL] A_d lost the unstable mode (max |eig|=%.4f) — check linearise_upright.m\n', max(abs(eig_Ad)));
    nfail = nfail + 1;
end

% ---------------------------------------------------------------------------
% Test 3 — smoke test of the three classical filters, including NaN coasting
% ---------------------------------------------------------------------------
x_probe = [0.0; 0.3; 0.4; 0.5];
if norm(plant_cartpole(0, x_probe, 1.0, p)) < 1e-12
    fprintf('plant_cartpole.m is still a stub -> using local_dynamics.\n');
    ref = @(t,x,u,pp) local_dynamics(x, u, pp);
else
    fprintf('plant_cartpole.m implemented -> validating against it.\n');
    ref = @plant_cartpole;
end

filters = { @est_lowpass,       'lowpass';
            @est_complementary, 'complementary';
            @est_kf,             'kf' };

for fi = 1:size(filters,1)
    est_fn = filters{fi,1};
    name   = filters{fi,2};
    [ok, msg] = run_filter_smoke(est_fn, ref, [0;0;deg2rad(5);0], p);
    if ok
        fprintf('[PASS] %s smoke (incl. NaN coasting): %s\n', name, msg);
        npass = npass + 1;
    else
        fprintf('[FAIL] %s smoke: %s\n', name, msg);
        nfail = nfail + 1;
    end
end

% ---------------------------------------------------------------------------
fprintf('--------------------------------------\n');
fprintf('Total: %d passed, %d failed\n', npass, nfail);
if nfail > 0
    error('test_classical_filters: %d test(s) failed.', nfail);
end
fprintf('ALL TESTS PASSED\n');

% ---------------------------------------------------------------------------
% Regenerate the report figures (Seção "Estimadores de Estado Clássicos").
% ---------------------------------------------------------------------------
generate_report_figures(p, ref);

end


% ===========================================================================
% Smoke-test helper: propagate truth (open loop, small perturbation about
% upright), feed noisy (occasionally dropped) measurements to the estimator,
% check numerics stay finite/bounded.
% ===========================================================================
function [ok, msg] = run_filter_smoke(est_fn, ref, x0_true, p)
    N       = 200;
    x_true  = x0_true;
    xhat    = p.x0hat;
    P       = p.P0;
    u       = 0;   % open loop is enough to exercise the filter numerically
    err_obs = zeros(N,1);

    for k = 1:N
        if mod(k,10) == 0
            z = [NaN; NaN];   % dropped measurement (S4 convention)
        else
            z = [x_true(1); x_true(3)] + [p.sigma_x; p.sigma_theta].*randn(2,1);
        end

        [xhat, P] = est_fn(xhat, P, z, u, p);

        if any(~isfinite(xhat)) || any(~isfinite(P(:)))
            ok = false; msg = sprintf('non-finite at step %d', k); return
        end

        err_obs(k) = norm(xhat([1 3]) - x_true([1 3]));
        x_true = rk4_step(@(t,x) ref(t,x,u,p), 0, x_true, p.dt);
    end

    obs_err = mean(err_obs(end-49:end));
    div_tol = 50;   % loose bound: open-loop pole falls, filter must not blow up faster than truth
    if obs_err > div_tol
        ok = false; msg = sprintf('diverged (obs err %.1f)', obs_err); return
    end
    ok  = true;
    msg = sprintf('obs err %.4f', obs_err);
end


% ===========================================================================
% generate_report_figures — builds the two figures referenced in main.tex
% (\ref{fig:s1-classicos-theta}, \ref{fig:classicos-rmse}) and saves them to
% results/.
%
% IMPORTANT: the closed loop below (u = -K_demo*x_true, K_demo taken from the
% report's already-published LQR design, Seção "Projeto do controlador LQR")
% exists ONLY to produce a representative near-upright trajectory to feed the
% three classical estimators. It is a self-contained demo, not a wiring of
% Pessoa A's controller_lqr.m/swingup.m (both still TODO at the time of
% writing) — nothing here is written to those files. Once Pessoa A's
% controller and plant are ready, main.m + run_scenario.m should be used
% instead to regenerate these figures from the real closed loop / real
% scenarios (S1, S3).
% ===========================================================================
function generate_report_figures(p, ref)

results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

K_demo = [-40.00, -42.14, 217.20, 48.46];   % report, Seção "Projeto do controlador LQR"
u_max  = 20;

x0_true = [0; 0; deg2rad(5); 0];            % same initial condition as S1
Tf      = 10;
N       = round(Tf / p.dt);
t       = (0:N-1)' * p.dt;

noise_levels = struct('name', {'nominal', 'S3 (10x)'}, 'mult', {1, 10});

rmse_theta = zeros(3, numel(noise_levels));   % rows: lowpass, comp, kf
rmse_x     = zeros(3, numel(noise_levels));
theta_true_nom = zeros(N,1);
theta_hat_nom  = zeros(N,3);   % columns: lowpass, comp, kf

for ni = 1:numel(noise_levels)
    mult = noise_levels(ni).mult;
    p_run   = p;
    p_run.R = diag([(mult*p.sigma_x)^2, (mult*p.sigma_theta)^2]);  % fair KF R for this noise level

    rng(0);   % same noise realisation across estimators/noise-levels for a fair comparison
    x_true = x0_true;
    xh_lp = p.x0hat; P_lp = p.P0;
    xh_cp = p.x0hat; P_cp = p.P0;
    xh_kf = p.x0hat; P_kf = p.P0;

    err_theta = zeros(N,3);
    err_x     = zeros(N,3);

    for k = 1:N
        z = [x_true(1); x_true(3)] + mult*[p.sigma_x; p.sigma_theta].*randn(2,1);

        [xh_lp, P_lp] = est_lowpass(xh_lp, P_lp, z, 0, p_run);
        [xh_cp, P_cp] = est_complementary(xh_cp, P_cp, z, 0, p_run);
        [xh_kf, P_kf] = est_kf(xh_kf, P_kf, z, 0, p_run);

        err_theta(k,:) = [xh_lp(3), xh_cp(3), xh_kf(3)] - x_true(3);
        err_x(k,:)     = [xh_lp(1), xh_cp(1), xh_kf(1)] - x_true(1);

        if ni == 1
            theta_true_nom(k) = x_true(3);
            theta_hat_nom(k,:) = [xh_lp(3), xh_cp(3), xh_kf(3)];
        end

        % Demo closed loop: full-state feedback on the TRUE state (not on the
        % estimate) — see the comment above this function.
        u = -K_demo * x_true;
        u = max(-u_max, min(u_max, u));
        x_true = rk4_step(@(tt,xx) ref(tt,xx,u,p), 0, x_true, p.dt);
    end

    rmse_theta(:,ni) = sqrt(mean(err_theta.^2, 1))';
    rmse_x(:,ni)     = sqrt(mean(err_x.^2, 1))';
end

names = {'Baseline', 'Complementar', 'KF linear'};
fprintf('\n--- RMSE (demo, near-upright, K publicado no relatorio) ---\n');
for ei = 1:3
    fprintf('%-14s | theta: nominal=%.4f rad (%.3f deg), S3=%.4f rad (%.3f deg) | x: nominal=%.4f m, S3=%.4f m\n', ...
        names{ei}, rmse_theta(ei,1), rad2deg(rmse_theta(ei,1)), rmse_theta(ei,2), rad2deg(rmse_theta(ei,2)), ...
        rmse_x(ei,1), rmse_x(ei,2));
end

% --- Figure 1: theta tracking (nominal noise) ---
fig1 = figure('Visible', 'off', 'Color', 'white');
plot(t, rad2deg(theta_true_nom), 'k-', 'LineWidth', 1.6); hold on
plot(t, rad2deg(theta_hat_nom(:,1)), '--', 'LineWidth', 1.1);
plot(t, rad2deg(theta_hat_nom(:,2)), '-.', 'LineWidth', 1.1);
plot(t, rad2deg(theta_hat_nom(:,3)), ':',  'LineWidth', 1.6);
grid on
ax1 = gca;
force_light_axes(ax1);
xlabel('Tempo (s)', 'Color', 'k'); ylabel('\theta (graus)', 'Color', 'k');
lg1 = legend({'Verdadeiro', 'Baseline', 'Complementar', 'KF linear'}, 'Location', 'best');
set(lg1, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', [0.4 0.4 0.4]);
title('Rastreamento de \theta — cenario tipo S1 (ruido nominal)', 'Color', 'k');
exportgraphics(fig1, fullfile(results_dir, 'S1_classicos_theta.png'), 'Resolution', 150, 'BackgroundColor', 'white');
close(fig1);

% --- Figure 2: RMSE comparison, nominal vs S3 (10x noise) ---
fig2 = figure('Visible', 'off', 'Color', 'white');
bar(rad2deg(rmse_theta));
ax2 = gca;
force_light_axes(ax2);
set(ax2, 'XTickLabel', names);
ylabel('RMSE de \theta (graus)', 'Color', 'k');
lg2 = legend({'S1 (nominal)', 'S3 (ruido 10x)'}, 'Location', 'best');
set(lg2, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', [0.4 0.4 0.4]);
title('RMSE de \theta por estimador classico — S1 vs S3', 'Color', 'k');
grid on
exportgraphics(fig2, fullfile(results_dir, 'S1_S3_classicos_rmse.png'), 'Resolution', 150, 'BackgroundColor', 'white');
close(fig2);

fprintf('Figures saved to %s\n', results_dir);

end


% ===========================================================================
% force_light_axes — MATLAB batch mode picks up a dark app theme by default
% (R2025b), which makes axes/legend text nearly invisible once exported onto
% a white background for the printed report. Force a plain light theme.
% ===========================================================================
function force_light_axes(ax)
    set(ax, 'Color', 'white', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.3);
end


% ===========================================================================
% local_dynamics — synthetic cart-pole EOM (KEEP IN SYNC with the equations in
% plant_cartpole.m's header / jacobian_f.m). Duplicated from
% tests/test_nonlinear_filters.m so this file works standalone; used only
% until Pessoa A implements the plant.
% ===========================================================================
function xdot = local_dynamics(x, u, p)
    M = p.M; m = p.m; L = p.L; g = p.g; b = p.b;
    xd  = x(2); th = x(3); thd = x(4);
    a = M + m;
    d = m*L^2/3 + m*L^2;      % I + m L^2
    k = m*L;
    c = cos(th); s = sin(th);
    det  = a*d - k^2*c^2;
    rhs1 = u - b*xd - k*thd^2*s;
    rhs2 = g*k*s;
    xddot  = (d*rhs1 + k*c*rhs2) / det;
    thddot = (a*rhs2 + k*c*rhs1) / det;
    xdot = [xd; xddot; thd; thddot];
end


% --- RK4 integrator step -----------------------------------------------------
function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
