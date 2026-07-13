% build_report_results.m — Aggregate results/*.mat into report-ready tables and
% figures (Fase 2, full harness, all estimators x all scenarios).
% Run from the repo root:
%   matlab -batch "addpath('src','tests'); build_report_results"
%
% Outputs (written to results/):
%   summary_table.csv          — RMSE (per state), NEES/NIS mean, t_converge, cpu_time
%                                 for every (scenario, estimator) pair
%   S2_swingup_theta.png       — theta(t) true vs all 5 estimators, swing-up (S2)
%   S1_S3_rmse_theta.png       — RMSE(theta) grouped bar, all estimators, S1 vs S3
%   S2_nees.png                — NEES(t) EKF vs UKF vs KF during swing-up (S2), with
%                                 chi-square 95% bounds (consistency check)
%   cpu_time.png               — mean cpu_time per estimator (log scale)
%   S0_ideal_traj.png          — theta(t)/x(t) de referência sem ruído de medição (S0)
%   S0_ideal_rmse.png          — RMSE(theta) por estimador em S0 (escala log), mostra
%                                 o piso de erro de baseline/complementar mesmo sem ruído

function build_report_results()

% Force a light theme (white background, black text/axes) regardless of the
% environment's default — figures are embedded in a printed LaTeX report.
set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesColor', 'white');
set(groot, 'DefaultAxesXColor', 'black');
set(groot, 'DefaultAxesYColor', 'black');
set(groot, 'DefaultAxesGridColor', [0.15 0.15 0.15]);
set(groot, 'DefaultTextColor', 'black');
set(groot, 'DefaultLegendTextColor', 'black');

results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');

est_names = {'lowpass','complementary','kf','ekf','ukf'};
est_labels = {'Baseline','Complementar','KF','EKF','UKF'};
scen_names = {'S1_stabilisation','S2_swingup','S3_high_noise', ...
              'S4_low_rate','S5_mistuned_QR','S6_bad_init'};

% --- Load everything into a struct array ------------------------------------
n_s = numel(scen_names); n_e = numel(est_names);
R = cell(n_s, n_e);
for si = 1:n_s
    for ei = 1:n_e
        fname = fullfile(results_dir, sprintf('%s_%s.mat', scen_names{si}, est_names{ei}));
        R{si,ei} = load(fname);
    end
end

% --- Summary CSV --------------------------------------------------------------
fid = fopen(fullfile(results_dir, 'summary_table.csv'), 'w');
fprintf(fid, 'scenario,estimator,rmse_x,rmse_xdot,rmse_theta_deg,rmse_thetadot,nees_mean,nis_mean,t_converge,cpu_time_ms\n');
for si = 1:n_s
    for ei = 1:n_e
        m = R{si,ei}.metrics;
        fprintf(fid, '%s,%s,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f\n', ...
            scen_names{si}, est_names{ei}, ...
            m.rmse(1), m.rmse(2), rad2deg(m.rmse(3)), m.rmse(4), ...
            m.nees_mean, nz(m.nis_mean), nz(m.t_converge), m.cpu_time*1000);
    end
end
fclose(fid);
fprintf('Wrote %s\n', fullfile(results_dir, 'summary_table.csv'));

% --- Print to console too (handy for copy-pasting into the .tex) ------------
fprintf('\n=== RMSE theta (deg) ===\n');
print_metric_grid(R, scen_names, est_labels, @(m) rad2deg(m.rmse(3)));
fprintf('\n=== RMSE x (m) ===\n');
print_metric_grid(R, scen_names, est_labels, @(m) m.rmse(1));
fprintf('\n=== RMSE theta_dot (rad/s) ===\n');
print_metric_grid(R, scen_names, est_labels, @(m) m.rmse(4));
fprintf('\n=== NEES mean (~4 if consistent) ===\n');
print_metric_grid(R, scen_names, est_labels, @(m) m.nees_mean);
fprintf('\n=== NIS mean (~2 if consistent; NaN = filter does not expose innovation) ===\n');
print_metric_grid(R, scen_names, est_labels, @(m) m.nis_mean);
fprintf('\n=== t_converge (s) ===\n');
print_metric_grid(R, scen_names, est_labels, @(m) m.t_converge);
fprintf('\n=== cpu_time (ms, total over run) ===\n');
print_metric_grid(R, scen_names, est_labels, @(m) m.cpu_time*1000);

% =============================================================================
% Figure 1 — S2 (swing-up) theta(t): true vs all 5 estimators
% =============================================================================
s2 = 2;
figure('Visible','off','Position',[100 100 900 500]); theme(gcf,'light');
hold on;
t = R{s2,1}.log.t;
plot(t, rad2deg(wrap_pi(R{s2,1}.log.x_true(:,3))), 'k', 'LineWidth', 2, 'DisplayName','Verdade');
colors = lines(n_e);
for ei = 1:n_e
    plot(t, rad2deg(wrap_pi(R{s2,ei}.log.xhat(:,3))), 'Color', colors(ei,:), ...
        'LineWidth', 1.1, 'DisplayName', est_labels{ei});
end
xlabel('t (s)'); ylabel('\theta (graus)');
title('Swing-up (S2): \theta estimado vs verdade — todos os estimadores');
legend('Location','southeast'); grid on;
saveas(gcf, fullfile(results_dir, 'S2_swingup_theta.png'));
close(gcf);

% =============================================================================
% Figure 2 — RMSE(theta) grouped bar, S1 vs S3, all 5 estimators
% =============================================================================
figure('Visible','off','Position',[100 100 800 500]); theme(gcf,'light');
rmse_s1 = arrayfun(@(ei) rad2deg(R{1,ei}.metrics.rmse(3)), 1:n_e);
rmse_s3 = arrayfun(@(ei) rad2deg(R{3,ei}.metrics.rmse(3)), 1:n_e);
bar([rmse_s1; rmse_s3]');
set(gca, 'XTickLabel', est_labels);
ylabel('RMSE \theta (graus)');
legend({'S1 (ruído nominal)','S3 (ruído 10\times)'}, 'Location','northwest');
title('RMSE de \theta por estimador — S1 vs S3');
grid on;
saveas(gcf, fullfile(results_dir, 'S1_S3_rmse_theta.png'));
close(gcf);

% =============================================================================
% Figure 3 — NEES(t) during swing-up (S2): KF vs EKF vs UKF, chi2 95% bounds
% =============================================================================
figure('Visible','off','Position',[100 100 900 500]); theme(gcf,'light');
hold on;
lo = chi2inv(0.025, 4); hi = chi2inv(0.975, 4);
kf_i = find(strcmp(est_names,'kf')); ekf_i = find(strcmp(est_names,'ekf')); ukf_i = find(strcmp(est_names,'ukf'));
plot(t, R{s2,kf_i}.metrics.nees, 'Color', colors(kf_i,:), 'DisplayName','KF');
plot(t, R{s2,ekf_i}.metrics.nees, 'Color', colors(ekf_i,:), 'DisplayName','EKF');
plot(t, R{s2,ukf_i}.metrics.nees, 'Color', colors(ukf_i,:), 'DisplayName','UKF');
yline(lo, '--k', 'HandleVisibility','off'); yline(hi, '--k', 'HandleVisibility','off');
set(gca,'YScale','log');
xlabel('t (s)'); ylabel('NEES');
title('NEES durante o swing-up (S2) — faixa tracejada = intervalo \chi^2_4 95%');
legend('Location','northeast'); grid on;
saveas(gcf, fullfile(results_dir, 'S2_nees.png'));
close(gcf);

% =============================================================================
% Figure 4 — mean cpu_time per estimator (log scale), averaged over scenarios
% =============================================================================
figure('Visible','off','Position',[700 100 700 500]); theme(gcf,'light');
cpu_mean = zeros(1,n_e);
for ei = 1:n_e
    vals = arrayfun(@(si) R{si,ei}.metrics.cpu_time, 1:n_s);
    cpu_mean(ei) = mean(vals) * 1000;
end
bar(cpu_mean);
set(gca, 'XTickLabel', est_labels, 'YScale','log');
ylabel('cpu\_time médio (ms, escala log)');
title('Custo computacional médio por estimador (média sobre S1..S6)');
grid on;
saveas(gcf, fullfile(results_dir, 'cpu_time.png'));
close(gcf);

% =============================================================================
% Figures 5/6 — S0 (comportamento ideal, sem ruído de medição)
% Carregado à parte de scen_names/R para não alterar índices/médias já usados
% acima (e já referenciados no relatório) pelas figuras/tabela de S1..S6.
% =============================================================================
R0 = struct();
for ei = 1:n_e
    fname = fullfile(results_dir, sprintf('S0_ideal_%s.mat', est_names{ei}));
    R0.(est_names{ei}) = load(fname);
end

% Trajetória de referência: verdade sob controle real (swing-up -> LQR) com
% sensores perfeitos. Com R~0 as 3 estimadores baseados em modelo (KF/EKF/UKF)
% ficam indistinguíveis da verdade; usa-se a verdade do UKF como referência.
t0 = R0.ukf.log.t;
figure('Visible','off','Position',[100 100 900 600]); theme(gcf,'light');
subplot(2,1,1);
plot(t0, rad2deg(wrap_pi(R0.ukf.log.x_true(:,3))), 'k', 'LineWidth', 1.6);
ylabel('\theta (graus)'); grid on;
title('Comportamento ideal sem ruído (S0): swing-up + estabilização');
subplot(2,1,2);
plot(t0, R0.ukf.log.x_true(:,1), 'k', 'LineWidth', 1.6);
xlabel('t (s)'); ylabel('x (m)'); grid on;
saveas(gcf, fullfile(results_dir, 'S0_ideal_traj.png'));
close(gcf);

% RMSE(theta) por estimador em S0: mesmo sem ruído de medição, baseline e
% complementar mantêm erro residual (limitação estrutural, não do ruído);
% KF/EKF/UKF colapsam para ~0.
figure('Visible','off','Position',[100 100 800 500]); theme(gcf,'light');
rmse_s0 = arrayfun(@(ei) rad2deg(R0.(est_names{ei}).metrics.rmse(3)), 1:n_e);
bar(max(rmse_s0, 1e-6));   % floor para o eixo log não estourar em zero exato
set(gca, 'XTickLabel', est_labels, 'YScale', 'log');
ylabel('RMSE \theta (graus, escala log)');
title('RMSE de \theta por estimador em S0 (sem ruído de medição)');
grid on;
saveas(gcf, fullfile(results_dir, 'S0_ideal_rmse.png'));
close(gcf);

% Acrescenta as linhas de S0 ao summary_table.csv (não mexe nas linhas S1..S6
% já escritas acima).
fid = fopen(fullfile(results_dir, 'summary_table.csv'), 'a');
for ei = 1:n_e
    m = R0.(est_names{ei}).metrics;
    fprintf(fid, '%s,%s,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f\n', ...
        'S0_ideal', est_names{ei}, ...
        m.rmse(1), m.rmse(2), rad2deg(m.rmse(3)), m.rmse(4), ...
        m.nees_mean, nz(m.nis_mean), nz(m.t_converge), m.cpu_time*1000);
end
fclose(fid);

fprintf('\nFigures written to %s\n', results_dir);

end

% --- helpers -------------------------------------------------------------
function v = nz(x)
    if isnan(x), v = -1; else, v = x; end
end

function w = wrap_pi(th)
    w = mod(th + pi, 2*pi) - pi;
end

function print_metric_grid(R, scen_names, est_labels, fn)
    n_s = numel(scen_names); n_e = numel(est_labels);
    fprintf('%-16s', '');
    for ei = 1:n_e, fprintf('%14s', est_labels{ei}); end
    fprintf('\n');
    for si = 1:n_s
        fprintf('%-16s', scen_names{si});
        for ei = 1:n_e
            fprintf('%14.4f', fn(R{si,ei}.metrics));
        end
        fprintf('\n');
    end
end
