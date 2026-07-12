% main.m — Entry point: runs all estimators on all scenarios and saves results.
%
% Usage (MATLAB command window):
%   run('src/main.m')
%
% Results are saved to results/<scenario>_<estimator>.mat

addpath(fileparts(mfilename('fullpath')));   % add src/ to path
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'scenarios'));

p = params();

estimators = {
    @est_lowpass,       'lowpass';
    @est_complementary, 'complementary';
    @est_kf,            'kf';
    @est_ekf,           'ekf';
    @est_ukf,           'ukf';
};

scenarios = {
    S1_stabilisation();
    S2_swingup();
    S3_high_noise();
    S4_low_rate();
    S5_mistuned_QR();
    S6_bad_init();
};

results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

for si = 1:numel(scenarios)
    sc = scenarios{si};
    for ei = 1:size(estimators, 1)
        est_fn   = estimators{ei,1};
        est_name = estimators{ei,2};

        fprintf('Running %s / %s ...\n', sc.name, est_name);
        log     = run_scenario(est_fn, sc, p);   % log.cpu_time is set inside
        metrics = compute_metrics(log, p);

        fname = fullfile(results_dir, sprintf('%s_%s.mat', sc.name, est_name));
        save(fname, 'log', 'metrics', 'sc', 'p');
    end
end

fprintf('Done. Results saved to %s\n', results_dir);
