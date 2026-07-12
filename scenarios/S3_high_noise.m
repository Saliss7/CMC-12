% S3 — High measurement noise (10× default sigma)
%
% s.sigma_x/s.sigma_theta drive the actual noise injected by sensor_model.m;
% s.R is the covariance the ESTIMATORS assume (run_scenario.m only reads
% scenario.R, not the sigma_* fields). The two must be kept in sync here so
% S3 tests robustness to a genuinely higher noise floor with correctly-tuned
% filters — a mistuned R is a separate, deliberate stress test (S5).

function s = S3_high_noise()
    p = params();
    s.name        = 'S3_high_noise';
    s.x0_true     = [0; 0; deg2rad(5); 0];
    s.Tf          = 10;
    s.sigma_x     = 10 * p.sigma_x;
    s.sigma_theta = 10 * p.sigma_theta;
    s.R           = diag([s.sigma_x^2, s.sigma_theta^2]);
end
