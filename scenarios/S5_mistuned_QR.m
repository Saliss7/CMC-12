% S5 — Mistuned Q/R (estimator uses wrong noise covariances)
% The sensor_model still uses nominal noise; estimators receive wrong p.Q / p.R.

function s = S5_mistuned_QR()
    p = params();
    s.name    = 'S5_mistuned_QR';
    s.x0_true = [0; 0; deg2rad(5); 0];
    s.Tf      = 10;
    % Override covariances seen by estimators (10× wrong)
    s.Q = p.Q * 10;
    s.R = p.R / 10;
end
