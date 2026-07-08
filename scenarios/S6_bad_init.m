% S6 — Poor initialisation (large initial covariance and wrong prior)

function s = S6_bad_init()
    s.name    = 'S6_bad_init';
    s.x0_true = [0; 0; deg2rad(5); 0];
    s.Tf      = 10;
    % Estimator starts far from truth with high uncertainty
    s.x0hat   = [0.5; 0; deg2rad(30); 0];
    s.P0      = 100 * diag([1, 1, 1, 1]);
end
