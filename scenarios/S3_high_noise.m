% S3 — High measurement noise (10× default sigma)

function s = S3_high_noise()
    p = params();
    s.name        = 'S3_high_noise';
    s.x0_true     = [0; 0; deg2rad(5); 0];
    s.Tf          = 10;
    s.sigma_x     = 10 * p.sigma_x;
    s.sigma_theta = 10 * p.sigma_theta;
end
