% S2 — Swing-up (strongly nonlinear regime)
% Pole starts hanging down; energy-based swing-up then hands off to LQR.

function s = S2_swingup()
    s.name       = 'S2_swingup';
    s.x0_true    = [0; 0; pi; 0];   % hanging down
    s.Tf         = 15;               % s (needs longer for swing-up)
end
