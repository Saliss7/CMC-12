% S1 — Stabilisation (near-linear regime)
% Pole starts near upright; only LQR phase is exercised.

function s = S1_stabilisation()
    s.name       = 'S1_stabilisation';
    s.x0_true    = [0; 0; deg2rad(5); 0];   % 5 deg off upright
    s.Tf         = 10;                        % s
    % no sensor overrides — use defaults from params
end
