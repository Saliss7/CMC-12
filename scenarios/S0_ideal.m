% S0 — Comportamento ideal (sem ruído de medição)
% Mesma condição inicial e horizonte do swing-up (S2), mas com sensores
% perfeitos (sigma_x = sigma_theta = 0, sem dropout/outlier). Serve como
% referência do comportamento-alvo do controlador (swing-up + LQR) antes de
% qualquer degradação de sensor ser introduzida.

function s = S0_ideal()
    s.name         = 'S0_ideal';
    s.x0_true      = [0; 0; pi; 0];   % hanging down
    s.Tf           = 15;               % s (needs longer for swing-up)
    s.sigma_x      = 0;
    s.sigma_theta  = 0;
    s.dropout_prob = 0;
    s.outlier_prob = 0;
    s.R            = diag([1e-10, 1e-10]);   % near-zero: estimators trust the (exact) measurement
end
