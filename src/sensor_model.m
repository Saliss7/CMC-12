% sensor_model.m — Simulated sensor: adds noise, dropout, and outliers.
% Owner: Matheus Felipe
%
% Usage:
%   z = sensor_model(x_true, p, scenario)
%
% Inputs
%   x_true   — true state [cart_pos; cart_vel; pole_angle; pole_ang_vel] (4×1)
%   p        — params struct from params.m
%   scenario — (optional) scenario struct with overrides:
%                .sigma_x, .sigma_theta   — encoder noise std (m, rad)
%                .dropout_prob            — P(measurement lost) per step   [S4]
%                .outlier_prob            — P(spurious spike) per channel  [S3+]
%                .outlier_scale           — outlier size, in units of sigma
%
% Output
%   z — noisy measurement [cart_pos; pole_angle] (2×1), OR [] when the
%       measurement is dropped (see dropout below). Callers must treat an empty
%       z as "no measurement this step" and perform a prediction-only update.

function z = sensor_model(x_true, p, scenario)

if nargin < 3
    scenario = struct();
end

sx        = getfield_default(scenario, 'sigma_x',       p.sigma_x);
stheta    = getfield_default(scenario, 'sigma_theta',   p.sigma_theta);
p_drop    = getfield_default(scenario, 'dropout_prob',  0);
p_out     = getfield_default(scenario, 'outlier_prob',  0);
out_scale = getfield_default(scenario, 'outlier_scale', 30);   % spike vs sigma

% --- Dropout (S4) ---------------------------------------------------------
% With probability p_drop the whole reading is lost. Signalled by an empty z;
% the harness/filters treat it as a prediction-only step (no measurement update).
if rand < p_drop
    z = [];
    return;
end

% --- Nominal reading: cart position and pole angle + Gaussian encoder noise
sigma = [sx; stheta];
z = [x_true(1); x_true(3)] + sigma .* randn(2,1);

% --- Outliers (S3+) -------------------------------------------------------
% Each channel independently gets a large spurious spike with probability p_out.
% The spike is out_scale times the nominal std, with random sign — an unlikely
% reading a robust filter should reject.
if p_out > 0
    hit = rand(2,1) < p_out;
    z   = z + hit .* (out_scale .* sigma .* randn(2,1));
end

end

% --- helper ------------------------------------------------------------------
function v = getfield_default(s, field, default)
    if isfield(s, field)
        v = s.(field);
    else
        v = default;
    end
end
