% sensor_model.m — Simulated sensor: adds noise, dropout, and outliers.
% Owner: Pessoa A
%
% Usage:
%   z = sensor_model(x_true, p, scenario)
%
% Inputs
%   x_true   — true state [cart_pos; cart_vel; pole_angle; pole_ang_vel] (4×1)
%   p        — params struct from params.m
%   scenario — (optional) scenario struct with overrides:
%                .sigma_x, .sigma_theta, .dropout_prob, .outlier_prob
%
% Output
%   z — noisy measurement [cart_pos; pole_angle]  (2×1)

function z = sensor_model(x_true, p, scenario)

if nargin < 3
    scenario = struct();
end

sx     = getfield_default(scenario, 'sigma_x',      p.sigma_x);
stheta = getfield_default(scenario, 'sigma_theta',  p.sigma_theta);

% TODO (Pessoa A): add dropout and outlier logic for scenarios S3/S4

z = [x_true(1); x_true(3)] + [sx; stheta] .* randn(2,1);

end

% --- helper ------------------------------------------------------------------
function v = getfield_default(s, field, default)
    if isfield(s, field)
        v = s.(field);
    else
        v = default;
    end
end
