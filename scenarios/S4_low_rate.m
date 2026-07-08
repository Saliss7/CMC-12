% S4 — Low sampling rate and random measurement dropout

function s = S4_low_rate()
    s.name         = 'S4_low_rate';
    s.x0_true      = [0; 0; deg2rad(5); 0];
    s.Tf           = 10;
    s.dropout_prob = 0.2;   % 20 % of measurements dropped
    % TODO (Pessoa A): implement dropout logic in sensor_model.m
end
