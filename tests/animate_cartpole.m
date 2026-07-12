% animate_cartpole.m — Animate cart-pole trajectories from results/*.mat.
% Owner: Pessoa C (viz auxiliar)
%
% Draws the cart+pole for the TRUE state and for one or more ESTIMATED
% states side by side, subsampled in time, and exports a GIF (and, if
% VideoWriter/MP4 codec is available, an MP4) to results/.
%
% Usage (from repo root, after src/main.m has produced results/*.mat):
%   addpath('src', 'tests');
%   animate_cartpole('S2_swingup', {'kf','ekf','ukf'})
%   animate_cartpole('S1_stabilisation', {'ekf'})
%
% Inputs
%   scenario_name — e.g. 'S2_swingup' (must match a results/<name>_*.mat set)
%   est_names     — cell array of estimator name suffixes to overlay, e.g.
%                   {'kf','ekf','ukf'} (subset of lowpass/complementary/kf/ekf/ukf)
%   opts.fps      — playback frames per second in the output file (default 30)
%   opts.speed    — simulation seconds per output second (default 1, i.e.
%                   real-time; use e.g. 0.25 to slow down 4x)
%   opts.out      — output basename without extension (default
%                   '<scenario_name>_anim')
%
% Output
%   results/<out>.gif  (always)
%   results/<out>.mp4  (if an MPEG-4 VideoWriter profile is available)

function animate_cartpole(scenario_name, est_names, opts)

if nargin < 2 || isempty(est_names), est_names = {'ekf', 'ukf'}; end
if nargin < 3, opts = struct(); end
fps   = getfield_default(opts, 'fps', 30);
speed = getfield_default(opts, 'speed', 1.0);
out   = getfield_default(opts, 'out', [scenario_name '_anim']);

% Force a light theme regardless of the environment's default (matches
% build_report_results.m so all exported media look consistent).
set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesColor', 'white');
set(groot, 'DefaultAxesXColor', 'black');
set(groot, 'DefaultAxesYColor', 'black');
set(groot, 'DefaultAxesGridColor', [0.15 0.15 0.15]);
set(groot, 'DefaultTextColor', 'black');

results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');

% --- Load truth + estimate PER ESTIMATOR ---------------------------------------
% Note: the simulation is closed-loop (u_k depends on xhat), so x_true differs
% across estimator runs — each panel gets its own truth, not a shared one.
n_e = numel(est_names);
x_true_all = cell(1, n_e);
xhat_all   = cell(1, n_e);
t = []; L = [];
for i = 1:n_e
    Si = load(fullfile(results_dir, sprintf('%s_%s.mat', scenario_name, est_names{i})));
    x_true_all{i} = Si.log.x_true;
    xhat_all{i}   = Si.log.xhat;
    if isempty(t), t = Si.log.t; L = Si.p.L; end
end
labels_map = struct('lowpass','Baseline', 'complementary','Complementar', ...
                     'kf','KF', 'ekf','EKF', 'ukf','UKF');
labels = cellfun(@(n) labels_map.(n), est_names, 'UniformOutput', false);

% --- Subsample to the target output frame rate --------------------------------
dt        = t(2) - t(1);
step      = max(1, round((speed / fps) / dt));
frame_idx = 1:step:numel(t);

% --- Figure layout: one panel per estimator, EACH WITH ITS OWN SCALE ---------
% A diverging filter (e.g. KF in swing-up) drives the cart far from the origin;
% sharing one axis across panels would make the well-behaved filters invisible.
n_panels = n_e;
fig = figure('Visible', 'off', 'Position', [100 100 380*n_panels+40 420]);
set(fig, 'Color', 'white');

axs = gobjects(1, n_panels);
cart_w = 0.2; cart_h = 0.1;
xlim_pad = zeros(1, n_panels);
for pnl = 1:n_panels
    span = max(abs([x_true_all{pnl}(:,1); xhat_all{pnl}(:,1)]));
    xlim_pad(pnl) = max(0.6, span * 1.15);
end

colors = lines(n_panels);

ylim_pad = max(0.6, 2.3 * L);

for pnl = 1:n_panels
    axs(pnl) = subplot(1, n_panels, pnl);
    hold(axs(pnl), 'on'); axis(axs(pnl), 'equal');
    xlim(axs(pnl), [-xlim_pad(pnl), xlim_pad(pnl)]); ylim(axs(pnl), [-ylim_pad, ylim_pad]);
    grid(axs(pnl), 'on');
    title(axs(pnl), sprintf('Verdade vs %s', labels{pnl}));
    xlabel(axs(pnl), 'x (m)');
    yline(axs(pnl), 0, 'k-');
end

frames(numel(frame_idx)) = struct('cdata', [], 'colormap', []);

for fi = 1:numel(frame_idx)
    k = frame_idx(fi);
    for pnl = 1:n_panels
        cla(axs(pnl));
        hold(axs(pnl), 'on'); axis(axs(pnl), 'equal');
        xlim(axs(pnl), [-xlim_pad(pnl), xlim_pad(pnl)]); ylim(axs(pnl), [-ylim_pad, ylim_pad]);
        grid(axs(pnl), 'on');
        yline(axs(pnl), 0, 'k-');

        draw_cartpole(axs(pnl), x_true_all{pnl}(k,1), x_true_all{pnl}(k,3), L, ...
                      cart_w, cart_h, [0.6 0.6 0.6], '--', 0.5);
        draw_cartpole(axs(pnl), xhat_all{pnl}(k,1), xhat_all{pnl}(k,3), L, ...
                      cart_w, cart_h, colors(pnl,:), '-', 2.0);

        title(axs(pnl), sprintf('%s — t = %.2f s', labels{pnl}, t(k)));
        xlabel(axs(pnl), 'x (m)');
    end
    drawnow;
    frames(fi) = getframe(fig);
end
close(fig);

% --- Write GIF -----------------------------------------------------------------
gif_path = fullfile(results_dir, [out '.gif']);
for fi = 1:numel(frames)
    [imind, cm] = rgb2ind(frames(fi).cdata, 256);
    if fi == 1
        imwrite(imind, cm, gif_path, 'gif', 'Loopcount', inf, 'DelayTime', 1/fps);
    else
        imwrite(imind, cm, gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end
fprintf('Wrote %s\n', gif_path);

% --- Write a video file (best-effort, profile availability is platform-dependent) ---
video_profiles = {'MPEG-4', '.mp4'; 'Motion JPEG AVI', '.avi'; 'Uncompressed AVI', '.avi'};
video_written = false;
for vi = 1:size(video_profiles, 1)
    try
        video_path = fullfile(results_dir, [out video_profiles{vi,2}]);
        vw = VideoWriter(video_path, video_profiles{vi,1});
        vw.FrameRate = fps;
        open(vw);
        for fi = 1:numel(frames)
            writeVideo(vw, frames(fi).cdata);
        end
        close(vw);
        fprintf('Wrote %s\n', video_path);
        video_written = true;
        break
    catch
        % try next profile
    end
end
if ~video_written
    fprintf('No video profile available on this system. GIF was written above.\n');
end

end


% ===========================================================================
% draw_cartpole — draw cart (rectangle), pole (line) and bob (marker) for a
% given (x, theta). theta = 0 is UP (project convention, see params.m).
% ===========================================================================
function draw_cartpole(ax, x, theta, L, cart_w, cart_h, color, ls, lw)
    % Cart body
    rectangle(ax, 'Position', [x-cart_w/2, -cart_h/2, cart_w, cart_h], ...
               'EdgeColor', color, 'LineStyle', ls, 'LineWidth', lw);
    % Pole: theta measured from upward vertical
    tip_x = x + 2*L*sin(theta);
    tip_y = 0 + 2*L*cos(theta);
    plot(ax, [x, tip_x], [0, tip_y], 'Color', color, 'LineStyle', ls, 'LineWidth', lw);
    plot(ax, tip_x, tip_y, 'o', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 6);
end

function v = getfield_default(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end
