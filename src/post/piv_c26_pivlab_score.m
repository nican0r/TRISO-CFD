function score = piv_c26_pivlab_score(export_dir, variant, opts)
%PIV_C26_PIVLAB_SCORE  C26-style scoring of a PIVlab run, in MATLAB with MATLAB defaults.
%
%   score = piv_c26_pivlab_score(export_dir, variant)
%
% Takes the PIVlab fields written by piv_c26_pivlab_run.m for both views
% (<export_dir>/<variant>/{spout,wall}/fields.mat) and computes what C26
% (Conry et al. 2026, NED 456:115009) reports, the way C26 did it (MATLAB /
% PIVlab), with MATLAB's defaults wherever C26 gives no detail:
%
%   calibration : PIVlab's (u, v) * calxy / time with calxy = pixel size [m/px]
%                 and time = frame interval [s] (export.json); y down, so +v
%                 is C26's "downward positive".
%   spout map   : time mean of the vectors over all pairs, pixel by pixel
%                 (C26 3.1), then the magnitude of the mean vector; the map is
%                 cut to C26's Fig. 4 window (2 in. wide, from the annulus
%                 surface; export.json "map") and put on the grid of the
%                 digitized experiment panel with interp2 (default linear).
%   metrics     : CC = corrcoef, SSI = ssim (MATLAB defaults), MSE = immse,
%                 each against the digitized C26 Fig. 4a experiment map, and
%                 the same three for the digitized C26 Fig. 4b model map so
%                 the paper's own model is scored through the same code.
%   wall        : PIVlab's area extraction (mean over the vectors in the ~4 mm
%                 box, C26 Fig. 1c, 'omitnan' as +extract/get_mean_of_selection)
%                 per pair, then the time mean and SD (C26 3.2, Fig. 6d).
%
% Output: <export_dir>/<variant>/score.json and map_fig4grid.csv (top row
% first, as data/digitized).  Returns the score struct.

arguments
    export_dir (1,1) string
    variant (1,1) string
    opts.digitized_dir (1,1) string = ""     % default: <repo>/data/digitized
end

meta = jsondecode(fileread(fullfile(export_dir, 'export.json')));
if strlength(opts.digitized_dir) == 0
    here = fileparts(mfilename('fullpath'));
    opts.digitized_dir = fullfile(here, '..', '..', 'data', 'digitized');
end
px = meta.pixel_m;            % [m/px]
dt = meta.frame_dt_s;         % [s]
cal = px / dt;                % PIVlab: calu = calv = calxy / time  [m/s per px/frame]

score = struct('run', meta.run, 'case', meta.case, 'flow_slpm', meta.flow_slpm, 'variant', char(variant), ...
    'pixel_m', px, 'frame_dt_s', dt, 'window_t0_s', meta.frame_times_s(1), 'window_t1_s', meta.frame_times_s(end));

%% ---- spout map -----------------------------------------------------------------
S = load(fullfile(export_dir, variant, 'spout', 'fields.mat'), 'x', 'y', 'u', 'v', 'stats');
score.n_pairs_spout = size(S.u, 3);
Um = mean(S.u, 3, 'omitnan') * cal;         % time mean of the vectors [m/s], x right
Vm = mean(S.v, 3, 'omitnan') * cal;         % y down (C26 downward positive)
mag = hypot(Um, Vm);                        % |mean vector|
mag(isnan(mag)) = 0;
xs = S.x(1, :, 1); ys = S.y(:, 1, 1)';     % PIVlab grid, 1-based pixel coordinates
ny = meta.views.spout.shape_px(1);
x_m = meta.views.spout.x0_m + (xs - 0.5) * px;            % pixel j centre
y_m = meta.views.spout.y0_m + (ny - ys + 0.5) * px;       % row from the bottom
% Fig. 4 window and the digitized experiment panel of the same flow
R_c = meta.map.R_c_m; y_map = meta.map.y_map_m; H = meta.map.H_m;
flow_panel = meta.map.panel_flow_slpm;
exp_map = read_digitized(fullfile(opts.digitized_dir, sprintf('c26_fig4_exp_%gslpm.csv', flow_panel)));
mod_map = read_digitized(fullfile(opts.digitized_dir, sprintf('c26_fig4_model_%gslpm.csv', flow_panel)));
[nrow, ncol] = size(exp_map);
xc = -R_c + (2 * R_c) * ((1:ncol) - 0.5) / ncol;          % panel cell centres
yc = y_map + H * ((nrow:-1:1) - 0.5) / nrow;              % top row first, as the CSV
[XI, YI] = meshgrid(xc, yc);
[XG, YG] = meshgrid(x_m, y_m);
map = interp2(XG, YG, mag, XI, YI);                       % MATLAB default: linear
map(isnan(map)) = 0;
score.map_mean_m_s = mean(map(:)); score.map_max_m_s = max(map(:));
score.map_grid = size(map);
score.spout_nan_frac_final = S.stats.nan_frac_final;
score.spout_rejected_frac = S.stats.rejected_frac;
score.metrics_vs_exp = metrics(map, exp_map);
score.metrics_c26_model_vs_exp = metrics(mod_map, exp_map);
score.metrics_vs_c26_model = metrics(map, mod_map);
if abs(flow_panel - meta.flow_slpm) > 1e-9
    score.metrics_note = sprintf('no C26 panel at %g slpm; scored against the %g slpm panel', meta.flow_slpm, flow_panel);
end
writematrix(map, fullfile(export_dir, variant, 'map_fig4grid.csv'));

%% ---- wall velocity ---------------------------------------------------------------
W = load(fullfile(export_dir, variant, 'wall', 'fields.mat'), 'x', 'y', 'v', 'stats');
box = meta.views.wall.pivlab_box_px;
xw = W.x(:, :, 1); yw = W.y(:, :, 1);
BW = xw >= box.x_min & xw <= box.x_max & yw >= box.y_min & yw <= box.y_max;
n_pairs_w = size(W.v, 3);
per_pair = zeros(n_pairs_w, 1);
for k = 1:n_pairs_w
    vk = W.v(:, :, k) * cal;
    per_pair(k) = extract.get_mean_of_selection(BW, vk);   % PIVlab area mean, downward +
end
score.n_pairs_wall = n_pairs_w;
score.wall_vectors_in_box = nnz(BW);
score.wall_piv_m_s = mean(per_pair, 'omitnan');
score.wall_piv_sd_m_s = std(per_pair, 'omitnan');
score.wall_piv_per_pair_m_s = per_pair';
score.wall_nan_frac_final = W.stats.nan_frac_final;

%% ---- settings used and output ------------------------------------------------
sj = jsondecode(fileread(fullfile(export_dir, variant, 'spout', 'settings.json')));
score.pivlab_version = pivlab_version(); score.matlab_version = sj.matlab_version;
score.pivlab_settings = rmfield(sj, intersect(fieldnames(sj), {'export_dir', 'view', 'out_dir', 'n_frames', 'n_pairs', 'date'}));
score.scoring = struct('time_mean', 'mean over pairs, omitnan, per grid point', ...
    'map', '|mean vector|, interp2 linear onto the digitized Fig. 4a panel grid', ...
    'CC', 'corrcoef', 'SSI', 'ssim (MATLAB defaults, reference = experiment)', 'MSE', 'immse [m^2/s^2]', ...
    'wall', 'extract.get_mean_of_selection of +v in the box per pair, then mean / std over pairs');
fid = fopen(fullfile(export_dir, variant, 'score.json'), 'w');
fwrite(fid, jsonencode(score, 'PrettyPrint', true)); fclose(fid);
fprintf('%s %s: map mean %.4f m/s; vs exp CC %.3f SSI %.3f MSE %.4f (C26 model vs exp: %.3f %.3f %.4f); wall %.3e +- %.1e m/s (%d vectors in box, %d pairs)\n', ...
    meta.run, variant, score.map_mean_m_s, score.metrics_vs_exp.CC, score.metrics_vs_exp.SSI, score.metrics_vs_exp.MSE_m2_s2, ...
    score.metrics_c26_model_vs_exp.CC, score.metrics_c26_model_vs_exp.SSI, score.metrics_c26_model_vs_exp.MSE_m2_s2, ...
    score.wall_piv_m_s, score.wall_piv_sd_m_s, nnz(BW), n_pairs_w);
end

function m = metrics(a, ref)
% CC (corrcoef), SSI (ssim, MATLAB defaults; ref = experiment), MSE (immse)
r = corrcoef(a(:), ref(:));
m = struct('CC', r(1, 2), 'SSI', ssim(a, ref), 'MSE_m2_s2', immse(a, ref));
end

function ver = pivlab_version()
% PIVlab_GUI.m line "version = '3.14';" of the PIVlab on the path
ver = 'unknown';
try
    txt = fileread(which('PIVlab_GUI'));
    m = regexp(txt, 'version\s*=\s*''(\d+\.\d+[\.\d]*)''', 'tokens', 'once');
    if ~isempty(m), ver = m{1}; end
catch
end
end

function a = read_digitized(path)
% data/digitized CSV: '#' comment lines, then comma rows, top row first [m/s]
a = readmatrix(path, 'CommentStyle', '#', 'FileType', 'text');
end
