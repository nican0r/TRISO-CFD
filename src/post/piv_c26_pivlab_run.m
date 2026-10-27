function out = piv_c26_pivlab_run(export_dir, view, out_dir, opts)
%PIV_C26_PIVLAB_RUN  PIVlab (GUI defaults) over the exported C26 frames of one view.
%
%   piv_c26_pivlab_run(export_dir, view, out_dir)
%   piv_c26_pivlab_run(..., binarize=1)
%
% Re-runs the project's synthetic PIV (src/post/piv_c26.py) with PIVlab, the
% program C26 (Conry et al. 2026, NED 456:115009) used, on the frames written
% by src/post/piv_c26_export_frames.py:
%
%   <export_dir>/<view>/frame_NNNN.tif   16-bit side-view frames, 2 ms apart
%
% Settings taken from the papers (C26 2.1, 2.2.2, 3.1; CES25 = Conry et al.
% 2025, CES 306:121197, 2.3, which C26 cites for the full PIV process):
%   - consecutive frames, time-resolved sequencing 1-2, 2-3, ... (CES25
%     "frame-by-frame sequence"; C26 "on a frame-by-frame basis");
%   - PIVlab's multipass FFT window deformation (CES25 "adaptive multipass");
%   - no calibration here: x, y in px and u, v in px/frame, converted with the
%     export's pixel size and frame interval by src.post.piv_c26_pivlab_compare;
%   - C26 Fig. 4d shows the simulated frames as PIVlab saw them: white
%     parcel-sized disks on black (the CES25 binarization applied to the
%     render); that is binarize=1, binarize_thresh=0, binarize_dilate_px=0,
%     image_class='uint8' (variant pivlab_c26img in piv_c26_pivlab_batch.sh);
%   - binarize=1 applies CES25's pre-processing ("binarization, background
%     thresholding and dilation", CES25 2.3) with a threshold and dilation
%     that CES25 does not give (opts.binarize_thresh, binarize_dilate_px).
% Every other setting is PIVlab's default (PIVlab_settings_default.mat of
% PIVlab 3.14 and the GUI code in +piv/DCC_and_DFT_analyze_all.m,
% +piv/CorrQuality.m, +validate/filtervectors.m):
%   pre-processing  : auto contrast stretch per image (stretchlim), CLAHE 64 px,
%                     no high-pass, no intensity capping, no denoise;
%   PIV             : FFT window deformation, pass 1 64 px step 32, pass 2
%                     32 px step 16, passes 3-4 off, Gauss 2x3-point sub-pixel,
%                     correlation robustness "standard" (linear image
%                     deformation, no repeated correlation, circular
%                     correlation), no repeated last pass, no auto-mask;
%   validation      : no velocity limits, std filter 8 sd, normalized median
%                     filter threshold 3 (eps 0.1 px, Westerweel & Scarano 2005),
%                     second-peak substitution for rejected vectors (GUI
%                     default), "interpolate missing data" on (inpaint_nans 4);
%   exclusions      : no ROI, no mask (the exported views are unmasked
%                     rectangles above the annulus surface / around the wall box).
%
% Output (out_dir): pair_NNNN.txt (x [px], y [px], u [px/frame], v [px/frame],
% typevector; y down, 1-based, PIVlab's convention), fields.mat (all fields,
% raw and validated), settings.json (every setting used).  Returns a struct
% with the settings and the fields.
%
% Run headless, e.g.
%   matlab -batch "addpath('src/post'); piv_c26_pivlab_run('results/pivlab_export/c02_apex_v232','wall','results/pivlab_export/c02_apex_v232/pivlab_default/wall')"

arguments
    export_dir (1,1) string
    view (1,1) string
    out_dir (1,1) string
    opts.binarize (1,1) double = 0
    opts.binarize_thresh (1,1) double = 0.9    % intensity 0-1 kept as "dot" (CES25 gives none; notes 2026-10-25 A/B)
    opts.binarize_dilate_px (1,1) double = 1   % dilation radius [px] (CES25 gives none)
    opts.image_class (1,1) string = "uint16"  % 'uint8': 8-bit frames, as a screenshot / video frame (C26 Fig. 4d)
    opts.clahe (1,1) double = 1               % PIVlab default on; 0 only for the diagnostic A/B
    opts.passes (1,1) double = 2              % PIVlab default 2 (64 -> 32 px); 3 adds int3
    opts.int3 (1,1) double = 32               % pass 3 window [px] (used when passes >= 3)
    opts.int4 (1,1) double = 32               % pass 4 window [px] (used when passes >= 4)
    opts.max_pairs (1,1) double = Inf
    opts.pivlab_root (1,1) string = fullfile(userpath, 'PIVlab')
    opts.quiet (1,1) logical = false
end

addpath(char(opts.pivlab_root));
if ~usejava('desktop')
    % headless: shadow the desktop-only COMMANDWINDOW that piv_FFTmulti calls
    % after a bad-data warning (see matlab_shims/commandwindow.m)
    addpath(fullfile(fileparts(mfilename('fullpath')), 'matlab_shims'));
end

%% ---- settings (PIVlab defaults unless noted) --------------------------------
s = struct();
s.source = 'PIVlab 3.14 PIVlab_settings_default.mat + GUI code; C26 / CES25 where noted';
s.sequencing = 'time-resolved 1-2, 2-3, ... (C26 3.1, CES25 2.3)';
% pre-processing
s.pre.binarize = opts.binarize;                 % CES25 2.3 (threshold / dilation assumed)
s.pre.binarize_thresh = opts.binarize_thresh;
s.pre.binarize_dilate_px = opts.binarize_dilate_px;
s.pre.image_class = char(opts.image_class);
s.pre.autolimit = 1;                            % "Auto contrast stretch" (stretchlim per image)
s.pre.clahe = opts.clahe;  s.pre.clahesize = 64;
s.pre.highp = 0;  s.pre.highpsize = 15;
s.pre.intenscap = 0;
s.pre.wienerwurst = 0; s.pre.wienerwurstsize = 3;
s.pre.roirect = [];                             % no ROI
% PIV
s.piv.algorithm = 'FFT window deformation (piv_FFTmulti)';
s.piv.interrogationarea = 64; s.piv.step = 32;  % pass 1
s.piv.passes = opts.passes; s.piv.int2 = 32;    % pass 2 (step = int2/2 inside piv_FFTmulti)
s.piv.int3 = opts.int3; s.piv.int4 = opts.int4;  % passes 3-4 (default off)
s.piv.subpixfinder = 1;                         % Gauss 2x3-point
s.piv.corr_quality = 'standard';                % CorrQuality 1 ->
s.piv.imdeform = '*linear'; s.piv.repeat = 0; s.piv.do_linear_correlation = 0;
s.piv.mask_auto = 0;
s.piv.repeat_last_pass = 0; s.piv.delta_diff_min = 0.025;
s.piv.limit_peak_search_area = 1;
s.piv.compute_uncertainty = 0;
s.piv.mask = 'none';
% validation
s.post.calu = 1; s.post.calv = 1;               % uncalibrated (px, px/frame)
s.post.valid_vel = [];                          % velocity limits off
s.post.do_stdev_check = 1; s.post.stdthresh = 8;
s.post.do_local_median = 1; s.post.neigh_thresh = 3; s.post.eps_px = 0.1;
s.post.second_peak_substitution = 1;            % filtervectors.m default
s.post.interpol_missing = 1;                    % inpaint_nans(..., 4)
s.post.image_based_filters = 0; s.post.correlation_filter = 0; s.post.notch_filter = 0;

%% ---- frames ------------------------------------------------------------------
files = dir(fullfile(export_dir, view, 'frame_*.tif'));
if isempty(files)
    error('piv_c26_pivlab_run:noFrames', 'no frame_*.tif in %s', fullfile(export_dir, view));
end
names = sort({files.name});
n_pairs = min(numel(names) - 1, opts.max_pairs);
if n_pairs < 1
    error('piv_c26_pivlab_run:tooFew', 'need at least two frames');
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
s.export_dir = char(export_dir); s.view = char(view); s.out_dir = char(out_dir);
s.n_frames = numel(names); s.n_pairs = n_pairs;
s.pivlab_version = pivlab_version(opts.pivlab_root);
s.matlab_version = version;
s.date = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

fid = fopen(fullfile(out_dir, 'settings.json'), 'w');
fwrite(fid, jsonencode(s, 'PrettyPrint', true)); fclose(fid);

%% ---- main loop ----------------------------------------------------------------
x = []; y = []; u_raw = []; v_raw = []; u = []; v = []; tv_raw = []; tv = []; corr = [];
n_rejected = zeros(n_pairs, 1); n_second_peak = zeros(n_pairs, 1); n_total = zeros(n_pairs, 1);
t_start = tic;
for k = 1:n_pairs
    A = load_frame(fullfile(export_dir, view, names{k}), s.pre);
    B = load_frame(fullfile(export_dir, view, names{k + 1}), s.pre);
    % Auto contrast stretch, per image, as DCC_and_DFT_analyze_all.m does
    stA = stretchlim(A); stB = stretchlim(B);
    A = preproc.PIVlab_preproc(in=A, roirect=s.pre.roirect, clahe=s.pre.clahe, clahesize=s.pre.clahesize, ...
        highp=s.pre.highp, highpsize=s.pre.highpsize, intenscap=s.pre.intenscap, ...
        wienerwurst=s.pre.wienerwurst, wienerwurstsize=s.pre.wienerwurstsize, ...
        minintens=stA(1), maxintens=stA(2));
    B = preproc.PIVlab_preproc(in=B, roirect=s.pre.roirect, clahe=s.pre.clahe, clahesize=s.pre.clahesize, ...
        highp=s.pre.highp, highpsize=s.pre.highpsize, intenscap=s.pre.intenscap, ...
        wienerwurst=s.pre.wienerwurst, wienerwurstsize=s.pre.wienerwurstsize, ...
        minintens=stB(1), maxintens=stB(2));
    [xk, yk, uk, vk, tvk, cmk, ~, ~, u2k, v2k] = piv.piv_FFTmulti( ...
        image1=A, image2=B, interrogationarea=s.piv.interrogationarea, step=s.piv.step, ...
        subpixfinder=s.piv.subpixfinder, mask_inpt=[], roi_inpt=s.pre.roirect, ...
        passes=s.piv.passes, int2=s.piv.int2, int3=s.piv.int3, int4=s.piv.int4, ...
        imdeform=s.piv.imdeform, repeat=s.piv.repeat, mask_auto=s.piv.mask_auto, ...
        do_linear_correlation=s.piv.do_linear_correlation, do_correlation_matrices=0, ...
        repeat_last_pass=s.piv.repeat_last_pass, delta_diff_min=s.piv.delta_diff_min, ...
        limit_peak_search_area=s.piv.limit_peak_search_area, compute_uncertainty=s.piv.compute_uncertainty);
    [uf, vf, tvf, nrej, nsub] = validate_like_gui(uk, vk, tvk, u2k, v2k, s.post);
    if k == 1
        [nr, nc] = size(uk);
        x = zeros(nr, nc, n_pairs); y = x; u_raw = x; v_raw = x; u = x; v = x; tv_raw = x; tv = x; corr = x;
    end
    x(:,:,k) = xk; y(:,:,k) = yk; u_raw(:,:,k) = uk; v_raw(:,:,k) = vk; tv_raw(:,:,k) = tvk;
    u(:,:,k) = uf; v(:,:,k) = vf; tv(:,:,k) = tvf; corr(:,:,k) = cmk;
    n_rejected(k) = nrej; n_second_peak(k) = nsub; n_total(k) = numel(uk);
    write_ascii(fullfile(out_dir, sprintf('pair_%04d.txt', k - 1)), xk, yk, uf, vf, tvf, names{k}, names{k + 1});
    if ~opts.quiet && (mod(k, 25) == 0 || k == n_pairs)
        fprintf('  %s %s: pair %d / %d (%.0f s)\n', export_dir, view, k, n_pairs, toc(t_start));
    end
end

stats = struct('n_pairs', n_pairs, 'grid', size(u, [1 2]), ...
    'rejected_frac', sum(n_rejected) / sum(n_total), 'second_peak_frac', sum(n_second_peak) / sum(n_total), ...
    'nan_frac_raw', mean(isnan(u_raw(:))), 'nan_frac_final', mean(isnan(u(:))), 'seconds', toc(t_start));
save(fullfile(out_dir, 'fields.mat'), 'x', 'y', 'u_raw', 'v_raw', 'u', 'v', 'tv_raw', 'tv', 'corr', 's', 'stats', '-v7.3');
fid = fopen(fullfile(out_dir, 'stats.json'), 'w'); fwrite(fid, jsonencode(stats, 'PrettyPrint', true)); fclose(fid);
out = struct('settings', s, 'stats', stats, 'x', x, 'y', y, 'u', u, 'v', v, 'typevector', tv);
if ~opts.quiet
    fprintf('done: %d pairs, grid %dx%d, rejected %.1f %%, second peak %.1f %%, NaN raw %.1f %% -> final %.1f %%, %.0f s\n', ...
        n_pairs, stats.grid(1), stats.grid(2), 100 * stats.rejected_frac, 100 * stats.second_peak_frac, ...
        100 * stats.nan_frac_raw, 100 * stats.nan_frac_final, stats.seconds);
end
end

%% ---- helpers -------------------------------------------------------------------
function img = load_frame(path, pre)
img = imread(path);
if size(img, 3) > 1, img = rgb2gray(img); end
if pre.binarize
    % CES25 2.3: binarize, threshold the background, dilate.  Keeps bright
    % "dots" (parcel centres) on black, as CES25's glints.
    % thresh 0 keeps every lit pixel: the full parcel silhouette, the white
    % disks of C26 Fig. 4d; thresh 0.9 keeps the bright centres (dots).
    d = im2double(img) > pre.binarize_thresh;
    if pre.binarize_dilate_px > 0
        d = imdilate(d, strel('disk', pre.binarize_dilate_px, 0));
    end
    img = uint16(d) * 65535;
end
if strcmp(pre.image_class, 'uint8')
    img = im2uint8(img);
end
end

function [u, v, tv, n_rejected, n_substituted] = validate_like_gui(u0, v0, tv0, u2, v2, p)
% +validate/filtervectors.m with the default GUI settings: PIVlab_postproc
% (velocity limits off, std filter, normalized median filter), second-peak
% substitution, mask restore, interpolate missing data.
[u, v] = postproc.PIVlab_postproc(u=u0, v=v0, calu=p.calu, calv=p.calv, valid_vel=p.valid_vel, ...
    do_stdev_check=p.do_stdev_check, stdthresh=p.stdthresh, ...
    do_local_median=p.do_local_median, neigh_thresh=p.neigh_thresh);
tv = tv0;
tv(isnan(u)) = 2; tv(isnan(v)) = 2;
n_rejected = nnz(tv == 2 & tv0 ~= 0);
n_substituted = 0;
if p.second_peak_substitution && ~isempty(u2) && ~isempty(v2)
    u2 = single(u2); v2 = single(v2);
    candidates = (tv == 2) & ~isnan(u2) & ~isnan(v2) & (tv0 ~= 0);
    if any(candidates(:))
        u_sc = misc.inpaint_nans(u, 4); v_sc = misc.inpaint_nans(v, 4);
        accepted = candidates;
        if p.do_local_median == 1
            b = 2;
            MedU = medfilt2(u_sc, [2*b+1, 2*b+1], 'symmetric'); MedV = medfilt2(v_sc, [2*b+1, 2*b+1], 'symmetric');
            ResU = medfilt2(abs(u_sc - MedU), [2*b+1, 2*b+1], 'symmetric');
            ResV = medfilt2(abs(v_sc - MedV), [2*b+1, 2*b+1], 'symmetric');
            NU = abs((u2 - MedU) ./ (ResU + p.eps_px)); NV = abs((v2 - MedV) ./ (ResV + p.eps_px));
            accepted = accepted & (sqrt(NU.^2 + NV.^2) <= p.neigh_thresh);
        end
        if p.do_stdev_check == 1
            mu = mean(u(:), 'omitnan'); mv = mean(v(:), 'omitnan');
            su = std(u(:), 'omitnan'); sv = std(v(:), 'omitnan');
            accepted = accepted & (u2 >= mu - p.stdthresh*su) & (u2 <= mu + p.stdthresh*su) ...
                                & (v2 >= mv - p.stdthresh*sv) & (v2 <= mv + p.stdthresh*sv);
        end
        u(accepted) = u2(accepted); v(accepted) = v2(accepted);
        tv(accepted) = 3;
        n_substituted = nnz(accepted);
    end
end
tv(tv0 == 0) = 0;
if p.interpol_missing
    try, u = misc.inpaint_nans(u, 4); catch, end
    try, v = misc.inpaint_nans(v, 4); catch, end
end
u = double(u); v = double(v);
end

function write_ascii(path, x, y, u, v, tv, nameA, nameB)
% PIVlab-style ASCII: one row per vector, y down, 1-based pixel coordinates,
% u, v in px/frame (uncalibrated).  Header lines are non-numeric and are
% skipped by src.post.piv_c26_pivlab_compare.
fid = fopen(path, 'w');
fprintf(fid, '%% PIVlab batch (piv_c26_pivlab_run.m), frames %s -> %s\n', nameA, nameB);
fprintf(fid, 'x [px],y [px],u [px/frame],v [px/frame],typevector\n');
fprintf(fid, '%.4f,%.4f,%.6f,%.6f,%d\n', [x(:)'; y(:)'; u(:)'; v(:)'; double(tv(:))']);
fclose(fid);
end

function ver = pivlab_version(root)
ver = 'unknown';
try
    txt = fileread(fullfile(root, 'PIVlab_GUI.m'));
    m = regexp(txt, 'version\s*=\s*''(\d+\.\d+[\.\d]*)''', 'tokens', 'once');   % PIVlab_GUI.m: version = '3.14';
    if ~isempty(m), ver = m{1}; end
catch
end
end
