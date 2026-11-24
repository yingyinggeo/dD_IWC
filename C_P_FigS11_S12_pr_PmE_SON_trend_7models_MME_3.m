%% ========================================================================
%  SON precipitation (pr) and P − E trends: 7 models + MME, 4 rows x 2 cols
%  Two figures:
%     Fig 1: pr       SON linear trend (mm day^-1 decade^-1)
%     Fig 2: P − E    SON linear trend (mm day^-1 decade^-1),  P − E = pr − ev
%  - Shading: trend; green dots: p < 0.05
%  - Land is masked in light grey (1° topo nearest point > 0, same as the MFC / ω figures); coastline drawn
%  - Boxes (WIO / EIO) are NOT drawn; box-mean trends are still computed and printed
%  - Colour map: precip4_diff_19lev with the two central white levels removed;
%    colour-bar ticks sit on level edges
%  Data: C:\Users\yingy\Downloads\Rice_mac\WISOMIP_Other_Variabiles\Total.<model>_Monthly.nc
%        (128 x 64; variables pr and ev, units mm/day)
%  Calculation domain: 60°S–60°N, 0–180°E; plotting domain: 20°S–40°N, 30–120°E
%  Panel aspect ratio is true to geography (1° lon = 1° lat)
%  Period: 2003–2020 if every model has SON 2020, otherwise 2003–2019
%  MME: ensemble mean of the yearly SON fields, then the trend
%% ========================================================================
clear; clc; close all;

%% ===== Settings =====
data_dir = 'C:\Users\yingy\Downloads\Rice_mac\WISOMIP_Other_Variabiles';
out_dir  = 'C:\Users\yingy\Downloads\Rice_mac';
models   = {'CAM5', 'CAM6', 'ECHAM', 'GISS', 'GSM', 'LMDZ', 'MIROC'};   % swap in NICAM if needed
n_models = numel(models);

yr_start = 2003;  yr_try_end = 2019;  yr_fallback_end = 2019;
years_all = (yr_start:yr_try_end)';
nyr_all   = numel(years_all);
son = [9 10 11];

calc_lon = [0 180];   calc_lat = [-60 60];
plot_lon = [30 120];  plot_lat = [-20 40];

clim_pr   = [];          % colour limit ± (mm day^-1 decade^-1); empty → from 98th percentile of |trend|,
clim_pme  = [];          %   rounded so that each of the 19 levels has a "nice" width
tick_lbl_step = 2;       % colour bar: 1 = label every level edge; 2 = every other edge (symmetric)
sig_level = 0.05;
fs        = 10;
dot_sz    = 3;

boxes = { ...
    [50 70],  [-15 5], [1.00 0.50 0.00];   % WIO, orange
    [70 100], [-20 0], [0.00 0.45 0.90]};  % EIO, blue
box_lw = 1.2;
draw_boxes = false;     % do not draw WIO / EIO boxes; box means are still computed and printed

%% ===== Read pr and ev for each model → yearly SON means =====
PR_all  = [];  PME_all = [];
have    = false(n_models, 1);
ok_last = false(n_models, 1);

for im = 1:n_models
    f = fullfile(data_dir, sprintf('Total.%s_Monthly.nc', models{im}));
    fprintf('\n===== %s =====\n', models{im});
    if ~exist(f, 'file'); warning('File not found: %s', f); continue; end

    try
        lon = double(ncread(f, 'lon'));
        lat = double(ncread(f, 'lat'));
        [yr, mo] = model_time(f);
        fprintf('Time: %d-%02d to %d-%02d (%d months)\n', yr(1), mo(1), yr(end), mo(end), numel(yr));
        P = double(ncread(f, 'pr'));   P(P <= -998) = NaN;     % lon x lat x time
        E = double(ncread(f, 'ev'));   E(E <= -998) = NaN;
    catch ME
        warning('%s read failed: %s', models{im}, ME.message);
        continue
    end

    % --- sanity checks: units and evaporation sign ---
    mP = mean(P(:), 'omitnan');  mE = mean(E(:), 'omitnan');
    if abs(mP) < 1e-2            % looks like kg m^-2 s^-1
        P = P * 86400;  fprintf('  pr looks like kg m-2 s-1 → converted to mm/day\n');
    end
    if abs(mE) < 1e-2
        E = E * 86400;  fprintf('  ev looks like kg m-2 s-1 → converted to mm/day\n');
    end
    if mean(E(:), 'omitnan') < 0  % evaporation should be positive upward
        E = -E;  fprintf('  ev mean < 0 → sign flipped (evaporation taken as positive)\n');
    end
    fprintf('Mean pr = %.2f, mean ev = %.2f mm/day\n', mean(P(:),'omitnan'), mean(E(:),'omitnan'));

    [lonS, ix] = sort(mod(lon, 360));
    [latS, iy] = sort(lat);
    mlon = lonS >= calc_lon(1) & lonS <= calc_lon(2);
    mlat = latS >= calc_lat(1) & latS <= calc_lat(2);
    lon_c = lonS(mlon);  lat_c = latS(mlat);
    P = permute(P(ix(mlon), iy(mlat), :), [2 1 3]);    % lat x lon x t
    E = permute(E(ix(mlon), iy(mlat), :), [2 1 3]);

    if isempty(PR_all)
        lon_g = lon_c;  lat_g = lat_c;
        PR_all  = NaN(numel(lat_g), numel(lon_g), nyr_all, n_models);
        PME_all = PR_all;
    elseif numel(lon_c) ~= numel(lon_g) || numel(lat_c) ~= numel(lat_g) || ...
           max(abs(lon_c - lon_g)) > 1e-3 || max(abs(lat_c - lat_g)) > 1e-3
        error('%s grid differs from the other models', models{im});
    end

    for k = 1:nyr_all
        kk = ismember(mo, son) & yr == years_all(k);
        if sum(kk) ~= 3; continue; end
        Ap = mean(P(:,:,kk), 3);
        Ae = mean(E(:,:,kk), 3);
        if all(isnan(Ap(:))); continue; end
        PR_all(:,:,k,im)  = Ap;
        PME_all(:,:,k,im) = Ap - Ae;
    end
    ok_last(im) = any(~isnan(reshape(PR_all(:,:,end,im), [], 1)));
    have(im) = true;
    fprintf('%s done; SON %d %s\n', models{im}, yr_try_end, ternary(ok_last(im), 'available', 'missing'));
    clear P E
end
if ~any(have); error('No usable model data'); end

%% ===== Period =====
if all(ok_last(have))
    yr_end = yr_try_end;
else
    yr_end = yr_fallback_end;
    fprintf('\nModels missing SON %d: %s → using %d–%d for all\n', yr_try_end, ...
        strjoin(models(have & ~ok_last), ', '), yr_start, yr_end);
end
reg_years = (yr_start:yr_end)';
keep = ismember(years_all, reg_years);
PR_all  = PR_all(:,:,keep,:);
PME_all = PME_all(:,:,keep,:);
fprintf('Period: %d–%d (%d yr), models = %d (%s)\n', yr_start, yr_end, numel(reg_years), ...
    sum(have), strjoin(models(have), ', '));

%% ===== Land mask (for box means) and coastline =====
load topo;
topo_lon = 0.5:359.5;  topo_lat = -89.5:89.5;
[TLON, TLAT] = meshgrid(topo_lon, topo_lat);
[GLON, GLAT] = meshgrid(lon_g, lat_g);
land_g = interp2(TLON, TLAT, topo, GLON, GLAT, 'nearest', 0) > 0;

topo_ext     = [topo, topo(:,1)];
topo_lon_ext = [topo_lon, 360.5];
buf = 5;
tlon_idx    = topo_lon_ext >= plot_lon(1)-buf & topo_lon_ext <= plot_lon(2)+buf;
tlat_idx    = topo_lat     >= plot_lat(1)-buf & topo_lat     <= plot_lat(2)+buf;
coast.lon   = topo_lon_ext(tlon_idx);
coast.lat   = topo_lat(tlat_idx);
coast.elev  = topo_ext(tlat_idx, tlon_idx);

%% ===== Trends and figures for both variables =====
load Precip4_diff_19lev.mat;               % variable name inside the file: precip4_diff_19lev
cmap_wet = double(precip4_diff_19lev);     % negative (drier) → positive (wetter); use flipud() if reversed
if max(cmap_wet(:)) > 1; cmap_wet = cmap_wet / 255; end
% remove the two white levels closest to the centre of the colour map
white_thr = 0.95;                          % RGB all ≥ this → treated as white
n0      = size(cmap_wet, 1);
i_white = find(all(cmap_wet >= white_thr, 2));
[~, ord] = sort(abs(i_white - (n0 + 1) / 2));
i_rm    = i_white(ord(1:min(2, numel(i_white))));
cmap_wet(i_rm, :) = [];
fprintf('Colour map: %d levels, removed %d central white level(s) (rows %s) → %d levels\n', ...
    n0, numel(i_rm), mat2str(sort(i_rm)'), size(cmap_wet, 1));
if numel(i_rm) ~= 2; warning('Did not find exactly 2 central white levels; check white_thr'); end

V(1).name = 'pr';   V(1).data = PR_all;   V(1).clim = clim_pr;
V(1).label = 'SON precipitation trend';   V(1).tag = 'pr';
V(2).name = 'P−E';  V(2).data = PME_all;  V(2).clim = clim_pme;
V(2).label = 'SON P − E trend';           V(2).tag = 'PmE';

panel_names = [models, {'MME'}];
n_panels    = n_models + 1;

for iv = 1:numel(V)
    R = struct([]);
    for ip = 1:n_panels
        if ip <= n_models
            if ~have(ip); R(ip).ok = false; continue; end
            M = V(iv).data(:,:,:,ip);
        else
            M = mean(V(iv).data(:,:,:,have), 4, 'omitnan');
        end
        [R(ip).tr, R(ip).p] = trend_decade(M, reg_years);
        R(ip).tr(land_g) = NaN;  R(ip).p(land_g) = NaN;   % land → NaN → shown in light grey
        R(ip).ok = true;

        % box-mean trends (ocean points, cos(lat) weighted)
        w = box_mean(M, lon_g, lat_g, boxes{1,1}, boxes{1,2}, land_g);
        e = box_mean(M, lon_g, lat_g, boxes{2,1}, boxes{2,2}, land_g);
        [R(ip).wio_tr, R(ip).wio_p] = ts_trend(w, reg_years, 10);
        [R(ip).eio_tr, R(ip).eio_p] = ts_trend(e, reg_years, 10);
    end

    fprintf('\n===== Box-mean %s trends (mm day^-1 decade^-1) =====\n', V(iv).name);
    fprintf('%-6s %10s %8s %10s %8s %12s\n', 'Panel', 'WIO', 'p', 'EIO', 'p', 'WIO-EIO');
    for ip = 1:n_panels
        if ~R(ip).ok; continue; end
        fprintf('%-6s %+10.4f %8.3f %+10.4f %8.3f %+12.4f\n', panel_names{ip}, ...
            R(ip).wio_tr, R(ip).wio_p, R(ip).eio_tr, R(ip).eio_p, R(ip).wio_tr - R(ip).eio_tr);
    end

    save(fullfile(out_dir, sprintf('%s_SON_trend_%d_%d_7models_MME.mat', V(iv).tag, yr_start, yr_end)), ...
        'R', 'panel_names', 'lon_g', 'lat_g', 'reg_years', 'land_g');

    fig = plot_panels(R, panel_names, lon_g, lat_g, plot_lon, plot_lat, V(iv).clim, ...
        cmap_wet, coast, boxes, box_lw, draw_boxes, sig_level, dot_sz, fs, tick_lbl_step, ...
        sprintf('%s, %d–%d (mm day^{-1} decade^{-1})', V(iv).label, yr_start, yr_end));

    base = fullfile(out_dir, sprintf('FigS_%s_SON_trend_%d_%d_7models_MME', V(iv).tag, yr_start, yr_end));
    exportgraphics(fig, [base '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
    exportgraphics(fig, [base '.tif'], 'Resolution', 600, 'BackgroundColor', 'white');
    fprintf('Saved: %s.png / .tif\n', base);
end

%% ===== Local functions =====
function fig = plot_panels(R, panel_names, lon_g, lat_g, plot_lon, plot_lat, clim_v, ...
    cmap, coast, boxes, box_lw, draw_boxes, sig_level, dot_sz, fs, tick_lbl_step, cb_label)
% 4 x 2 panels, true geographic aspect ratio, shared horizontal colour bar
n_panels = numel(panel_names);
mlon_p = lon_g >= plot_lon(1) & lon_g <= plot_lon(2);
mlat_p = lat_g >= plot_lat(1) & lat_g <= plot_lat(2);
lon_p  = lon_g(mlon_p);  lat_p = lat_g(mlat_p);

vals = [];
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    R(ip).trp = R(ip).tr(mlat_p, mlon_p);
    R(ip).pp  = R(ip).p (mlat_p, mlon_p);
    vals = [vals; abs(R(ip).trp(isfinite(R(ip).trp)))]; %#ok<AGROW>
end
% Colour levels: choose a "nice" level width so that the colour-bar ticks sit on level edges
n_lev = size(cmap, 1);
if isempty(clim_v)
    lev_w  = nice_ceil(2 * prctile(vals, 98) / n_lev);   % level width (1, 2, 5 × 10^k)
    clim_v = n_lev * lev_w / 2;
else
    lev_w  = 2 * clim_v / n_lev;
end
cb_edges = linspace(-clim_v, clim_v, n_lev + 1);
k_edge   = round(cb_edges / lev_w * 2);                 % half-level index of each edge
k_min    = min(abs(k_edge));
show_lb  = mod((abs(k_edge) - k_min) / 2, tick_lbl_step) == 0;
cb_lbls  = arrayfun(@(x) sprintf('%g', round(x, 6)), cb_edges, 'UniformOutput', false);
cb_lbls(~show_lb) = {''};
fprintf('Colour limit = ±%g (%d levels, level width %g)\n', clim_v, n_lev, lev_w);

cm2px    = 96/2.54;
fig_w_cm = 18;
pw_cm  = 0.32*18;
ph_cm  = pw_cm * diff(plot_lat) / diff(plot_lon);      % true lon/lat proportion
row_dy = ph_cm + 0.94;
nrow   = 4;
y_bot  = 2.6;
y_cm   = y_bot + (nrow-1:-1:0) * row_dy;
fig_h_cm = y_cm(1) + ph_cm + 0.7;

fig = figure('Color','w', ...
             'Position',[50 50 round(fig_w_cm*cm2px) round(fig_h_cm*cm2px)], ...
             'Renderer','opengl', 'GraphicsSmoothing','off', ...
             'PaperUnits','centimeters', 'PaperSize',[fig_w_cm, fig_h_cm], ...
             'PaperPositionMode','auto');
set(fig, 'DefaultAxesFontSize', fs, 'DefaultTextFontSize', fs);

pw = pw_cm / fig_w_cm;   ph = ph_cm / fig_h_cm;
xs = [0.07 0.46];
[LONs, LATs] = meshgrid(lon_p, lat_p);
silver_rgb = [0.85 0.85 0.85];

letters = 'abcdefgh';
ax_last = [];
for ip = 1:n_panels
    ir = ceil(ip / 2);  ic = mod(ip-1, 2) + 1;
    pos = [xs(ic), y_cm(ir)/fig_h_cm, pw, ph];
    ax  = axes('Position', pos);  hold(ax, 'on');

    if R(ip).ok
        dat = R(ip).trp;
        imagesc(ax, lon_p, lat_p, dat, 'AlphaData', ~isnan(dat));
        sg = R(ip).pp < sig_level & ~isnan(dat);
        plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize', dot_sz, ...
             'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerEdgeColor', 'none');
        ttl = sprintf('(%c) %s', letters(ip), panel_names{ip});
    else
        ttl = sprintf('(%c) %s (no data)', letters(ip), panel_names{ip});
    end

    set(ax, 'YDir', 'normal', 'Color', silver_rgb);
    colormap(ax, cmap);
    caxis(ax, [-clim_v clim_v]);
    contour(ax, coast.lon, coast.lat, coast.elev, [0 0], 'k-', 'LineWidth', 0.5);
    if draw_boxes
        for ib = 1:size(boxes, 1)
            bx = boxes{ib, 1};  by = boxes{ib, 2};
            plot(ax, bx([1 2 2 1 1]), by([1 1 2 2 1]), '-', 'Color', boxes{ib, 3}, 'LineWidth', box_lw);
        end
    end
    xlim(ax, plot_lon);  ylim(ax, plot_lat);
    daspect(ax, [1 1 1]);

    if ir == nrow
        set(ax, 'XTick', 30:30:120, 'XTickLabel', {'30°E','60°E','90°E','120°E'});
    else
        set(ax, 'XTick', 30:30:120, 'XTickLabel', {});
    end
    if ic == 1
        set(ax, 'YTick', -20:20:40, 'YTickLabel', {'20°S','0°','20°N','40°N'});
    else
        set(ax, 'YTick', -20:20:40, 'YTickLabel', {});
    end
    set(ax, 'FontSize', fs, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
    ax.Position = pos;

    annotation(fig, 'textbox', [pos(1), pos(2)+ph+0.003, pw, 0.015], ...
               'String', ttl, 'FontSize', fs, 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
               'EdgeColor', 'none', 'FitBoxToText', 'off', 'Interpreter', 'tex');
    ax.Position = pos;
    ax_last = ax;
end

cb = colorbar(ax_last, 'southoutside');
cb.Position = [xs(1)+0.1*pw, (y_bot - 1.3)/fig_h_cm, xs(2)+pw-xs(1)-0.2*pw, 0.3/fig_h_cm];
cb.Limits     = [-clim_v clim_v];
cb.Ticks      = cb_edges;          % ticks on every level edge
cb.TickLabels = cb_lbls;           % labels on a symmetric subset
cb.Label.String   = cb_label;
cb.FontSize       = fs;
cb.Label.FontSize = fs;
cb.TickDirection  = 'out';
ax_last.Position  = [xs(2), y_cm(end)/fig_h_cm, pw, ph];

set(findall(fig, 'Type', 'axes'), 'TickLength', [0.025 0.035]);
end

function [yr, mo] = model_time(file)
% "months since YYYY-M-D": rebuild as consecutive months (guards against bad/NaN values)
info  = ncinfo(file);
names = {info.Variables.Name};
if any(strcmp(names, 'time')); tn = 'time'; else; tn = 'time_counter'; end
t_raw = double(ncread(file, tn));
units = ncreadatt(file, tn, 'units');
tok = regexp(units, '(\w+)\s+since\s+(\d{4})-(\d{1,2})', 'tokens', 'once');
if isempty(tok) || ~contains(lower(tok{1}), 'month')
    error('Time units are not "months since": %s', units);
end
y0 = str2double(tok{2});  m0 = str2double(tok{3});
nt = numel(t_raw);
t1 = t_raw(find(isfinite(t_raw), 1));
if isempty(t1); t1 = 0; end
tm = round(t1) + (0:nt-1)';
bad = find(~isfinite(t_raw(:)) | abs(t_raw(:) - tm) > 0.5);
if ~isempty(bad)
    fprintf('  %d bad time values; rebuilt as consecutive months\n', numel(bad));
end
yr = y0 + floor((m0 - 1 + tm) / 12);
mo = mod(m0 - 1 + tm, 12) + 1;
end

function [slope, pval] = trend_decade(Y, years)
% grid-point OLS trend (per decade) + two-sided t test; NaN where any year is missing
[d1, d2, n] = size(Y);
x   = years(:) - mean(years);
Sxx = sum(x.^2);
Y2  = reshape(Y, d1*d2, n)';
ok  = ~any(isnan(Y2), 1);
yb  = mean(Y2, 1);
b   = (x' * (Y2 - yb)) / Sxx;
res = Y2 - yb - x * b;
s2  = sum(res.^2, 1) / (n - 2);
t   = b ./ sqrt(s2 / Sxx);
p   = 2 * (1 - tcdf(abs(t), n - 2));
b(~ok) = NaN;  p(~ok) = NaN;
slope = reshape(b * 10, d1, d2);
pval  = reshape(p, d1, d2);
end

function ts = box_mean(M, lon, lat, blon, blat, land)
mask = (lat(:) >= blat(1) & lat(:) <= blat(2)) & (lon(:)' >= blon(1) & lon(:)' <= blon(2)) & ~land;
W = repmat(cosd(lat(:)), 1, numel(lon)) .* mask;
ts = NaN(size(M, 3), 1);
for k = 1:size(M, 3)
    A = M(:,:,k);
    ok = W > 0 & isfinite(A);
    if any(ok(:)); ts(k) = sum(A(ok) .* W(ok)) / sum(W(ok)); end
end
end

function [b, p] = ts_trend(y, years, per)
x = years(:);  y = y(:);
ok = isfinite(y);  x = x(ok);  y = y(ok);  n = numel(y);
if n < 4; b = NaN; p = NaN; return; end
x = x - mean(x);
Sxx = sum(x.^2);
bb  = sum(x .* (y - mean(y))) / Sxx;
res = y - mean(y) - bb * x;
t   = bb / sqrt(sum(res.^2) / (n - 2) / Sxx);
p   = 2 * (1 - tcdf(abs(t), n - 2));
b   = bb * per;
end

function y = nice_ceil(x)
if ~isfinite(x) || x <= 0; y = 1; return; end
e = 10^floor(log10(x));  m = x / e;
if m <= 1; y = e; elseif m <= 2; y = 2*e; elseif m <= 5; y = 5*e; else; y = 10*e; end
end

function out = ternary(c, a, b)
if c; out = a; else; out = b; end
end
