%% ========================================================================
%  WisoMIP：比湿 q 的 SON 趋势（或气候态），7 个模式 + MME，4 行 2 列
%   两张图：
%     图 1：~850 hPa（取最接近 850 hPa 的模式气压层）
%     图 2：400–600 hPa 层平均（层内所有模式气压层取平均）
%   - 逐年 SON 平均 → OLS 趋势（g kg^-1 decade^-1），双侧 t 检验 p<0.05 → 绿点
%   - MME：先对各模式逐年 SON 场取集合平均，再算趋势
%   - plot_mode = 'clim' 时画 SON 气候态（g kg^-1），不打点
%   imagesc 填色；陆地屏蔽为浅灰；不画 WIO / EIO 框（框内平均仍打印）；
%   20°S–40°N，30–120°E，真实长宽比（1° 经度 = 1° 纬度）
%   色标：precip4_diff_19lev 去掉中间两个白色色阶（同 P / P−E 趋势图），刻度放在色阶分界上
%   数据：C:\Users\yingy\Downloads\Rice_mac\WISOMIP_Other_Variabiles\Total.<model>_Monthly.nc
%         (128 x 64 网格，17 个气压层 p (mb)；q 单位 kg/kg)
%   时段：优先 2003–2020；若有模式缺 2020 年 SON，则统一用 2003–2019
%% ========================================================================
clear; clc; close all;

%% ===== 基本设置 =====
data_dir = 'C:\Users\yingy\Downloads\Rice_mac\WISOMIP_Other_Variabiles';
out_dir  = 'C:\Users\yingy\Downloads\Rice_mac';
models   = {'CAM5', 'CAM6', 'ECHAM', 'GISS', 'GSM', 'LMDZ', 'MIROC'};   % 可换成 NICAM
n_models = numel(models);

plot_mode = 'trend';                     % 'trend' = SON 趋势；'clim' = SON 气候态
p_low     = 850;                         % 低层：取最接近的模式气压层
p_layer   = [400 600];                   % 中层：层内所有模式气压层取平均
yr_start = 2003;  yr_try_end = 2019;  yr_fallback_end = 2019;
years_all = (yr_start:yr_try_end)';
nyr_all   = numel(years_all);
son = [9 10 11];

plot_lon = [30 120];  plot_lat = [-20 40];

% 色标：[] = 自动（趋势：按 |值| 98% 分位选整齐的色阶宽度；气候态：2%–98% 分位）
clim_low = [];                           % 例如 0.5 → ±0.5 g kg^-1 decade^-1
clim_mid = [];                           % 例如 0.2
tick_lbl_step = 2;                       % 色标：1 = 每个色阶分界都标数字；2 = 隔一个标一个（对称）

% colormap：precip4_diff_19lev，去掉最中间的两个白色色阶
load Precip4_diff_19lev.mat;             % 文件内变量名：precip4_diff_19lev
cmap_div = double(precip4_diff_19lev);   % 负值（变干）→ 正值（变湿）；颜色反了就用 flipud()
if max(cmap_div(:)) > 1; cmap_div = cmap_div / 255; end
white_thr = 0.95;                        % RGB 三个分量都 ≥ 此值视为白色
n0      = size(cmap_div, 1);
i_white = find(all(cmap_div >= white_thr, 2));
[~, ord] = sort(abs(i_white - (n0 + 1) / 2));
i_rm    = i_white(ord(1:min(2, numel(i_white))));
cmap_div(i_rm, :) = [];
fprintf('色标：原 %d 个色阶，去掉中间白色 %d 个（第 %s 行），剩余 %d 个\n', ...
    n0, numel(i_rm), mat2str(sort(i_rm)'), size(cmap_div, 1));
if numel(i_rm) ~= 2; warning('中间白色色阶不是 2 个，请检查 white_thr'); end
sig_level = 0.05;
fs        = 10;
dot_sz    = 3;

boxes = { ...
    [50 70],  [-15 5], [1.00 0.50 0.00];   % WIO，橙色
    [70 100], [-20 0], [0.00 0.45 0.90]};  % EIO，蓝色
box_lw     = 1.2;
draw_boxes = false;                      % 不画 WIO / EIO 框；框内平均仍计算并打印

%% ===== 逐模式：读 q → 逐年 SON 平均（单位 g/kg）=====
Q_low = [];  Q_mid = [];
have    = false(n_models, 1);
ok_last = false(n_models, 1);
lev_txt_low = cell(n_models, 1);

for im = 1:n_models
    f = fullfile(data_dir, sprintf('Total.%s_Monthly.nc', models{im}));
    fprintf('\n===== %s =====\n', models{im});
    if ~exist(f, 'file'); warning('找不到 %s', f); continue; end

    try
        lon = double(ncread(f, 'lon'));
        lat = double(ncread(f, 'lat'));
        p   = double(ncread(f, 'p'));
        [~, i850] = min(abs(p - p_low));
        il = find(p >= p_layer(1) & p <= p_layer(2));
        if isempty(il); error('没有 %d–%d hPa 的气压层', p_layer(1), p_layer(2)); end
        fprintf('低层 p = %g hPa；中层 p = %s hPa\n', p(i850), mat2str(p(il)'));
        lev_txt_low{im} = sprintf('%g', p(i850));

        [yr, mo] = model_time(f);
        fprintf('时间：%d-%02d 到 %d-%02d（%d 个月）\n', yr(1), mo(1), yr(end), mo(end), numel(yr));

        qL = read_layer(f, 'q', i850) * 1000;     % lon x lat x time，g/kg
        qM = read_layer(f, 'q', il)   * 1000;
    catch ME
        warning('%s 读取失败：%s', models{im}, ME.message);
        continue
    end
    fprintf('q850 范围 [%.2f, %.2f]，q400–600 范围 [%.2f, %.2f] g/kg\n', ...
        min(qL(:)), max(qL(:)), min(qM(:)), max(qM(:)));

    [lonS, ix] = sort(mod(lon, 360));
    [latS, iy] = sort(lat);
    qL = permute(qL(ix, iy, :), [2 1 3]);       % lat(南→北) x lon x t
    qM = permute(qM(ix, iy, :), [2 1 3]);

    if isempty(Q_low)
        lon_g = lonS;  lat_g = latS;
        Q_low = NaN(numel(lat_g), numel(lon_g), nyr_all, n_models);
        Q_mid = Q_low;
    elseif numel(lonS) ~= numel(lon_g) || numel(latS) ~= numel(lat_g) || ...
           max(abs(lonS - lon_g)) > 1e-3 || max(abs(latS - lat_g)) > 1e-3
        error('%s 的网格与其他模式不一致', models{im});
    end

    for k = 1:nyr_all
        kk = ismember(mo, son) & yr == years_all(k);
        if sum(kk) ~= 3; continue; end
        Q_low(:,:,k,im) = mean(qL(:,:,kk), 3);
        Q_mid(:,:,k,im) = mean(qM(:,:,kk), 3);
    end
    ok_last(im) = any(~isnan(reshape(Q_mid(:,:,end,im), [], 1)));
    have(im) = true;
    clear qL qM
end
if ~any(have); error('没有可用的模式数据'); end

%% ===== 确定时段 =====
if all(ok_last(have)); yr_end = yr_try_end; else; yr_end = yr_fallback_end; end
reg_years = (yr_start:yr_end)';
kk = ismember(years_all, reg_years);
Q_low = Q_low(:,:,kk,:);  Q_mid = Q_mid(:,:,kk,:);
fprintf('\n时段：%d–%d（%d 年），模式数 = %d (%s)\n', yr_start, yr_end, numel(reg_years), ...
    sum(have), strjoin(models(have), ', '));

%% ===== 陆地 mask 与海岸线 =====
load topo;
topo_lon = 0.5:359.5;  topo_lat = -89.5:89.5;
[TLON, TLAT] = meshgrid(topo_lon, topo_lat);
[GLON, GLAT] = meshgrid(lon_g, lat_g);
land_g = interp2(TLON, TLAT, topo, GLON, GLAT, 'nearest', 0) > 0;

topo_ext     = [topo, topo(:,1)];
topo_lon_ext = [topo_lon, 360.5];
buf = 5;
tlon_idx = topo_lon_ext >= plot_lon(1)-buf & topo_lon_ext <= plot_lon(2)+buf;
tlat_idx = topo_lat     >= plot_lat(1)-buf & topo_lat     <= plot_lat(2)+buf;
TOPO.lon  = topo_lon_ext(tlon_idx);
TOPO.lat  = topo_lat(tlat_idx);
TOPO.elev = topo_ext(tlat_idx, tlon_idx);

%% ===== 公共画图参数 =====
S.models = models;  S.have = have;  S.n_models = n_models;
S.reg_years = reg_years;  S.yr_start = yr_start;  S.yr_end = yr_end;
S.lon_g = lon_g;  S.lat_g = lat_g;  S.land_g = land_g;
S.plot_lon = plot_lon;  S.plot_lat = plot_lat;
S.boxes = boxes;  S.box_lw = box_lw;  S.draw_boxes = draw_boxes;
S.sig_level = sig_level;  S.fs = fs;  S.dot_sz = dot_sz;
S.mode = plot_mode;  S.out_dir = out_dir;  S.topo = TOPO;
S.cmap = cmap_div;  S.tick_lbl_step = tick_lbl_step;

%% ===== 图 1：~850 hPa =====
S.lev_name = sprintf('%d hPa', p_low);
S.tag      = sprintf('q%d', p_low);
S.clim     = clim_low;
make_fig(Q_low, S);

%% ===== 图 2：400–600 hPa =====
S.lev_name = sprintf('%d–%d hPa', p_layer(1), p_layer(2));
S.tag      = sprintf('q%d_%d', p_layer(1), p_layer(2));
S.clim     = clim_mid;
make_fig(Q_mid, S);

%% ===== 子函数 =====
function make_fig(Q, S)
% Q: lat x lon x year x model（g/kg）
n_panels    = S.n_models + 1;
panel_names = [S.models, {'MME'}];
is_tr = strcmpi(S.mode, 'trend');

% ---- 计算每个面板 ----
R = struct([]);
for ip = 1:n_panels
    if ip <= S.n_models
        if ~S.have(ip); R(ip).ok = false; continue; end
        M = Q(:,:,:,ip);
    else
        M = mean(Q(:,:,:,S.have), 4, 'omitnan');
    end
    if is_tr
        [v, pv] = trend_decade(M, S.reg_years);
    else
        v = mean(M, 3, 'omitnan');  pv = NaN(size(v));
    end
    v(S.land_g) = NaN;  pv(S.land_g) = NaN;
    R(ip).v = v;  R(ip).p = pv;  R(ip).ok = true;

    w = box_mean(M, S.lon_g, S.lat_g, S.boxes{1,1}, S.boxes{1,2}, S.land_g);
    e = box_mean(M, S.lon_g, S.lat_g, S.boxes{2,1}, S.boxes{2,2}, S.land_g);
    if is_tr
        [R(ip).wio, R(ip).wio_p] = ts_trend(w, S.reg_years, 10);
        [R(ip).eio, R(ip).eio_p] = ts_trend(e, S.reg_years, 10);
    else
        R(ip).wio = mean(w, 'omitnan');  R(ip).wio_p = NaN;
        R(ip).eio = mean(e, 'omitnan');  R(ip).eio_p = NaN;
    end
end

if is_tr; unit = 'g kg^{-1} decade^{-1}'; else; unit = 'g kg^{-1}'; end
fprintf('\n===== %s SON q %s：区域平均 =====\n', S.lev_name, S.mode);
fprintf('%-6s %10s %8s %10s %8s\n', 'Panel', 'WIO', 'p', 'EIO', 'p');
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    fprintf('%-6s %+10.3f %8.3f %+10.3f %8.3f\n', panel_names{ip}, ...
        R(ip).wio, R(ip).wio_p, R(ip).eio, R(ip).eio_p);
end

base = fullfile(S.out_dir, sprintf('FigS_%s_SON_%s_%d_%d_7models_MME', ...
    S.tag, S.mode, S.yr_start, S.yr_end));
lon_g = S.lon_g; lat_g = S.lat_g; land_g = S.land_g; reg_years = S.reg_years; %#ok<NASGU>
save([base '.mat'], 'R', 'panel_names', 'lon_g', 'lat_g', 'reg_years', 'land_g');

% ---- 裁剪、色标 ----
mlon_p = S.lon_g >= S.plot_lon(1) & S.lon_g <= S.plot_lon(2);
mlat_p = S.lat_g >= S.plot_lat(1) & S.lat_g <= S.plot_lat(2);
lon_p  = S.lon_g(mlon_p);  lat_p = S.lat_g(mlat_p);
vals = [];
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    P(ip).v = R(ip).v(mlat_p, mlon_p); %#ok<AGROW>
    P(ip).p = R(ip).p(mlat_p, mlon_p);
    vals = [vals; P(ip).v(isfinite(P(ip).v))]; %#ok<AGROW>
end

c_lbls = {};
if is_tr
    % precip4_diff_19lev（去掉中间白色）：选整齐的色阶宽度，使刻度落在色阶分界上
    cmap  = S.cmap;
    n_lev = size(cmap, 1);
    cl = S.clim;
    if isempty(cl)
        lev_w = nice_ceil(2 * prctile(abs(vals), 98) / n_lev);
        cl    = n_lev * lev_w / 2;
    else
        lev_w = 2 * cl / n_lev;
    end
    c_lim   = [-cl cl];
    c_ticks = linspace(-cl, cl, n_lev + 1);                % 每个色阶分界
    k_edge  = round(c_ticks / lev_w * 2);
    k_min   = min(abs(k_edge));
    show_lb = mod((abs(k_edge) - k_min) / 2, S.tick_lbl_step) == 0;
    c_lbls  = arrayfun(@(x) sprintf('%g', round(x, 6)), c_ticks, 'UniformOutput', false);
    c_lbls(~show_lb) = {''};
    fprintf('色阶数 %d，色阶宽度 %g\n', n_lev, lev_w);
else
    if isempty(S.clim)
        c_lim = [floor(prctile(vals, 2)), ceil(prctile(vals, 98))];
    else
        c_lim = S.clim;
    end
    c_ticks = [];
    nC   = size(S.cmap, 1);
    cmap = S.cmap(ceil(nC/2)+1:end, :);                    % 气候态：只用湿端一半（浅 → 深）
end
fprintf('色标 = [%g, %g] %s\n', c_lim(1), c_lim(2), unit);

% ---- 布局：4 行 2 列 ----
silver_rgb = [0.85 0.85 0.85];
fs = S.fs;
cm2px    = 96/2.54;
fig_w_cm = 18;
pw_cm  = 0.32*18;
ph_cm  = pw_cm * diff(S.plot_lat) / diff(S.plot_lon);
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
land_p   = S.land_g(mlat_p, mlon_p);
land_img = repmat(reshape(silver_rgb,1,1,3), [size(land_p) 1]);

letters = 'abcdefgh';
ax_last = [];
for ip = 1:n_panels
    ir = ceil(ip / 2);  ic = mod(ip-1, 2) + 1;
    pos = [xs(ic), y_cm(ir)/fig_h_cm, pw, ph];
    ax  = axes('Position', pos);  hold(ax, 'on');

    if R(ip).ok
        dat = P(ip).v;
        imagesc(ax, lon_p, lat_p, land_img);
        imagesc(ax, lon_p, lat_p, dat, 'AlphaData', ~isnan(dat));
        if is_tr
            sg = P(ip).p < S.sig_level & ~isnan(dat);
            plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize', S.dot_sz, ...
                 'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerEdgeColor', 'none');
        end
        ttl = sprintf('(%c) %s', letters(ip), panel_names{ip});
    else
        ttl = sprintf('(%c) %s (no data)', letters(ip), panel_names{ip});
    end

    set(ax, 'YDir', 'normal', 'Color', silver_rgb);
    colormap(ax, cmap);
    caxis(ax, c_lim);
    contour(ax, S.topo.lon, S.topo.lat, S.topo.elev, [0 0], 'k-', 'LineWidth', 0.5);
    if S.draw_boxes
        for ib = 1:size(S.boxes, 1)
            bx = S.boxes{ib, 1};  by = S.boxes{ib, 2};
            plot(ax, bx([1 2 2 1 1]), by([1 1 2 2 1]), '-', 'Color', S.boxes{ib, 3}, 'LineWidth', S.box_lw);
        end
    end
    xlim(ax, S.plot_lon);  ylim(ax, S.plot_lat);
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

% ---- 底部横向色标 ----
cb = colorbar(ax_last, 'southoutside');
cb.Position = [xs(1)+0.1*pw, (y_bot - 1.3)/fig_h_cm, xs(2)+pw-xs(1)-0.2*pw, 0.3/fig_h_cm];
cb.Limits = c_lim;
if ~isempty(c_ticks); cb.Ticks = c_ticks; end
if ~isempty(c_lbls);  cb.TickLabels = c_lbls; end
if is_tr
    cb.Label.String = sprintf('SON %s specific humidity trend, %d–%d (%s)', ...
                              S.lev_name, S.yr_start, S.yr_end, unit);
else
    cb.Label.String = sprintf('SON %s specific humidity, %d–%d (%s)', ...
                              S.lev_name, S.yr_start, S.yr_end, unit);
end
cb.FontSize       = fs;
cb.Label.FontSize = fs;
cb.TickDirection  = 'out';
cb.TickLength     = 0.015;
ax_last.Position  = [xs(2), y_cm(end)/fig_h_cm, pw, ph];

set(findall(fig, 'Type', 'axes'), 'TickLength', [0.025 0.035]);

exportgraphics(fig, [base '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig, [base '.tif'], 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('Saved: %s.png / .tif\n', base);
end

function F = read_layer(file, var, il)
% 读若干气压层并取层平均，返回 lon x lat x time；-999 → NaN
k0 = min(il);  nk = max(il) - k0 + 1;
F = double(ncread(file, var, [1 1 k0 1], [Inf Inf nk Inf]));
F = F(:,:,il - k0 + 1,:);
F(F <= -998) = NaN;
F = squeeze(mean(F, 3, 'omitnan'));
end

function [yr, mo] = model_time(file)
info  = ncinfo(file);
names = {info.Variables.Name};
if any(strcmp(names, 'time')); tn = 'time'; else; tn = 'time_counter'; end
t_raw = double(ncread(file, tn));
units = ncreadatt(file, tn, 'units');
tok = regexp(units, '(\w+)\s+since\s+(\d{4})-(\d{1,2})', 'tokens', 'once');
if isempty(tok) || ~contains(lower(tok{1}), 'month')
    error('时间单位不是 months since：%s', units);
end
y0 = str2double(tok{2});  m0 = str2double(tok{3});
nt = numel(t_raw);
t1 = t_raw(find(isfinite(t_raw), 1));
if isempty(t1); t1 = 0; end
tm = round(t1) + (0:nt-1)';
yr = y0 + floor((m0 - 1 + tm) / 12);
mo = mod(m0 - 1 + tm, 12) + 1;
end

function [slope, pval] = trend_decade(Y, years)
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
