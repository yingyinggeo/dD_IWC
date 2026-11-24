%% ========================================================================
%  SON ω 趋势（400–600 hPa 平均）：7 个模式 + MME，4 行 2 列
%   参照 Fig2 (a) 与 MFC 多模式图的版式：
%     - 填色：400–600 hPa 平均 ω 的 SON 线性趋势
%     - ω 趋势显著 (p<0.05) → 绿点
%     - 陆地（1° topo 最近点 > 0）屏蔽为浅灰
%   数据：C:\Users\yingy\Downloads\Rice_mac\omega\omega_Total.<model>_Monthly.nc
%         (128 x 64 网格，17 个气压层，单位 mb；变量 omega，lon x lat x p x time)
%   单位：若原始 ω 量级像 Pa s^-1，自动 ×864 转为 hPa day^-1（同 ERA5 处理）
%   计算范围：60°S–60°N, 0–180°E；画图范围：40°S–40°N, 30–120°E
%   时段：优先 2003–2020；若有模式缺 2020 年 SON，则统一用 2003–2019
%   MME：先对各模式逐年 SON 场取集合平均，再算趋势
%% ========================================================================
clear; clc; close all;

%% ===== 基本设置 =====
data_dir = 'C:\Users\yingy\Downloads\Rice_mac\omega';
models   = {'CAM5', 'CAM6', 'ECHAM', 'GISS', 'GSM', 'LMDZ', 'MIROC'};   % 可换成 NICAM
n_models = numel(models);

p_band   = [400 600];                    % hPa，垂直平均范围
yr_start = 2003;  yr_try_end = 2020;  yr_fallback_end = 2019;
years_all = (yr_start:yr_try_end)';
nyr_all   = numel(years_all);
son = [9 10 11];

calc_lon = [0 180];   calc_lat = [-60 60];
plot_lon = [30 120];  plot_lat = [-40 40];

trend_per = 1;            % 趋势单位：1 = per year（同 Fig2 a）；10 = per decade
clim_w    = [];           % 色标 ±；留空则按所有面板 |趋势| 的 98% 分位自动取整
sig_level = 0.05;
fs        = 10;

% 区域框（Fig. 1c）：{经度范围, 纬度范围, 颜色}
boxes = { ...
    [50 70],  [-15 5], [1.00 0.50 0.00];   % western Indian Ocean (WIO)，橙色
    [70 100], [-20 0], [0.00 0.45 0.90]};  % eastern Indian Ocean (EIO)，蓝色
box_lw = 1.2;

%% ===== 逐模式读取 → 400–600 hPa 平均 ω → 逐年 SON 平均 =====
W_all = [];
have    = false(n_models, 1);
ok_last = false(n_models, 1);
conv    = NaN(n_models, 1);

for im = 1:n_models
    f = fullfile(data_dir, sprintf('omega_Total.%s_Monthly.nc', models{im}));
    fprintf('\n===== %s =====\n', models{im});
    if ~exist(f, 'file'); warning('找不到 %s', f); continue; end

    try
        lon = double(ncread(f, 'lon'));
        lat = double(ncread(f, 'lat'));
        p   = double(ncread(f, 'p'));
        kk  = find(p >= p_band(1) & p <= p_band(2));
        if isempty(kk); error('没有 %d–%d hPa 的层', p_band); end
        fprintf('垂直平均的层：%s hPa\n', mat2str(p(kk)'));

        [yr, mo] = model_time(f);
        fprintf('时间：%d-%02d 到 %d-%02d（%d 个月）\n', yr(1), mo(1), yr(end), mo(end), numel(yr));

        % 只读 SON、所需年份的月份，减少内存
        it = find(ismember(mo, son) & ismember(yr, years_all));
        W  = read_band(f, 'omega', kk, it);           % lon x lat x numel(it)
        yr_it = yr(it);
    catch ME
        warning('%s 读取失败：%s', models{im}, ME.message);
        continue
    end

    % 单位判断：Pa s^-1 的典型量级 ~0.01–0.1；hPa day^-1 ~1–100
    typ = median(abs(W(isfinite(W))));
    if typ < 1
        conv(im) = 864;   unit_msg = 'Pa s^-1 → ×864 转为 hPa day^-1';
    else
        conv(im) = 1;     unit_msg = '视为已是 hPa day^-1';
    end
    W = W * conv(im);
    fprintf('|ω| 中位数 = %.4g（原始），%s\n', typ, unit_msg);

    % 整理为 lat(南→北) x lon(0→360) x t，裁剪到计算范围
    [lonS, ix] = sort(mod(lon, 360));
    [latS, iy] = sort(lat);
    mlon = lonS >= calc_lon(1) & lonS <= calc_lon(2);
    mlat = latS >= calc_lat(1) & latS <= calc_lat(2);
    lon_c = lonS(mlon);  lat_c = latS(mlat);
    W = permute(W(ix(mlon), iy(mlat), :), [2 1 3]);

    if isempty(W_all)
        lon_g = lon_c;  lat_g = lat_c;
        W_all = NaN(numel(lat_g), numel(lon_g), nyr_all, n_models);
    elseif numel(lon_c) ~= numel(lon_g) || numel(lat_c) ~= numel(lat_g) || ...
           max(abs(lon_c - lon_g)) > 1e-3 || max(abs(lat_c - lat_g)) > 1e-3
        error('%s 的网格与其他模式不一致', models{im});
    end

    for k = 1:nyr_all
        sel = yr_it == years_all(k);
        if sum(sel) ~= 3; continue; end
        A = mean(W(:,:,sel), 3);
        if all(isnan(A(:))); continue; end
        W_all(:,:,k,im) = A;
    end
    ok_last(im) = any(~isnan(reshape(W_all(:,:,end,im), [], 1)));
    have(im) = true;
    fprintf('%s 完成；%d 年 SON %s\n', models{im}, yr_try_end, ternary(ok_last(im), '可用', '缺失'));
    clear W
end

if ~any(have); error('没有可用的模式数据'); end
if numel(unique(conv(have))) > 1
    warning('各模式 ω 的单位判断不一致：%s，请检查原始单位', mat2str(conv(have)'));
end

%% ===== 确定时段 =====
if all(ok_last(have))
    yr_end = yr_try_end;
else
    yr_end = yr_fallback_end;
    fprintf('\n以下模式缺 %d 年 SON：%s → 统一使用 %d–%d\n', yr_try_end, ...
        strjoin(models(have & ~ok_last), ', '), yr_start, yr_end);
end
reg_years = (yr_start:yr_end)';
W_all = W_all(:,:,ismember(years_all, reg_years),:);
fprintf('时段：%d–%d（%d 年），模式数 = %d (%s)\n', yr_start, yr_end, numel(reg_years), ...
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
tlon_idx    = topo_lon_ext >= plot_lon(1)-buf & topo_lon_ext <= plot_lon(2)+buf;
tlat_idx    = topo_lat     >= plot_lat(1)-buf & topo_lat     <= plot_lat(2)+buf;
topo_lon_v  = topo_lon_ext(tlon_idx);
topo_lat_v  = topo_lat(tlat_idx);
topo_elev_v = topo_ext(tlat_idx, tlon_idx);

%% ===== 趋势：7 个模式 + MME =====
n_panels    = n_models + 1;
panel_names = [models, {'MME'}];
R = struct([]);
for ip = 1:n_panels
    if ip <= n_models
        if ~have(ip); R(ip).ok = false; continue; end
        M = W_all(:,:,:,ip);
    else
        M = mean(W_all(:,:,:,have), 4, 'omitnan');
    end
    [tr, pv] = trend_per_unit(M, reg_years, trend_per);
    tr(land_g) = NaN;  pv(land_g) = NaN;
    R(ip).tr = tr;  R(ip).p = pv;  R(ip).ok = true;
end

save(fullfile('C:\Users\yingy\Downloads\Rice_mac', ...
    sprintf('omega%d_%d_SON_trend_%d_%d_7models_MME.mat', p_band, yr_start, yr_end)), ...
    'R', 'panel_names', 'lon_g', 'lat_g', 'reg_years', 'land_g', 'trend_per', 'p_band');

%% ===== 偶极子诊断：WIO（橙框）vs EIO（蓝框）区域平均 =====
%  区域平均只用海洋格点，cos(lat) 加权；对每个模式和 MME：
%   1) 两个框各自区域平均 ω 的 SON 趋势及显著性
%   2) 偶极子指数 D = WIO − EIO 的趋势及显著性
%   3) 去趋势后 WIO 与 EIO 年际相关（负相关 → 年际尺度上反向变化）
%  判断：两框趋势符号相反 → 趋势呈偶极子型；
%        两框都显著 → 强；只有一个显著或 D 显著 → 中等；都不显著 → 弱
nyr = numel(reg_years);
nb  = size(boxes, 1);
TS  = NaN(nyr, nb, n_panels);
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    if ip <= n_models
        M = W_all(:,:,:,ip);
    else
        M = mean(W_all(:,:,:,have), 4, 'omitnan');
    end
    for ib = 1:nb
        TS(:,ib,ip) = box_mean(M, lon_g, lat_g, boxes{ib,1}, boxes{ib,2}, land_g);
    end
end

unit_tr = sprintf('hPa day^-1 %s^-1', ternary(trend_per==1, 'yr', 'decade'));
Name = {}; WIO_tr = []; WIO_p = []; EIO_tr = []; EIO_p = [];
D_tr = []; D_p = []; r_det = []; r_p = []; Opposite = {}; Dipole = {};
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    w = TS(:,1,ip);  e = TS(:,2,ip);
    [bw, pw_] = ts_trend(w, reg_years, trend_per);
    [be, pe ] = ts_trend(e, reg_years, trend_per);
    [bd, pd ] = ts_trend(w - e, reg_years, trend_per);
    wd = detrend(w);  ed = detrend(e);
    [rr, rp] = corr(wd, ed);

    opp = sign(bw) ~= sign(be);
    if ~opp
        lab = 'No (same sign)';
    elseif pw_ < sig_level && pe < sig_level
        lab = 'Yes - strong (both sig.)';
    elseif pw_ < sig_level || pe < sig_level || pd < sig_level
        lab = 'Yes - moderate';
    else
        lab = 'Yes - weak (not sig.)';
    end

    Name{end+1,1} = panel_names{ip}; %#ok<SAGROW>
    WIO_tr(end+1,1) = bw;  WIO_p(end+1,1) = pw_;   %#ok<SAGROW>
    EIO_tr(end+1,1) = be;  EIO_p(end+1,1) = pe;    %#ok<SAGROW>
    D_tr(end+1,1)   = bd;  D_p(end+1,1)   = pd;    %#ok<SAGROW>
    r_det(end+1,1)  = rr;  r_p(end+1,1)   = rp;    %#ok<SAGROW>
    Opposite{end+1,1} = ternary(opp, 'Yes', 'No'); %#ok<SAGROW>
    Dipole{end+1,1}   = lab;                       %#ok<SAGROW>
end
T = table(Name, WIO_tr, WIO_p, EIO_tr, EIO_p, D_tr, D_p, r_det, r_p, Opposite, Dipole);
fprintf('\n===== 偶极子诊断（SON %d–%d hPa ω，%d–%d；趋势单位 %s）=====\n', ...
    p_band(1), p_band(2), yr_start, yr_end, unit_tr);
fprintf('ω > 0 为下沉，ω < 0 为上升；趋势 > 0 表示下沉加强 / 上升减弱\n');
disp(T);
n_dip = sum(strcmp(T.Opposite(1:end-1), 'Yes'));
fprintf('%d/%d 个模式的两框趋势符号相反；MME：%s\n', n_dip, height(T)-1, T.Dipole{end});
writetable(T, fullfile('C:\Users\yingy\Downloads\Rice_mac', ...
    sprintf('omega%d_%d_SON_dipole_WIO_EIO_%d_%d.csv', p_band, yr_start, yr_end)));

% 柱状图：各模式与 MME 的 WIO / EIO 区域平均趋势（* 表示 p<0.05）
figD = figure('Color','w', 'Units','centimeters', 'Position',[2 2 18 7]);
axD = axes(figD, 'Position', [0.09 0.2 0.88 0.7]);  hold(axD, 'on');
hb = bar(axD, [T.WIO_tr T.EIO_tr], 'grouped', 'EdgeColor', 'none');
hb(1).FaceColor = boxes{1,3};  hb(2).FaceColor = boxes{2,3};
yl = max(abs([T.WIO_tr; T.EIO_tr])) * 1.3;
for ib = 1:2
    xe = hb(ib).XEndPoints;
    ye = hb(ib).YData;
    pp = ternary(ib == 1, T.WIO_p, T.EIO_p);
    for j = 1:numel(xe)
        if pp(j) < sig_level
            text(axD, xe(j), ye(j) + sign(ye(j))*0.06*yl, '*', 'FontSize', 12, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end
plot(axD, [0.4 height(T)+0.6], [0 0], 'k-', 'LineWidth', 0.5);
set(axD, 'XTick', 1:height(T), 'XTickLabel', T.Name, 'FontSize', fs, ...
    'TickDir', 'out', 'Box', 'off', 'XLim', [0.4 height(T)+0.6], 'YLim', [-yl yl]);
ylabel(axD, sprintf('\\omega trend (hPa day^{-1} %s^{-1})', ternary(trend_per==1, 'yr', 'decade')));
legend(hb, {'WIO (50–70°E, 15°S–5°N)', 'EIO (70–100°E, 20°S–0°)'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', fs);
exportgraphics(figD, fullfile('C:\Users\yingy\Downloads\Rice_mac', ...
    sprintf('FigS_omega%d_%d_SON_dipole_bar_%d_%d.png', p_band, yr_start, yr_end)), ...
    'Resolution', 600, 'BackgroundColor', 'white');

%% ===== 裁剪到画图范围、确定色标 =====
mlon_p = lon_g >= plot_lon(1) & lon_g <= plot_lon(2);
mlat_p = lat_g >= plot_lat(1) & lat_g <= plot_lat(2);
lon_p  = lon_g(mlon_p);  lat_p = lat_g(mlat_p);
land_p = land_g(mlat_p, mlon_p);
vals = [];
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    P(ip).tr = R(ip).tr(mlat_p, mlon_p); %#ok<SAGROW>
    P(ip).p  = R(ip).p (mlat_p, mlon_p);
    vals = [vals; abs(P(ip).tr(isfinite(P(ip).tr)))]; %#ok<AGROW>
end
if isempty(clim_w)
    clim_w = nice_ceil(prctile(vals, 98));
end
fprintf('色标 = ±%g hPa day^-1 %s^-1\n', clim_w, ternary(trend_per==1, 'yr', 'decade'));

%% ===== 画图：4 行 2 列 =====
load 'BlueWhiteOrangeRed.mat';
silver_rgb = [0.85 0.85 0.85];

cm2px    = 96/2.54;
fig_w_cm = 18;
pw_cm = 0.32*18;   ph_cm = 4.23;   row_dy = 5.17;
nrow  = 4;
y_bot = 2.6;
y_cm  = y_bot + (nrow-1:-1:0) * row_dy;
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
land_img = repmat(reshape(silver_rgb,1,1,3), [size(land_p) 1]);

letters = 'abcdefgh';
ax_last = [];
for ip = 1:n_panels
    ir = ceil(ip / 2);  ic = mod(ip-1, 2) + 1;
    pos = [xs(ic), y_cm(ir)/fig_h_cm, pw, ph];
    ax  = axes('Position', pos);  hold(ax, 'on');

    if R(ip).ok
        dat = P(ip).tr;
        imagesc(ax, lon_p, lat_p, land_img);
        imagesc(ax, lon_p, lat_p, dat, 'AlphaData', ~isnan(dat));
        sg = P(ip).p < sig_level & ~isnan(dat);
        plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize', 3, ...
             'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerEdgeColor', 'none');
        ttl = sprintf('(%c) %s', letters(ip), panel_names{ip});
    else
        ttl = sprintf('(%c) %s (no data)', letters(ip), panel_names{ip});
    end

    set(ax, 'YDir', 'normal', 'Color', silver_rgb);
    colormap(ax, BlueWhiteOrangeRed);
    caxis(ax, [-clim_w clim_w]);
    contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth', 0.5);

    % 区域框（Fig. 1c）：WIO 橙色，EIO 蓝色
    for ib = 1:size(boxes, 1)
        bx = boxes{ib, 1};  by = boxes{ib, 2};
        plot(ax, bx([1 2 2 1 1]), by([1 1 2 2 1]), '-', ...
             'Color', boxes{ib, 3}, 'LineWidth', box_lw);
    end
    xlim(ax, plot_lon);  ylim(ax, plot_lat);

    if ir == nrow
        set(ax, 'XTick', 30:30:120, 'XTickLabel', {'30°E','60°E','90°E','120°E'});
    else
        set(ax, 'XTick', 30:30:120, 'XTickLabel', {});
    end
    if ic == 1
        set(ax, 'YTick', -40:20:40, 'YTickLabel', {'40°S','20°S','0°','20°N','40°N'});
    else
        set(ax, 'YTick', -40:20:40, 'YTickLabel', {});
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

%% ===== 底部横向色标 =====
cb = colorbar(ax_last, 'southoutside');
cb.Position = [xs(1)+0.1*pw, (y_bot - 1.3)/fig_h_cm, xs(2)+pw-xs(1)-0.2*pw, 0.3/fig_h_cm];
cb.Ticks = clim_w * [-1 -0.5 0 0.5 1];
cb.Label.String = sprintf('SON %d–%d-hPa \\omega trend, %d–%d (hPa day^{-1} %s^{-1})', ...
    p_band(1), p_band(2), yr_start, yr_end, ternary(trend_per==1, 'yr', 'decade'));
cb.FontSize       = fs;
cb.Label.FontSize = fs;
cb.TickDirection  = 'out';
ax_last.Position  = [xs(2), y_cm(end)/fig_h_cm, pw, ph];

set(findall(fig, 'Type', 'axes'), 'TickLength', [0.025 0.035]);

%% ===== 保存 =====
out_dir = 'C:\Users\yingy\Downloads\Rice_mac';
base = fullfile(out_dir, sprintf('FigS_omega%d_%d_SON_trend_%d_%d_7models_MME', p_band, yr_start, yr_end));
exportgraphics(fig, [base '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig, [base '.tif'], 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('Saved: %s.png / .tif\n', base);

%% ===== 子函数 =====
function W = read_band(file, var, kk, it)
% 读 kk 指定的连续气压层与 it 指定的时次，垂直平均（忽略缺测），返回 lon x lat x numel(it)
k0 = min(kk);  nk = max(kk) - k0 + 1;
W = NaN(0);
for j = 1:numel(it)
    A = ncread(file, var, [1 1 k0 it(j)], [Inf Inf nk 1]);   % lon x lat x nk
    A = double(A);
    A(A <= -998) = NaN;
    A = A(:,:,kk - k0 + 1);
    if j == 1; W = NaN(size(A,1), size(A,2), numel(it)); end
    W(:,:,j) = mean(A, 3, 'omitnan');
end
end

function [yr, mo] = model_time(file)
% months since YYYY-M-D：按逐月连续重建（防止时间轴有异常值或 NaN）
t_raw = double(ncread(file, 'time'));
units = ncreadatt(file, 'time', 'units');
tok = regexp(units, '(\w+)\s+since\s+(\d{4})-(\d{1,2})', 'tokens', 'once');
if isempty(tok) || ~contains(lower(tok{1}), 'month')
    error('时间单位不是 months since：%s', units);
end
y0 = str2double(tok{2});  m0 = str2double(tok{3});
nt = numel(t_raw);
t1 = t_raw(find(isfinite(t_raw), 1));
if isempty(t1); t1 = 0; end
tm = round(t1) + (0:nt-1)';
bad = find(~isfinite(t_raw(:)) | abs(t_raw(:) - tm) > 0.5);
if ~isempty(bad)
    fprintf('  时间轴有 %d 个异常值，已按逐月连续重建\n', numel(bad));
end
yr = y0 + floor((m0 - 1 + tm) / 12);
mo = mod(m0 - 1 + tm, 12) + 1;
end

function [slope, pval] = trend_per_unit(Y, years, per)
% 逐格点 OLS 趋势（per 年为单位，per=1 为每年，per=10 为每十年）+ 双侧 t 检验
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
slope = reshape(b * per, d1, d2);
pval  = reshape(p, d1, d2);
end

function ts = box_mean(M, lon, lat, blon, blat, land)
% M: lat x lon x nyr；海洋格点、cos(lat) 加权的区域平均，返回 nyr x 1
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
% 单条时间序列的 OLS 趋势（per 年为单位）+ 双侧 t 检验
x = years(:) - mean(years);  y = y(:);
ok = isfinite(y);  x = x(ok);  y = y(ok);  n = numel(y);
x = x - mean(x);
if n < 4; b = NaN; p = NaN; return; end
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
