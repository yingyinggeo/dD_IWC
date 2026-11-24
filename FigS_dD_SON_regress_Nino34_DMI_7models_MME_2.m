%% ========================================================================
%  SON δD 距平（400–600 hPa 平均）对 ENSO / IOD 指数的回归：7 个模式 + MME
%  三幅图，每幅 4 行 2 列（版式同 ω 趋势多模式图）：
%    图 1：SON δD 距平 对 D(0)JF(1) Niño3.4 回归
%    图 2：SON δD 距平 对 SON DMI 回归
%    图 3：去除 Niño3.4 的 SON δD 距平 对 去除 Niño3.4 的 SON DMI 回归（偏回归）
%          —— δD 与 DMI 各自先对 D(0)JF(1) Niño3.4 做线性回归，取残差
%  δD：各模式水汽 δD（Total.<model>_Monthly.nc 中的 dD，17 个气压层），400–600 hPa 平均
%  指数：每个模式用各自的 ts（WISOMIP_Other_Variabiles）计算；
%        相对研究时段的逐月气候态求距平，海洋格点 cos(lat) 加权区域平均
%    Niño3.4 = SSTA(170°W–120°W, 5°S–5°N)，取当年 12 月 + 次年 1、2 月
%    DMI     = SSTA(50–70°E, 10°S–10°N) − SSTA(90–110°E, 10°S–0°)，取当年 SON
%  MME：δD 和 ts 距平场先做集合平均，再由 MME 的 ts 算指数并回归
%  回归：指数标准化，系数单位 ‰ per SD；双侧 t 检验 p<0.05 → 绿点
%        （图 3 多估计一个参数，自由度 n−3）
%  只画图 3（图 1、2 仍会计算并保存到 .mat；想画可改 figs_to_plot）
%  计算范围：60°S–60°N, 0–180°E；画图范围：20°S–40°N, 30–120°E
%        （面板高度按纬度跨度缩放，保持真实长宽比；不画 WIO/EIO 框）
%  时段：优先 2003–2020（需要 2021 年 1–2 月 ts）；不够则统一用 2003–2019
%% ========================================================================
clear; clc; close all;

%% ===== 基本设置 =====
ts_dir    = 'C:\Users\yingy\Downloads\Rice_mac\WISOMIP_Other_Variabiles';
out_dir   = 'C:\Users\yingy\Downloads\Rice_mac';
models    = {'CAM5', 'CAM6', 'ECHAM', 'GISS', 'GSM', 'LMDZ', 'MIROC'};
n_models  = numel(models);

p_band   = [400 600];                    % hPa，垂直平均范围
yr_start = 2003;  yr_try_end = 2020;  yr_fallback_end = 2019;
years_all = (yr_start:yr_try_end)';
nyr_all   = numel(years_all);
son = [9 10 11];

detrend_flag = true;      % 回归前对 δD、指数去线性趋势（同模式 SSTA 回归）；false 则不去
calc_lon = [0 180];   calc_lat = [-60 60];
plot_lon = [30 120];  plot_lat = [-20 40];      % 画图范围 20°S–40°N
plot_lat_ref = [-40 40];   % 原画图纬度范围：用于保持长宽比，以及自动色标取值（与原图一致）
figs_to_plot = 3;          % 只画图 3；改为 1:3 则三幅都画

clim_reg  = [];           % 色标 ±（‰ per SD）；留空则每幅图按 98% 分位自动取整
sig_level = 0.05;
fs        = 10;
dot_sz    = 3;

% 指数区域（经度 0–360）
box_n34 = {[190 240], [-5 5]};
box_wio = {[50 70],   [-10 10]};
box_eio = {[90 110],  [-10 0]};

% 图上的区域框（Fig. 1c）
boxes = { ...
    [50 70],  [-15 5], [1.00 0.50 0.00];   % WIO，橙色
    [70 100], [-20 0], [0.00 0.45 0.90]};  % EIO，蓝色
box_lw = 1.2;

load topo;
topo_lon = 0.5:359.5;  topo_lat = -89.5:89.5;
[TLON, TLAT] = meshgrid([topo_lon 360.5], topo_lat);
topo_ext = [topo topo(:,1)];

%% ===== 1. 逐模式读取 δD → 400–600 hPa 平均 → 逐年 SON 平均 =====
W_all = [];  have_w = false(n_models,1);  okw_last = false(n_models,1);
for im = 1:n_models
    f = fullfile(ts_dir, sprintf('Total.%s_Monthly.nc', models{im}));
    fprintf('\n===== %s：dD =====\n', models{im});
    if ~exist(f, 'file'); warning('找不到 %s', f); continue; end
    try
        lon = double(ncread(f, 'lon'));  lat = double(ncread(f, 'lat'));
        p   = double(ncread(f, 'p'));
        kk  = find(p >= p_band(1) & p <= p_band(2));
        if isempty(kk); error('没有 %d–%d hPa 的层', p_band); end
        [yr, mo] = model_time(f);
        it = find(ismember(mo, son) & ismember(yr, years_all));
        W  = read_band(f, 'dD', kk, it);
        yr_it = yr(it);
    catch ME
        warning('%s dD 读取失败：%s', models{im}, ME.message);  continue
    end
    fprintf('垂直平均层 %s hPa；δD 范围 [%.1f, %.1f] ‰\n', mat2str(p(kk)'), ...
        min(W(:), [], 'omitnan'), max(W(:), [], 'omitnan'));

    [lonS, ix] = sort(mod(lon, 360));  [latS, iy] = sort(lat);
    mlon = lonS >= calc_lon(1) & lonS <= calc_lon(2);
    mlat = latS >= calc_lat(1) & latS <= calc_lat(2);
    lon_c = lonS(mlon);  lat_c = latS(mlat);
    W = permute(W(ix(mlon), iy(mlat), :), [2 1 3]);
    if isempty(W_all)
        lon_g = lon_c;  lat_g = lat_c;
        W_all = NaN(numel(lat_g), numel(lon_g), nyr_all, n_models);
    elseif numel(lon_c) ~= numel(lon_g) || numel(lat_c) ~= numel(lat_g) || ...
           max(abs(lon_c - lon_g)) > 1e-3 || max(abs(lat_c - lat_g)) > 1e-3
        error('%s 的 dD 网格与其他模式不一致', models{im});
    end
    for k = 1:nyr_all
        sel = yr_it == years_all(k);
        if sum(sel) ~= 3; continue; end
        A = mean(W(:,:,sel), 3);
        if ~all(isnan(A(:))); W_all(:,:,k,im) = A; end
    end
    okw_last(im) = any(isfinite(reshape(W_all(:,:,end,im), [], 1)));
    have_w(im) = true;
end

%% ===== 2. 逐模式读取 ts → 逐月 SSTA → SON / D(0)JF(1) 场 =====
S_son = [];  S_djf = [];  have_t = false(n_models,1);  okt_last = false(n_models,1);
for im = 1:n_models
    f = fullfile(ts_dir, sprintf('Total.%s_Monthly.nc', models{im}));
    fprintf('\n===== %s：ts =====\n', models{im});
    if ~exist(f, 'file'); warning('找不到 %s', f); continue; end
    try
        lon = double(ncread(f, 'lon'));  lat = double(ncread(f, 'lat'));
        [yr, mo] = model_time(f);
        X = double(ncread(f, 'ts'));  X(X <= -998) = NaN;
    catch ME
        warning('%s ts 读取失败：%s', models{im}, ME.message);  continue
    end
    [lonS, ix] = sort(mod(lon, 360));  [latS, iy] = sort(lat);
    X = permute(X(ix, iy, :), [2 1 3]);                  % lat x lon(0–360) x t
    if isempty(S_son)
        lon_t = lonS;  lat_t = latS;
        S_son = NaN(numel(lat_t), numel(lon_t), nyr_all, n_models);  S_djf = S_son;
    elseif numel(lonS) ~= numel(lon_t) || numel(latS) ~= numel(lat_t) || ...
           max(abs(lonS - lon_t)) > 1e-3 || max(abs(latS - lat_t)) > 1e-3
        error('%s 的 ts 网格与其他模式不一致', models{im});
    end
    % 逐月气候态：先用 2003–2020，时段确定后若改为 2019 也影响很小
    inper = yr >= yr_start & yr <= yr_try_end;
    for m = 1:12
        km = mo == m;
        X(:,:,km) = X(:,:,km) - mean(X(:,:,km & inper), 3, 'omitnan');
    end
    for k = 1:nyr_all
        y  = years_all(k);
        ks = yr == y & ismember(mo, son);
        kd = (yr == y & mo == 12) | (yr == y+1 & mo <= 2);
        if sum(ks) == 3; S_son(:,:,k,im) = mean(X(:,:,ks), 3); end
        if sum(kd) == 3; S_djf(:,:,k,im) = mean(X(:,:,kd), 3); end
    end
    okt_last(im) = any(isfinite(reshape(S_djf(:,:,end,im), [], 1)));
    have_t(im) = true;
    clear X
end

have = have_w & have_t;
if ~any(have); error('没有同时具备 dD 和 ts 的模式'); end

%% ===== 3. 时段 =====
if all(okw_last(have) & okt_last(have)); yr_end = yr_try_end; else; yr_end = yr_fallback_end; end
iyr = ismember(years_all, yr_start:yr_end);
reg_years = years_all(iyr);
W_all = W_all(:,:,iyr,:);  S_son = S_son(:,:,iyr,:);  S_djf = S_djf(:,:,iyr,:);
fprintf('\n时段：%d–%d（%d 年）；模式：%s\n', yr_start, yr_end, numel(reg_years), ...
    strjoin(models(have), ', '));

%% ===== 4. 陆地 mask =====
[GLON, GLAT] = meshgrid(lon_g, lat_g);
land_g = interp2(TLON, TLAT, topo_ext, GLON, GLAT, 'nearest', 0) > 0;     % δD 网格
[GLT, GLA] = meshgrid(lon_t, lat_t);
land_t = interp2(TLON, TLAT, topo_ext, GLT, GLA, 'nearest', 0) > 0;       % ts 网格（算指数用）

%% ===== 5. 指数与回归：7 个模式 + MME =====
n_panels    = n_models + 1;
panel_names = [models, {'MME'}];
RG = struct([]);
fprintf('\n%-6s %10s %10s %12s %8s %14s\n', 'Panel', 'sd(N34)', 'sd(DMI)', 'r(N34,DMI)', 'p', 'DMI保留方差');
for ip = 1:n_panels
    if ip <= n_models
        if ~have(ip); RG(ip).ok = false; continue; end
        Wf = W_all(:,:,:,ip);  Ss = S_son(:,:,:,ip);  Sd = S_djf(:,:,:,ip);
    else
        Wf = mean(W_all(:,:,:,have), 4, 'omitnan');
        Ss = mean(S_son(:,:,:,have), 4, 'omitnan');
        Sd = mean(S_djf(:,:,:,have), 4, 'omitnan');
    end
    n34 = box_mean(Sd, lon_t, lat_t, box_n34{:}, land_t);
    dmi = box_mean(Ss, lon_t, lat_t, box_wio{:}, land_t) - ...
          box_mean(Ss, lon_t, lat_t, box_eio{:}, land_t);
    if any(isnan([n34; dmi])); error('%s 的指数有缺测年份', panel_names{ip}); end
    sd_raw = [std(n34) std(dmi)];

    Wa = Wf - mean(Wf, 3, 'omitnan');                    % SON δD 距平
    if detrend_flag
        n34 = detrend(n34);  dmi = detrend(dmi);  Wa = detrend_field(Wa);
    end
    [rr, rp] = corr(n34, dmi);

    % 图 1、2
    [b1, p1] = regress_field(Wa, n34, 0);
    [b2, p2] = regress_field(Wa, dmi, 0);
    % 图 3：δD 与 DMI 都去除 Niño3.4
    xn    = n34 - mean(n34);
    dmi_r = (dmi - mean(dmi)) - xn * (xn' * (dmi - mean(dmi))) / (xn' * xn);
    Wa_r  = remove_index(Wa, xn);
    [b3, p3] = regress_field(Wa_r, dmi_r, 1);

    B = {b1, b2, b3};  P = {p1, p2, p3};
    for k = 1:3; B{k}(land_g) = NaN;  P{k}(land_g) = NaN; end
    RG(ip).b = B;  RG(ip).p = P;  RG(ip).ok = true;
    RG(ip).n34 = n34;  RG(ip).dmi = dmi;  RG(ip).dmi_r = dmi_r;

    fprintf('%-6s %10.3f %10.3f %12.2f %8.3f %13.0f%%\n', panel_names{ip}, sd_raw, rr, rp, ...
        100 * var(dmi_r) / var(dmi));
end

% 两个框内的区域平均回归系数（便于判断偶极子）
reg_names = {'N34', 'DMI', 'DMI_noN34'};
fprintf('\n===== 框内平均回归系数（‰ per SD）=====\n');
fprintf('%-6s', 'Panel');
for k = 1:3; fprintf(' %11s %11s', [reg_names{k} '_WIO'], [reg_names{k} '_EIO']); end
fprintf('\n');
for ip = 1:n_panels
    if ~RG(ip).ok; continue; end
    fprintf('%-6s', panel_names{ip});
    for k = 1:3
        w = box_mean(RG(ip).b{k}, lon_g, lat_g, boxes{1,1}, boxes{1,2}, land_g);
        e = box_mean(RG(ip).b{k}, lon_g, lat_g, boxes{2,1}, boxes{2,2}, land_g);
        fprintf(' %+11.3f %+11.3f', w, e);
    end
    fprintf('\n');
end

save(fullfile(out_dir, sprintf('dD%d_%d_SON_regress_Nino34_DMI_%d_%d_7models_MME.mat', ...
    p_band, yr_start, yr_end)), 'RG', 'panel_names', 'lon_g', 'lat_g', 'reg_years', ...
    'land_g', 'detrend_flag', 'p_band', 'box_n34', 'box_wio', 'box_eio');

%% ===== 6. 画图：三幅，各 4 行 2 列 =====
load 'BlueWhiteOrangeRed.mat';
silver_rgb = [0.85 0.85 0.85];
mlon_p = lon_g >= plot_lon(1) & lon_g <= plot_lon(2);
mlat_p = lat_g >= plot_lat(1) & lat_g <= plot_lat(2);
lon_p  = lon_g(mlon_p);  lat_p = lat_g(mlat_p);
land_p = land_g(mlat_p, mlon_p);
[LONs, LATs] = meshgrid(lon_p, lat_p);
land_img = repmat(reshape(silver_rgb,1,1,3), [size(land_p) 1]);

buf = 5;
topo_lon_ext = [topo_lon, 360.5];
tl = topo_lon_ext >= plot_lon(1)-buf & topo_lon_ext <= plot_lon(2)+buf;
tt = topo_lat     >= plot_lat(1)-buf & topo_lat     <= plot_lat(2)+buf;
topo_lon_v = topo_lon_ext(tl);  topo_lat_v = topo_lat(tt);  topo_elev_v = topo_ext(tt, tl);

cb_labs = { ...
    'SON \deltaD anomaly regressed on D(0)JF(1) Niño3.4', ...
    'SON \deltaD anomaly regressed on SON DMI', ...
    'Niño3.4-removed SON \deltaD anomaly regressed on Niño3.4-removed SON DMI'};
file_tags = {'Nino34DJF', 'DMI', 'DMI_noNino34'};

for k = figs_to_plot
    cl = clim_reg;
    if isempty(cl)
        % 自动色标仍按原 40°S–40°N 范围取 98% 分位，保证与原图色标一致
        mlat_ref = lat_g >= plot_lat_ref(1) & lat_g <= plot_lat_ref(2);
        vals = [];
        for ip = 1:n_panels
            if ~RG(ip).ok; continue; end
            A = RG(ip).b{k}(mlat_ref, mlon_p);
            vals = [vals; abs(A(isfinite(A)))]; %#ok<AGROW>
        end
        cl = nice_ceil(prctile(vals, 98));
    end
    fprintf('图 %d 色标 = ±%g ‰ per SD\n', k, cl);

    cm2px = 96/2.54;  fig_w_cm = 18;
    % 保持长宽比：面板宽度不变，高度按纬度跨度缩放（60°/80° × 4.23 ≈ 3.17 cm）
    ph_ref = 4.23;  gap_cm = 5.17 - ph_ref;
    pw_cm  = 0.32*18;
    ph_cm  = ph_ref * diff(plot_lat) / diff(plot_lat_ref);
    row_dy = ph_cm + gap_cm;  nrow = 4;  y_bot = 2.6;
    y_cm  = y_bot + (nrow-1:-1:0) * row_dy;
    fig_h_cm = y_cm(1) + ph_cm + 0.7;
    fig = figure('Color','w', 'Position',[50 50 round(fig_w_cm*cm2px) round(fig_h_cm*cm2px)], ...
                 'Renderer','opengl', 'GraphicsSmoothing','off', ...
                 'PaperUnits','centimeters', 'PaperSize',[fig_w_cm, fig_h_cm], 'PaperPositionMode','auto');
    set(fig, 'DefaultAxesFontSize', fs, 'DefaultTextFontSize', fs);
    pw = pw_cm / fig_w_cm;  ph = ph_cm / fig_h_cm;  xs = [0.07 0.46];

    letters = 'abcdefgh';  ax_last = [];
    for ip = 1:n_panels
        ir = ceil(ip/2);  ic = mod(ip-1,2) + 1;
        pos = [xs(ic), y_cm(ir)/fig_h_cm, pw, ph];
        ax = axes('Position', pos);  hold(ax, 'on');
        if RG(ip).ok
            dat = RG(ip).b{k}(mlat_p, mlon_p);
            pv  = RG(ip).p{k}(mlat_p, mlon_p);
            imagesc(ax, lon_p, lat_p, land_img);
            imagesc(ax, lon_p, lat_p, dat, 'AlphaData', ~isnan(dat));
            sg = pv < sig_level & ~isnan(dat);
            plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize', dot_sz, ...
                 'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerEdgeColor', 'none');
            ttl = sprintf('(%c) %s', letters(ip), panel_names{ip});
        else
            ttl = sprintf('(%c) %s (no data)', letters(ip), panel_names{ip});
        end
        set(ax, 'YDir','normal', 'Color', silver_rgb);
        colormap(ax, BlueWhiteOrangeRed);  caxis(ax, [-cl cl]);
        contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth', 0.5);
        % 不画 WIO / EIO 方框（框内平均回归系数仍会打印）
        xlim(ax, plot_lon);  ylim(ax, plot_lat);
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
        set(ax, 'FontSize', fs, 'TickDir','out', 'Box','on', 'Layer','top');
        ax.Position = pos;
        annotation(fig, 'textbox', [pos(1), pos(2)+ph+0.003, pw, 0.015], 'String', ttl, ...
                   'FontSize', fs, 'FontWeight','bold', 'HorizontalAlignment','left', ...
                   'VerticalAlignment','bottom', 'EdgeColor','none', 'FitBoxToText','off');
        ax.Position = pos;  ax_last = ax;
    end

    cb = colorbar(ax_last, 'southoutside');
    cb.Position = [xs(1)+0.1*pw, (y_bot - 1.3)/fig_h_cm, xs(2)+pw-xs(1)-0.2*pw, 0.3/fig_h_cm];
    cb.Ticks = cl * [-1 -0.5 0 0.5 1];
    cb.Label.String = sprintf('%s, %d–%d (‰ per SD)', cb_labs{k}, yr_start, yr_end);
    cb.FontSize = fs;  cb.Label.FontSize = fs;  cb.TickDirection = 'out';
    ax_last.Position = [xs(2), y_cm(end)/fig_h_cm, pw, ph];
    set(findall(fig, 'Type','axes'), 'TickLength', [0.025 0.035]);

    base = fullfile(out_dir, sprintf('FigS_dD%d_%d_SON_regress_%s_%d_%d_7models_MME', ...
        p_band, file_tags{k}, yr_start, yr_end));
    exportgraphics(fig, [base '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
    exportgraphics(fig, [base '.tif'], 'Resolution', 600, 'BackgroundColor', 'white');
    fprintf('Saved: %s.png / .tif\n', base);
end

%% ===== 子函数 =====
function W = read_band(file, var, kk, it)
% 读 kk 层与 it 时次，垂直平均（忽略缺测），返回 lon x lat x numel(it)
k0 = min(kk);  nk = max(kk) - k0 + 1;
for j = 1:numel(it)
    A = double(ncread(file, var, [1 1 k0 it(j)], [Inf Inf nk 1]));
    A(A <= -998) = NaN;
    A = A(:,:,kk - k0 + 1);
    if j == 1; W = NaN(size(A,1), size(A,2), numel(it)); end
    W(:,:,j) = mean(A, 3, 'omitnan');
end
end

function [yr, mo] = model_time(file)
% months since YYYY-M-D：按逐月连续重建（防止时间轴有异常值或 NaN）
info  = ncinfo(file);  names = {info.Variables.Name};
if any(strcmp(names, 'time')); tn = 'time'; else; tn = 'time_counter'; end
t_raw = double(ncread(file, tn));
units = ncreadatt(file, tn, 'units');
tok = regexp(units, '(\w+)\s+since\s+(\d{4})-(\d{1,2})', 'tokens', 'once');
if isempty(tok) || ~contains(lower(tok{1}), 'month')
    error('时间单位不是 months since：%s', units);
end
y0 = str2double(tok{2});  m0 = str2double(tok{3});
nt = numel(t_raw);
t1 = t_raw(find(isfinite(t_raw), 1));  if isempty(t1); t1 = 0; end
tm = round(t1) + (0:nt-1)';
bad = find(~isfinite(t_raw(:)) | abs(t_raw(:) - tm) > 0.5);
if ~isempty(bad); fprintf('  时间轴有 %d 个异常值，已按逐月连续重建\n', numel(bad)); end
yr = y0 + floor((m0 - 1 + tm) / 12);
mo = mod(m0 - 1 + tm, 12) + 1;
end

function [b, p] = regress_field(Y, x, k_extra)
% Y: lat x lon x n；x 标准化后回归 → 每 1 个标准差的变化；k_extra 为额外估计的参数个数
x = (x(:) - mean(x)) / std(x);
[d1, d2, n] = size(Y);
Y2 = reshape(Y, d1*d2, n)';
ok = ~any(isnan(Y2), 1);
Y2 = Y2 - mean(Y2, 1);
Sxx = sum(x.^2);
bb  = (x' * Y2) / Sxx;
res = Y2 - x * bb;
dof = n - 2 - k_extra;
t   = bb ./ sqrt(sum(res.^2, 1) / dof / Sxx);
pp  = 2 * (1 - tcdf(abs(t), dof));
bb(~ok) = NaN;  pp(~ok) = NaN;
b = reshape(bb, d1, d2);  p = reshape(pp, d1, d2);
end

function Fr = remove_index(F, xn)
% 逐格点线性去除指数 xn（n x 1，已去均值）的影响
[d1, d2, n] = size(F);
Y2 = reshape(F, d1*d2, n)';
Y2 = Y2 - mean(Y2, 1);
Y2 = Y2 - xn * ((xn' * Y2) / (xn' * xn));
Fr = reshape(Y2', d1, d2, n);
end

function Yd = detrend_field(Y)
[d1, d2, n] = size(Y);
Y2 = reshape(Y, d1*d2, n)';
ok = ~any(isnan(Y2), 1);
Y2(:, ok) = detrend(Y2(:, ok));
Yd = reshape(Y2', d1, d2, n);
end

function ts = box_mean(M, lon, lat, blon, blat, land)
% M: lat x lon x n；海洋格点、cos(lat) 加权区域平均，返回 n x 1
mask = (lat(:) >= blat(1) & lat(:) <= blat(2)) & (lon(:)' >= blon(1) & lon(:)' <= blon(2)) & ~land;
W = repmat(cosd(lat(:)), 1, numel(lon)) .* mask;
ts = NaN(size(M, 3), 1);
for k = 1:size(M, 3)
    A = M(:,:,k);
    ok = W > 0 & isfinite(A);
    if any(ok(:)); ts(k) = sum(A(ok) .* W(ok)) / sum(W(ok)); end
end
end

function y = nice_ceil(x)
if ~isfinite(x) || x <= 0; y = 1; return; end
e = 10^floor(log10(x));  m = x / e;
if m <= 1; y = e; elseif m <= 2; y = 2*e; elseif m <= 5; y = 5*e; else; y = 10*e; end
end

function out = ternary(c, a, b)
if c; out = a; else; out = b; end
end
