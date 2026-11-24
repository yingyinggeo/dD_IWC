%% ========================================================================
%  ERA5 SON ω 距平（400–600 hPa 平均）对 IOD 指数的回归：只画原图的 (c)，1 个面板
%   去除 ENSO 的 SON ω 距平 对 去除 ENSO 的 SON DMI 回归
%   （(a)(b) 仍会计算并保存到 .mat，但不画图）
%       —— ω 与 DMI 各自先对 D(0)JF(1) Niño3.4 做线性回归，取残差（自由度 n−3）
%   指数由 HadISST 计算（相对研究时段逐月气候态的距平，cos(lat) 加权区域平均）：
%     Niño3.4 = SSTA(170°W–120°W, 5°S–5°N)，取当年 12 月 + 次年 1、2 月
%     DMI     = SSTA(50–70°E, 10°S–10°N) − SSTA(90–110°E, 10°S–0°)，取当年 SON
%   ω：ERA5（同 Fig. 2），×864 转为 hPa day^-1，1° → 5° 块平均，
%      相对研究时段逐月气候态求距平，400–600 hPa 平均，再做 SON 平均
%   回归：指数标准化，系数单位 hPa day^-1 per SD；双侧 t 检验 p<0.05 → 绿点
%   版式：地图范围 30–120°E, 40°S–40°N（同 Fig. 2），陆地浅灰
%% ========================================================================
clear; clc; close all;

%% ===== 基础设置 =====
lon_range = [30 120];
lat_range = [-40 40];
yr_start  = 2003;
yr_end    = 2020;            % Fig. 2 的 ω 趋势用的是 2019，如需一致改为 2019
reg_years = (yr_start:yr_end)';
n_reg     = numel(reg_years);
son       = [9 10 11];
exclude_years = [];          % 不参与回归的年份，例如 [2003]

detrend_flag = false;        % 同观测 δD 回归：不去趋势；改为 true 则对 ω 和指数都去线性趋势
clim_reg     = [];           % 色标 ±（hPa day^-1 per SD）；留空则三个面板统一按 98% 分位自动取整
show_boxes   = true;         % 是否画 WIO / EIO 框（Fig. 1c）
fs           = 9;

nc_w    = 'C:\Users\yingy\Downloads\Rice_mac\ERA5\era5_omega_90N90S_1deg_2000_2024.nc';
nc_sst  = 'C:\Users\yingy\Downloads\Rice_mac\HadISST_sst.nc';
out_dir = 'C:\Users\yingy\Downloads\Rice_mac';

boxes = { ...
    [50 70],  [-15 5], [1.00 0.50 0.00];   % WIO，橙色
    [70 100], [-20 0], [0.00 0.45 0.90]};  % EIO，蓝色

%% ===== 1. HadISST → 逐月 SSTA → Niño3.4 (DJF) 与 DMI (SON) =====
fprintf('读取 HadISST...\n');
sst_lat = double(ncread(nc_sst, 'latitude'));
sst_lon = double(ncread(nc_sst, 'longitude'));
sst_t   = ncread(nc_sst, 'time');
sst     = ncread(nc_sst, 'sst');
sst_dates  = datetime(1870,1,1) + days(double(sst_t));
sst_years  = year(sst_dates);
sst_months = month(sst_dates);
sst = permute(sst, [2 1 3]);
if sst_lat(1) > sst_lat(end); sst_lat = flipud(sst_lat); sst = flipud(sst); end
sst_lon = mod(sst_lon, 360);
[sst_lon, si] = sort(sst_lon);  sst = sst(:, si, :);
sst(sst < -2) = NaN;

keep = sst_years >= yr_start - 1 & sst_years <= yr_end + 1;
sst = double(sst(:,:,keep));  sst_years = sst_years(keep);  sst_months = sst_months(keep);

inper = sst_years >= yr_start & sst_years <= yr_end;
for m = 1:12
    km = sst_months == m;
    sst(:,:,km) = sst(:,:,km) - mean(sst(:,:,km & inper), 3, 'omitnan');
end

n34_m = box_mean_sst(sst, sst_lon, sst_lat, [190 240], [-5 5]);
dmi_m = box_mean_sst(sst, sst_lon, sst_lat, [50 70],   [-10 10]) - ...
        box_mean_sst(sst, sst_lon, sst_lat, [90 110],  [-10 0]);
clear sst

NINO_djf = NaN(n_reg, 1);   DMI_son = NaN(n_reg, 1);
for iy = 1:n_reg
    y  = reg_years(iy);
    kd = (sst_years == y & sst_months == 12) | (sst_years == y+1 & sst_months <= 2);
    ks =  sst_years == y & ismember(sst_months, son);
    if sum(kd) == 3; NINO_djf(iy) = mean(n34_m(kd)); end
    if sum(ks) == 3; DMI_son(iy)  = mean(dmi_m(ks)); end
end
if any(isnan(NINO_djf))
    fprintf('  缺 DJF 的年份：%s（HadISST 未覆盖到次年 1–2 月，这些年份不参与回归）\n', ...
        mat2str(reg_years(isnan(NINO_djf))'));
end

%% ===== 2. ERA5 ω → 400–600 hPa → 5° → 逐月距平 → SON 平均 =====
fprintf('读取 ERA5 omega（只读 400–600 hPa）...\n');
lat  = double(ncread(nc_w, 'latitude'));
lon  = double(ncread(nc_w, 'longitude'));
plev = double(ncread(nc_w, 'pressure_level'));
vt   = ncread(nc_w, 'valid_time');
dates  = datetime(1970,1,1,0,0,0) + seconds(double(vt));
years  = year(dates);  months = month(dates);

kk = find(plev >= 400 & plev <= 600);
fprintf('ω 垂直平均层：%s hPa\n', mat2str(plev(kk)'));
k0 = min(kk);  nk = max(kk) - k0 + 1;
w = double(ncread(nc_w, 'w', [1 1 k0 1], [Inf Inf nk Inf]));   % lon x lat x nk x t
w = w(:,:,kk - k0 + 1,:);
w = squeeze(mean(w, 3, 'omitnan')) * 864;                      % lon x lat x t，hPa day^-1
w = permute(w, [2 1 3]);                                         % lat x lon x t
if lat(1) > lat(end); lat = flipud(lat); w = flipud(w); end
lon = mod(lon, 360);  [lon, li] = sort(lon);  w = w(:, li, :);

% 1° → 5° 块平均（同 Fig. 2）
factor = 5;
[nlat0, nlon0, ntime] = size(w);
nla = floor(nlat0/factor);  nlo = floor(nlon0/factor);
lat = mean(reshape(lat(1:nla*factor), factor, nla), 1)';
lon = mean(reshape(lon(1:nlo*factor), factor, nlo), 1)';
w   = w(1:nla*factor, 1:nlo*factor, :);
w   = squeeze(mean(mean(reshape(w, factor, nla, factor, nlo, ntime), 1, 'omitnan'), 3, 'omitnan'));

la = lat >= lat_range(1) & lat <= lat_range(2);
lo = lon >= lon_range(1) & lon <= lon_range(2);
lat_sub = lat(la);  lon_sub = lon(lo);
w = w(la, lo, :);

inper = years >= yr_start & years <= yr_end;
for m = 1:12
    km = months == m;
    w(:,:,km) = w(:,:,km) - mean(w(:,:,km & inper), 3, 'omitnan');
end
w_son = NaN(numel(lat_sub), numel(lon_sub), n_reg);
for iy = 1:n_reg
    ks = years == reg_years(iy) & ismember(months, son);
    if sum(ks) == 3; w_son(:,:,iy) = mean(w(:,:,ks), 3); end
end
clear w

%% ===== 3. 陆地 mask 与 topo =====
load topo;
topo_lon = 0.5:359.5;  topo_lat = -89.5:89.5;
[TLON, TLAT] = meshgrid(topo_lon, topo_lat);
[OLON, OLAT] = meshgrid(lon_sub, lat_sub);
land_mask = interp2(TLON, TLAT, topo, OLON, OLAT, 'nearest', 0) > 0;
topo_ext = [topo, topo(:,1)];  topo_lon_ext = [topo_lon, 360.5];
buf = 5;
tl = topo_lon_ext >= lon_range(1)-buf & topo_lon_ext <= lon_range(2)+buf;
tt = topo_lat     >= lat_range(1)-buf & topo_lat     <= lat_range(2)+buf;
topo_lon_p = topo_lon_ext(tl);  topo_lat_p = topo_lat(tt);  topo_elev_p = topo_ext(tt, tl);

%% ===== 4. 回归 =====
ok_idx = isfinite(NINO_djf) & isfinite(DMI_son) & ~ismember(reg_years, exclude_years) & ...
         squeeze(any(any(isfinite(w_son), 1), 2));
fprintf('参与回归的年份：%s（共 %d 年）\n', mat2str(reg_years(ok_idx)'), sum(ok_idx));
Y  = w_son(:,:,ok_idx);
xN = NINO_djf(ok_idx);
xD = DMI_son(ok_idx);
if detrend_flag
    xN = detrend(xN);  xD = detrend(xD);  Y = detrend_field(Y);
end
[r_nd, p_nd] = corr(xN, xD);
fprintf('r(D(0)JF(1) Niño3.4, SON DMI) = %.2f (p = %.3f)\n', r_nd, p_nd);

[b_n34, p_n34] = regress_field(Y, xN, 0);
[b_dmi, p_dmi] = regress_field(Y, xD, 0);

% (c) ω 与 DMI 都去除 D(0)JF(1) Niño3.4
xn   = xN - mean(xN);
d    = xD - mean(xD);
xD_r = d - xn * (xn' * d) / (xn' * xn);
Y_r  = remove_index(Y, xn);
fprintf('去除 Niño3.4 后 DMI 保留方差 = %.0f%%\n', 100 * var(xD_r) / var(d));
[b_dmr, p_dmr] = regress_field(Y_r, xD_r, 1);

B = {b_n34, b_dmi, b_dmr};  P = {p_n34, p_dmi, p_dmr};
for k = 1:3; B{k}(land_mask) = NaN;  P{k}(land_mask) = NaN; end

% 两个框内的平均回归系数（ω>0 下沉）
fprintf('\n框内平均回归系数（hPa day^-1 per SD；ω>0 为下沉）：\n');
tags = {'(a) Niño3.4', '(b) DMI', '(c) DMI, ENSO removed'};
for k = 1:3
    wv = box_mean_field(B{k}, lon_sub, lat_sub, boxes{1,1}, boxes{1,2}, land_mask);
    ev = box_mean_field(B{k}, lon_sub, lat_sub, boxes{2,1}, boxes{2,2}, land_mask);
    fprintf('  %-22s WIO = %+7.3f   EIO = %+7.3f\n', tags{k}, wv, ev);
end

save(fullfile(out_dir, 'omega_ERA5_SON_regress_Nino34DJF_DMI.mat'), 'B', 'P', 'lon_sub', 'lat_sub', ...
    'land_mask', 'reg_years', 'NINO_djf', 'DMI_son', 'ok_idx', 'detrend_flag');

%% ===== 5. 画图：只画 (c) 去除 ENSO 的 ω 对去除 ENSO 的 DMI 回归，1 个面板 =====
load BlueWhiteOrangeRed.mat;
silver_rgb = [0.85 0.85 0.85];
data = B{3};  pv = P{3};
cl = clim_reg;
if isempty(cl); cl = nice_ceil(prctile(abs(data(isfinite(data))), 98)); end
fprintf('色标 = ±%g hPa day^-1 per SD\n', cl);

ttl = 'ENSO-removed SON \omega (400–600 hPa) regression on ENSO-removed SON DMI';

cm2px    = 96/2.54;
fig_w_cm = 18;
fig_w    = round(fig_w_cm * cm2px);
aspect   = diff(lon_range) / diff(lat_range);
pw_norm    = 0.55;
panel_h_px = round(pw_norm * fig_w / aspect);
mt_px = round(1.0 * cm2px);  mb_px = round(1.0 * cm2px);
fig_h = mt_px + panel_h_px + mb_px;
ph_norm = panel_h_px / fig_h;
ml_norm = (1 - pw_norm) / 2 - 0.03;
pos = [ml_norm, mb_px/fig_h, pw_norm, ph_norm];

fig = figure('Color','w', 'Position',[50 50 fig_w fig_h], 'Renderer','painters', ...
    'PaperUnits','centimeters', 'PaperSize',[fig_w_cm, fig_h*2.54/cm2px]);
[LONs, LATs] = meshgrid(lon_sub, lat_sub);
land_img = repmat(reshape(silver_rgb,1,1,3), [size(land_mask) 1]);

ax = axes('Position', pos);  hold(ax, 'on');
imagesc(ax, lon_sub, lat_sub, land_img);
imagesc(ax, lon_sub, lat_sub, data, 'AlphaData', ~isnan(data));
set(ax, 'YDir','normal', 'Color', silver_rgb);
colormap(ax, BlueWhiteOrangeRed);
clim(ax, [-cl cl]);
contour(ax, topo_lon_p, topo_lat_p, topo_elev_p, [0 0], 'k-', 'LineWidth', 0.5);
sg = pv < 0.05 & ~isnan(data);
plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none', 'HandleVisibility','off');
if show_boxes
    for ib = 1:size(boxes, 1)
        bx = boxes{ib,1};  by = boxes{ib,2};
        plot(ax, bx([1 2 2 1 1]), by([1 1 2 2 1]), '-', 'Color', boxes{ib,3}, 'LineWidth', 1.2);
    end
end
xlim(ax, lon_range);  ylim(ax, lat_range);
daspect(ax, [1 1 1]);
set(ax, 'XTick', 30:30:120, 'XTickLabel', {'30°E','60°E','90°E','120°E'}, ...
    'YTick', -40:20:40, 'YTickLabel', {'40°S','20°S','0°','20°N','40°N'}, ...
    'FontSize', fs, 'TickDir','out', 'Box','on', 'Layer','top');
cb = colorbar(ax, 'Location','eastoutside');
cb.Position = [pos(1)+pw_norm+0.025, pos(2), 0.02, ph_norm];
cb.Ticks = cl * [-1 -0.5 0 0.5 1];
cb.Label.String = 'hPa day^{-1} / \sigma_{DMI}';
cb.FontSize = fs;  cb.Label.FontSize = fs;  cb.TickDirection = 'out';
ax.Position = pos;
title(ax, ttl, 'FontSize', fs, 'FontWeight','bold', 'HorizontalAlignment','left', ...
    'Units','normalized', 'Position',[0 1.02 0], 'Interpreter','tex');

base = fullfile(out_dir, sprintf('FigS_omega_ERA5_SON_regress_DMI_ENSOremoved_%d_%d', yr_start, yr_end));
exportgraphics(fig, [base '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig, [base '.tif'], 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('Saved: %s.png / .tif\n', base);

%% ===== 子函数 =====
function ts = box_mean_sst(A, lon, lat, blon, blat)
% A: lat x lon x t；cos(lat) 加权区域平均（NaN = 陆地/海冰，自动跳过）
la = lat >= blat(1) & lat <= blat(2);
lo = lon >= blon(1) & lon <= blon(2);
sub = A(la, lo, :);
w = repmat(cosd(lat(la)), 1, sum(lo));
ts = NaN(size(A, 3), 1);
for k = 1:size(A, 3)
    s = sub(:,:,k);  ok = isfinite(s);
    if any(ok(:)); ts(k) = sum(s(ok) .* w(ok)) / sum(w(ok)); end
end
end

function v = box_mean_field(A, lon, lat, blon, blat, land)
% 2D 场在框内的海洋格点、cos(lat) 加权平均
mask = (lat(:) >= blat(1) & lat(:) <= blat(2)) & (lon(:)' >= blon(1) & lon(:)' <= blon(2)) & ~land;
W = repmat(cosd(lat(:)), 1, numel(lon)) .* mask;
ok = W > 0 & isfinite(A);
v = sum(A(ok) .* W(ok)) / sum(W(ok));
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

function y = nice_ceil(x)
if ~isfinite(x) || x <= 0; y = 1; return; end
e = 10^floor(log10(x));  m = x / e;
if m <= 1; y = e; elseif m <= 2; y = 2*e; elseif m <= 5; y = 5*e; else; y = 10*e; end
end
