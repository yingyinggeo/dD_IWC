%% ========================================================================
%  原始 AIRS δD 距平（SON，400–600 hPa）对 ENSO / IOD 指数的回归
%   (a) SON δD 距平 对 D(0)JF(1) Niño3.4 回归
%   (b) SON δD 距平 对 SON DMI 回归
%   (c) SON δD 距平 对 SON DMI 回归（δD 与 DMI 先线性去除 D(0)JF(1) Niño3.4）
%   指数由 HadISST 计算（相对 2003–2020 逐月气候态的距平，海洋格点 cos(lat) 加权）：
%     Niño3.4 = SSTA(170°W–120°W, 5°S–5°N)，取当年 12 月 + 次年 1、2 月
%     DMI     = SSTA(50–70°E, 10°S–10°N) − SSTA(90–110°E, 10°S–0°)，取当年 SON
%   δD：AIRS 原始 δD（非 δD_dyn），相对 2003–2020 逐月气候态的距平，
%       取与 Fig. 3 相同的 350–650 hPa 层（619、511、422 hPa）平均，再做 SON 平均
%   回归：指数标准化，系数单位 ‰ per SD；双侧 t 检验 p<0.05 打绿点
%   版式沿用 Fig. 3 的地图面板（30–120°E, 20°S–40°N）
%% ========================================================================
clear; clc; close all;

%% ===== 基础设置 =====
lon_range = [30 180];
lat_range = [-20 40];
yr_start  = 2003;
yr_end    = 2020;
reg_years = (yr_start:yr_end)';
n_reg     = numel(reg_years);
son       = [9 10 11];
exclude_years = [2003];   % 不参与回归的年份（距平的气候态基准期仍为 yr_start–yr_end）

detrend_flag = false;     % 原 Fig. 3 回归未去趋势；改为 true 则对 δD 和指数都去线性趋势
clim_reg     = 10;        % 色标 ±（‰ per SD），同 Fig. 3
fs           = 9;

nc_sst  = 'C:\Users\yingy\Downloads\Rice_mac\HadISST_sst.nc';
nc_airs = 'C:\Users\yingy\Downloads\Rice_mac\Isotope-IOD\AIRS_RS_monthly_5deg.nc';
out_dir = 'C:\Users\yingy\Downloads\Rice_mac';

%% ===== 1. HadISST → 逐月 SSTA → Niño3.4 (DJF) 与 DMI (SON) =====
fprintf('读取 HadISST...\n');
sst_lat = double(ncread(nc_sst, 'latitude'));
sst_lon = double(ncread(nc_sst, 'longitude'));
sst_t   = ncread(nc_sst, 'time');
sst     = ncread(nc_sst, 'sst');
sst_dates  = datetime(1870,1,1) + days(double(sst_t));
sst_years  = year(sst_dates);
sst_months = month(sst_dates);
sst = permute(sst, [2 1 3]);                         % lat x lon x t
if sst_lat(1) > sst_lat(end); sst_lat = flipud(sst_lat); sst = flipud(sst); end
sst_lon = mod(sst_lon, 360);
[sst_lon, si] = sort(sst_lon);  sst = sst(:, si, :);
sst(sst < -2) = NaN;                                 % 海冰/缺测

% 只保留 2002-12 到 2021-02 附近的月份，减少内存
keep = sst_years >= yr_start - 1 & sst_years <= yr_end + 1;
sst = double(sst(:,:,keep));  sst_years = sst_years(keep);  sst_months = sst_months(keep);

% 逐月气候态（2003–2020）与距平
inper = sst_years >= yr_start & sst_years <= yr_end;
ssta  = NaN(size(sst));
for m = 1:12
    km = sst_months == m;
    ssta(:,:,km) = sst(:,:,km) - mean(sst(:,:,km & inper), 3, 'omitnan');
end
clear sst

n34_m = box_mean_sst(ssta, sst_lon, sst_lat, [190 240], [-5 5]);
wio_m = box_mean_sst(ssta, sst_lon, sst_lat, [50 70],   [-10 10]);
eio_m = box_mean_sst(ssta, sst_lon, sst_lat, [90 110],  [-10 0]);
dmi_m = wio_m - eio_m;
clear ssta

NINO_djf = NaN(n_reg, 1);   DMI_son = NaN(n_reg, 1);
for iy = 1:n_reg
    y  = reg_years(iy);
    kd = (sst_years == y & sst_months == 12) | (sst_years == y+1 & sst_months <= 2);
    ks =  sst_years == y & ismember(sst_months, son);
    if sum(kd) == 3; NINO_djf(iy) = mean(n34_m(kd)); end
    if sum(ks) == 3; DMI_son(iy)  = mean(dmi_m(ks)); end
end
fprintf('Niño3.4 (D(0)JF(1)) 可用年份：%d/%d；DMI (SON)：%d/%d\n', ...
    sum(isfinite(NINO_djf)), n_reg, sum(isfinite(DMI_son)), n_reg);
if any(isnan(NINO_djf))
    fprintf('  缺 DJF 的年份：%s（HadISST 文件可能未覆盖到次年 1–2 月，这些年份不参与回归）\n', ...
        mat2str(reg_years(isnan(NINO_djf))'));
end
ok_idx = isfinite(NINO_djf) & isfinite(DMI_son) & ~ismember(reg_years, exclude_years);
fprintf('参与回归的年份：%d–%d，去掉 %s，共 %d 年\n', yr_start, yr_end, ...
    mat2str(exclude_years), sum(ok_idx));
[r_nd, p_nd] = corr(NINO_djf(ok_idx), DMI_son(ok_idx));
fprintf('r(D(0)JF(1) Niño3.4, SON DMI) = %.2f (p = %.3f)\n', r_nd, p_nd);

%% ===== 2. AIRS 原始 δD → 逐月距平 → 400–600 hPa → SON 平均 =====
fprintf('读取 AIRS dD...\n');
lat = double(ncread(nc_airs, 'lat'));
lon = double(ncread(nc_airs, 'lon'));
dD  = ncread(nc_airs, 'dD');
dD  = permute(dD, [2 1 3 4]);
dD  = flipud(dD);  lat = flipud(lat);
lon = mod(lon, 360);
[lon, li] = sort(lon);  dD = dD(:, li, :, :);
pressure_levels = [1017.2227, 1000, 908.51398, 825.40198, 749.89301, ...
                   681.29102, 618.96600, 510.89801, 421.69800, 348.06900, ...
                   287.29800, 237.13699, 177.82899, 133.35201, 74.989601, ...
                   28.729900, 0.10000];
lev_mid = find(pressure_levels >= 350 & pressure_levels <= 650);
fprintf('δD 垂直平均层：%s hPa\n', mat2str(round(pressure_levels(lev_mid))));

t_airs    = ncread(nc_airs, 'time');
dates     = datetime(2000,1,1) + days(double(t_airs));
years_dD  = year(dates);
months_dD = month(dates);

li = lon >= lon_range(1) & lon <= lon_range(2);
la = lat >= lat_range(1) & lat <= lat_range(2);
lon_sub = lon(li);  lat_sub = lat(la);
nlat = numel(lat_sub);  nlon = numel(lon_sub);
dD_mid = squeeze(mean(double(dD(la, li, lev_mid, :)), 3, 'omitnan'));   % lat x lon x t
clear dD

inper = years_dD >= yr_start & years_dD <= yr_end;
dD_anom = NaN(size(dD_mid));
for m = 1:12
    km = months_dD == m;
    dD_anom(:,:,km) = dD_mid(:,:,km) - mean(dD_mid(:,:,km & inper), 3, 'omitnan');
end

dD_son = NaN(nlat, nlon, n_reg);
for iy = 1:n_reg
    ks = years_dD == reg_years(iy) & ismember(months_dD, son);
    if any(ks); dD_son(:,:,iy) = mean(dD_anom(:,:,ks), 3, 'omitnan'); end
end

%% ===== 3. 陆地 mask 与 topo =====
load topo;
topo_lon_v = 0.5:359.5;  topo_lat_v = -89.5:89.5;
[TLON, TLAT] = meshgrid(topo_lon_v, topo_lat_v);
[OLON, OLAT] = meshgrid(lon_sub, lat_sub);
land_mask = interp2(TLON, TLAT, topo, OLON, OLAT, 'nearest', 0) > 0;
topo_pack.lon  = [topo_lon_v, 360.5];
topo_pack.lat  = topo_lat_v;
topo_pack.elev = [topo, topo(:,1)];

%% ===== 4. 回归 =====
Y  = dD_son(:,:,ok_idx);
xN = NINO_djf(ok_idx);
xD = DMI_son(ok_idx);
if detrend_flag
    xN = detrend(xN);  xD = detrend(xD);  Y = detrend_field(Y);
end

[b_n34, p_n34] = regress_field(Y, xN, 0);
[b_dmi, p_dmi] = regress_field(Y, xD, 0);

% (c) 先从 δD 和 DMI 中线性去除 D(0)JF(1) Niño3.4，再回归
xn  = xN - mean(xN);
d   = xD - mean(xD);
xD_r = d - xn * (xn' * d) / (xn' * xn);
Y_r  = remove_index(Y, xn);
fprintf('去除 Niño3.4 后 DMI 保留方差 = %.0f%%\n', 100 * var(xD_r) / var(d));
[b_dmr, p_dmr] = regress_field(Y_r, xD_r, 1);

B = {b_n34, b_dmi, b_dmr};  P = {p_n34, p_dmi, p_dmr};
for k = 1:3; B{k}(land_mask) = NaN;  P{k}(land_mask) = NaN; end

save(fullfile(out_dir, 'dD_SON_regress_Nino34DJF_DMI.mat'), 'B', 'P', 'lon_sub', 'lat_sub', ...
    'land_mask', 'reg_years', 'NINO_djf', 'DMI_son', 'ok_idx', 'detrend_flag');

%% ===== 5. 画图：3 行 1 列，沿用 Fig. 3 地图面板版式 =====
load BlueWhiteOrangeRed.mat;
titles = {
    '(a) SON \deltaD regression on D(0)JF(1) Niño3.4 (400–600 hPa)';
    '(b) SON \deltaD regression on SON DMI (400–600 hPa)';
    '(c) SON \deltaD regression on SON DMI, Niño3.4 removed (400–600 hPa)'};
cb_labels = {'‰ / \sigma_{Niño3.4}', '‰ / \sigma_{DMI}', '‰ / \sigma_{DMI}'};

cm2px    = 96/2.54;
fig_w_cm = 18;
fig_w    = round(fig_w_cm * cm2px);
aspect   = diff(lon_range) / diff(lat_range);
pw_norm    = 0.70;                 % 经度范围扩大到 180°E，面板加宽以保持高度
panel_h_px = round(pw_norm * fig_w / aspect);
mt_px = round(1.2 * cm2px);  mb_px = round(1.2 * cm2px);  vgap_px = round(1.5 * cm2px);
fig_h = mt_px + 3*panel_h_px + 2*vgap_px + mb_px;
ph_norm = panel_h_px / fig_h;
ml_norm = 0.10;                    % 左边距（右侧留给色标）
pos_y = 1 - mt_px/fig_h - ph_norm - (0:2) * (ph_norm + vgap_px/fig_h);

fig = figure('Color','w', 'Position',[50 50 fig_w fig_h], 'Renderer','painters', ...
    'PaperUnits','centimeters', 'PaperSize',[fig_w_cm, fig_h*2.54/cm2px]);

buf = 5;
tl = topo_pack.lon >= lon_range(1)-buf & topo_pack.lon <= lon_range(2)+buf;
tt = topo_pack.lat >= lat_range(1)-buf & topo_pack.lat <= lat_range(2)+buf;
topo_lon_p = topo_pack.lon(tl);  topo_lat_p = topo_pack.lat(tt);  topo_elev_p = topo_pack.elev(tt, tl);
[LONs, LATs] = meshgrid(lon_sub, lat_sub);

for k = 1:3
    pos = [ml_norm, pos_y(k), pw_norm, ph_norm];
    ax = axes('Position', pos);  hold(ax, 'on');
    data = B{k};
    imagesc(ax, lon_sub, lat_sub, data, 'AlphaData', ~isnan(data));
    set(ax, 'YDir','normal', 'Color',[0.85 0.85 0.85]);
    colormap(ax, BlueWhiteOrangeRed);
    clim(ax, [-clim_reg clim_reg]);
    contour(ax, topo_lon_p, topo_lat_p, topo_elev_p, [0 0], 'k-', 'LineWidth', 0.5);
    sg = P{k} < 0.05 & ~isnan(data);
    plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize', 2.5, ...
        'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none', 'HandleVisibility','off');
    xlim(ax, lon_range);  ylim(ax, lat_range);
    if k == 3
        set(ax, 'XTick', 30:30:180, 'XTickLabel', {'30°E','60°E','90°E','120°E','150°E','180°'});
    else
        set(ax, 'XTick', 30:30:180, 'XTickLabel', {});
    end
    set(ax, 'YTick', -20:20:40, 'YTickLabel', {'20°S','0°','20°N','40°N'}, ...
        'FontSize', fs, 'TickDir','out', 'Box','off', 'Layer','top');
    cb = colorbar(ax, 'Location','eastoutside');
    cb.Position = [ml_norm+pw_norm+0.025, pos_y(k), 0.02, ph_norm];
    cb.Label.String = cb_labels{k};
    cb.FontSize = fs;  cb.Label.FontSize = fs;  cb.TickDirection = 'out';
    ax.Position = pos;
    xl = ax.XLim;  yl = ax.YLim;
    plot(ax, [xl(2) xl(2)], yl, 'k-', 'LineWidth', 0.5, 'HandleVisibility','off');
    plot(ax, xl, [yl(2) yl(2)], 'k-', 'LineWidth', 0.5, 'HandleVisibility','off');
    title(ax, titles{k}, 'FontSize', fs, 'FontWeight','bold', 'HorizontalAlignment','left', ...
        'Units','normalized', 'Position',[0 1.02 0], 'Interpreter','tex');
end

exportgraphics(fig, fullfile(out_dir, 'FigS_dD_SON_regress_Nino34DJF_DMI.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(out_dir, 'FigS_dD_SON_regress_Nino34DJF_DMI.tif'), 'Resolution', 600);
fprintf('Saved: FigS_dD_SON_regress_Nino34DJF_DMI.png / .tif\n');

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

function [b, p] = regress_field(Y, x, k_extra)
% Y: lat x lon x n；x 标准化后回归 → 每 1 个标准差的变化；k_extra 为额外估计的参数个数
x = (x - mean(x)) / std(x);
[d1, d2, n] = size(Y);
b = NaN(d1, d2);  p = NaN(d1, d2);
for i = 1:d1
    for j = 1:d2
        y = squeeze(Y(i,j,:));
        v = isfinite(y);
        nv = sum(v);
        if nv < 5 + k_extra; continue; end
        xv = x(v) - mean(x(v));  yv = y(v) - mean(y(v));
        Sxx = sum(xv.^2);
        bb  = sum(xv .* yv) / Sxx;
        res = yv - bb * xv;
        dof = nv - 2 - k_extra;
        se  = sqrt(sum(res.^2) / dof / Sxx);
        b(i,j) = bb;
        p(i,j) = 2 * (1 - tcdf(abs(bb / se), dof));
    end
end
end

function Fr = remove_index(F, xn)
% 逐格点线性去除指数 xn（已去均值）的影响
[d1, d2, n] = size(F);
Fr = NaN(d1, d2, n);
for i = 1:d1
    for j = 1:d2
        y = squeeze(F(i,j,:));
        v = isfinite(y);
        if sum(v) < 4; continue; end
        yv = y(v) - mean(y(v));  xv = xn(v) - mean(xn(v));
        r  = NaN(n, 1);
        r(v) = yv - xv * (xv' * yv) / (xv' * xv);
        Fr(i,j,:) = r;
    end
end
end

function Yd = detrend_field(Y)
[d1, d2, n] = size(Y);
Yd = NaN(d1, d2, n);
for i = 1:d1
    for j = 1:d2
        y = squeeze(Y(i,j,:));
        if all(isfinite(y)); Yd(i,j,:) = detrend(y); end
    end
end
end
