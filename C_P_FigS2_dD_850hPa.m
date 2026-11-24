% ============================================================
%  FigS2：AIRS 原始 dD 趋势空间图（不做 Rayleigh 分解），850–700 hPa 平均
%  使用 AIRS 层：825.4 + 749.9 hPa
%  画 MAM + JJA + SON + DJF 共 4 个 panel（2 行 × 2 列）
%  显著性水平 p < 0.1，字号统一 10
% ============================================================
clear; clc; close all;

sig_level = 0.1;   % 显著性水平
fs        = 10;    % 统一字号

%% ===== 0. 读取 AIRS dD =====
nc_file = 'C:\Users\yingy\Downloads\Rice_mac\Isotope-IOD\AIRS_RS_monthly_5deg.nc';
lat      = ncread(nc_file, 'lat');
lon      = ncread(nc_file, 'lon');
dD       = ncread(nc_file, 'dD');
dD_per   = permute(dD, [2 1 3 4]);
dD_fn    = flipud(dD_per);
lat_fn   = flipud(lat);
lon_fn   = lon;
clear lat lon dD dD_per;
lat = lat_fn; lon = lon_fn; dD = dD_fn;
clear lat_fn lon_fn dD_fn;

% 经度转 [0, 360]
lon_360 = lon;
lon_360(lon_360 < 0) = lon_360(lon_360 < 0) + 360;
[lon_sorted, sort_idx] = sort(lon_360);
dD  = dD(:, sort_idx, :, :);
lon = lon_sorted;
[~, ~, nlev, ntime] = size(dD);

pressure_levels = [1017.2227, 1000, 908.51398, 825.40198, 749.89301, ...
                   681.29102, 618.96600, 510.89801, 421.69800, 348.06900, ...
                   287.29800, 237.13699, 177.82899, 133.35201, 74.989601, ...
                   28.729900, 0.10000];
dates  = datetime(2002, 9, 1) + calmonths(0:ntime-1);
years  = year(dates);
months = month(dates);

%% 区域: 印度洋
lon_range = [30 120];
lat_range = [-40 40];
lon_idx = lon >= lon_range(1) & lon <= lon_range(2);
lat_idx = lat >= lat_range(1) & lat <= lat_range(2);
lon_sub = lon(lon_idx);
lat_sub = lat(lat_idx);
nlat    = sum(lat_idx);
nlon    = sum(lon_idx);
dD_sub  = dD(lat_idx, lon_idx, :, :);
clear dD;
fprintf('AIRS dD 加载: %d lat × %d lon × %d lev × %d time\n', nlat, nlon, nlev, ntime);

%% ===== 1. 计算 dD 趋势：4 季节 =====
% 第 3 列 cross_year：DJF 跨年（Dec 算当年，Jan/Feb 算次年）
seasons = {
    'MAM',    [3 4 5],  false;
    'JJA',    [6 7 8],  false;
    'SON',    [9 10 11],false;
    'DJF',    [12 1 2], true;
};
yr_start = 2003;
yr_end   = 2020;
all_years = (yr_start:yr_end)';
n_years   = numel(all_years);

% 层选择（850–700 hPa）：AIRS 825.4 + 749.9 hPa
lev_mid = find(pressure_levels >= 700 & pressure_levels <= 850);
fprintf('使用 AIRS 层: %s hPa\n', mat2str(round(pressure_levels(lev_mid),1)));

% topo + 陆地 mask
load topo;
topo_lon     = 0.5:359.5;
topo_lat     = -89.5:89.5;
topo_ext     = [topo, topo(:,1)];
topo_lon_ext = [topo_lon, 360.5];
[TLON, TLAT] = meshgrid(topo_lon, topo_lat);
[OLON, OLAT] = meshgrid(lon_sub, lat_sub);
land_mask = interp2(TLON, TLAT, topo, OLON, OLAT, 'nearest', 0) > 0;
topo_pack.lon  = topo_lon_ext;
topo_pack.lat  = topo_lat;
topo_pack.elev = topo_ext;
load BlueWhiteOrangeRed.mat;

% 存结果
n_panels  = size(seasons, 1);
trend_mid = NaN(nlat, nlon, n_panels);
pval_mid  = NaN(nlat, nlon, n_panels);

fprintf('\n=== 计算 dD 趋势 (4 季节) ===\n');
for is = 1:n_panels
    s_name   = seasons{is, 1};
    s_months = seasons{is, 2};
    s_cross  = seasons{is, 3};
    fprintf('  时段: %s\n', s_name);
    mid_samples = NaN(nlat, nlon, n_years);
    for iy = 1:n_years
        yr = all_years(iy);
        mid_samples(:,:,iy) = extract_season_mean(dD_sub, years, months, ...
            yr, s_months, s_cross, lev_mid);
    end
    [trend_mid(:,:,is), pval_mid(:,:,is)] = calc_trend_map(mid_samples, all_years, land_mask);
end
fprintf('dD 趋势计算完成\n');

%% ===== 2. 画图 =====
clim_val = 1.5;   % 色标范围（‰/yr）— 视实际情况调

% 低层 850–700 hPa
fig_mid = make_trend_figure(trend_mid, pval_mid, seasons, ...
    lon_sub, lat_sub, lon_range, lat_range, ...
    BlueWhiteOrangeRed, topo_pack, land_mask, clim_val, sig_level, fs, ...
    'δD Trend (850–700 hPa, 2003–2020)');
exportgraphics(fig_mid, 'FigS2.png', 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig_mid, 'FigS2.tif', 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('Saved: FigS2.png / FigS2.tif（当前 MATLAB 工作目录：%s）\n', pwd);

disp('=== All dD trend figures saved. ===');

%% ============================================================
%% Local functions
%% ============================================================
function data_2d = extract_season_mean(field, yrs, mos, ...
                                        yr, season_months, cross_year, lev_idx)
% DJF (cross_year = true) 时, Dec 取当年, Jan/Feb 取次年
    [nlat, nlon, ~, ~] = size(field);
    n_mo = numel(season_months);
    season_3d = NaN(nlat, nlon, n_mo);
    for im = 1:n_mo
        mo = season_months(im);
        if cross_year && mo == 12
            target_yr = yr;
        elseif cross_year && mo < 6
            target_yr = yr + 1;
        else
            target_yr = yr;
        end
        idx = find((yrs == target_yr) & (mos == mo));
        if isempty(idx); continue; end
        season_3d(:,:,im) = mean(field(:,:,lev_idx,idx(1)), 3, 'omitnan');
    end
    data_2d = mean(season_3d, 3, 'omitnan');
end

%% ============================================================
function draw_box(ax, lon_lim, lat_lim, col)
% 画矩形区域框
    plot(ax, [lon_lim(1) lon_lim(2) lon_lim(2) lon_lim(1) lon_lim(1)], ...
             [lat_lim(1) lat_lim(1) lat_lim(2) lat_lim(2) lat_lim(1)], ...
             '-', 'Color', col, 'LineWidth', 1.8, 'HandleVisibility', 'off');
end

%% ============================================================
function [trend_map, pval_map] = calc_trend_map(samples, all_years, land_mask)
    [nlat, nlon, ~] = size(samples);
    trend_map = NaN(nlat, nlon);
    pval_map  = NaN(nlat, nlon);
    x = double(all_years);
    for ilat = 1:nlat
        for ilon = 1:nlon
            if land_mask(ilat, ilon); continue; end
            y = squeeze(samples(ilat, ilon, :));
            valid = ~isnan(y);
            if sum(valid) < 4; continue; end
            x_v = x(valid);
            y_v = y(valid);
            n   = sum(valid);
            p     = polyfit(x_v, y_v, 1);
            yhat  = polyval(p, x_v);
            resid = y_v - yhat;
            s2    = sum(resid.^2) / (n - 2);
            Sxx   = sum((x_v - mean(x_v)).^2);
            se    = sqrt(s2 / Sxx);
            t_val = p(1) / se;
            trend_map(ilat, ilon) = p(1);
            pval_map(ilat, ilon)  = 2 * (1 - tcdf(abs(t_val), n - 2));
        end
    end
end

%% ============================================================
function fig = make_trend_figure(trend_4d, pval_4d, seasons, ...
                                  lon_sub, lat_sub, lon_range, lat_range, ...
                                  cmap, topo_pack, land_mask, clim_val, ...
                                  sig_level, fs, sup_title)
    n_panels = size(seasons, 1);
    nrow = 2; ncol = 2;
    panel_labels = {'(a)','(b)','(c)','(d)'};

    lon_span = lon_range(2) - lon_range(1);
    lat_span = lat_range(2) - lat_range(1);
    aspect   = lon_span / lat_span;

    panel_h_px = 180;
    panel_w_px = round(panel_h_px * aspect);
    gap_x_px = 55;
    gap_y_px = 30;
    mt_px  = 60;
    mb_px  = 95;
    ml_px  = 60;
    mr_px  = 30;
    fig_w = ml_px + ncol*panel_w_px + (ncol-1)*gap_x_px + mr_px;
    fig_h = mt_px + nrow*panel_h_px + (nrow-1)*gap_y_px + mb_px;
    fig = figure('Color', 'w', 'Position', [50 50 fig_w fig_h]);

    ml   = ml_px / fig_w;
    mt   = mt_px / fig_h;
    gx   = gap_x_px / fig_w;
    gy   = gap_y_px / fig_h;
    pw_n = panel_w_px / fig_w;
    ph_n = panel_h_px / fig_h;

    buf = 5;
    tlon_idx = topo_pack.lon >= lon_range(1)-buf & topo_pack.lon <= lon_range(2)+buf;
    tlat_idx = topo_pack.lat >= lat_range(1)-buf & topo_pack.lat <= lat_range(2)+buf;
    topo_lon_v  = topo_pack.lon(tlon_idx);
    topo_lat_v  = topo_pack.lat(tlat_idx);
    topo_elev_v = topo_pack.elev(tlat_idx, tlon_idx);

    last_ax = [];
    for ip = 1:n_panels
        s_name = seasons{ip, 1};
        ir = ceil(ip / ncol);
        ic = mod(ip-1, ncol) + 1;
        pos_x = ml + (ic-1) * (pw_n + gx);
        pos_y = 1 - mt - ir * ph_n - (ir-1) * gy;
        ax = axes('Position', [pos_x, pos_y, pw_n, ph_n]);
        hold on;

        data_2d = trend_4d(:,:,ip);
        pval_2d = pval_4d(:,:,ip);
        data_2d(land_mask) = NaN;
        imagesc(ax, lon_sub, lat_sub, data_2d, 'AlphaData', ~isnan(data_2d));
        set(ax, 'YDir', 'normal', 'Color', [0.80 0.82 0.84]);   % 冷调银灰
        colormap(ax, cmap);
        clim(ax, [-clim_val clim_val]);
        contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], ...
            'k-', 'LineWidth', 0.5);
        set(ax, 'Clipping', 'on');

        % 区域框（仅图 c）：橙色 (50–70°E, 15°S–5°N)，蓝色 (70–100°E, 20°S–0°)
        if ip == 3
            draw_box(ax, [50 70],  [-15 5], [0.93 0.49 0.13]);
            draw_box(ax, [70 100], [-20 0], [0.13 0.40 0.75]);
        end

        % 显著性点（p < sig_level）
        for r = 1:numel(lat_sub)
            for c = 1:numel(lon_sub)
                if land_mask(r,c); continue; end
                if isnan(pval_2d(r,c)); continue; end
                if pval_2d(r,c) < sig_level
                    plot(ax, lon_sub(c), lat_sub(r), 'o', ...
                        'MarkerSize', 3, 'MarkerFaceColor', [0.1 0.7 0.1], ...
                        'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
                end
            end
        end

        xlim(ax, lon_range);
        ylim(ax, lat_range);
        set(ax, 'DataAspectRatio', [1 1 1], 'FontSize', fs);

        % x 轴标签只在最下一行
        if ir == nrow
            set(ax, 'XTick', 30:30:120, ...
                'XTickLabel', {'30°E','60°E','90°E','120°E'});
        else
            set(ax, 'XTick', 30:30:120, 'XTickLabel', {});
        end
        % y 轴标签只在最左一列
        if ic == 1
            set(ax, 'YTick', -40:20:40, ...
                'YTickLabel', {'40°S','20°S','0°','20°N','40°N'});
        else
            set(ax, 'YTick', -40:20:40, 'YTickLabel', {});
        end
        set(ax, 'TickDir', 'out', 'TickLength', [0.025 0.025], ...
            'Box', 'off', 'Layer', 'top');

        xl = ax.XLim; yl = ax.YLim;
        plot(ax, [xl(2) xl(2)], yl, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        plot(ax, xl, [yl(2) yl(2)], 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');

        title(ax, sprintf('%s %s', panel_labels{ip}, s_name), ...
            'FontSize', fs, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0, 1.02, 0]);
        last_ax = ax;
    end

    % 共用 colorbar（横跨两列）
    cb = colorbar(last_ax, 'Position', [ml, 0.12, ncol*pw_n + gx, 0.025], ...
        'Orientation', 'horizontal');
    cb.Label.String   = 'δD trend (‰ yr^{-1})';
    cb.FontSize       = fs;
    cb.Label.FontSize = fs;
    cb.TickDirection  = 'out';

    annotation(fig, 'textbox', [0, 1 - mt, 1, mt], ...
        'String', sup_title, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', fs, 'FontWeight', 'bold', ...
        'Interpreter', 'tex');
end
