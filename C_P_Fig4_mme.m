%% ================================================================
%% SON dD trend over the Indian Ocean - 7 models + MME (4行2列)
%%   数据: WISOMIP_AK 下 regrid5deg_<model>.nc（5°×5°）
%%   使用原始 dD（不做 Rayleigh 分解）
%% ================================================================
clear; clc; close all;

%% ===== Configuration =====
data_dir = 'C:\Users\yingy\Downloads\Rice_mac\WISOMIP_AK\';
models   = {'CAM5', 'CAM6', 'ECHAM', 'GISS', 'GSM', 'LMDZ', 'MIROC'};
n_models = numel(models);

out_dir = 'C:\Users\yingy\Downloads\Rice_mac\ERA5';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

var_dD = 'dD';                  % dD 变量名（如不同请修改）

lon_range = [30  120];
lat_range = [-20  40];

yr_start    = 2003;
yr_end      = 2019;
press_range = [400 600];        % hPa
clim_val    = 2;
sig_level   = 0.1;              % 显著性水平
fs          = 10;               % 统一字号

season_months = [9 10 11];
season_cross  = false;

all_years = (yr_start:yr_end)';
n_years   = numel(all_years);

silver = [0.80 0.82 0.84];      % 冷调银灰

%% ===== Loop: read each model, compute SON yearly dD =====
samples_all = [];
lon_t = []; lat_t = [];

for im = 1:n_models
    model_name = models{im};
    nc_file = fullfile(data_dir, sprintf('regrid5deg_%s.nc', model_name));
    fprintf('\n========== %s ==========\n', model_name);

    if ~exist(nc_file, 'file')
        fprintf('  *** 文件不存在，跳过: %s\n', nc_file);
        continue;
    end

    info  = ncinfo(nc_file);
    vnames = {info.Variables.Name};

    % ---- 自动识别坐标变量名 ----
    lon_name  = pick_name(vnames, {'lon','longitude'});
    lat_name  = pick_name(vnames, {'lat','latitude'});
    lev_name  = pick_name(vnames, {'p','lev','plev','lev_p','levels','presnivs','pressure','level'});
    time_name = pick_name(vnames, {'time','time_counter'});

    lon = double(ncread(nc_file, lon_name));
    lat = double(ncread(nc_file, lat_name));
    p   = double(ncread(nc_file, lev_name));
    if max(p) > 2000; p = p / 100; end          % Pa -> hPa
    dD  = double(ncread(nc_file, var_dD));      % (lon, lat, lev, time)
    dD(dD <= -998) = NaN;

    ntime = size(dD, 4);
    try
        time_raw   = ncread(nc_file, time_name);
        time_units = ncreadatt(nc_file, time_name, 'units');
        dates      = parse_nc_time(time_raw, time_units);
    catch ME
        warning('%s: time 读取失败 (%s)，使用默认起始 1979-01。', ...
            model_name, ME.message);
        dates = datetime(1979,1,1) + calmonths(0:ntime-1);
        dates = dates(:);
    end
    yrs = year(dates);
    mos = month(dates);

    % 经度转 [0, 360]，纬度升序
    lon(lon < 0) = lon(lon < 0) + 360;
    [lon, ilon_sort] = sort(lon);
    dD = permute(dD(ilon_sort, :, :, :), [2 1 3 4]);   % -> (lat, lon, lev, time)
    if lat(1) > lat(end)
        lat = flipud(lat);
        dD  = flipud(dD);
    end

    lon_idx = lon >= lon_range(1) & lon <= lon_range(2);
    lat_idx = lat >= lat_range(1) & lat <= lat_range(2);
    lon_sub = lon(lon_idx);
    lat_sub = lat(lat_idx);
    dD_sub  = dD(lat_idx, lon_idx, :, :);
    clear dD;

    lev_mid = find(p >= press_range(1) & p <= press_range(2));
    fprintf('  选用层 (hPa): %s\n', mat2str(round(p(lev_mid)', 1)));

    samples = NaN(numel(lat_sub), numel(lon_sub), n_years);
    for iy = 1:n_years
        samples(:,:,iy) = extract_season_mean(dD_sub, yrs, mos, ...
            all_years(iy), season_months, season_cross, lev_mid);
    end
    clear dD_sub;

    % 以第一个模式的网格为目标网格；其他模式若网格不同则插值
    if isempty(lon_t)
        lon_t = lon_sub; lat_t = lat_sub;
        samples_all = NaN(numel(lat_t), numel(lon_t), n_years, n_models);
    end
    if ~isequal(lon_sub, lon_t) || ~isequal(lat_sub, lat_t)
        fprintf('  网格与第一个模式不同，插值到目标网格\n');
        [OLON, OLAT] = meshgrid(lon_sub, lat_sub);
        [TLON, TLAT] = meshgrid(lon_t, lat_t);
        for iy = 1:n_years
            samples_all(:,:,iy,im) = interp2(OLON, OLAT, samples(:,:,iy), ...
                TLON, TLAT, 'linear');
        end
    else
        samples_all(:,:,:,im) = samples;
    end
end

nlat_t = numel(lat_t);
nlon_t = numel(lon_t);

%% ===== Land mask & topo =====
load topo;
topo_lon     = 0.5:359.5;
topo_lat     = -89.5:89.5;
[TLON_topo, TLAT_topo] = meshgrid(topo_lon, topo_lat);
[TLON_t, TLAT_t]       = meshgrid(lon_t, lat_t);
land_mask_t = interp2(TLON_topo, TLAT_topo, topo, ...
                      TLON_t, TLAT_t, 'nearest', 0) > 0;

topo_ext     = [topo, topo(:,1)];
topo_lon_ext = [topo_lon, 360.5];
buf = 5;
tlon_idx    = topo_lon_ext >= lon_range(1)-buf & topo_lon_ext <= lon_range(2)+buf;
tlat_idx    = topo_lat     >= lat_range(1)-buf & topo_lat     <= lat_range(2)+buf;
topo_lon_v  = topo_lon_ext(tlon_idx);
topo_lat_v  = topo_lat(tlat_idx);
topo_elev_v = topo_ext(tlat_idx, tlon_idx);

load BlueWhiteOrangeRed.mat;

%% ===== Trends: 7 models + MME =====
fprintf('\n========== Computing trends ==========\n');
n_panels  = n_models + 1;
trend_all = NaN(nlat_t, nlon_t, n_panels);
pval_all  = NaN(nlat_t, nlon_t, n_panels);

for im = 1:n_models
    [trend_all(:,:,im), pval_all(:,:,im)] = ...
        calc_trend_map(samples_all(:,:,:,im), all_years);
end

mme_avg = mean(samples_all, 4, 'omitnan');
[trend_all(:,:,end), pval_all(:,:,end)] = calc_trend_map(mme_avg, all_years);

panel_names = [models, {'MME'}];
fprintf('Trend 计算完成\n');

%% ===== Dipole 判定 =====
% 定义：橙色框区域平均 δD 趋势 > 0，且蓝色框区域平均 δD 趋势 < 0
% 区域平均：cos(lat) 面积加权，仅海洋格点，格点中心严格位于框内
box_W = struct('lon', [50 70],  'lat', [-15 5]);   % 橙色框
box_E = struct('lon', [70 100], 'lat', [-20 0]);   % 蓝色框

trend_W    = NaN(n_panels, 1);
trend_E    = NaN(n_panels, 1);
trend_diff = NaN(n_panels, 1);
pval_diff  = NaN(n_panels, 1);
for ip = 1:n_panels
    if ip <= n_models
        smp = samples_all(:,:,:,ip);
    else
        smp = mme_avg;
    end
    % 每年 SON 的区域平均及其差值序列（即 IODDI）
    W_ts = NaN(n_years, 1);
    E_ts = NaN(n_years, 1);
    for iy = 1:n_years
        W_ts(iy) = box_mean(smp(:,:,iy), lon_t, lat_t, box_W, land_mask_t);
        E_ts(iy) = box_mean(smp(:,:,iy), lon_t, lat_t, box_E, land_mask_t);
    end
    trend_W(ip) = lin_trend(all_years, W_ts);
    trend_E(ip) = lin_trend(all_years, E_ts);
    [trend_diff(ip), pval_diff(ip)] = lin_trend(all_years, W_ts - E_ts);
end
is_dipole = (trend_W > 0) & (trend_E < 0);

% 打印各框实际选中的格点中心，便于核对
fprintf('\n橙色框选中 lon: %s\n', mat2str(lon_t(lon_t > box_W.lon(1) & lon_t < box_W.lon(2))'));
fprintf('橙色框选中 lat: %s\n',   mat2str(lat_t(lat_t > box_W.lat(1) & lat_t < box_W.lat(2))'));
fprintf('蓝色框选中 lon: %s\n',   mat2str(lon_t(lon_t > box_E.lon(1) & lon_t < box_E.lon(2))'));
fprintf('蓝色框选中 lat: %s\n',   mat2str(lat_t(lat_t > box_E.lat(1) & lat_t < box_E.lat(2))'));

dipole_table = table(panel_names(:), round(trend_W, 3), round(trend_E, 3), ...
    round(trend_diff, 3), round(pval_diff, 3), is_dipole, ...
    'VariableNames', {'Model', 'Trend_orange', 'Trend_blue', ...
                      'Trend_orange_minus_blue', 'p_diff', 'Dipole'});
fprintf('\n===== Dipole 判定 (SON, 橙 > 0 且 蓝 < 0；趋势单位 ‰/yr) =====\n');
disp(dipole_table);
writetable(dipole_table, fullfile(out_dir, 'dipole_check_SON_7models_MME.csv'));

%% ================================================================
%% Figure: 7 models + MME, 4行2列
%% ================================================================
panel_labels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
nrow = 4; ncol = 2;

aspect = (lon_range(2) - lon_range(1)) / (lat_range(2) - lat_range(1));

panel_h_px = 110;
panel_w_px = round(panel_h_px * aspect);
gap_x_px   = 55;
gap_y_px   = 30;
mt_px      = 60;
mb_px      = 95;
ml_px      = 70;
mr_px      = 25;

fig_w = ml_px + ncol*panel_w_px + (ncol-1)*gap_x_px + mr_px;
fig_h = mt_px + nrow*panel_h_px + (nrow-1)*gap_y_px + mb_px;

fig1 = figure('Color','w','Position',[50 50 fig_w fig_h]);
set(fig1, 'DefaultAxesFontSize', fs, 'DefaultTextFontSize', fs);

ml_n = ml_px / fig_w;
mt_n = mt_px / fig_h;
gx_n = gap_x_px / fig_w;
gy_n = gap_y_px / fig_h;
pw_n = panel_w_px / fig_w;
ph_n = panel_h_px / fig_h;

last_ax = [];
for ip = 1:n_panels
    ir = ceil(ip / ncol);
    ic = mod(ip-1, ncol) + 1;
    pos_x = ml_n + (ic-1) * (pw_n + gx_n);
    pos_y = 1 - mt_n - ir * ph_n - (ir-1) * gy_n;

    ax = axes('Position', [pos_x, pos_y, pw_n, ph_n]);
    hold on;

    data_2d = trend_all(:,:,ip);
    pval_2d = pval_all(:,:,ip);
    data_2d(land_mask_t) = NaN;

    h = imagesc(ax, lon_t, lat_t, data_2d);
    h.AlphaData = ~isnan(data_2d);
    set(ax, 'YDir','normal', 'Color', silver);
    colormap(ax, BlueWhiteOrangeRed);
    clim(ax, [-clim_val clim_val]);

    contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], ...
        'k-', 'LineWidth', 0.5);
    set(ax, 'Clipping','on');

    % 区域框：橙色 (50–70°E, 15°S–5°N)，蓝色 (70–100°E, 20°S–0°)
    draw_box(ax, [50 70],  [-15 5], [0.93 0.49 0.13]);
    draw_box(ax, [70 100], [-20 0], [0.13 0.40 0.75]);

    % 显著性点（仅海洋，p < sig_level）
    for r = 1:nlat_t
        for c = 1:nlon_t
            if land_mask_t(r,c) || isnan(pval_2d(r,c)); continue; end
            if pval_2d(r,c) < sig_level
                plot(ax, lon_t(c), lat_t(r), 'o', ...
                    'MarkerSize', 2.5, 'MarkerFaceColor', [0.1 0.7 0.1], ...
                    'MarkerEdgeColor', 'none', 'HandleVisibility','off');
            end
        end
    end

    xlim(ax, lon_range); ylim(ax, lat_range);
    set(ax, 'DataAspectRatio', [1 1 1], 'FontSize', fs);

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
    set(ax, 'TickDir','out', 'TickLength', [0.025 0.025], ...
        'Box','off', 'Layer','top');

    xl = ax.XLim; yl = ax.YLim;
    plot(ax, [xl(2) xl(2)], yl, 'k-', 'LineWidth', 0.5, 'HandleVisibility','off');
    plot(ax, xl, [yl(2) yl(2)], 'k-', 'LineWidth', 0.5, 'HandleVisibility','off');

    title(ax, sprintf('%s %s', panel_labels{ip}, panel_names{ip}), ...
        'FontSize', fs, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'Units', 'normalized', 'Position', [0, 1.02, 0]);

    last_ax = ax;
end

% 共用 colorbar
cb1 = colorbar(last_ax, 'Position', [ml_n, 0.04, ncol*pw_n + gx_n, 0.013], ...
    'Orientation', 'horizontal');
cb1.Label.String   = 'δD trend (‰ yr^{-1})';
cb1.FontSize       = fs;
cb1.Label.FontSize = fs;
cb1.TickDirection  = 'out';

% 总标题
annotation(fig1, 'textbox', [0, 1 - mt_n, 1, mt_n], ...
    'String', sprintf('SON δD Trend (400–600 hPa, %d–%d)', yr_start, yr_end), ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', fs, 'FontWeight', 'bold', 'Interpreter', 'tex');

fname1 = fullfile(out_dir, 'dD_trend_SON_7models_MME.png');
exportgraphics(fig1, fname1, 'Resolution', 600);
fprintf('图已保存: %s\n', fname1);

%% ================================================================
%% Helper functions
%% ================================================================
function draw_box(ax, lon_lim, lat_lim, col)
% 画矩形区域框
    plot(ax, [lon_lim(1) lon_lim(2) lon_lim(2) lon_lim(1) lon_lim(1)], ...
             [lat_lim(1) lat_lim(1) lat_lim(2) lat_lim(2) lat_lim(1)], ...
             '-', 'Color', col, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

function val = box_mean(map2d, lon_v, lat_v, box, land_mask)
% 区域平均：cos(lat) 面积加权，仅海洋格点，忽略 NaN
% 只选格点中心严格位于框内的格点（不含落在边界线上的格点）
    ilon = lon_v > box.lon(1) & lon_v < box.lon(2);
    ilat = lat_v > box.lat(1) & lat_v < box.lat(2);
    sub  = map2d(ilat, ilon);
    msk  = land_mask(ilat, ilon);
    w    = repmat(cosd(lat_v(ilat)), 1, sum(ilon));
    valid = ~isnan(sub) & ~msk;
    if ~any(valid(:)); val = NaN; return; end
    val = sum(sub(valid) .* w(valid)) / sum(w(valid));
end

function [slope, pval] = lin_trend(x, y)
% OLS 线性趋势及 t 检验 p 值
    x = double(x(:)); y = double(y(:));
    v = ~isnan(y); x = x(v); y = y(v);
    n = numel(y);
    if n < 4; slope = NaN; pval = NaN; return; end
    pp    = polyfit(x, y, 1);
    resid = y - polyval(pp, x);
    se    = sqrt(sum(resid.^2) / (n - 2) / sum((x - mean(x)).^2));
    slope = pp(1);
    pval  = 2 * (1 - tcdf(abs(slope / se), n - 2));
end

function name = pick_name(vnames, candidates)
    name = '';
    for k = 1:numel(candidates)
        if any(strcmp(vnames, candidates{k}))
            name = candidates{k}; return;
        end
    end
    error('找不到变量，候选名: %s', strjoin(candidates, ', '));
end

function dates = parse_nc_time(time_raw, units_str)
    tokens = regexp(units_str, ...
        '(\w+)\s+since\s+(\d{4})-(\d{1,2})-(\d{1,2})', 'tokens', 'once');
    if isempty(tokens)
        error('Cannot parse time units: %s', units_str);
    end
    unit_str = lower(tokens{1});
    ref = datetime(str2double(tokens{2}), str2double(tokens{3}), str2double(tokens{4}));
    switch unit_str
        case {'day','days'};       dates = ref + days(time_raw);
        case {'hour','hours'};     dates = ref + hours(time_raw);
        case {'minute','minutes'}; dates = ref + minutes(time_raw);
        case {'second','seconds'}; dates = ref + seconds(time_raw);
        case {'month','months'};   dates = ref + calmonths(round(time_raw));
        otherwise; error('Unknown time unit: %s', unit_str);
    end
    dates = dates(:);
end

function data_2d = extract_season_mean(field, yrs, mos, ...
                                       yr, season_months, cross_year, lev_idx)
    [nlat, nlon, ~, ~] = size(field);
    n_mo = numel(season_months);
    season_3d = NaN(nlat, nlon, n_mo);
    for k = 1:n_mo
        mo = season_months(k);
        if cross_year && mo == 12
            target_yr = yr;
        elseif cross_year && mo < 6
            target_yr = yr + 1;
        else
            target_yr = yr;
        end
        idx = find((yrs == target_yr) & (mos == mo));
        if isempty(idx); continue; end
        season_3d(:,:,k) = mean(field(:,:,lev_idx,idx(1)), 3, 'omitnan');
    end
    data_2d = mean(season_3d, 3, 'omitnan');
end

function [trend_map, pval_map] = calc_trend_map(samples, all_years)
    [nlat, nlon, ~] = size(samples);
    trend_map = NaN(nlat, nlon);
    pval_map  = NaN(nlat, nlon);
    x = double(all_years);
    for ilat = 1:nlat
        for ilon = 1:nlon
            y = squeeze(samples(ilat, ilon, :));
            valid = ~isnan(y);
            if sum(valid) < 4; continue; end
            x_v = x(valid);
            y_v = y(valid);
            n   = sum(valid);
            pp    = polyfit(x_v, y_v, 1);
            resid = y_v - polyval(pp, x_v);
            se    = sqrt(sum(resid.^2) / (n - 2) / sum((x_v - mean(x_v)).^2));
            t_val = pp(1) / se;
            trend_map(ilat, ilon) = pp(1);
            pval_map(ilat, ilon)  = 2 * (1 - tcdf(abs(t_val), n - 2));
        end
    end
end
