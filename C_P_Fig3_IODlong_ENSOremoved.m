%% ============================================================
%  综合图：3行1列
% ============================================================
%% ===== 基础设置 =====
lon_range = [30 120];
lat_range = [-20 40];
yr_start  = 2003;
yr_end    = 2020;
reg_years = (yr_start:yr_end)';
n_reg     = numel(reg_years);
s_months = [9 10 11];
s_cross  = false;
%% ===== 1. 读取 HadISST =====
fprintf('读取 HadISST...\n');
nc_sst       = 'C:\Users\yingy\Downloads\Rice_mac\HadISST_sst.nc';
sst_lat      = ncread(nc_sst, 'latitude');
sst_lon      = ncread(nc_sst, 'longitude');
sst_time_raw = ncread(nc_sst, 'time');
sst_raw      = ncread(nc_sst, 'sst');
sst_dates  = datetime(1870,1,1) + days(double(sst_time_raw));
sst_years  = year(sst_dates);
sst_months = month(sst_dates);
sst_raw = permute(sst_raw, [2 1 3]);
if sst_lat(1) > sst_lat(end)
    sst_lat = flipud(sst_lat);
    sst_raw = flipud(sst_raw);
end
sst_lon(sst_lon < 0) = sst_lon(sst_lon < 0) + 360;
[sst_lon, sidx] = sort(sst_lon);
sst_raw = sst_raw(:, sidx, :);
sst_raw(sst_raw < -2) = NaN;
lat_idx_sst = sst_lat >= lat_range(1) & sst_lat <= lat_range(2);
lon_idx_sst = sst_lon >= lon_range(1) & sst_lon <= lon_range(2);
sst_lat_sub = sst_lat(lat_idx_sst);
sst_lon_sub = sst_lon(lon_idx_sst);
sst_sub     = sst_raw(lat_idx_sst, lon_idx_sst, :);
[nlat_sst, nlon_sst, ~] = size(sst_sub);
sst_son = NaN(nlat_sst, nlon_sst, n_reg);
for iy = 1:n_reg
    yr  = reg_years(iy);
    tmp = NaN(nlat_sst, nlon_sst, 3);
    for im = 1:3
        mo  = s_months(im);
        idx = find(sst_years == yr & sst_months == mo, 1);
        if ~isempty(idx); tmp(:,:,im) = sst_sub(:,:,idx); end
    end
    sst_son(:,:,iy) = mean(tmp, 3, 'omitnan');
end
t_vec     = reg_years - mean(reg_years);
X_reg     = [ones(n_reg,1), t_vec];
sst_trend = NaN(nlat_sst, nlon_sst);
sst_pval  = NaN(nlat_sst, nlon_sst);
for ilat = 1:nlat_sst
    for ilon = 1:nlon_sst
        y = squeeze(sst_son(ilat,ilon,:));
        valid = ~isnan(y);
        if sum(valid) < 6; continue; end
        y_v = y(valid); X_v = X_reg(valid,:);
        n_v = sum(valid); df = n_v-2;
        b   = X_v \ y_v;
        resid = y_v - X_v*b;
        s2    = sum(resid.^2)/df;
        se_b  = sqrt(s2 * inv(X_v'*X_v));
        t_val = b(2)/se_b(2,2);
        sst_trend(ilat,ilon) = b(2)*10;
        sst_pval(ilat,ilon)  = 2*(1-tcdf(abs(t_val),df));
    end
end
fprintf('SST trend 完成\n');
lon_5_sst = (32.5:5:117.5)';
lat_5_sst = (-17.5:5:37.5)';
nlon5_sst = numel(lon_5_sst);
nlat5_sst = numel(lat_5_sst);
sst_trend_5 = NaN(nlat5_sst, nlon5_sst);
sst_pval_5  = NaN(nlat5_sst, nlon5_sst);
for ilat = 1:nlat5_sst
    lat_mask = sst_lat_sub >= lat_5_sst(ilat)-2.5 & sst_lat_sub < lat_5_sst(ilat)+2.5;
    for ilon = 1:nlon5_sst
        lon_mask = sst_lon_sub >= lon_5_sst(ilon)-2.5 & sst_lon_sub < lon_5_sst(ilon)+2.5;
        blk_t = sst_trend(lat_mask, lon_mask);
        blk_p = sst_pval(lat_mask, lon_mask);
        sst_trend_5(ilat,ilon) = mean(blk_t(:),'omitnan');
        sst_pval_5(ilat,ilon)  = mean(blk_p(:),'omitnan');
    end
end
fprintf('SST 5° 插值完成\n');
%% ===== 2. 读取 AIRS dD =====
fprintf('读取 AIRS dD...\n');
nc_file = 'C:\Users\yingy\Downloads\Rice_mac\Isotope-IOD\AIRS_RS_monthly_5deg.nc';
lat     = ncread(nc_file, 'lat');
lon     = ncread(nc_file, 'lon');
dD      = ncread(nc_file, 'dD');
dD      = permute(dD, [2 1 3 4]);
dD      = flipud(dD);
lat     = flipud(lat);
lon_360 = lon;
lon_360(lon_360 < 0) = lon_360(lon_360 < 0) + 360;
[lon_sorted, sort_idx] = sort(lon_360);
dD  = dD(:, sort_idx, :, :);
lon = lon_sorted;
[~, ~, ~, ntime] = size(dD);
pressure_levels = [1017.2227, 1000, 908.51398, 825.40198, 749.89301, ...
                   681.29102, 618.96600, 510.89801, 421.69800, 348.06900, ...
                   287.29800, 237.13699, 177.82899, 133.35201, 74.989601, ...
                   28.729900, 0.10000];
time_raw  = ncread(nc_file, 'time');
dates     = datetime(2000,1,1) + days(double(time_raw));
years_dD  = year(dates);
months_dD = month(dates);
lon_idx_dD = lon >= lon_range(1) & lon <= lon_range(2);
lat_idx_dD = lat >= lat_range(1) & lat <= lat_range(2);
lon_sub    = lon(lon_idx_dD);
lat_sub    = lat(lat_idx_dD);
nlat       = sum(lat_idx_dD);
nlon       = sum(lon_idx_dD);
dD_sub     = dD(lat_idx_dD, lon_idx_dD, :, :);
%% ===== 3. 原始 AIRS δD 距平（不做 Rayleigh 分离，不需要 ERA5）=====
fprintf('计算原始 δD 月距平...\n');
nlev_airs   = numel(pressure_levels);
climatology = NaN(nlat, nlon, nlev_airs, 12);
for mo = 1:12
    idx = (months_dD == mo);
    climatology(:,:,:,mo) = mean(dD_sub(:,:,:,idx), 4, 'omitnan');
end
dD_anom = NaN(size(dD_sub));
for t = 1:ntime
    dD_anom(:,:,:,t) = dD_sub(:,:,:,t) - climatology(:,:,:,months_dD(t));
end
fprintf('δD 距平完成\n');
lev_mid = find(pressure_levels >= 350 & pressure_levels <= 650);
%% ===== 4. topo + land mask =====
load topo;
topo_lon_v   = 0.5:359.5;
topo_lat_v   = -89.5:89.5;
topo_ext     = [topo, topo(:,1)];
topo_lon_ext = [topo_lon_v, 360.5];
[TLON, TLAT] = meshgrid(topo_lon_v, topo_lat_v);
[OLON, OLAT] = meshgrid(lon_sub, lat_sub);
land_mask = interp2(TLON, TLAT, topo, OLON, OLAT, 'nearest', 0) > 0;
[OLON_5, OLAT_5] = meshgrid(lon_5_sst, lat_5_sst);
land_mask_sst5 = interp2(TLON, TLAT, topo, OLON_5, OLAT_5, 'nearest', 0) > 0;
topo_pack.lon  = topo_lon_ext;
topo_pack.lat  = topo_lat_v;
topo_pack.elev = topo_ext;
%% ===== 西印度洋区域 SST trend 统计（45–80°E, 20°S–20°N）=====
wio_lon_range = [45 80];
wio_lat_range = [-20 20];
wio_lon_idx5 = lon_5_sst >= wio_lon_range(1) & lon_5_sst <= wio_lon_range(2);
wio_lat_idx5 = lat_5_sst >= wio_lat_range(1) & lat_5_sst <= wio_lat_range(2);
wio_trend5 = sst_trend_5(wio_lat_idx5, wio_lon_idx5);
wio_land5  = land_mask_sst5(wio_lat_idx5, wio_lon_idx5);
wio_lat5_v = lat_5_sst(wio_lat_idx5);
wio_trend5(wio_land5) = NaN;
w  = cosd(wio_lat5_v);
Wm = repmat(w(:), 1, size(wio_trend5,2));
Wm(isnan(wio_trend5)) = NaN;
wio_mean5 = sum(wio_trend5 .* Wm, 'all', 'omitnan') / sum(Wm, 'all', 'omitnan');
wio_min5  = min(wio_trend5, [], 'all', 'omitnan');
wio_max5  = max(wio_trend5, [], 'all', 'omitnan');
fprintf('\n===== 西印度洋 (45–80°E, 20°S–20°N) SST trend（5°格点，面积加权）=====\n');
fprintf('  区域平均: %.3f °C/10yr\n', wio_mean5);
fprintf('  格点最小值: %.3f °C/10yr\n', wio_min5);
fprintf('  格点最大值: %.3f °C/10yr\n', wio_max5);
wio_lon_idx1 = sst_lon_sub >= wio_lon_range(1) & sst_lon_sub <= wio_lon_range(2);
wio_lat_idx1 = sst_lat_sub >= wio_lat_range(1) & sst_lat_sub <= wio_lat_range(2);
wio_trend1 = sst_trend(wio_lat_idx1, wio_lon_idx1);
wio_lat1_v = sst_lat_sub(wio_lat_idx1);
w1  = cosd(wio_lat1_v);
Wm1 = repmat(w1(:), 1, size(wio_trend1,2));
Wm1(isnan(wio_trend1)) = NaN;
wio_mean1 = sum(wio_trend1 .* Wm1, 'all', 'omitnan') / sum(Wm1, 'all', 'omitnan');
wio_min1  = min(wio_trend1, [], 'all', 'omitnan');
wio_max1  = max(wio_trend1, [], 'all', 'omitnan');
fprintf('\n===== 西印度洋 (45–80°E, 20°S–20°N) SST trend（原始1°格点，面积加权）=====\n');
fprintf('  区域平均: %.3f °C/10yr\n', wio_mean1);
fprintf('  格点最小值: %.3f °C/10yr\n', wio_min1);
fprintf('  格点最大值: %.3f °C/10yr\n', wio_max1);
%% ===== 5. 读取 IOD index =====
fprintf('读取 IOD index...\n');
iod_data    = load('C:\Users\yingy\Downloads\Rice_mac\Isotope-IOD\IOD_2000_2020CE_Hadisst.mat');
iod_vars    = fieldnames(iod_data);
iod_annual  = double(iod_data.(iod_vars{1})); iod_annual = iod_annual(:);
iod_yr_axis = (2000:2020)';
IOD_reg = NaN(n_reg,1);
for iy = 1:n_reg
    idx = find(iod_yr_axis == reg_years(iy));
    if ~isempty(idx); IOD_reg(iy) = iod_annual(idx); end
end
IOD_norm = (IOD_reg - mean(IOD_reg,'omitnan')) / std(IOD_reg,'omitnan');

%% ===== 5b. 图c 用：由 HadISST 直接计算的长序列 SON DMI（不去除 ENSO）=====
%   DMI = WTIO (50–70°E, 10°S–10°N) − SETIO (90–110°E, 10°S–0°)，SON 平均
%   时段可在下面修改；结束年份自动截到 HadISST 最后一个完整 SON
iod_long_start = 1980;
iod_long_end   = 2024;
fprintf('由 HadISST 计算长序列 SON DMI（不去除 ENSO）...\n');
wtio_lat = sst_lat >= -10 & sst_lat <= 10;  wtio_lon = sst_lon >= 50 & sst_lon <= 70;
setio_lat = sst_lat >= -10 & sst_lat <= 0;  setio_lon = sst_lon >= 90 & sst_lon <= 110;
wtio_ts  = area_mean_ts(sst_raw(wtio_lat,  wtio_lon,  :), sst_lat(wtio_lat));
setio_ts = area_mean_ts(sst_raw(setio_lat, setio_lon, :), sst_lat(setio_lat));
dmi_mon  = wtio_ts - setio_ts;
last_son_yr = max(sst_years(sst_months == 11));
iod_long_end = min(iod_long_end, last_son_yr);
iod_years = (iod_long_start:iod_long_end)';
IOD_long  = NaN(numel(iod_years),1);
for iy = 1:numel(iod_years)
    sel = sst_years == iod_years(iy) & ismember(sst_months, s_months);
    if sum(sel) == 3; IOD_long(iy) = mean(dmi_mon(sel), 'omitnan'); end
end
% 与 dD 重叠年份（2003–2020）上的原始 DMI，用于图c的相关系数
IOD_raw_reg = NaN(n_reg,1);
[tf, loc] = ismember(reg_years, iod_years);
IOD_raw_reg(tf) = IOD_long(loc(tf));
fprintf('长序列 DMI：%d–%d（%d 年）\n', iod_years(1), iod_years(end), numel(iod_years));
%% ===== 6. 计算 Nino3.4 index =====
fprintf('计算 Nino3.4 index...\n');
nino_lat = [-5 5]; nino_lon = [190 240];
nino_lat_idx = sst_lat >= nino_lat(1) & sst_lat <= nino_lat(2);
nino_lon_idx = sst_lon >= nino_lon(1) & sst_lon <= nino_lon(2);
sst_nino     = sst_raw(nino_lat_idx, nino_lon_idx, :);
nino_ts      = squeeze(mean(mean(sst_nino,1,'omitnan'),2,'omitnan'));
NINO_reg = NaN(n_reg,1);
for iy = 1:n_reg
    idx_sel = (sst_years == reg_years(iy)) & ismember(sst_months, s_months);
    if any(idx_sel)
        NINO_reg(iy) = mean(nino_ts(idx_sel),'omitnan');
    end
end
NINO_norm = (NINO_reg - mean(NINO_reg,'omitnan')) / std(NINO_reg,'omitnan');
%% ===== 6b. 去除 ENSO 信号（线性回归去除 Niño3.4）=====
%   enso_season = 'DJF'：用 D(0)JF(1) Niño3.4（ENSO 成熟期，与前文偏回归一致）
%   enso_season = 'SON'：用同期 SON Niño3.4
enso_season = 'DJF';
if strcmp(enso_season, 'DJF')
    % DJF 需要下一年的 1–2 月：长序列结束年份截到有完整 D(0)JF(1) 的年份
    last_djf_yr = max(sst_years(sst_months == 2)) - 1;
    keep_long   = iod_years <= last_djf_yr;
    iod_years   = iod_years(keep_long);
    IOD_long    = IOD_long(keep_long);
end
nino_long = season_index(nino_ts, sst_years, sst_months, iod_years, enso_season);
nino_reg  = season_index(nino_ts, sst_years, sst_months, reg_years, enso_season);
% 长序列 IOD：在其全时段上去除 ENSO（用于图b）
IOD_long_res = regress_out(IOD_long, nino_long);
% 2003–2020 IOD：在同一时段上去除 ENSO（用于图a回归和图b相关）
IOD_reg_res  = regress_out(IOD_raw_reg, nino_reg);
IOD_norm     = (IOD_reg_res - mean(IOD_reg_res,'omitnan')) / std(IOD_reg_res,'omitnan');
fprintf('ENSO 去除：%s Niño3.4；长序列 IOD %d–%d\n', enso_season, iod_years(1), iod_years(end));
%% ===== 7. SON 原始 δD 距平年均 =====
fprintf('计算 SON δD...\n');
dD_son = NaN(nlat, nlon, n_reg);
for iy = 1:n_reg
    yr = reg_years(iy);
    dD_son(:,:,iy) = extract_season_mean_fn(dD_anom, years_dD, months_dD, ...
        yr, s_months, s_cross, lev_mid);
end
% 逐格点去除 ENSO（2003–2020 时段）
for ilat = 1:nlat
    for ilon = 1:nlon
        dD_son(ilat,ilon,:) = reshape(regress_out(squeeze(dD_son(ilat,ilon,:)), nino_reg), 1, 1, []);
    end
end
fprintf('δD 已逐格点去除 ENSO\n');
%% ===== 8. 逐格点回归 =====
fprintf('计算回归图...\n');
beta_iod  = NaN(nlat, nlon);
pval_iod  = NaN(nlat, nlon);
beta_nino = NaN(nlat, nlon);
pval_nino = NaN(nlat, nlon);
for ilat = 1:nlat
    for ilon = 1:nlon
        if land_mask(ilat,ilon); continue; end
        y = squeeze(dD_son(ilat,ilon,:));
        valid = ~isnan(y) & ~isnan(IOD_norm) & ~isnan(NINO_norm);
        if sum(valid) < 5; continue; end
        y_v = y(valid); n_v = sum(valid);
        x_v = IOD_norm(valid);
        % 两者都已去除 ENSO → 等价于偏回归，自由度 n−3
        pp    = polyfit(x_v, y_v, 1);
        yhat  = polyval(pp, x_v);
        resid = y_v - yhat;
        s2    = sum(resid.^2)/(n_v-3);
        Sxx   = sum((x_v-mean(x_v)).^2);
        se    = sqrt(s2/Sxx);
        t_val = pp(1)/se;
        beta_iod(ilat,ilon) = pp(1);
        pval_iod(ilat,ilon) = 2*(1-tcdf(abs(t_val),n_v-3));
        x_v   = NINO_norm(valid);
        pp    = polyfit(x_v, y_v, 1);
        yhat  = polyval(pp, x_v);
        resid = y_v - yhat;
        s2    = sum(resid.^2)/(n_v-2);
        Sxx   = sum((x_v-mean(x_v)).^2);
        se    = sqrt(s2/Sxx);
        t_val = pp(1)/se;
        beta_nino(ilat,ilon) = pp(1);
        pval_nino(ilat,ilon) = 2*(1-tcdf(abs(t_val),n_v-2));
    end
end
fprintf('回归计算完成\n');
%% ===== 9. 时间序列 =====
fprintf('计算时间序列...\n');
r1_lat_idx = lat_sub >= -15 & lat_sub <= 5;
r1_lon_idx = lon_sub >=  50 & lon_sub <= 70;
r2_lat_idx = lat_sub >= -15 & lat_sub <= 5;
r2_lon_idx = lon_sub >=  70 & lon_sub <= 100;
dD_diff = NaN(n_reg,1);
for iy = 1:n_reg
    tmp = dD_son(:,:,iy);
    r1  = mean(tmp(r1_lat_idx, r1_lon_idx),'all','omitnan');
    r2  = mean(tmp(r2_lat_idx, r2_lon_idx),'all','omitnan');
    dD_diff(iy) = r1 - r2;
end

% 保存图b曲线数据（原始 δD，已去 ENSO）
save('Fig3_curve_data_rawdD_ENSOremoved.mat', 'reg_years', 'IOD_reg_res', 'nino_reg', 'dD_diff', 'iod_years', 'IOD_long_res', 'enso_season');

% 图b：去 ENSO 后的 DMI 与 δDDI 在 2003–2020 的相关（偏相关，自由度 n−3）
valid_ts = ~isnan(dD_diff) & ~isnan(IOD_reg_res);
n_ts     = sum(valid_ts);
r_val    = corr(IOD_reg_res(valid_ts), dD_diff(valid_ts), 'Type','Pearson');
t_ts     = r_val * sqrt((n_ts-3) / (1-r_val^2));
p_val_ts = 2*(1-tcdf(abs(t_ts), n_ts-3));
fprintf('时间序列偏相关（去 ENSO 的 DMI vs δDDI, %d–%d, n=%d）: r=%.3f, p=%.4f\n', ...
    reg_years(1), reg_years(end), n_ts, r_val, p_val_ts);
%% ===== 图c trend 输出 =====
x_fit_iod = double(iod_years);
x_fit     = double(reg_years);
IOD_norm_plot   = (IOD_long_res - mean(IOD_long_res,'omitnan')) / std(IOD_long_res,'omitnan');
dD_norm_plot = (dD_diff - mean(dD_diff,'omitnan')) / std(dD_diff,'omitnan');
valid_iod = ~isnan(IOD_norm_plot);
p_iod  = polyfit(x_fit_iod(valid_iod), IOD_norm_plot(valid_iod), 1);
n_iod  = sum(valid_iod);
yhat   = polyval(p_iod, x_fit_iod(valid_iod));
resid  = IOD_norm_plot(valid_iod) - yhat;
s2     = sum(resid.^2)/(n_iod-2);
Sxx    = sum((x_fit_iod(valid_iod)-mean(x_fit_iod(valid_iod))).^2);
se     = sqrt(s2/Sxx);
t_iod  = p_iod(1)/se;
pv_iod = 2*(1-tcdf(abs(t_iod), n_iod-2));
valid_dd = ~isnan(dD_norm_plot);
p_dd   = polyfit(x_fit(valid_dd), dD_norm_plot(valid_dd), 1);
n_dd   = sum(valid_dd);
yhat   = polyval(p_dd, x_fit(valid_dd));
resid  = dD_norm_plot(valid_dd) - yhat;
s2     = sum(resid.^2)/(n_dd-2);
Sxx    = sum((x_fit(valid_dd)-mean(x_fit(valid_dd))).^2);
se     = sqrt(s2/Sxx);
t_dd   = p_dd(1)/se;
pv_dd  = 2*(1-tcdf(abs(t_dd), n_dd-2));
fprintf('\n===== 图c 趋势分析 =====\n');
fprintf('IOD index (ENSO removed, normalized, %d–%d): trend=%.4f/yr, p=%.4f\n', iod_years(1), iod_years(end), p_iod(1), pv_iod);
fprintf('δDDI (ENSO removed, normalized, %d–%d): trend=%.4f/yr, p=%.4f\n', reg_years(1), reg_years(end), p_dd(1), pv_dd);
%% ===== 10. colormap =====
load BlueWhiteOrangeRed.mat;
%% ===== 11. 画图 =====
fig = make_2x2_figure(sst_trend_5, sst_pval_5, lon_5_sst, lat_5_sst, land_mask_sst5, ...
    beta_iod, pval_iod, beta_nino, pval_nino, ...
    lon_sub, lat_sub, land_mask, ...
    IOD_long_res, iod_years, dD_diff, reg_years, r_val, p_val_ts, ...
    lon_range, lat_range, BlueWhiteOrangeRed, topo_pack);
print(fig, 'Fig.3_new_IODlong_rawdD_ENSOremoved.tif', '-dpng', '-r600', '-painters');
fig.Visible = 'on';
disp('=== 图保存完成 ===');
%% ============================================================
function fig = make_2x2_figure(sst_trend, sst_pval, sst_lon, sst_lat, land_mask_sst, ...
    beta_iod, pval_iod, beta_nino, pval_nino, ...
    lon_sub, lat_sub, land_mask, ...
    IOD_long, iod_years, dD_diff, reg_years, r_val, p_val_ts, ...
    lon_range, lat_range, cmap, topo_pack)
    fs = 9;
    cm2px    = 96/2.54;
    fig_w_cm = 18;
    fig_w    = round(fig_w_cm * cm2px);
    lon_span = lon_range(2) - lon_range(1);
    lat_span = lat_range(2) - lat_range(1);
    aspect   = lon_span / lat_span;
    pw_norm    = 0.42;
    panel_w_px = round(pw_norm * fig_w);
    panel_h_px = round(panel_w_px / aspect);
    ts_h_px = 180;
    mt_cm     = 1.2;
    mb_cm     = 1.8;
    vgap12_cm = 1.5;
    vgap23_cm = 1.5;
    mt_px     = round(mt_cm     * cm2px);
    mb_px     = round(mb_cm     * cm2px);
    vgap12_px = round(vgap12_cm * cm2px);
    vgap23_px = round(vgap23_cm * cm2px);
    % 三幅：(a) SST 趋势，(b) δD 对 IOD 回归，(c) 时间序列
    fig_h = mt_px + 2*panel_h_px + vgap12_px + vgap23_px + ts_h_px + mb_px;
    ph_norm   = panel_h_px / fig_h;
    ts_norm   = ts_h_px    / fig_h;
    mt_norm   = mt_px      / fig_h;
    vg12_norm = vgap12_px  / fig_h;
    vg23_norm = vgap23_px  / fig_h;
    ml_norm = (1 - pw_norm) / 2;
    pos_y_R1 = 1 - mt_norm - ph_norm;
    pos_y_R2 = pos_y_R1 - vg12_norm - ph_norm;
    pos_y_R3 = pos_y_R2 - vg23_norm - ts_norm;
    map_pos = @(pos_y) [ml_norm, pos_y, pw_norm, ph_norm];
    fig = figure('Color','w','Position',[50 50 fig_w fig_h], ...
        'Renderer','painters', ...
        'PaperUnits','centimeters', ...
        'PaperSize',[fig_w_cm, fig_h*2.54/cm2px]);
    buf = 5;
    tlon_idx    = topo_pack.lon >= lon_range(1)-buf & topo_pack.lon <= lon_range(2)+buf;
    tlat_idx    = topo_pack.lat >= lat_range(1)-buf & topo_pack.lat <= lat_range(2)+buf;
    topo_lon_v  = topo_pack.lon(tlon_idx);
    topo_lat_v  = topo_pack.lat(tlat_idx);
    topo_elev_v = topo_pack.elev(tlat_idx, tlon_idx);
    clim_sst = 0.5;
    clim_reg = 10;
    %% (a) SST Trend
    ax1 = axes('Position', map_pos(pos_y_R1));
    hold on;
    data1 = sst_trend; data1(land_mask_sst) = NaN;
    imagesc(ax1, sst_lon, sst_lat, data1, 'AlphaData', ~isnan(data1));
    set(ax1, 'YDir','normal', 'Color',[0.85 0.85 0.85]);
    colormap(ax1, cmap);
    clim(ax1, [-clim_sst clim_sst]);
    contour(ax1, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth',0.5);
    set(ax1, 'Clipping','on');
    for r = 1:numel(sst_lat)
        for c = 1:numel(sst_lon)
            if land_mask_sst(r,c); continue; end
            if isnan(sst_pval(r,c)); continue; end
            if sst_pval(r,c) < 0.05
                plot(ax1, sst_lon(c), sst_lat(r), 'o', 'MarkerSize',3, ...
                    'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none', ...
                    'HandleVisibility','off');
            end
        end
    end
    xlim(ax1, lon_range); ylim(ax1, lat_range);
    set(ax1, 'XTick',30:30:120, 'XTickLabel',{}, 'FontSize',fs);
    set(ax1, 'YTick',-20:20:40, ...
        'YTickLabel',{'20°S','0°','20°N','40°N'}, ...
        'FontSize',fs, 'TickDir','out', 'Box','off', 'Layer','top');
    cb1 = colorbar(ax1, 'Location','eastoutside');
    cb1.Position      = [ml_norm+pw_norm+0.025, pos_y_R1, 0.02, ph_norm];
    cb1.Label.String  = '°C / 10yr';
    cb1.Ticks         = -0.5:0.25:0.5;
    cb1.FontSize      = fs;
    cb1.Label.FontSize= fs;
    cb1.TickDirection = 'out';
    ax1.Position = map_pos(pos_y_R1);
    xl=ax1.XLim; yl=ax1.YLim;
    plot(ax1,[xl(2) xl(2)],yl,'k-','LineWidth',0.5,'HandleVisibility','off');
    plot(ax1,xl,[yl(2) yl(2)],'k-','LineWidth',0.5,'HandleVisibility','off');
    title(ax1,'(a) SST Trend (°C/10yr, 2003–2020)', ...
        'FontSize',fs,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized','Position',[0,1.02,0]);
    %% (b) δD regression on IOD
    ax2 = axes('Position', map_pos(pos_y_R2));
    hold on;
    data2 = beta_iod; data2(land_mask) = NaN;
    imagesc(ax2, lon_sub, lat_sub, data2, 'AlphaData', ~isnan(data2));
    set(ax2, 'YDir','normal', 'Color',[0.85 0.85 0.85]);
    colormap(ax2, cmap);
    clim(ax2, [-clim_reg clim_reg]);
    contour(ax2, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth',0.5);
    set(ax2, 'Clipping','on');
    for r = 1:numel(lat_sub)
        for c = 1:numel(lon_sub)
            if land_mask(r,c); continue; end
            if isnan(pval_iod(r,c)); continue; end
            if pval_iod(r,c) < 0.05
                plot(ax2, lon_sub(c), lat_sub(r), 'o', 'MarkerSize',2.5, ...
                    'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none', ...
                    'HandleVisibility','off');
            end
        end
    end
    xlim(ax2, lon_range); ylim(ax2, lat_range);
    set(ax2, 'XTick',30:30:120, ...
        'XTickLabel',{'30°E','60°E','90°E','120°E'}, 'FontSize',fs);
    set(ax2, 'YTick',-20:20:40, ...
        'YTickLabel',{'20°S','0°','20°N','40°N'}, ...
        'FontSize',fs, 'TickDir','out', 'Box','off', 'Layer','top');
    cb2 = colorbar(ax2, 'Location','eastoutside');
    cb2.Position      = [ml_norm+pw_norm+0.025, pos_y_R2, 0.02, ph_norm];
    cb2.Label.String  = '‰ / \sigma_{IOD}';
    cb2.FontSize      = fs;
    cb2.Label.FontSize= fs;
    cb2.TickDirection = 'out';
    ax2.Position = map_pos(pos_y_R2);
    xl=ax2.XLim; yl=ax2.YLim;
    plot(ax2,[xl(2) xl(2)],yl,'k-','LineWidth',0.5,'HandleVisibility','off');
    plot(ax2,xl,[yl(2) yl(2)],'k-','LineWidth',0.5,'HandleVisibility','off');
    title(ax2,'(b) \deltaD regression on IOD (400–600 hPa, ENSO removed)', ...
        'FontSize',fs,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized','Position',[0,1.02,0], ...
        'Interpreter','tex');
    %% (c) 时间序列
    ax_ts = axes('Position', [ml_norm, pos_y_R3, pw_norm, ts_norm]);
    hold on;
    % IOD：长序列（已去 ENSO），在其全时段上标准化；δDDI：2003–2020（已去 ENSO），在其时段上标准化
    IOD_norm_plot   = (IOD_long   - mean(IOD_long,'omitnan'))   / std(IOD_long,'omitnan');
    dD_norm_plot = (dD_diff - mean(dD_diff,'omitnan')) / std(dD_diff,'omitnan');
    plot(ax_ts, iod_years, IOD_norm_plot, 'b-', ...
        'LineWidth',1.2, 'DisplayName','IOD index');
    plot(ax_ts, reg_years, dD_norm_plot, 'r-', ...
        'LineWidth',1.2, 'DisplayName','\deltaDDI');
    x_iod     = double(iod_years);
    valid_iod = ~isnan(IOD_norm_plot);
    p_iod     = polyfit(x_iod(valid_iod), IOD_norm_plot(valid_iod), 1);
    plot(ax_ts, iod_years, polyval(p_iod, x_iod), 'b--', 'LineWidth',1.5, 'HandleVisibility','off');
    x_fit     = double(reg_years);
    valid_dd  = ~isnan(dD_norm_plot);
    p_dd      = polyfit(x_fit(valid_dd), dD_norm_plot(valid_dd), 1);
    plot(ax_ts, reg_years, polyval(p_dd, x_fit), 'r--', 'LineWidth',1.5, 'HandleVisibility','off');
    yline(ax_ts, 0, 'k--', 'LineWidth',0.5, 'HandleVisibility','off');
    xlim(ax_ts, [iod_years(1)-0.5, iod_years(end)+0.5]);
    set(ax_ts, 'XTick', ceil(iod_years(1)/5)*5 : 5 : iod_years(end), ...
        'FontSize',fs, 'TickDir','out', 'Box','off');
    xlabel(ax_ts, 'Year', 'FontSize',fs);
    ylabel(ax_ts, 'Normalized anomaly', 'FontSize',fs);
    lg = legend(ax_ts, 'FontSize',fs, 'Box','off', 'Interpreter','tex');
    lg.Position(1) = ml_norm + pw_norm*0.05;
    title(ax_ts, '(c) IOD vs \deltaDDI (ENSO removed)', ...
        'FontSize',fs,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized','Position',[0,1.02,0], ...
        'Interpreter','tex');
    set(fig, 'PaperPositionMode','auto');
    text(ax_ts, 0.02, 0.08, sprintf('r = %.2f, p = %.3f (%d–%d)', r_val, p_val_ts, reg_years(1), reg_years(end)), ...
        'Units','normalized', 'FontSize',fs, 'FontWeight','bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');
end
%% ============================================================
function fig = make_nino_figure(beta_nino, pval_nino, lon_sub, lat_sub, land_mask, ...
    lon_range, lat_range, cmap, topo_pack)
    fs = 9;
    cm2px    = 96/2.54;
    fig_w_cm = 10;
    fig_w    = round(fig_w_cm * cm2px);
    lon_span = lon_range(2) - lon_range(1);
    lat_span = lat_range(2) - lat_range(1);
    aspect   = lon_span / lat_span;
    pw_norm    = 0.60;
    panel_w_px = round(pw_norm * fig_w);
    panel_h_px = round(panel_w_px / aspect);
    mt_cm = 0.9; mb_cm = 1.2;
    mt_px = round(mt_cm * cm2px);
    mb_px = round(mb_cm * cm2px);
    fig_h = mt_px + panel_h_px + mb_px;
    ph_norm = panel_h_px / fig_h;
    mt_norm = mt_px      / fig_h;
    ml_norm = (1 - pw_norm) / 2;
    pos_y   = 1 - mt_norm - ph_norm;
    fig = figure('Color','w','Position',[50 50 fig_w fig_h], ...
        'Renderer','painters', ...
        'PaperUnits','centimeters', ...
        'PaperSize',[fig_w_cm, fig_h*2.54/cm2px]);
    buf = 5;
    tlon_idx    = topo_pack.lon >= lon_range(1)-buf & topo_pack.lon <= lon_range(2)+buf;
    tlat_idx    = topo_pack.lat >= lat_range(1)-buf & topo_pack.lat <= lat_range(2)+buf;
    topo_lon_v  = topo_pack.lon(tlon_idx);
    topo_lat_v  = topo_pack.lat(tlat_idx);
    topo_elev_v = topo_pack.elev(tlat_idx, tlon_idx);
    clim_reg = 10;
    ax = axes('Position', [ml_norm, pos_y, pw_norm, ph_norm]);
    hold on;
    data = beta_nino; data(land_mask) = NaN;
    imagesc(ax, lon_sub, lat_sub, data, 'AlphaData', ~isnan(data));
    set(ax, 'YDir','normal', 'Color',[0.85 0.85 0.85]);
    colormap(ax, cmap);
    clim(ax, [-clim_reg clim_reg]);
    contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth',0.5);
    set(ax, 'Clipping','on');
    for r = 1:numel(lat_sub)
        for c = 1:numel(lon_sub)
            if land_mask(r,c); continue; end
            if isnan(pval_nino(r,c)); continue; end
            if pval_nino(r,c) < 0.05
                plot(ax, lon_sub(c), lat_sub(r), 'o', 'MarkerSize',2.5, ...
                    'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none', ...
                    'HandleVisibility','off');
            end
        end
    end
    xlim(ax, lon_range); ylim(ax, lat_range);
    set(ax, 'XTick',30:30:120, ...
        'XTickLabel',{'30°E','60°E','90°E','120°E'}, 'FontSize',fs);
    set(ax, 'YTick',-20:20:40, ...
        'YTickLabel',{'20°S','0°','20°N','40°N'}, ...
        'FontSize',fs, 'TickDir','out', 'Box','off', 'Layer','top');
    cb = colorbar(ax,'Location','eastoutside');
    cb.Label.String  = '‰ / \sigma_{Nino3.4}';
    cb.FontSize      = fs;
    cb.Label.FontSize= fs;
    cb.TickDirection = 'out';
    ax.Position = [ml_norm, pos_y, pw_norm, ph_norm];
    xl=ax.XLim; yl=ax.YLim;
    plot(ax,[xl(2) xl(2)],yl,'k-','LineWidth',0.5,'HandleVisibility','off');
    plot(ax,xl,[yl(2) yl(2)],'k-','LineWidth',0.5,'HandleVisibility','off');
    title(ax,'SON δD_{dyn} regression on Nino3.4 (400–600 hPa)', ...
        'FontSize',fs,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized','Position',[0,1.02,0], ...
        'Interpreter','tex');
    set(fig, 'PaperPositionMode','auto');
end
%% ============================================================
function idx = season_index(ts, yrs, mos, years_out, season)
% 逐年季节平均：'SON' = 当年 9–11 月；'DJF' = 当年 12 月 + 次年 1–2 月
    idx = NaN(numel(years_out),1);
    for k = 1:numel(years_out)
        y = years_out(k);
        if strcmp(season, 'DJF')
            sel = (yrs == y & mos == 12) | (yrs == y+1 & (mos == 1 | mos == 2));
        else
            sel = yrs == y & ismember(mos, [9 10 11]);
        end
        if sum(sel) == 3; idx(k) = mean(ts(sel), 'omitnan'); end
    end
end
%% ============================================================
function r = regress_out(y, x)
% 线性回归去除 x 的信号：r = y − [a + b·x]，保留 y 的均值；缺测处为 NaN
    y = y(:); x = x(:);
    r  = NaN(size(y));
    ok = isfinite(y) & isfinite(x);
    if sum(ok) < 3; return; end
    X  = [ones(sum(ok),1), x(ok) - mean(x(ok))];
    b  = X \ y(ok);
    r(ok) = y(ok) - X(:,2)*b(2);
end
%% ============================================================
function ts = area_mean_ts(F, lat_v)
% cos(lat) 加权区域平均；F: lat x lon x time
    w  = repmat(cosd(lat_v(:)), 1, size(F,2));
    ts = NaN(size(F,3),1);
    for t = 1:size(F,3)
        A  = F(:,:,t);
        ok = isfinite(A);
        if any(ok(:)); ts(t) = sum(A(ok).*w(ok)) / sum(w(ok)); end
    end
end
%% ============================================================
function dD_prof = rayleigh_column(T_K, p_hPa, dD0, aD, bD, cD, eps_w)
    n = numel(T_K);
    dD_prof = nan(n,1);
    start = find(~isnan(T_K),1,'first');
    if isempty(start); return; end
    dD_prof(start) = dD0;
    T_C = T_K - 273.15;
    es  = 6.112 * exp(17.67.*T_C./(T_C+243.5));
    qs  = eps_w.*es./(p_hPa-(1-eps_w).*es);
    for k = (start+1):n
        prev = find(~isnan(dD_prof(1:k-1)),1,'last');
        if isempty(prev); continue; end
        if isnan(T_K(k))||isnan(qs(k))||qs(k)<=0||qs(prev)<=0; continue; end
        F = qs(k)/qs(prev);
        if F<=0||F>=1; dD_prof(k)=dD_prof(prev); continue; end
        T_avg = 0.5*(T_K(k)+T_K(prev));
        alpha = exp(aD./T_avg.^2+bD./T_avg+cD);
        ratio = F.^(alpha-1);
        dD_prof(k) = (1+dD_prof(prev)/1000)*ratio*1000-1000;
    end
end
%% ============================================================
function data_2d = extract_season_mean_fn(dD_anom, years, months, ...
    yr, season_months, cross_year, lev_idx)
    [nlat, nlon, ~, ~] = size(dD_anom);
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
        idx = find((years == target_yr) & (months == mo));
        if isempty(idx); continue; end
        season_3d(:,:,im) = mean(dD_anom(:,:,lev_idx,idx(1)), 3, 'omitnan');
    end
    data_2d = mean(season_3d, 3, 'omitnan');
end
