%% ========================================================================
%  SON ω 与 850-hPa 水汽通量辐合：3 行 2 列（左列 = 趋势，右列 = 气候态）
%   (a) ω trend, 400–600 hPa           (b) ω climatology, 400–600 hPa
%   (c) ω trend, 15S–5N 经度-气压剖面   (d) ω climatology, 15S–5N 剖面
%   (e) 850-hPa MFC trend + (qu, qv)    (f) 850-hPa MFC climatology + (qu, qv)
%   色标：Precip4_diff_19lev（全部 6 个面板）
%   地图范围统一 30–120°E, 40°S–40°N；陆地统一为浅灰 [0.85 0.85 0.85]
%   (e)(f) 时段与 ω trend 一致（yr_start–yr_end）
%% ========================================================================
%% ===== 读取 ERA5 omega 数据 =====
nc_file = 'C:\Users\yingy\Downloads\Rice_mac\ERA5\era5_omega_90N90S_1deg_2000_2024.nc';

lat        = ncread(nc_file, 'latitude');
lon        = ncread(nc_file, 'longitude');
plev       = ncread(nc_file, 'pressure_level');
valid_time = ncread(nc_file, 'valid_time');
w          = ncread(nc_file, 'w');

dates  = datetime(1970,1,1,0,0,0) + seconds(double(valid_time));
years  = year(dates);
months = month(dates);
ntime  = numel(valid_time);

w = permute(w, [2 1 3 4]);
if lat(1) > lat(end)
    lat = flipud(lat);
    w   = flipud(w);
end
w = w * 864;

fprintf('omega数据覆盖: %d-%02d 到 %d-%02d\n', ...
    years(1), months(1), years(end), months(end));

%% ===== 空间分辨率 1° -> 5° 块平均 =====
factor = 5;
[nlat0, nlon0, nlev, ~] = size(w);
nlat_new = floor(nlat0 / factor);
nlon_new = floor(nlon0 / factor);
lat = lat(1:nlat_new*factor);
lon = lon(1:nlon_new*factor);
w   = w(1:nlat_new*factor, 1:nlon_new*factor, :, :);
lat = mean(reshape(lat, factor, nlat_new), 1)';
lon = mean(reshape(lon, factor, nlon_new), 1)';
w_r = reshape(w, factor, nlat_new, factor, nlon_new, nlev, ntime);
w   = squeeze(mean(mean(w_r, 1, 'omitnan'), 3, 'omitnan'));
clear w_r;
[nlat_full, nlon_full, nlev_full, ntime] = size(w);

%% ===== pressure levels =====
lev_mid  = find(plev >= 400 & plev <= 600);
lev_all  = 1:numel(plev);

%% ===== 区域范围 =====
lon_range = [30 120];
lat_range = [-40 40];

lon_idx = lon >= lon_range(1) & lon <= lon_range(2);
lat_idx = lat >= lat_range(1) & lat <= lat_range(2);
lon_sub = lon(lon_idx);
lat_sub = lat(lat_idx);
nlat    = sum(lat_idx);
nlon    = sum(lon_idx);

lat_xsec_idx_sub = lat_sub >= -15 & lat_sub <= 5;

w_sub = w(lat_idx, lon_idx, :, :);
clear w

%% ===== 计算月气候态 =====
fprintf('计算 omega 月气候态...\n');
climatology_full = NaN(nlat, nlon, nlev_full, 12);
for mo = 1:12
    idx = (months == mo);
    climatology_full(:,:,:,mo) = mean(w_sub(:,:,:,idx), 4, 'omitnan');
end
fprintf('月气候态计算完成！\n');

%% ===== 计算距平 =====
fprintf('计算距平...\n');
w_anom = NaN(size(w_sub));
for t = 1:ntime
    mo = months(t);
    w_anom(:,:,:,t) = w_sub(:,:,:,t) - climatology_full(:,:,:,mo);
end
clear w_sub
fprintf('距平计算完成！\n');

son_months = [9 10 11];

%% ===== Climatology数据 =====
% 左下：400-600hPa 水平分布
clim_son_mid = zeros(nlat, nlon);
for mo = son_months
    clim_son_mid = clim_son_mid + ...
        mean(climatology_full(:,:,lev_mid,mo), 3, 'omitnan') / numel(son_months);
end

% 右下：15S-5N 经度-气压剖面
clim_son_xsec = zeros(nlon, nlev_full);
for mo = son_months
    tmp = climatology_full(lat_xsec_idx_sub, :, :, mo);
    tmp_mean = squeeze(mean(tmp, 1, 'omitnan'));
    clim_son_xsec = clim_son_xsec + tmp_mean / numel(son_months);
end
clim_son_xsec = clim_son_xsec';   % [nlev nlon]
clear climatology_full

%% ===== Trend数据 =====
fprintf('计算 SON trend...\n');

% 提取SON年份列表
yr_start = 2003; yr_end = 2019;
reg_years = (yr_start:yr_end)';
n_reg = numel(reg_years);

% 提取每年SON平均
son_mid_annual  = NaN(nlat, nlon, n_reg);           % 400-600hPa [nlat nlon nyear]
son_xsec_annual = NaN(nlon, nlev_full, n_reg);      % 15S-5N xsec [nlon nlev nyear]

for iy = 1:n_reg
    yr = reg_years(iy);
    % 400-600hPa 水平
    tmp3d = zeros(nlat, nlon);
    for mo = son_months
        idx = find(years == yr & months == mo);
        if isempty(idx); continue; end
        tmp3d = tmp3d + mean(w_anom(:,:,lev_mid,idx(1)), 3, 'omitnan') / numel(son_months);
    end
    son_mid_annual(:,:,iy) = tmp3d;

    % 15S-5N 所有层
    tmp_xsec = zeros(nlon, nlev_full);
    for mo = son_months
        idx = find(years == yr & months == mo);
        if isempty(idx); continue; end
        tmp = squeeze(w_anom(lat_xsec_idx_sub,:,:,idx(1)));   % [nlat_band nlon nlev]
        tmp_mean = squeeze(mean(tmp, 1, 'omitnan'));            % [nlon nlev]
        tmp_xsec = tmp_xsec + tmp_mean / numel(son_months);
    end
    son_xsec_annual(:,:,iy) = tmp_xsec;
end
clear w_anom

% 逐格点线性trend
t_vec = (reg_years - mean(reg_years)) / std(reg_years);   % 标准化时间轴
X = [ones(n_reg,1), t_vec];

% 左上：400-600hPa trend [nlat nlon]
trend_mid  = NaN(nlat, nlon);
pval_mid   = NaN(nlat, nlon);
for ilat = 1:nlat
    for ilon = 1:nlon
        y = squeeze(son_mid_annual(ilat,ilon,:));
        valid = ~isnan(y);
        if sum(valid) < 4; continue; end
        y_v = y(valid); X_v = X(valid,:);
        b = X_v \ y_v;
        resid = y_v - X_v*b;
        df = sum(valid) - 2;
        s2 = sum(resid.^2)/df;
        se = sqrt(s2 * inv(X_v'*X_v));
        t_stat = b(2)/se(2,2);
        trend_mid(ilat,ilon)  = b(2) / std(reg_years);   % per year
        pval_mid(ilat,ilon)   = 2*(1-tcdf(abs(t_stat),df));
    end
end

% 右上：15S-5N xsec trend [nlev nlon]
trend_xsec = NaN(nlev_full, nlon);
for ilev = 1:nlev_full
    for ilon = 1:nlon
        y = squeeze(son_xsec_annual(ilon,ilev,:));
        valid = ~isnan(y);
        if sum(valid) < 4; continue; end
        y_v = y(valid); X_v = X(valid,:);
        b = X_v \ y_v;
        resid = y_v - X_v*b;
        df = sum(valid) - 2;
        s2 = sum(resid.^2)/df;
        se = sqrt(s2 * inv(X_v'*X_v));
        trend_xsec(ilev,ilon) = b(2) / std(reg_years);
    end
end
fprintf('Trend计算完成！\n');

%% ===== 850 hPa 水汽通量 (qu, qv) 与辐合 MFC = -∇·(qV) =====
%  算法与 SON MFC 气候态/趋势代码一致：逐月计算 → 气候平均 / 逐年 SON 平均后求趋势
fprintf('计算 850-hPa 水汽通量辐合...\n');
era_dir = 'C:\Users\yingy\Downloads\Rice_mac\ERA5';
ufile_m = fullfile(era_dir, 'era5_uwind_90N90S_1deg_2000_2024.nc');
vfile_m = fullfile(era_dir, 'era5_vwind_90N90S_1deg_2000_2024.nc');
qfile_m = fullfile(era_dir, 'era5_qs_90N90S_1deg_2000_2024.nc');

lon_m = double(ncread(ufile_m,'longitude'));
lat_m = double(ncread(ufile_m,'latitude'));
tm    = double(ncread(ufile_m,'valid_time'));
dm    = datetime(1970,1,1) + seconds(tm);
yr_m  = year(dm);  mo_m = month(dm);
it_m  = find(ismember(yr_m, reg_years) & ismember(mo_m, son_months));   % 只用 SON
yr_it = yr_m(it_m);
fprintf('MFC 使用月份数 = %d\n', numel(it_m));

u850 = read850(ufile_m, 'u', it_m);
v850 = read850(vfile_m, 'v', it_m);
q850 = read850(qfile_m, 'q', it_m);
QU = q850 .* u850;   QV = q850 .* v850;                  % kg kg^-1 m s^-1
clear u850 v850 q850

Re   = 6.371e6;
phiM = reshape(deg2rad(lat_m), 1, []);                   % 1 x nlat
dlon_m = abs(median(diff(lon_m)));
dQU  = (circshift(QU,-1,1) - circshift(QU,1,1)) / deg2rad(2*dlon_m);   % 经度周期
QVc  = QV .* cos(phiM);
dQV  = NaN(size(QV));
dQV(:,2:end-1,:) = (QVc(:,3:end,:) - QVc(:,1:end-2,:)) ./ (phiM(3:end) - phiM(1:end-2));
MFC  = -(dQU + dQV) ./ (Re * cos(phiM)) * 1e6;          % 10^-6 kg kg^-1 s^-1
clear dQU dQV QVc

% 气候态
MFC_clim = mean(MFC,3,'omitnan');
QU_clim  = mean(QU, 3,'omitnan') * 1e3;                  % g kg^-1 m s^-1
QV_clim  = mean(QV, 3,'omitnan') * 1e3;

% 逐年 SON 平均 → 趋势 (per decade) + 双侧 t 检验
nyr = numel(reg_years);
MFC_y = NaN([size(MFC,1) size(MFC,2) nyr]);  QU_y = MFC_y;  QV_y = MFC_y;
for iy = 1:nyr
    k = yr_it == reg_years(iy);
    MFC_y(:,:,iy) = mean(MFC(:,:,k),3,'omitnan');
    QU_y(:,:,iy)  = mean(QU(:,:,k), 3,'omitnan') * 1e3;
    QV_y(:,:,iy)  = mean(QV(:,:,k), 3,'omitnan') * 1e3;
end
clear MFC QU QV
[MFC_tr, MFC_p] = trend_decade(MFC_y, reg_years);
[QU_tr,  QU_p ] = trend_decade(QU_y,  reg_years);
[QV_tr,  QV_p ] = trend_decade(QV_y,  reg_years);
clear MFC_y QU_y QV_y

% → lat(南→北) x lon(0→360)，截取 30–120E, 40S–40N
[lonM, ix] = sort(mod(lon_m,360));
[latM, iy] = sort(lat_m);
fix2 = @(A) A(ix, iy)';
mlon = lonM >= lon_range(1) & lonM <= lon_range(2);
mlat = latM >= lat_range(1) & latM <= lat_range(2);
lon_mfc = lonM(mlon);  lat_mfc = latM(mlat);
cut = @(A) subsref(fix2(A), struct('type','()','subs',{{mlat, mlon}}));

MFC_clim = cut(MFC_clim);  QU_clim = cut(QU_clim);  QV_clim = cut(QV_clim);
MFC_tr   = cut(MFC_tr);    MFC_p   = cut(MFC_p);
QU_tr    = cut(QU_tr);     QU_p    = cut(QU_p);
QV_tr    = cut(QV_tr);     QV_p    = cut(QV_p);
fprintf('MFC 计算完成\n');

%% ===== topo + 陆地mask =====
load topo;
topo_lon     = 0.5:359.5;
topo_lat     = -89.5:89.5;
topo_ext     = [topo, topo(:,1)];
topo_lon_ext = [topo_lon, 360.5];
[TLON, TLAT] = meshgrid(topo_lon, topo_lat);
[OLON, OLAT] = meshgrid(lon_sub, lat_sub);
land_mask = interp2(TLON, TLAT, topo, OLON, OLAT, 'nearest', 0) > 0;

% MFC 1° 格点的陆地 mask（最近 1° 地形格点，同 MFC 原代码）
[OLONm, OLATm] = meshgrid(lon_mfc, lat_mfc);
land_mfc = interp2(TLON, TLAT, topo, OLONm, OLATm, 'nearest', 0) > 0;
MFC_clim(land_mfc) = NaN;  QU_clim(land_mfc) = NaN;  QV_clim(land_mfc) = NaN;
MFC_tr(land_mfc)   = NaN;  QU_tr(land_mfc)   = NaN;  QV_tr(land_mfc)   = NaN;
MFC_p(land_mfc)    = NaN;  QU_p(land_mfc)    = NaN;  QV_p(land_mfc)    = NaN;

topo_pack.lon  = topo_lon_ext;
topo_pack.lat  = topo_lat;
topo_pack.elev = topo_ext;


% topo 裁剪
buf = 5;
tlon_idx    = topo_pack.lon >= lon_range(1)-buf & topo_pack.lon <= lon_range(2)+buf;
tlat_idx    = topo_pack.lat >= lat_range(1)-buf & topo_pack.lat <= lat_range(2)+buf;
topo_lon_v  = topo_pack.lon(tlon_idx);
topo_lat_v  = topo_pack.lat(tlat_idx);
topo_elev_v = topo_pack.elev(tlat_idx, tlon_idx);

%% ===== 色标范围 =====
clim_trend = max(abs([nanmax(trend_mid(:)),  nanmin(trend_mid(:)), ...
                      nanmax(trend_xsec(:)), nanmin(trend_xsec(:))]));
clim_trend = ceil(clim_trend * 10) / 10;   % 取整到0.1

clim_clim  = max(abs([nanmax(clim_son_mid(:)),  nanmin(clim_son_mid(:)), ...
                      nanmax(clim_son_xsec(:)), nanmin(clim_son_xsec(:))]));
clim_clim  = ceil(clim_clim / 5) * 5;

fprintf('Trend色标: ±%.2f hPa/day/yr\n', clim_trend);
fprintf('Clim色标:  ±%.1f hPa/day\n',    clim_clim);


%% ===== 色标：Precip4_diff_19lev =====
% ω 面板 (a)–(d)：沿用原色标 BlueWhiteOrangeRed
Sw     = load('BlueWhiteOrangeRed.mat');
cmap_w = Sw.BlueWhiteOrangeRed;
% MFC 面板 (e)(f)：Precip4_diff_19lev（文件名 P 大写，变量名 precip4_diff_19lev 小写）
Sm     = load('Precip4_diff_19lev.mat');
cmap_m = Sm.precip4_diff_19lev;

%% ===== 画图：3行2列（左列 = 趋势，右列 = 气候态）=====
cm2px    = 96/2.54;
fig_w_cm = 18;
pw_cm = 5.4;   ph_cm = 4.23;   row_dy = 5.17;  % 行间距 = 面板高 + 标题/刻度
y3_cm = 1.0;
y2_cm = y3_cm + row_dy;
y1_cm = y2_cm + row_dy;
fig_h_cm = y1_cm + ph_cm + 0.7;

fig_w = round(fig_w_cm * cm2px);
fig_h = round(fig_h_cm * cm2px);
fig = figure('Color','w', ...
             'Position',[50 50 fig_w fig_h], ...
             'Renderer','opengl', ...
             'GraphicsSmoothing','off', ...
             'PaperUnits','centimeters', ...
             'PaperSize',[fig_w_cm, fig_h_cm], ...
             'PaperPositionMode','auto');

%% ---- 布局参数（归一化）----
pw      = pw_cm / fig_w_cm;
ph      = ph_cm / fig_h_cm;
left_x  = 0.07;
right_x = 0.565;                       % 左列色标 + 右列纵坐标留出空间
row_y   = [y1_cm y2_cm y3_cm] / fig_h_cm;
cb_w    = 0.015;  cb_gap = 0.012;
cbL_x   = left_x  + pw + cb_gap;       % 左列（趋势）色标
cbR_x   = right_x + pw + cb_gap;       % 右列（气候态）色标
ttl_w   = pw + 0.1;
silver_rgb = [0.85 0.85 0.85];
ylat_lab = {'40°S','20°S','0°','20°N','40°N'};

%% ---- (a) ω trend, 400–600 hPa ----
pos = [left_x, row_y(1), pw, ph];
ax  = axes('Position', pos);
draw_wmap(ax, lon_sub, lat_sub, trend_mid, land_mask, pval_mid, silver_rgb, ...
          cmap_w, clim_trend, topo_lon_v, topo_lat_v, topo_elev_v, lon_range, lat_range);
ax.Position = pos;
add_title(fig, [left_x, row_y(1)+ph+0.005, ttl_w, 0.03], '(a) SON \omega Trend, 400–600 hPa');

%% ---- (b) ω climatology, 400–600 hPa ----
pos = [right_x, row_y(1), pw, ph];
ax  = axes('Position', pos);
draw_wmap(ax, lon_sub, lat_sub, clim_son_mid, land_mask, [], silver_rgb, ...
          cmap_w, clim_clim, topo_lon_v, topo_lat_v, topo_elev_v, lon_range, lat_range);
ax.Position = pos;
add_title(fig, [right_x, row_y(1)+ph+0.005, ttl_w, 0.03], '(b) SON \omega Climatology, 400–600 hPa');

%% ---- (c) ω trend, 15S–5N section ----
pos = [left_x, row_y(2), pw, ph];
ax  = axes('Position', pos);
ti  = clim_trend / 12;
draw_xsec(ax, lon_sub, double(plev), trend_xsec, cmap_w, clim_trend, ...
          [-clim_trend-1, -clim_trend:ti/2:clim_trend, clim_trend+1], ...   % 原图填色级数
          ti:ti:clim_trend, -clim_trend:ti:-ti, lon_range);
ax.Position = pos;
add_title(fig, [left_x, row_y(2)+ph+0.005, ttl_w, 0.03], '(c) SON \omega Trend, 15S–5N Mean');

%% ---- (d) ω climatology, 15S–5N section ----
pos = [right_x, row_y(2), pw, ph];
ax  = axes('Position', pos);
ci  = 6.25;
draw_xsec(ax, lon_sub, double(plev), clim_son_xsec, cmap_w, clim_clim, ...
          [-clim_clim-1, -clim_clim:ci/2:clim_clim, clim_clim+1], ...      % 原图填色级数
          ci:ci:clim_clim, -40:ci:-ci, lon_range);
ax.Position = pos;
add_title(fig, [right_x, row_y(2)+ph+0.005, ttl_w, 0.03], '(d) SON \omega Climatology, 15S–5N Mean');

%% ---- 第 1–2 行色标：左列 ω trend，右列 ω climatology（各跨两行）----
cb_h  = row_y(1) + ph - row_y(2);        % 两行总高度
cbL_s = 0.75;                            % (a)(c) 色标长度 = 原来的 0.75 倍，上下居中
cbL_h = cb_h * cbL_s;
cbL_y = row_y(2) + (cb_h - cbL_h)/2;
add_cbar([cbL_x, cbL_y, cb_w, cbL_h], cmap_w, clim_trend, [], ...
         '\omega trend (hPa day^{-1} yr^{-1})');
add_cbar([cbR_x, cbL_y, cb_w, cbL_h], cmap_w, clim_clim, [], ...     % (b)(d) 同样缩短为 0.75 倍并居中
         '\omega climatology (hPa day^{-1})');

%% ---- (e) MFC trend（左） / (f) MFC climatology（右）----
vec_skip  = 6;                         % 矢量间隔（1° 格点 → 6°）
sig_skip  = 2;                         % 显著性打点间隔 2°
clim_ef   = 0.03;                      % (e)(f) 共用色标 ±0.03
blank_f   = 0.04;                      % 气候态中 |MFC| > 0.04 的格点画成白色
iq = 1:vec_skip:numel(lon_mfc);  jq = 1:vec_skip:numel(lat_mfc);
[LONq, LATq] = meshgrid(lon_mfc(iq), lat_mfc(jq));
is = 1:sig_skip:numel(lon_mfc);  js = 1:sig_skip:numel(lat_mfc);
[LONs, LATs] = meshgrid(lon_mfc(is), lat_mfc(js));
land_img_m = repmat(reshape(silver_rgb,1,1,3), [size(land_mfc) 1]);

% (e)(f) 各自色标的标签（两行，避免超出面板高度）
cbar_label_e = {'MFC trend (10^{-6} kg kg^{-1}', 's^{-1} decade^{-1})'};
cbar_label_f = {'MFC climatology', '(10^{-6} kg kg^{-1} s^{-1})'};

for k = 1:2
    if k == 1                          % (e) 趋势，左列
        xpos = left_x;   dat = MFC_tr;   U = QU_tr;   V = QV_tr;
        Usig = (QU_p < 0.05) | (QV_p < 0.05);         % qu 或 qv 显著
        ttl  = sprintf('(e) SON 850-hPa MFC Trend (%d–%d)', yr_start, yr_end);
        ref_unit = 'g kg^{-1} m s^{-1} decade^{-1}';
    else                               % (f) 气候态，右列
        xpos = right_x;  dat = MFC_clim; U = QU_clim; V = QV_clim;
        Usig = [];
        ttl  = sprintf('(f) SON 850-hPa MFC Climatology (%d–%d)', yr_start, yr_end);
        ref_unit = 'g kg^{-1} m s^{-1}';
    end
    pos = [xpos, row_y(3), pw, ph];
    ax  = axes('Position', pos);  hold(ax,'on');

    bg = land_img_m;
    if k == 2
        out = abs(dat) > blank_f & ~land_mfc;        % 超出 ±0.04 → 白色
        dat(out) = NaN;
        for c3 = 1:3
            tmp = bg(:,:,c3);  tmp(out) = 1;  bg(:,:,c3) = tmp;
        end
    end
    imagesc(ax, lon_mfc, lat_mfc, bg);
    imagesc(ax, lon_mfc, lat_mfc, dat, 'AlphaData', ~isnan(dat));
    set(ax, 'YDir','normal', 'Color',silver_rgb);
    colormap(ax, cmap_m);
    caxis(ax, [-clim_ef clim_ef]);

    % 显著性打点（仅趋势图）
    if k == 1
        sg = MFC_p(js,is) < 0.05 & ~isnan(dat(js,is));
        plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize',1.5, ...
             'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none');
    end

    % 矢量：统一手动比例（参考箭头 = 90% 分位，长度 ref_len 度）
    Uq = U(jq,iq);  Vq = V(jq,iq);
    mag = hypot(Uq(:), Vq(:));  mag = mag(isfinite(mag));
    ref = nice_ceil(prctile(mag, 90));
    ref_len = 5;  sc = ref_len / ref;
    if isempty(Usig)
        quiver(ax, LONq, LATq, Uq*sc, Vq*sc, 0, 'k', 'LineWidth',1.1, 'MaxHeadSize',0.6);
    else
        S  = Usig(jq,iq);
        Un = Uq;  Vn = Vq;  Un(S) = NaN;  Vn(S) = NaN;          % 不显著：黑
        Us = Uq;  Vs = Vq;  Us(~S) = NaN; Vs(~S) = NaN;         % 显著：紫
        quiver(ax, LONq, LATq, Un*sc, Vn*sc, 0, 'k', 'LineWidth',1.1, 'MaxHeadSize',0.6);
        quiver(ax, LONq, LATq, Us*sc, Vs*sc, 0, 'Color',[0.75 0 0.95], ...
               'LineWidth',1.5, 'MaxHeadSize',0.6);
    end

    contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth',0.5);
    xlim(ax, lon_range); ylim(ax, lat_range);
    set(ax, 'XTick',30:30:120, 'XTickLabel',{'30°E','60°E','90°E','120°E'}, ...
            'YTick',-40:20:40, 'YTickLabel',ylat_lab, ...
            'FontSize',10, 'TickDir','out', 'Box','on', 'Layer','top');
    ax.Position = pos;

    draw_ref_arrow(ax, lon_range, lat_range, ref_len, sprintf('%g', ref));
    add_title(fig, [xpos, row_y(3)+ph+0.005, ttl_w, 0.02], ttl);
    fprintf('(%c) 参考箭头 = %g %s\n', 'e'+k-1, ref, ref_unit);
end

%% ---- 第 3 行色标：(e)(f) 各一个，高度与面板一致 ----
tk_ef = clim_ef*[-1 -0.5 0 0.5 1];
add_cbar([cbL_x, row_y(3), cb_w, ph], cmap_m, clim_ef, tk_ef, cbar_label_e);   % (e) 右侧
add_cbar([cbR_x, row_y(3), cb_w, ph], cmap_m, clim_ef, tk_ef, cbar_label_f);   % (f) 右侧

%% ===== 刻度线加长（所有坐标区）=====
set(findall(fig,'Type','axes'), 'TickLength',[0.025 0.035]);

%% ===== 保存 =================================================
out_dir = 'C:\Users\yingy\Downloads\Rice_mac';
if ~exist(out_dir,'dir'); mkdir(out_dir); end
fname = fullfile(out_dir, 'Fig2_omega_MFC_SON_6panels_v2.tif');
exportgraphics(fig, fname, 'Resolution',600, 'BackgroundColor','white');
fprintf('Saved: %s\n', fname);
pngname = fullfile(out_dir, 'Fig2_omega_MFC_SON_6panels_v2.png');
exportgraphics(fig, pngname, 'Resolution',600, 'BackgroundColor','white');
fprintf('Saved: %s\n', pngname);

%% ===== 子函数 ===============================================
function F = read850(file, var, it)
% 读 850 hPa 单层，只保留 it 指定的月份，返回 lon x lat x time
lev = double(ncread(file,'pressure_level'));
[~, il] = min(abs(lev - 850));
if abs(lev(il) - 850) > 0.1; error('%s 中没有 850 hPa', file); end
F = double(ncread(file, var, [1 1 il 1], [Inf Inf 1 Inf]));
F = reshape(F(:,:,1,it), size(F,1), size(F,2), numel(it));
end

function [slope, pval] = trend_decade(Y, years)
% 逐格点 OLS 趋势（per decade）+ 双侧 t 检验；有缺测年份的格点为 NaN
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

function y = nice_ceil(x)
if ~isfinite(x) || x <= 0; y = 1; return; end
e = 10^floor(log10(x));  m = x / e;
if m <= 1; y = e; elseif m <= 2; y = 2*e; elseif m <= 5; y = 5*e; else; y = 10*e; end
end

function draw_ref_arrow(ax, lon_range, lat_range, ref_len, label)
% 面板右上角白底框内画参考箭头和数值
y0 = lat_range(2) - 6;
ht = text(ax, lon_range(2) - 1.5, y0, label, 'FontSize',10, ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle');
ext = ht.Extent;  delete(ht);
x0 = ext(1) - ref_len - 1.5;
yb = y0 - ext(4)/2 - 1.5;
rectangle(ax, 'Position',[x0-1.5, yb, lon_range(2)-(x0-1.5), lat_range(2)-yb], ...
          'FaceColor','w', 'EdgeColor','k', 'LineWidth',0.5);
quiver(ax, x0, y0, ref_len, 0, 0, 'k', 'LineWidth',1.1, 'MaxHeadSize',0.8);
text(ax, lon_range(2) - 1.5, y0, label, 'FontSize',10, ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle');
end

function draw_wmap(ax, lon, lat, dat, land_mask, pval, silver_rgb, cmap, cl, ...
                   tlon, tlat, telev, lon_range, lat_range)
% ω 400–600 hPa 平面图：陆地浅灰，海洋 imagesc，可选显著性打点
hold(ax,'on');
dat(land_mask) = NaN;
land_img = repmat(reshape(silver_rgb,1,1,3), [size(land_mask) 1]);
imagesc(ax, lon, lat, land_img);
imagesc(ax, lon, lat, dat, 'AlphaData', ~isnan(dat));
set(ax, 'YDir','normal', 'Color',silver_rgb);
colormap(ax, cmap);
caxis(ax, [-cl cl]);
contour(ax, tlon, tlat, telev, [0 0], 'k-', 'LineWidth',0.5);
set(ax, 'Clipping','on');
if ~isempty(pval)
    [LO, LA] = meshgrid(lon, lat);
    sg = ~land_mask & pval < 0.05;
    plot(ax, LO(sg), LA(sg), 'o', 'MarkerSize',1.5, ...
         'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none');
end
xlim(ax, lon_range); ylim(ax, lat_range);
set(ax, 'XTick',30:30:120, 'XTickLabel',{}, ...
        'YTick',-40:20:40, 'YTickLabel',{'40°S','20°S','0°','20°N','40°N'}, ...
        'FontSize',10, 'TickDir','out', 'Box','on', 'Layer','top');
end

function draw_xsec(ax, lon, plev, X, cmap, cl, lev_fill, lev_pos, lev_neg, lon_range)
% 15S–5N 经度-气压剖面：lev_fill = 填色级数，lev_pos/lev_neg = 实线/虚线等值线
hold(ax,'on');
[~, hf] = contourf(ax, lon, plev, X, 'LevelList', lev_fill);
hf.LineColor = 'none';
colormap(ax, cmap);
caxis(ax, [-cl cl]);
contour(ax, lon, plev, X, lev_pos, 'k-',  'LineWidth',0.8);
contour(ax, lon, plev, X, lev_neg, 'k--', 'LineWidth',0.8);
set(ax, 'YDir','reverse');
xlim(ax, lon_range); ylim(ax, [200 1000]);
set(ax, 'XTick',30:30:120, 'XTickLabel',{}, ...
        'YTick',[200 300 400 500 600 700 850 1000], ...
        'YTickLabel',{'200','300','400','500','600','700','850','1000'}, ...
        'FontSize',10, 'TickDir','out', 'Box','on', 'Layer','top');
ylabel(ax, 'hPa', 'FontSize',10);
end

function add_title(fig, pos, str)
annotation(fig, 'textbox', pos, 'String', str, ...
           'FontSize',10, 'FontWeight','bold', ...
           'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
           'EdgeColor','none', 'FitBoxToText','off', 'Interpreter','tex');
end

function add_cbar(pos, cmap, cl, ticks, lab)
% 竖向色标（挂在一个隐藏坐标区上）
axc = axes('Position', pos, 'Visible','off');
colormap(axc, cmap);
caxis(axc, [-cl cl]);
cb = colorbar(axc, 'eastoutside');
cb.Position = pos;
if ~isempty(ticks); cb.Ticks = ticks; end
cb.Label.String   = lab;
cb.FontSize       = 10;
cb.Label.FontSize = 10;
cb.TickDirection  = 'out';
cb.AxisLocation   = 'out';                % 刻度和标签都在色标右侧
end
