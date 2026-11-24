%% ========================================================================
%  2 列 4 行组合图（0–120°E, 20°S–40°N，保持真实经纬度长宽比）
%    左列 = 气候态 / 右列 = 趋势（g、h 两个都是趋势）；每个面板右侧一个 colorbar
%
%  (a) ERA5 SON 850-hPa 水汽通量辐合 MFC = −∇·(qV) 气候态 + (qu, qv) 矢量
%                                                      10^-6 kg kg^-1 s^-1
%  (b) ERA5 SON 850-hPa MFC 趋势 + (qu, qv) 趋势矢量    10^-6 kg kg^-1 s^-1 decade^-1
%      矢量：qu 或 qv 趋势显著 (p<0.05) 为紫色，否则黑色
%  (c) ERA5 SON 850-hPa 比湿 q 气候态                   g kg^-1
%  (d) ERA5 SON 850-hPa 比湿 q 趋势                     g kg^-1 decade^-1
%  (e) AIRS SON δD 气候态, ~850 hPa (825 + 908 hPa)      ‰
%  (f) 动力项 −V′·∇δD_clim, ~850 hPa                     ‰ day^-1 decade^-1
%      V′ = ERA5 800–900 hPa 风场的 SON 趋势（插值到 AIRS 5° 格点）
%      δD_clim = (e) 的 SON 气候态；绿点 = u 或 v 趋势显著 (p<0.05)
%  (g) ERA5 SON 对流降水 (CP) 趋势                     mm day^-1 decade^-1
%  (h) ERA5 SON P−E 趋势                                mm day^-1 decade^-1
%  所有面板同为 SON 2003–2019
%
%  版面: 画布宽 18 cm，面板宽 ≈ 6.6 cm、高 ≈ 3.3 cm，画布高 ≈ 17.7 cm；
%        colorbar 都在右侧，每个 5 个刻度标签。
%  陆地: 全部 8 个面板屏蔽为 NaN，显示为浅灰 [0.85 0.85 0.85]，只画海岸线。
%  显著性: 趋势图绿点 p < 0.05（双侧 t 检验）
%  色标: (a)(b)(c)(d)(g)(h) Precip4_diff_19lev（(c) 只用其湿端一半）；
%        (e)(f) BlueWhiteOrangeRed
%  原 (c)(d) 总平流趋势、δD 通量矢量、400–600 hPa q 的代码已整块注释（未删除）
%% ========================================================================
clear; clc; close all;

%% ========================================================================
% 0. 设置
% ========================================================================
airs_file = 'C:\Users\yingy\Downloads\Isotope-IOD\AIRS_RS_monthly_5deg.nc';
era_dir   = 'C:\Users\yingy\Downloads\Rice_mac\ERA5';
ufile     = fullfile(era_dir,'era5_uwind_90N90S_1deg_2000_2024.nc');
vfile     = fullfile(era_dir,'era5_vwind_90N90S_1deg_2000_2024.nc');
qfile     = fullfile(era_dir,'era5_qs_90N90S_1deg_2000_2024.nc');
ep_file   = 'C:\Users\yingy\Downloads\ERA5_monthly_ep_cp_lp_tp.nc';
out_dir   = 'C:\Users\yingy\Downloads\Rice_mac';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

season_month = [9 10 11];
season_name  = 'SON';
trend_years  = [2003 2019];   % 全部 6 个面板，沿用 Fig.5 代码

lon_range = [0 120];
lat_range = [-20 40];       % 20°S–40°N

p_mid_airs = [400 600];       % AIRS 510.9 + 421.7 hPa
p_mid_era  = [400 600];ncd
p_low_airs = [800 910];       % AIRS 825.4 + 908.5 hPa
p_low_era  = [800 900];       % ERA5 800–900 hPa 内可用层（q 和风都用这一层）
p_850      = [850 850];       % 850 hPa 单层（q、u、v，MFC 与比湿面板）
vec_spacing_m = 6;            % (a)(b) 水汽通量矢量间隔（度）
sig_deg_mfc = 2.5;            % MFC 显著性打点间隔（度）
elev_max   = 1000;            % m，~850 hPa 面板屏蔽格box平均海拔超过此值的格点

days_per_month = 30;          % 沿用 Fig.5 代码的换算，见第 5 节说明
sig_deg_ep  = 2.5;            % CP、P−E 显著性打点间隔（度）
sig_deg_adv = 5;              % 平流趋势显著性打点间隔（度，AIRS 5° 格点即每格一个）
vec_spacing = 10;             % (c)(d) 风矢量间隔（度）
sig_deg_q   = 5;              % q 显著性打点间隔（度）
ALPHA       = 0.05;

pressure_levels = [1017.2227, 1000, 908.51398, 825.40198, 749.89301, ...
                   681.29102, 618.96600, 510.89801, 421.69800, 348.06900, ...
                   287.29800, 237.13699, 177.82899, 133.35201, 74.989601, ...
                   28.729900, 0.10000];

% 目标月份（同一年的 3 个月相邻排列）
[mm,yy]  = ndgrid(season_month, trend_years(1):trend_years(2));
key_tr   = yy(:)*100 + mm(:);
yrs_tr   = (trend_years(1):trend_years(2))';
key_adv  = key_tr;            % 平流与其他面板用同一组月份

load topo;                                     % 180 x 360, 1°, 海拔 (m)
topo_lon = 0.5:359.5;
topo_lat = -89.5:89.5;

%% ========================================================================
% 1. AIRS δD
% ========================================================================
lat_a = double(ncread(airs_file,'lat'));
lon_a = double(ncread(airs_file,'lon'));
dD    = double(ncread(airs_file,'dD'));        % lon x lat x lev x time
dD    = permute(dD,[2 1 3 4]);                 % lat x lon x lev x time
[lat_a,ia] = sort(lat_a(:));   dD = dD(ia,:,:,:);
lon_a      = mod(lon_a(:),360);
[lon_a,io] = sort(lon_a);      dD = dD(:,io,:,:);
[nlat_a,nlon_a,nlev_a,nt_a] = size(dD);
if nlev_a ~= numel(pressure_levels)
    error('AIRS 层数 (%d) 与 pressure_levels (%d) 不一致', nlev_a, numel(pressure_levels));
end
dates_a = datetime(2002,9,1) + calmonths(0:nt_a-1);
key_a   = year(dates_a(:))*100 + month(dates_a(:));

[ok,it_a] = ismember(key_adv, key_a);
if any(~ok); error('AIRS 缺少 %d 个目标月份', sum(~ok)); end

D_mid = airs_layer(dD, pressure_levels, p_mid_airs, it_a);
D_low = airs_layer(dD, pressure_levels, p_low_airs, it_a);
clear dD

dlon_a = median(diff(lon_a));
dlat_a = median(diff(lat_a));

%% ========================================================================
% 2. ERA5 风场 → AIRS 格点（平流用）
% ========================================================================
lon_e0 = double(ncread(ufile,'longitude'));
lat_e0 = double(ncread(ufile,'latitude'));
lev_e  = double(ncread(ufile,'pressure_level'));
fprintf('ERA5 可用层: %s hPa\n', mat2str(lev_e(:)'));

[U_mid, lat_e, lon_e] = read_era_layer(ufile,'u',lev_e,p_mid_era,key_adv,lon_e0,lat_e0);
V_mid                 = read_era_layer(vfile,'v',lev_e,p_mid_era,key_adv,lon_e0,lat_e0);
U_low                 = read_era_layer(ufile,'u',lev_e,p_low_era,key_adv,lon_e0,lat_e0);
V_low                 = read_era_layer(vfile,'v',lev_e,p_low_era,key_adv,lon_e0,lat_e0);

Wlon = box_weights(lon_e, lon_a, dlon_a, true);
Wlat = box_weights(lat_e, lat_a, dlat_a, false);
U_mid = regrid(U_mid, Wlat, Wlon);  V_mid = regrid(V_mid, Wlat, Wlon);
U_low = regrid(U_low, Wlat, Wlon);  V_low = regrid(V_low, Wlat, Wlon);

%% ========================================================================
% 3. AIRS δD 气候态 (e) 与平流动力项 (f)，~850 hPa
% ========================================================================
elev_a = box_weights(topo_lat, lat_a, dlat_a, false) * double(topo) * ...
         box_weights(topo_lon, lon_a, dlon_a, true)';
% (c) ~850 hPa 不再屏蔽高地形：按要求陆地全部显示。注意青藏高原、东非高原等处
%     ~850 hPa 位于地面以下，那里的平流趋势是外推值，解读时请忽略。
% mask_a = repmat(elev_a > elev_max, 1, 1, size(D_low,3));
% D_low(mask_a) = NaN;  U_low(mask_a) = NaN;  V_low(mask_a) = NaN;

nmon = numel(season_month);

% ---- 原 (c)(d)：总平流趋势与 δD 通量矢量，暂不使用，整块注释（未删除）----
%{
[trend_am, pval_am] = adv_trend(D_mid, U_mid, V_mid, lat_a, lon_a, nmon, yrs_tr);
[trend_al, pval_al] = adv_trend(D_low, U_low, V_low, lat_a, lon_a, nmon, yrs_tr);

% (c)(d) 叠加的矢量：δD 通量 (δD·u, δD·v) 的 SON 平均（‰ m s^-1，AIRS 格点）
%   逐月相乘后再平均，即 mean(δD·u)，与 dD_flux_convergence 脚本 (a) 的矢量一致
%   注意 δD < 0，矢量方向与风向相反
Um_c = mean(D_mid .* U_mid, 3, 'omitnan');  Vm_c = mean(D_mid .* V_mid, 3, 'omitnan');
Ul_c = mean(D_low .* U_low, 3, 'omitnan');  Vl_c = mean(D_low .* V_low, 3, 'omitnan');
%}

% ---- (e) δD 气候态（SON 2003–2019，~850 hPa）----
dD_clim = mean(D_low, 3, 'omitnan');                        % lat x lon（AIRS 5°）
dD_clim(sum(~isnan(D_low),3) < 0.5*size(D_low,3)) = NaN;    % 有效月份不足一半的格点设为 NaN

% ---- (f) 动力项：风场趋势 × 气候态 δD 梯度 = −V′·∇δD_clim ----
[u_tr, pu] = trend_decade(son_mean(U_low, nmon), yrs_tr);   % m s^-1 decade^-1
[v_tr, pv] = trend_decade(son_mean(V_low, nmon), yrs_tr);
dyn_term   = dyn_adv(dD_clim, u_tr, v_tr, lat_a, lon_a);    % ‰ day^-1 decade^-1
pval_dyn   = min(pu, pv);                                   % u 或 v 趋势显著即打点
fprintf('动力项范围: %.2f 到 %.2f ‰ day^-1 decade^-1\n', min(dyn_term(:)), max(dyn_term(:)));
clear D_mid D_low U_mid V_mid U_low V_low

%% ========================================================================
% 4. ERA5 850 hPa：比湿 q 与水汽通量辐合 MFC = −∇·(qV)
%    （MFC 算法同 Fig. 2 原 (e)(f)：逐月计算 → 逐年 SON 平均 → 气候态 / 趋势）
% ========================================================================
lon_q0 = double(ncread(qfile,'longitude'));
lat_q0 = double(ncread(qfile,'latitude'));
lev_q  = double(ncread(qfile,'pressure_level'));

% ---- 原 400–600 hPa 与 800–900 hPa q 趋势：暂不使用，整块注释（未删除）----
%{
[q_mid, lat_q, lon_q] = read_era_layer(qfile,'q',lev_q,p_mid_era,key_tr,lon_q0,lat_q0);
q_low                 = read_era_layer(qfile,'q',lev_q,p_low_era,key_tr,lon_q0,lat_q0);
%}

% 全球场（MFC 的经向差分需要周期经度）；lat(南→北) x lon(0→360) x time
[q850, lat_q, lon_q] = read_era_layer(qfile,'q',lev_q,p_850,key_tr,lon_q0,lat_q0);
u850 = read_era_layer(ufile,'u',lev_e,p_850,key_tr,lon_e0,lat_e0);
v850 = read_era_layer(vfile,'v',lev_e,p_850,key_tr,lon_e0,lat_e0);
if ~isequal(size(u850), size(q850)) || ~isequal(size(v850), size(q850))
    error('u、v、q 的网格不一致，无法直接计算 qV');
end
QU = q850 .* u850;   QV = q850 .* v850;                  % kg kg^-1 m s^-1
clear u850 v850

Re   = 6.371e6;
phiq = deg2rad(lat_q(:));
lamq = deg2rad(lon_q(:)');
cq   = cos(phiq);
dlq  = median(diff(lon_q));
perq = abs(mod(lon_q(end) - lon_q(1) + dlq, 360)) < 1e-3;
MFC  = -(zonal_deriv(QU, lamq, perq) + merid_deriv(QV .* cq, phiq)) ./ (Re * cq) * 1e6;  % 10^-6 kg kg^-1 s^-1

% 截取区域
iqlon = lon_q >= lon_range(1) & lon_q <= lon_range(2);
iqlat = lat_q >= lat_range(1) & lat_q <= lat_range(2);
lon_qs = lon_q(iqlon);  lat_qs = lat_q(iqlat);
MFC  = MFC(iqlat,iqlon,:);
QU   = QU(iqlat,iqlon,:)   * 1e3;                         % g kg^-1 m s^-1
QV   = QV(iqlat,iqlon,:)   * 1e3;
q850 = q850(iqlat,iqlon,:) * 1e3;                         % g kg^-1

% 逐年 SON 平均 → 气候态 / 趋势 (per decade) + 双侧 t 检验
MFC_y = son_mean(MFC,  nmon);   QU_y = son_mean(QU, nmon);   QV_y = son_mean(QV, nmon);
q_y   = son_mean(q850, nmon);
clear MFC QU QV q850

MFC_clim = mean(MFC_y, 3, 'omitnan');
QU_clim  = mean(QU_y,  3, 'omitnan');   QV_clim = mean(QV_y, 3, 'omitnan');
q_clim   = mean(q_y,   3, 'omitnan');
[MFC_tr, MFC_p] = trend_decade(MFC_y, yrs_tr);
[QU_tr,  QU_p ] = trend_decade(QU_y,  yrs_tr);
[QV_tr,  QV_p ] = trend_decade(QV_y,  yrs_tr);
[trend_q, pval_q] = trend_decade(q_y, yrs_tr);
clear MFC_y QU_y QV_y q_y

%% ========================================================================
% 5. ERA5 CP 与 P−E 趋势
%
%  单位换算沿用 Fig.5 代码：× 1000 / days_per_month。
%  注意：CDS 的 "monthly averaged reanalysis" 里，累积量 (tp, cp, e) 是
%  "日累积量的月平均"，单位已是 m day^-1，只需 × 1000 就是 mm day^-1，
%  再除以 30 会小 30 倍。脚本下面会打印区域平均降水，请核对：
%  热带印度洋 SON 气候态降水应在几 mm/day 量级。若打印值只有 0.1–0.3，
%  说明应把 days_per_month 改成 1。
% ========================================================================
lat_ep0 = double(ncread(ep_file,'latitude'));
lon_ep0 = double(ncread(ep_file,'longitude'));
vt      = double(ncread(ep_file,'valid_time'));
dates_ep = datetime(1970,1,1) + seconds(vt);
key_ep   = year(dates_ep(:))*100 + month(dates_ep(:));
[ok,it_ep] = ismember(key_tr, key_ep);
if any(~ok); error('P−E 文件缺少 %d 个目标月份', sum(~ok)); end

jb = find(lat_ep0 >= lat_range(1) & lat_ep0 <= lat_range(2));
j1 = min(jb); nj = max(jb) - j1 + 1;
nT = numel(key_tr);
tp = NaN(numel(lon_ep0), nj, nT); ev = tp; cp = tp;
for t = 1:nT
    tp(:,:,t) = double(ncread(ep_file,'tp',[1 j1 it_ep(t)],[Inf nj 1]));
    ev(:,:,t) = double(ncread(ep_file,'e', [1 j1 it_ep(t)],[Inf nj 1]));
    cp(:,:,t) = double(ncread(ep_file,'cp',[1 j1 it_ep(t)],[Inf nj 1]));
end
tp_mm =  tp * 1000 / days_per_month;
e_mm  = -ev * 1000 / days_per_month;          % ERA5 蒸发为负
cp_mm =  cp * 1000 / days_per_month;
PE_mm = tp_mm - e_mm;
clear tp ev cp

[PE_mm, lat_ep, lon_ep] = to_latlon_sorted(PE_mm, lat_ep0(j1:j1+nj-1), lon_ep0);
[cp_mm, ~, ~]           = to_latlon_sorted(cp_mm, lat_ep0(j1:j1+nj-1), lon_ep0);
[tp_mm, ~, ~]           = to_latlon_sorted(tp_mm, lat_ep0(j1:j1+nj-1), lon_ep0);

ielon = lon_ep >= lon_range(1) & lon_ep <= lon_range(2);
lon_eps = lon_ep(ielon);  lat_eps = lat_ep;
PE_mm = PE_mm(:,ielon,:);  cp_mm = cp_mm(:,ielon,:);  tp_mm = tp_mm(:,ielon,:);

chk = tp_mm(lat_eps >= -10 & lat_eps <= 10, lon_eps >= 50 & lon_eps <= 100, :);
fprintf('\n[单位检查] 10S–10N, 50–100E SON 平均总降水 = %.2f mm/day\n', mean(chk(:),'omitnan'));

[trend_cp, pval_cp] = trend_decade(son_mean(cp_mm, numel(season_month)), yrs_tr);
[trend_pe, pval_pe] = trend_decade(son_mean(PE_mm, numel(season_month)), yrs_tr);
clear PE_mm cp_mm tp_mm

%% ========================================================================
% 6. 画图准备
% ========================================================================
if exist('BlueWhiteOrangeRed.mat','file')
    C = load('BlueWhiteOrangeRed.mat'); fn = fieldnames(C); cmap = C.(fn{1});
else
    warning('BlueWhiteOrangeRed.mat not found. Using turbo.'); cmap = turbo(256);
end
% (a)(b)(c)(d)(g)(h)：Precip4_diff_19lev（文件名 P 大写，变量名大小写以文件内为准）
C2 = load('Precip4_diff_19lev.mat'); fn2 = fieldnames(C2); cmap_p = C2.(fn2{1});
% (e)(f) δD 气候态与动力项用 BlueWhiteOrangeRed

% 海岸线
topo_ext     = [topo topo(:,1)];
topo_lon_ext = [topo_lon 360.5];
tl = topo_lon_ext >= lon_range(1)-5 & topo_lon_ext <= lon_range(2)+5;
tb = topo_lat     >= lat_range(1)-5 & topo_lat     <= lat_range(2)+5;
coast.lon = topo_lon_ext(tl); coast.lat = topo_lat(tb); coast.elev = topo_ext(tb,tl);

% AIRS 子区域
ialon = lon_a >= lon_range(1) & lon_a <= lon_range(2);
ialat = lat_a >= lat_range(1) & lat_a <= lat_range(2);
lon_as = lon_a(ialon);  lat_as = lat_a(ialat);
dD_clim  = dD_clim(ialat,ialon);
dyn_term = dyn_term(ialat,ialon);  pval_dyn = pval_dyn(ialat,ialon);

% 陆地全部屏蔽（显示为 NaN 的浅灰色）
%   ≥1° 格点：格box内陆地占比 > 50% 视为陆地（ERA5 1°、AIRS 5°）
%   <1° 格点：取最近的 1° 地形格点（ERA5 0.25° 的 CP、P−E）
lm_ep = land_mask(lon_eps, lat_eps, topo, topo_lon, topo_lat);
lm_q  = land_mask(lon_qs,  lat_qs,  topo, topo_lon, topo_lat);
lm_a  = land_mask(lon_as,  lat_as,  topo, topo_lon, topo_lat);
trend_cp(lm_ep) = NaN;  trend_pe(lm_ep) = NaN;
MFC_clim(lm_q)  = NaN;  MFC_tr(lm_q) = NaN;  MFC_p(lm_q) = NaN;
QU_clim(lm_q)   = NaN;  QV_clim(lm_q) = NaN;
QU_tr(lm_q)     = NaN;  QV_tr(lm_q)   = NaN;
q_clim(lm_q)    = NaN;  trend_q(lm_q) = NaN;  pval_q(lm_q) = NaN;
dD_clim(lm_a)   = NaN;  dyn_term(lm_a) = NaN; pval_dyn(lm_a) = NaN;

% (a)(b) 水汽通量矢量：按间隔抽样
vstep = max(1, round(vec_spacing_m / median(diff(lon_qs))));
iv = 1:vstep:numel(lon_qs);  jv = 1:vstep:numel(lat_qs);
[Xm, Ym] = meshgrid(lon_qs(iv), lat_qs(jv));
Usig = (QU_p < ALPHA) | (QV_p < ALPHA);                    % qu 或 qv 趋势显著
vlab = @(r) sprintf('%g', r);                              % 参考箭头数值（单位写在图注里）
vec_mfc_clim = struct('X',Xm,'Y',Ym,'U',QU_clim(jv,iv),'V',QV_clim(jv,iv), ...
                      'label',vlab,'ref',[],'sig',[]);
vec_mfc_tr   = struct('X',Xm,'Y',Ym,'U',QU_tr(jv,iv),'V',QV_tr(jv,iv), ...
                      'label',vlab,'ref',[],'sig',Usig(jv,iv));

% 色标范围
c_mfc = 0.03;                                              % (a)(b) 同 Fig. 2：±0.03
cr_dD = [floor(prctile(dD_clim(:),2)/10)*10, ceil(prctile(dD_clim(:),98)/10)*10]; % (e) ‰
c_q   = sym_lim(trend_q);
c_dyn = 1.5;                                               % (f) 固定色标 ±1.5 ‰ day^-1 decade^-1
c_cp  = sym_lim(trend_cp);
c_pe  = sym_lim(trend_pe);
fprintf('色标: δD 气候态 [%g %g]，q 趋势 ±%g，动力项 ±%g\n', cr_dD, c_q, c_dyn);

% 显著性点
sig_mfc = sig_points(MFC_p,    MFC_tr,   lon_qs,  lat_qs,  sig_deg_mfc, ALPHA);
sig_q   = sig_points(pval_q,   trend_q,  lon_qs,  lat_qs,  sig_deg_q,   ALPHA);
sig_dyn = sig_points(pval_dyn, dyn_term, lon_as,  lat_as,  sig_deg_adv, ALPHA);
sig_cp  = sig_points(pval_cp,  trend_cp, lon_eps, lat_eps, sig_deg_ep,  ALPHA);
sig_pe  = sig_points(pval_pe,  trend_pe, lon_eps, lat_eps, sig_deg_ep,  ALPHA);

% (c) 气候态比湿用 Precip4_diff_19lev 的湿端一半，去掉第一行白色
cmap_q = cmap_p(ceil(size(cmap_p,1)/2)+1:end, :);
nq     = size(cmap_q, 1);                                  % 颜色数（19 级色标 → 9 色）
% 色标范围 = 色块数 × 整数间隔，使刻度正好落在颜色变化处
q_lo = prctile(q_clim(:), 2);  q_hi = prctile(q_clim(:), 98);
for q_step = [0.5 1 2 2.5 5]                              % 取能覆盖数据范围的最小整齐间隔
    lo_try = floor(q_lo/q_step)*q_step;
    if lo_try + nq*q_step >= q_hi; break; end
end
q_lo = lo_try;
cr_q = [q_lo, q_lo + nq*q_step];                           % g/kg
tk_q = cr_q(1):q_step:cr_q(2);                             % nq+1 个刻度 = 每个颜色边界
fprintf('q 气候态色标: %g–%g g/kg，间隔 %g\n', cr_q, q_step);
%% ========================================================================
% 7. 版面：宽 18 cm，2 列 4 行，每个面板右侧跟 colorbar，保持长宽比
%    横向: 左边距(纬度标签) + [面板 + 色标区] + 列间距 + [面板 + 色标区] + 右边距
% ========================================================================
nrow   = 4;
fig_w  = 18;
ml     = 1.3;                                       % 左列纬度刻度标签
cb_sp  = 2.15;                                      % 色标间隔+色标+刻度+单位
col_g  = 0.5;                                       % 列间距（右列不标纬度）
mr     = 0.2;
pw     = (fig_w - ml - 2*cb_sp - col_g - mr) / 2;   % ≈ 6.6 cm
ph     = pw * diff(lat_range) / diff(lon_range);    % ≈ 3.3 cm
gap_t  = 0.65; gap_b = 0.2;  mt = 0.2;  mb = 1.1;
fig_h  = mt + nrow*(gap_t + ph) + (nrow-1)*gap_b + mb;   % ≈ 17.7 cm
xcol   = [ml, ml + pw + cb_sp + col_g];
yrow   = mb + (nrow-1:-1:0) * (ph + gap_b + gap_t);

fig = figure('Color','w','Units','centimeters','Position',[2 1 fig_w fig_h], ...
    'Renderer','painters','PaperUnits','centimeters', ...
    'PaperSize',[fig_w fig_h],'PaperPosition',[0 0 fig_w fig_h]);
fprintf('画布 = %.1f x %.1f cm，面板 = %.2f x %.2f cm\n', fig_w, fig_h, pw, ph);

ttr  = sprintf('%s %d–%d', season_name, trend_years(1), trend_years(2));
u_tr = 'mm day^{-1} decade^{-1}';
u_mf = '10^{-6} kg kg^{-1} s^{-1}';

% 列: lon, lat, Z, 色标范围 [lo hi], colormap, 标题, 色标单位, 显著性点, 矢量, 画法, 刻度（[] = 默认 5 个）
P = {
 lon_qs,  lat_qs,  MFC_clim, [-c_mfc c_mfc], cmap_p, ['(a) 850-hPa MFC climatology (' ttr ')'], u_mf,                       [],      vec_mfc_clim, 'image', []
 lon_qs,  lat_qs,  MFC_tr,   [-c_mfc c_mfc], cmap_p, '(b) 850-hPa MFC trend',                    [u_mf ' decade^{-1}'],      sig_mfc, vec_mfc_tr,   'image', []
 lon_qs,  lat_qs,  q_clim,   cr_q,           cmap_q, '(c) q climatology, 850 hPa',               'g kg^{-1}',                [],      [],           'image', tk_q
 lon_qs,  lat_qs,  trend_q,  [-c_q c_q],     cmap_p, '(d) q trend, 850 hPa',                     'g kg^{-1} decade^{-1}',    sig_q,   [],           'image', []
 lon_as,  lat_as,  dD_clim,  cr_dD,          cmap,   '(e) \deltaD climatology, ~850 hPa',        '‰',                        [],      [],           'image', []
 lon_as,  lat_as,  dyn_term, [-c_dyn c_dyn], cmap,   '(f) -\bfV\rm''\cdot\nabla\deltaD_{clim}, ~850 hPa', '‰ day^{-1} decade^{-1}', sig_dyn, [],   'image', []
 lon_eps, lat_eps, trend_cp, [-c_cp c_cp],   cmap_p, '(g) Convective precipitation trend',       u_tr,                       sig_cp,  [],           'image', []
 lon_eps, lat_eps, trend_pe, [-c_pe c_pe],   cmap_p, '(h) P−E trend',                            u_tr,                       sig_pe,  [],           'image', []
};

for k = 1:size(P,1)
    ir = ceil(k/2);  ic = 2 - mod(k,2);            % (a)(b) 第 1 行，(c)(d) 第 2 行 …
    draw_panel([xcol(ic) yrow(ir) pw ph], P{k,1}, P{k,2}, P{k,3}, P{k,4}, P{k,5}, ...
        P{k,6}, P{k,7}, ir == nrow, ic == 1, coast, P{k,8}, P{k,9}, lon_range, lat_range, P{k,10}, P{k,11});
end

%% ========================================================================
% 8. 保存
% ========================================================================
figfile = fullfile(out_dir, 'Fig_8panel_MFC_q_dD_CP_PE.tif');
exportgraphics(fig, figfile, 'Resolution',600, 'BackgroundColor','white');
fprintf('Saved: %s\n', figfile);

pngfile = fullfile(out_dir, 'Fig_8panel_MFC_q_dD_CP_PE.png');
exportgraphics(fig, pngfile, 'Resolution',600, 'BackgroundColor','white');
fprintf('Saved: %s\n', pngfile);
%% ========================================================================
% 子函数
% ========================================================================
function D = airs_layer(dD, plev, prange, it)
il = plev >= prange(1) & plev <= prange(2);
if ~any(il); error('AIRS 在 %d–%d hPa 没有层', prange(1), prange(2)); end
fprintf('AIRS %d–%d hPa 使用层: %s\n', prange(1), prange(2), mat2str(round(plev(il),1)));
[nl, nn, ~, nt] = size(dD);
D = reshape(mean(dD(:,:,il,:), 3, 'omitnan'), nl, nn, nt);
D = D(:,:,it);
end

function [F, lat_s, lon_s] = read_era_layer(file, var, lev, prange, keys, lon0, lat0)
% 读取 ERA5 某层范围内各等压面的平均，只保留 keys 对应月份；
% 返回 lat(南→北) x lon(0→360) x time
il = find(lev >= prange(1) & lev <= prange(2));
if isempty(il)
    error('%s: 没有 %d–%d hPa 的层。可用层: %s', file, prange(1), prange(2), mat2str(lev(:)'));
end
tt = double(ncread(file,'valid_time'));
dd = datetime(1970,1,1) + seconds(tt);
kk = year(dd(:))*100 + month(dd(:));
[ok, it] = ismember(keys, kk);
if any(~ok); error('%s: 缺少 %d 个目标月份', file, sum(~ok)); end
nlon = numel(lon0); nlat = numel(lat0); nT = numel(it);
F = zeros(nlon, nlat, nT);
for k = 1:numel(il)
    tmp = double(ncread(file, var, [1 1 il(k) 1], [Inf Inf 1 Inf]));
    F = F + reshape(tmp(:,:,1,it), nlon, nlat, nT);
end
F = F / numel(il);
fprintf('ERA5 %s %d–%d hPa 使用层: %s\n', var, prange(1), prange(2), mat2str(lev(il)'));
[F, lat_s, lon_s] = to_latlon_sorted(F, lat0, lon0);
end

function [F, lat_s, lon_s] = to_latlon_sorted(F, lat0, lon0)
% lon x lat x t → lat(南→北) x lon(0→360) x t
F = permute(F, [2 1 3]);
[lat_s, j] = sort(lat0(:));      F = F(j,:,:);
lon_s = mod(lon0(:), 360);
[lon_s, i] = sort(lon_s);        F = F(:,i,:);
end

function W = box_weights(xs, xt, d, periodic)
xs = xs(:)'; xt = xt(:);
dx = xs - xt;
if periodic; dx = mod(dx + 180, 360) - 180; end
tol = 1e-6;
W = double(abs(dx) < d/2 - tol) + 0.5 * double(abs(abs(dx) - d/2) <= tol);
s = sum(W, 2);
W = W ./ s;
W(s == 0, :) = NaN;
end

function G = regrid(F, Wlat, Wlon)
G = NaN(size(Wlat,1), size(Wlon,1), size(F,3));
for t = 1:size(F,3)
    G(:,:,t) = Wlat * F(:,:,t) * Wlon';
end
end

function [slope, pval] = adv_trend(D, U, V, lat, lon, nmon, years)
% −V·∇δD 逐月计算（‰ day^-1）→ 逐年季节平均 → 趋势（‰ day^-1 decade^-1）
% 某一年三个月全部缺测时该年为 NaN，trend_decade 会把该格点整体设为 NaN
Re  = 6.371e6;
phi = deg2rad(lat(:));
lam = deg2rad(lon(:)');
cphi = cos(phi); cphi(abs(lat) > 89.9) = NaN;
dl = median(diff(lon));
periodic = abs(mod(lon(end) - lon(1) + dl, 360)) < 1e-3;
gx = zonal_deriv(D, lam, periodic) ./ (Re * cphi);
gy = merid_deriv(D, phi) ./ Re;
ADV = -(U.*gx + V.*gy) * 86400;
[slope, pval] = trend_decade(son_mean(ADV, nmon), years);
end

function dF = zonal_deriv(F, lam, periodic)
n  = size(F, 2);
dF = NaN(size(F));
if periodic
    ip = [2:n 1]; im = [n 1:n-1];
    dl = mod(lam(ip) - lam(im), 2*pi);
    dF = (F(:,ip,:) - F(:,im,:)) ./ dl;
else
    dl = lam(3:n) - lam(1:n-2);
    dF(:,2:n-1,:) = (F(:,3:n,:) - F(:,1:n-2,:)) ./ dl;
end
end

function dF = merid_deriv(F, phi)
n  = size(F, 1);
dF = NaN(size(F));
dp = phi(3:n) - phi(1:n-2);
dF(2:n-1,:,:) = (F(3:n,:,:) - F(1:n-2,:,:)) ./ dp;
end

function M = clim_mean(X, nT, min_frac)
M = mean(X, 3, 'omitnan');
M(sum(~isnan(X), 3) < min_frac * nT) = NaN;
end

function Y = son_mean(F, nmon)
% 每年 nmon 个月相邻排列 → 逐年季节平均
[d1, d2, nt] = size(F);
Y = reshape(mean(reshape(F, d1, d2, nmon, nt/nmon), 3, 'omitnan'), d1, d2, nt/nmon);
end

function [slope, pval] = trend_decade(Y, years)
% 逐格点 OLS 趋势 (per decade) + 双侧 t 检验；有缺测的格点为 NaN
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
b(~ok) = NaN; p(~ok) = NaN;
slope = reshape(b * 10, d1, d2);
pval  = reshape(p, d1, d2);
end

function S = sig_points(pval, Z, lon, lat, sp_deg, alpha)
step = max(1, round(sp_deg / median(diff(lon))));
r = 1:step:numel(lat); c = 1:step:numel(lon);
[X, Y] = meshgrid(lon(c), lat(r));
s = pval(r,c) < alpha & ~isnan(Z(r,c));
S.X = X(s); S.Y = Y(s);
end

function c = sym_lim(Z)
z = abs(Z(:)); z = z(isfinite(z));
c = nice_ceil(prctile(z, 98));
end

function y = nice_ceil(x)
if ~isfinite(x) || x <= 0; y = 1; return; end
e = 10^floor(log10(x));
m = x / e;
if m <= 1; y = e; elseif m <= 2; y = 2*e; elseif m <= 5; y = 5*e; else; y = 10*e; end
end

function lm = land_mask(lon, lat, topo, topo_lon, topo_lat)
d = median(diff(lon));
if d >= 1
    frac = box_weights(topo_lat, lat, d, false) * double(topo > 0) * ...
           box_weights(topo_lon, lon, d, true)';
    lm = frac > 0.5;
else
    [X, Y] = meshgrid(mod(lon,360), lat);
    tl = [topo_lon(end)-360, topo_lon, topo_lon(1)+360];
    tt = double([topo(:,end), topo, topo(:,1)]);
    lm = interp2(tl, topo_lat, tt, X, Y, 'nearest', 0) > 0;
end
end

function A = dyn_adv(Dc, u, v, lat, lon)
% 动力项 −V′·∇δD_clim（‰ day^-1 decade^-1）
% Dc: 气候态 δD (lat x lon)；u, v: 风场趋势 (m s^-1 decade^-1)
Re   = 6.371e6;
phi  = deg2rad(lat(:));
lam  = deg2rad(lon(:)');
cphi = cos(phi); cphi(abs(lat(:)) > 89.9) = NaN;
dl   = median(diff(lon));
periodic = abs(mod(lon(end) - lon(1) + dl, 360)) < 1e-3;
gx = zonal_deriv(Dc, lam, periodic) ./ (Re * cphi);
gy = merid_deriv(Dc, phi) ./ Re;
A  = -(u.*gx + v.*gy) * 86400;
end

function draw_panel(pos, lon, lat, Z, cr, cmap, ttl, cblabel, show_x, show_y, coast, sig, vec, lon_range, lat_range, mode, ticks)
% cr = 色标范围 [lo hi]（趋势图为对称范围，气候态可以不对称）
fs = 10;                                            % 全图统一字号
ax = axes('Units','centimeters','Position',pos);
ax.Toolbar.Visible = 'off';                          % 屏幕上不显示坐标区工具栏
hold(ax,'on');
if strcmp(mode, 'contour')
    % 超出色标的值截断到边界（逐个替换，NaN 保持不变）
    Zc = Z;
    Zc(Zc > cr(2)) = cr(2);
    Zc(Zc < cr(1)) = cr(1);
    lonp = [lon_range(1); lon(:); lon_range(2)];
    latp = [lat_range(1); lat(:); lat_range(2)];
    Zp = Zc([1 1:end end], [1 1:end end]);
    contourf(ax, lonp, latp, Zp, linspace(cr(1), cr(2), 61), 'LineColor','none');
else
    imagesc(ax, lon, lat, Z, 'AlphaData', ~isnan(Z));
end
set(ax,'YDir','normal','Color',[0.85 0.85 0.85]);   % 陆地 NaN 显示为浅灰
colormap(ax, cmap);
clim(ax, cr);
if ~isempty(vec)
    % 手动比例：参考值取矢量长度的 90% 分位（取整到 1/2/5×10^n），
    % 参考箭头画成 ref_len 度长，所有箭头按同一比例缩放
    ref_len = 5;                                    % 参考箭头长度（经度°）
    mag = hypot(vec.U(:), vec.V(:));  mag = mag(isfinite(mag));
    if isfield(vec,'ref') && ~isempty(vec.ref)
        ref = vec.ref;                              % 手动指定
    else
        ref = nice_ceil(prctile(mag, 90));          % 自动
    end
    sc = ref_len / ref;
    if isfield(vec,'sig') && ~isempty(vec.sig)
        S  = vec.sig;
        Un = vec.U;  Vn = vec.V;  Un(S)  = NaN;  Vn(S)  = NaN;    % 不显著：黑
        Us = vec.U;  Vs = vec.V;  Us(~S) = NaN;  Vs(~S) = NaN;    % 显著：紫
        quiver(ax, vec.X, vec.Y, Un*sc, Vn*sc, 0, 'k', 'LineWidth',1.1, 'MaxHeadSize',0.6);
        quiver(ax, vec.X, vec.Y, Us*sc, Vs*sc, 0, 'Color',[0.75 0 0.95], ...
               'LineWidth',1.5, 'MaxHeadSize',0.6);
    else
        quiver(ax, vec.X, vec.Y, vec.U*sc, vec.V*sc, 0, 'k', 'LineWidth',1.1, 'MaxHeadSize',0.6);
    end
end
contour(ax, coast.lon, coast.lat, coast.elev, [0 0], 'k-', 'LineWidth',0.5);
if ~isempty(sig) && ~isempty(sig.X)
    plot(ax, sig.X, sig.Y, 'o', 'MarkerSize',2, ...
        'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none');
end
xlim(ax, lon_range); ylim(ax, lat_range);
daspect(ax, [1 1 1]);
set(ax, 'XTick',0:30:120, 'YTick',-20:20:40, 'XTickLabelRotation',0, ...
    'YTickLabel',{'20°S','0°','20°N','40°N'}, ...
    'FontSize',fs, 'TickDir','out', 'Box','off', 'Layer','top');
if show_x
    set(ax,'XTickLabel',{'0°','30°E','60°E','90°E','120°E'});
else
    set(ax,'XTickLabel',{});
end
if ~show_y
    set(ax,'YTickLabel',{});
end
if ~isempty(vec)
    % 参考箭头 + 数值：面板内右上角，白底黑框
    y0 = lat_range(2) - 5;
    ht = text(ax, lon_range(2) - 2, y0, vec.label(ref), 'FontSize',fs, ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', 'Interpreter','tex');
    ext = ht.Extent;  delete(ht);                   % 先量文字宽度（数据单位）
    x0 = ext(1) - ref_len - 2;
    rectangle(ax, 'Position',[x0-2, y0-ext(4)/2-1.5, lon_range(2)-(x0-2), lat_range(2)-(y0-ext(4)/2-1.5)], ...
        'FaceColor','w', 'EdgeColor','k', 'LineWidth',0.5);
    quiver(ax, x0, y0, ref_len, 0, 0, 'k', 'LineWidth',1.1, 'MaxHeadSize',0.8);
    text(ax, lon_range(2) - 2, y0, vec.label(ref), 'FontSize',fs, ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', 'Interpreter','tex');
    fprintf('%s：参考箭头 = %g\n', ttl, ref);
end
xl = ax.XLim; yl = ax.YLim;
plot(ax,[xl(2) xl(2)],yl,'k-','LineWidth',0.5);
plot(ax,xl,[yl(2) yl(2)],'k-','LineWidth',0.5);
title(ax, ttl, 'FontSize',fs, 'FontWeight','bold', 'HorizontalAlignment','left', ...
    'Units','normalized', 'Position',[0 1.03 0], 'Interpreter','tex');

if nargin < 17 || isempty(ticks)
    ticks = linspace(cr(1), cr(2), 5);
end
labs = arrayfun(@(v) sprintf('%g', round(v,6)), ticks, 'UniformOutput', false);
if numel(ticks) > 6                                  % 刻度多时隔一个标数字，刻度线仍在每个颜色边界
    labs(2:2:end) = {''};
end
cb = colorbar(ax);
cb.Units = 'centimeters';
cb.Position = [pos(1)+pos(3)+0.2, pos(2), 0.25, pos(4)];
cb.Ticks = ticks;
cb.TickLabels = labs;
cb.FontSize = fs;
cb.TickDirection = 'out';
cb.TickLength = 0.04;
cb.Label.String = cblabel;
cb.Label.FontSize = fs;
ax.Position = pos;
end
