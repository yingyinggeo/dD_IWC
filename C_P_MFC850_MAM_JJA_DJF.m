%% ========================================================================
%  850-hPa 水汽通量辐合 MFC = -∇·(qV)：春 / 夏 / 冬，3 行 2 列
%   左列：气候态 + (qu, qv) 矢量
%   右列：趋势 (per decade) + (qu, qv) 趋势矢量（qu 或 qv 显著 → 紫色）
%   (a) MAM climatology   (b) MAM trend
%   (c) JJA climatology   (d) JJA trend
%   (e) DJF climatology   (f) DJF trend
%   算法、配色、布局沿用 Fig2 的 (e)(f)；30–120°E, 40°S–40°N；陆地浅灰
%   DJF 年份按 1、2 月所在年份记：DJF 2003 = Dec 2002 + Jan 2003 + Feb 2003
%% ========================================================================

%% ===== 基本设置 =====
yr_start = 2003;  yr_end = 2019;
reg_years = (yr_start:yr_end)';
nyr = numel(reg_years);

lon_range = [30 120];
lat_range = [-40 40];

seasons = struct( ...
    'name',   {'MAM', 'JJA', 'DJF'}, ...
    'months', {[3 4 5], [6 7 8], [12 1 2]});
nsea = numel(seasons);

% 色标（每行左右共用一个色标）；气候态中 |MFC| > blank_clim 的格点画成白色
% 若某季节气候态大片变白，可单独调大该季节的 blank_clim 或 clim_lim
clim_lim   = [0.03 0.03 0.03];    % MAM, JJA, DJF 色标 ±（10^-6 kg kg^-1 s^-1）
blank_clim = [0.04 0.04 0.04];    % 气候态白色阈值；设为 Inf 则不变白

%% ===== 读取 850 hPa u, v, q =====
fprintf('读取 850-hPa u, v, q ...\n');
era_dir = 'C:\Users\yingy\Downloads\Rice_mac\ERA5';
ufile_m = fullfile(era_dir, 'era5_uwind_90N90S_1deg_2000_2024.nc');
vfile_m = fullfile(era_dir, 'era5_vwind_90N90S_1deg_2000_2024.nc');
qfile_m = fullfile(era_dir, 'era5_qs_90N90S_1deg_2000_2024.nc');

lon_m = double(ncread(ufile_m,'longitude'));
lat_m = double(ncread(ufile_m,'latitude'));
tm    = double(ncread(ufile_m,'valid_time'));
dm    = datetime(1970,1,1) + seconds(tm);
yr_m  = year(dm);  mo_m = month(dm);

% 只读需要的月份：Dec(yr_start-1) 到 Nov(yr_end)
t_first = datetime(yr_start-1,12,1);
t_last  = datetime(yr_end,11,30);
it_m = find(dm >= t_first & dm <= t_last);
yr_it = yr_m(it_m);  mo_it = mo_m(it_m);
fprintf('读取月份数 = %d\n', numel(it_m));

u850 = read850(ufile_m, 'u', it_m);
v850 = read850(vfile_m, 'v', it_m);
q850 = read850(qfile_m, 'q', it_m);
QU = q850 .* u850;   QV = q850 .* v850;                  % kg kg^-1 m s^-1
clear u850 v850 q850

%% ===== MFC = -∇·(qV)（球坐标中央差分）=====
Re   = 6.371e6;
phiM = reshape(deg2rad(lat_m), 1, []);
dlon_m = abs(median(diff(lon_m)));
dQU  = (circshift(QU,-1,1) - circshift(QU,1,1)) / deg2rad(2*dlon_m);
QVc  = QV .* cos(phiM);
dQV  = NaN(size(QV));
dQV(:,2:end-1,:) = (QVc(:,3:end,:) - QVc(:,1:end-2,:)) ./ (phiM(3:end) - phiM(1:end-2));
MFC  = -(dQU + dQV) ./ (Re * cos(phiM)) * 1e6;          % 10^-6 kg kg^-1 s^-1
clear dQU dQV QVc

%% ===== 网格整理、裁剪 =====
[lonM, ix] = sort(mod(lon_m,360));
[latM, iy] = sort(lat_m);
mlon = lonM >= lon_range(1) & lonM <= lon_range(2);
mlat = latM >= lat_range(1) & latM <= lat_range(2);
lon_mfc = lonM(mlon);  lat_mfc = latM(mlat);
ixc = ix(mlon);  iyc = iy(mlat);
% lon x lat x t → lat x lon x t，并裁剪
MFC = permute(MFC(ixc, iyc, :), [2 1 3]);
QU  = permute(QU (ixc, iyc, :), [2 1 3]) * 1e3;          % g kg^-1 m s^-1
QV  = permute(QV (ixc, iyc, :), [2 1 3]) * 1e3;

%% ===== 陆地 mask =====
load topo;
topo_lon     = 0.5:359.5;
topo_lat     = -89.5:89.5;
[TLON, TLAT] = meshgrid(topo_lon, topo_lat);
[OLONm, OLATm] = meshgrid(lon_mfc, lat_mfc);
land_mfc = interp2(TLON, TLAT, topo, OLONm, OLATm, 'nearest', 0) > 0;

buf = 5;
topo_lon_ext = [topo_lon, 360.5];
topo_ext     = [topo, topo(:,1)];
tlon_idx    = topo_lon_ext >= lon_range(1)-buf & topo_lon_ext <= lon_range(2)+buf;
tlat_idx    = topo_lat     >= lat_range(1)-buf & topo_lat     <= lat_range(2)+buf;
topo_lon_v  = topo_lon_ext(tlon_idx);
topo_lat_v  = topo_lat(tlat_idx);
topo_elev_v = topo_ext(tlat_idx, tlon_idx);

%% ===== 各季节：气候态 + 逐年季节平均 → 趋势 =====
R = struct([]);
for s = 1:nsea
    mos = seasons(s).months;
    sy  = yr_it;                                   % 季节年份
    if strcmp(seasons(s).name, 'DJF')
        sy(mo_it == 12) = sy(mo_it == 12) + 1;     % 12 月归入下一年的 DJF
    end
    use = ismember(mo_it, mos) & ismember(sy, reg_years);

    % 逐年季节平均（缺月份的年份记为 NaN）
    MFC_y = NaN([size(MFC,1) size(MFC,2) nyr]);  QU_y = MFC_y;  QV_y = MFC_y;
    for k = 1:nyr
        kk = use & sy == reg_years(k);
        if sum(kk) ~= numel(mos)
            warning('%s %d 只有 %d 个月', seasons(s).name, reg_years(k), sum(kk));
            continue
        end
        MFC_y(:,:,k) = mean(MFC(:,:,kk), 3);
        QU_y(:,:,k)  = mean(QU(:,:,kk),  3);
        QV_y(:,:,k)  = mean(QV(:,:,kk),  3);
    end

    % 气候态 = 逐年季节平均的平均
    R(s).clim    = mean(MFC_y, 3, 'omitnan');
    R(s).qu_clim = mean(QU_y,  3, 'omitnan');
    R(s).qv_clim = mean(QV_y,  3, 'omitnan');

    [R(s).tr,    R(s).p   ] = trend_decade(MFC_y, reg_years);
    [R(s).qu_tr, R(s).qu_p] = trend_decade(QU_y,  reg_years);
    [R(s).qv_tr, R(s).qv_p] = trend_decade(QV_y,  reg_years);

    % 陆地 → NaN
    f = {'clim','qu_clim','qv_clim','tr','p','qu_tr','qu_p','qv_tr','qv_p'};
    for j = 1:numel(f)
        A = R(s).(f{j});  A(land_mfc) = NaN;  R(s).(f{j}) = A;
    end
    fprintf('%s 完成\n', seasons(s).name);
end
clear MFC QU QV

%% ===== 画图：3 行 2 列 =====
load 'BlueWhiteOrangeRed.mat';
silver_rgb = [0.85 0.85 0.85];

cm2px    = 96/2.54;
fig_w_cm = 18;
pw_cm = 0.32*18;   ph_cm = 4.23;   row_dy = 5.17;
y3_cm = 1.0;
y_cm  = [y3_cm + 2*row_dy, y3_cm + row_dy, y3_cm];     % 第 1、2、3 行底边
fig_h_cm = y_cm(1) + ph_cm + 0.7;

fig = figure('Color','w', ...
             'Position',[50 50 round(fig_w_cm*cm2px) round(fig_h_cm*cm2px)], ...
             'Renderer','opengl', 'GraphicsSmoothing','off', ...
             'PaperUnits','centimeters', 'PaperSize',[fig_w_cm, fig_h_cm], ...
             'PaperPositionMode','auto');

pw      = pw_cm / fig_w_cm;
ph      = ph_cm / fig_h_cm;
left_x  = 0.07;
right_x = 0.46;
cb_x    = 0.82;   cb_w = 0.02;

vec_skip = 6;   sig_skip = 2;   ref_len = 5;
iq = 1:vec_skip:numel(lon_mfc);  jq = 1:vec_skip:numel(lat_mfc);
[LONq, LATq] = meshgrid(lon_mfc(iq), lat_mfc(jq));
is = 1:sig_skip:numel(lon_mfc);  js = 1:sig_skip:numel(lat_mfc);
[LONs, LATs] = meshgrid(lon_mfc(is), lat_mfc(js));
land_img_m = repmat(reshape(silver_rgb,1,1,3), [size(land_mfc) 1]);

cbar_label = {'Climatology (10^{-6} kg kg^{-1} s^{-1})', ...
              'Trend (10^{-6} kg kg^{-1} s^{-1} decade^{-1})'};
letters = 'abcdef';

for s = 1:nsea
    row_y = y_cm(s) / fig_h_cm;
    cl    = clim_lim(s);
    for k = 1:2
        ip = 2*(s-1) + k;                          % 面板序号 a–f
        if k == 1
            xpos = left_x;   dat = R(s).clim;
            U = R(s).qu_clim;  V = R(s).qv_clim;  Usig = [];
            ttl = sprintf('(%c) %s 850-hPa MFC Climatology (%d–%d)', ...
                          letters(ip), seasons(s).name, yr_start, yr_end);
            ref_unit = 'g kg^{-1} m s^{-1}';
        else
            xpos = right_x;  dat = R(s).tr;
            U = R(s).qu_tr;    V = R(s).qv_tr;
            Usig = (R(s).qu_p < 0.05) | (R(s).qv_p < 0.05);
            ttl = sprintf('(%c) %s 850-hPa MFC Trend (%d–%d)', ...
                          letters(ip), seasons(s).name, yr_start, yr_end);
            ref_unit = 'g kg^{-1} m s^{-1} decade^{-1}';
        end
        pos = [xpos, row_y, pw, ph];
        ax  = axes('Position', pos);  hold(ax,'on');

        bg = land_img_m;
        if k == 1
            out = abs(dat) > blank_clim(s) & ~land_mfc;   % 超出阈值 → 白色
            dat(out) = NaN;
            for c3 = 1:3
                tmp = bg(:,:,c3);  tmp(out) = 1;  bg(:,:,c3) = tmp;
            end
        end
        imagesc(ax, lon_mfc, lat_mfc, bg);
        imagesc(ax, lon_mfc, lat_mfc, dat, 'AlphaData', ~isnan(dat));
        set(ax, 'YDir','normal', 'Color',silver_rgb);
        colormap(ax, BlueWhiteOrangeRed);
        caxis(ax, [-cl cl]);

        % 显著性打点（仅趋势）
        if k == 2
            sg = R(s).p(js,is) < 0.05 & ~isnan(dat(js,is));
            plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize',1.5, ...
                 'MarkerFaceColor',[0.1 0.7 0.1], 'MarkerEdgeColor','none');
        end

        % 矢量：参考箭头 = 90% 分位（取整），统一手动比例
        Uq = U(jq,iq);  Vq = V(jq,iq);
        mag = hypot(Uq(:), Vq(:));  mag = mag(isfinite(mag));
        ref = nice_ceil(prctile(mag, 90));
        sc  = ref_len / ref;
        if isempty(Usig)
            quiver(ax, LONq, LATq, Uq*sc, Vq*sc, 0, 'k', 'LineWidth',0.6, 'MaxHeadSize',0.5);
        else
            S  = Usig(jq,iq);
            Un = Uq;  Vn = Vq;  Un(S) = NaN;  Vn(S) = NaN;
            Us = Uq;  Vs = Vq;  Us(~S) = NaN; Vs(~S) = NaN;
            quiver(ax, LONq, LATq, Un*sc, Vn*sc, 0, 'k', 'LineWidth',0.6, 'MaxHeadSize',0.5);
            quiver(ax, LONq, LATq, Us*sc, Vs*sc, 0, 'Color',[0.75 0 0.95], ...
                   'LineWidth',0.9, 'MaxHeadSize',0.5);
        end

        contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth',0.5);
        xlim(ax, lon_range); ylim(ax, lat_range);

        % x 刻度标签只在最后一行；y 刻度标签只在左列
        if s == nsea
            set(ax, 'XTick',30:30:120, 'XTickLabel',{'30°E','60°E','90°E','120°E'});
        else
            set(ax, 'XTick',30:30:120, 'XTickLabel',{});
        end
        if k == 1
            set(ax, 'YTick',-40:20:40, 'YTickLabel',{'40°S','20°S','0°','20°N','40°N'});
        else
            set(ax, 'YTick',-40:20:40, 'YTickLabel',{});
        end
        set(ax, 'FontSize',10, 'TickDir','out', 'Box','on', 'Layer','top');
        ax.Position = pos;

        draw_ref_arrow(ax, lon_range, lat_range, ref_len, sprintf('%g', ref));

        annotation(fig, 'textbox', [xpos, row_y+ph+0.005, pw+0.1, 0.02], ...
                   'String', ttl, 'FontSize',10, 'FontWeight','bold', ...
                   'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
                   'EdgeColor','none', 'FitBoxToText','off', 'Interpreter','tex');

        % 每行一个竖向色标（放在趋势图右侧），两行标签
        if k == 2
            cb = colorbar(ax, 'eastoutside');
            cb.Position = [cb_x+cb_w+0.005, row_y, 0.018, ph];
            cb.Ticks = cl * [-1 -0.5 0 0.5 1];
            cb.Label.String   = cbar_label;
            cb.FontSize       = 10;
            cb.Label.FontSize = 10;
            cb.TickDirection  = 'out';
        end
        ax.Position = pos;
        fprintf('(%c) 参考箭头 = %g %s\n', letters(ip), ref, ref_unit);
    end
end

set(findall(fig,'Type','axes'), 'TickLength',[0.025 0.035]);

%% ===== 保存 =====
out_dir = 'C:\Users\yingy\Downloads\Rice_mac';
if ~exist(out_dir,'dir'); mkdir(out_dir); end
fname = fullfile(out_dir, 'FigS_MFC850_MAM_JJA_DJF.tif');
exportgraphics(fig, fname, 'Resolution',600, 'BackgroundColor','white');
fprintf('Saved: %s\n', fname);
pngname = fullfile(out_dir, 'FigS_MFC850_MAM_JJA_DJF.png');
exportgraphics(fig, pngname, 'Resolution',600, 'BackgroundColor','white');
fprintf('Saved: %s\n', pngname);

%% ===== 子函数 =====
function F = read850(file, var, it)
% 读 850 hPa 单层，只保留 it 指定的时次，返回 lon x lat x time
lev = double(ncread(file,'pressure_level'));
[~, il] = min(abs(lev - 850));
if abs(lev(il) - 850) > 0.1; error('%s 中没有 850 hPa', file); end
t0 = min(it);  nt = max(it) - t0 + 1;                     % 只读需要的时间段
F = double(ncread(file, var, [1 1 il t0], [Inf Inf 1 nt]));
F = reshape(F(:,:,1,it - t0 + 1), size(F,1), size(F,2), numel(it));
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
quiver(ax, x0, y0, ref_len, 0, 0, 'k', 'LineWidth',0.6, 'MaxHeadSize',0.8);
text(ax, lon_range(2) - 1.5, y0, label, 'FontSize',10, ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle');
end
