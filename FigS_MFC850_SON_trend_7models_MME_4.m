%% ========================================================================
%  SON 850-hPa 水汽通量辐合 MFC = -∇·(qV) 趋势：7 个模式 + MME，4 行 2 列
%   参照 Fig2 (f)：MFC 趋势填色 + (qu, qv) 趋势矢量
%     - MFC 趋势显著 (p<0.05) → 绿点
%     - qu 或 qv 趋势显著     → 紫色箭头；不显著 → 黑色箭头
%     - 陆地（1° topo 最近点 > 0）屏蔽为浅灰
%   数据：C:\Users\yingy\Downloads\Rice_mac\WISOMIP_Other_Variabiles\Total.<model>_Monthly.nc
%         (128 x 64 网格，17 个气压层，单位 mb；变量 u, v, q)
%   计算范围：60°S–60°N, 0–180°E（MFC 先在全球网格上逐月计算，再裁剪）
%   画图范围：20°S–40°N, 30–120°E（纬向比例与原 40°S–40°N 版本相同）
%   时段：优先 2003–2020；若有模式缺 2020 年 SON，则统一用 2003–2019
%   MME：先对各模式逐年 SON 场取集合平均，再算趋势
%   色标：Precip4_diff_19lev 去掉最中间的白色色阶，色标范围 ±18（10^-9 单位），
%         刻度放在色阶分界上，保证刻度与色块对齐
%   WIO/EIO 方框不再画出，但偶极子判断仍按方框区域计算
%% ========================================================================
clear; clc; close all;

%% ===== 基本设置 =====
data_dir = 'C:\Users\yingy\Downloads\Rice_mac\WISOMIP_Other_Variabiles';
models   = {'CAM5', 'CAM6', 'ECHAM', 'GISS', 'GSM', 'LMDZ', 'MIROC'};   % 可换成 NICAM
n_models = numel(models);

p_target = 850;                          % hPa
yr_start = 2003;  yr_try_end = 2020;  yr_fallback_end = 2019;
years_all = (yr_start:yr_try_end)';
nyr_all   = numel(years_all);
son = [9 10 11];

calc_lon = [0 180];   calc_lat = [-60 60];      % 计算范围
plot_lon = [30 120];  plot_lat = [-20 40];      % 画图范围
plot_lat_ref = [-40 40];                        % 原画图纬度范围，用于保持长宽比

% ---- 色标设置（色阶分界与刻度对齐）----
%  数据单位：10^-6 kg kg^-1 s^-1 decade^-1
%  clim_tr ：色标范围 ±0.018，即标签上的 ±18（10^-9 单位）
%  去掉最中间的白色色阶后，色阶宽度 = 2 × clim_tr / 色阶数（自动计算）
%  tick_scale：刻度标签的换算系数。1e3 表示标签以 10^-9 为单位显示
%            （数据本身不变，只是标签更简洁）
%  tick_lbl_step：1 = 每个分界都标数字；2 = 隔一个标一个（对称）
clim_tr       = 0.018;
white_thr     = 0.95;                    % RGB 三个分量都 ≥ 此值视为白色
tick_scale    = 1e3;
tick_lbl_unit = '10^{-9}';
tick_lbl_step = 1;

sig_level = 0.05;
fs        = 10;
dot_sz    = 3;                           % 绿点大小（原 1.5）

% 区域框（Fig. 1c）：{经度范围, 纬度范围, 颜色}
boxes = { ...
    [50 70],  [-15 5], [1.00 0.50 0.00];   % western Indian Ocean (WIO)，橙色
    [70 100], [-20 0], [0.00 0.45 0.90]};  % eastern Indian Ocean (EIO)，蓝色
box_lw = 1.2;
show_boxes = false;                      % 不画 WIO/EIO 方框（偶极子判断仍使用方框区域）
no_dipole_tag_panels = 1:5;              % 面板 (a)–(e) 标题不加 "(dipole)"
show_dipole_tag = true;                  % 满足偶极子条件的面板，标题后加 "dipole"

%% ===== 逐模式读取 → 逐月 MFC, qu, qv → 逐年 SON 平均 =====
MFC_all = [];  QU_all = [];  QV_all = [];
have    = false(n_models, 1);
ok_last = false(n_models, 1);            % 是否有完整的 2020 年 SON

for im = 1:n_models
    f = fullfile(data_dir, sprintf('Total.%s_Monthly.nc', models{im}));
    fprintf('\n===== %s =====\n', models{im});
    if ~exist(f, 'file'); warning('找不到 %s', f); continue; end

    try
        lon = double(ncread(f, 'lon'));
        lat = double(ncread(f, 'lat'));
        p   = double(ncread(f, 'p'));
        [~, il] = min(abs(p - p_target));
        fprintf('使用气压层 p = %g hPa\n', p(il));

        [yr, mo] = model_time(f);
        fprintf('时间：%d-%02d 到 %d-%02d（%d 个月）\n', yr(1), mo(1), yr(end), mo(end), numel(yr));

        u = read_level(f, 'u', il);      % lon x lat x time
        v = read_level(f, 'v', il);
        q = read_level(f, 'q', il);
    catch ME
        warning('%s 读取失败：%s', models{im}, ME.message);
        continue
    end

    QU = q .* u;   QV = q .* v;                          % kg kg^-1 m s^-1
    clear u v q
    MFC = calc_mfc(QU, QV, lon, lat);                    % 10^-6 kg kg^-1 s^-1（全球网格上算）

    % 整理为 lat(南→北) x lon(0→360) x t，并裁剪到计算范围
    [lonS, ix] = sort(mod(lon, 360));
    [latS, iy] = sort(lat);
    mlon = lonS >= calc_lon(1) & lonS <= calc_lon(2);
    mlat = latS >= calc_lat(1) & latS <= calc_lat(2);
    lon_c = lonS(mlon);  lat_c = latS(mlat);
    MFC = permute(MFC(ix(mlon), iy(mlat), :), [2 1 3]);
    QU  = permute(QU (ix(mlon), iy(mlat), :), [2 1 3]) * 1e3;   % g kg^-1 m s^-1
    QV  = permute(QV (ix(mlon), iy(mlat), :), [2 1 3]) * 1e3;

    if isempty(MFC_all)
        lon_g = lon_c;  lat_g = lat_c;
        sz = [numel(lat_g) numel(lon_g) nyr_all n_models];
        MFC_all = NaN(sz);  QU_all = NaN(sz);  QV_all = NaN(sz);
    elseif numel(lon_c) ~= numel(lon_g) || numel(lat_c) ~= numel(lat_g) || ...
           max(abs(lon_c - lon_g)) > 1e-3 || max(abs(lat_c - lat_g)) > 1e-3
        error('%s 的网格与其他模式不一致', models{im});
    end

    % 逐年 SON 平均
    for k = 1:nyr_all
        kk = ismember(mo, son) & yr == years_all(k);
        if sum(kk) ~= 3; continue; end
        A = mean(MFC(:,:,kk), 3);
        if all(isnan(A(:))); continue; end
        MFC_all(:,:,k,im) = A;
        QU_all(:,:,k,im)  = mean(QU(:,:,kk), 3);
        QV_all(:,:,k,im)  = mean(QV(:,:,kk), 3);
    end
    ok_last(im) = any(~isnan(reshape(MFC_all(:,:,end,im), [], 1)));
    have(im) = true;
    fprintf('%s 完成；%d 年 SON %s\n', models{im}, yr_try_end, ...
        ternary(ok_last(im), '可用', '缺失'));
    clear MFC QU QV
end

if ~any(have); error('没有可用的模式数据'); end

%% ===== 确定时段 =====
if all(ok_last(have))
    yr_end = yr_try_end;
else
    yr_end = yr_fallback_end;
    fprintf('\n以下模式缺 %d 年 SON：%s → 统一使用 %d–%d\n', yr_try_end, ...
        strjoin(models(have & ~ok_last), ', '), yr_start, yr_end);
end
reg_years = (yr_start:yr_end)';
iyr = ismember(years_all, reg_years);
MFC_all = MFC_all(:,:,iyr,:);  QU_all = QU_all(:,:,iyr,:);  QV_all = QV_all(:,:,iyr,:);
fprintf('时段：%d–%d（%d 年），模式数 = %d (%s)\n', yr_start, yr_end, numel(reg_years), ...
    sum(have), strjoin(models(have), ', '));

%% ===== 陆地 mask（同 Fig2 f：1° topo 最近点 > 0）=====
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

%% ===== 趋势（计算范围内）：7 个模式 + MME =====
n_panels    = n_models + 1;
panel_names = [models, {'MME'}];
R = struct([]);
for ip = 1:n_panels
    if ip <= n_models
        if ~have(ip); R(ip).ok = false; continue; end
        M = MFC_all(:,:,:,ip);  U = QU_all(:,:,:,ip);  V = QV_all(:,:,:,ip);
    else
        M = mean(MFC_all(:,:,:,have), 4, 'omitnan');
        U = mean(QU_all(:,:,:,have),  4, 'omitnan');
        V = mean(QV_all(:,:,:,have),  4, 'omitnan');
    end
    [R(ip).tr,    R(ip).p   ] = trend_decade(M, reg_years);
    [R(ip).qu_tr, R(ip).qu_p] = trend_decade(U, reg_years);
    [R(ip).qv_tr, R(ip).qv_p] = trend_decade(V, reg_years);
    fl = {'tr','p','qu_tr','qu_p','qv_tr','qv_p'};
    for j = 1:numel(fl)
        A = R(ip).(fl{j});  A(land_g) = NaN;  R(ip).(fl{j}) = A;
    end
    R(ip).ok = true;
end

% 保存计算范围内的结果
save(fullfile('C:\Users\yingy\Downloads\Rice_mac', ...
    sprintf('MFC850_SON_trend_%d_%d_7models_MME.mat', yr_start, yr_end)), ...
    'R', 'panel_names', 'lon_g', 'lat_g', 'reg_years', 'land_g', 'calc_lon', 'calc_lat');

%% ===== 偶极子判断：WIO（橙框）区域平均 > 0 且 EIO（蓝框）区域平均 < 0 =====
%  区域平均只用海洋格点，cos(lat) 加权，对逐年 SON MFC 序列求趋势（per decade）
%  （与对趋势图做区域平均等价），并给出各自及 D = WIO − EIO 的显著性
nyr = numel(reg_years);
Name = {}; WIO_tr = []; WIO_p = []; EIO_tr = []; EIO_p = []; D_tr = []; D_p = [];
WIO_pos = {}; EIO_neg = {}; Dipole = {};
for ip = 1:n_panels
    R(ip).dipole = false;
    if ~R(ip).ok; continue; end
    if ip <= n_models
        M = MFC_all(:,:,:,ip);
    else
        M = mean(MFC_all(:,:,:,have), 4, 'omitnan');
    end
    w = box_mean(M, lon_g, lat_g, boxes{1,1}, boxes{1,2}, land_g);
    e = box_mean(M, lon_g, lat_g, boxes{2,1}, boxes{2,2}, land_g);
    [bw, pw_] = ts_trend(w, reg_years, 10);
    [be, pe ] = ts_trend(e, reg_years, 10);
    [bd, pd ] = ts_trend(w - e, reg_years, 10);
    dip = bw > 0 && be < 0;
    R(ip).dipole = dip;
    R(ip).box = struct('WIO_tr', bw, 'WIO_p', pw_, 'EIO_tr', be, 'EIO_p', pe, 'D_tr', bd, 'D_p', pd);

    Name{end+1,1}    = panel_names{ip};                  %#ok<SAGROW>
    WIO_tr(end+1,1)  = bw;  WIO_p(end+1,1) = pw_;         %#ok<SAGROW>
    EIO_tr(end+1,1)  = be;  EIO_p(end+1,1) = pe;          %#ok<SAGROW>
    D_tr(end+1,1)    = bd;  D_p(end+1,1)   = pd;          %#ok<SAGROW>
    WIO_pos{end+1,1} = ternary(bw > 0, 'Yes', 'No');     %#ok<SAGROW>
    EIO_neg{end+1,1} = ternary(be < 0, 'Yes', 'No');     %#ok<SAGROW>
    Dipole{end+1,1}  = ternary(dip, 'Yes', 'No');        %#ok<SAGROW>
end
T = table(Name, WIO_tr, WIO_p, WIO_pos, EIO_tr, EIO_p, EIO_neg, D_tr, D_p, Dipole);
fprintf('\n===== 偶极子判断：WIO 平均趋势 > 0 且 EIO 平均趋势 < 0 =====\n');
fprintf('SON 850-hPa MFC 趋势，%d–%d，单位 10^-6 kg kg^-1 s^-1 decade^-1（海洋格点、cos(lat) 加权）\n', ...
    yr_start, yr_end);
disp(T);
n_dip = sum(strcmp(T.Dipole(1:end-1), 'Yes'));
fprintf('满足偶极子条件的模式：%d/%d（%s）；MME：%s\n', n_dip, height(T)-1, ...
    strjoin(T.Name(strcmp(T.Dipole(1:end-1), 'Yes')), ', '), T.Dipole{end});
writetable(T, fullfile('C:\Users\yingy\Downloads\Rice_mac', ...
    sprintf('MFC850_SON_dipole_WIO_EIO_%d_%d.csv', yr_start, yr_end)));

%% ===== 裁剪到画图范围 =====
mlon_p = lon_g >= plot_lon(1) & lon_g <= plot_lon(2);
mlat_p = lat_g >= plot_lat(1) & lat_g <= plot_lat(2);
lon_p  = lon_g(mlon_p);  lat_p = lat_g(mlat_p);
land_p = land_g(mlat_p, mlon_p);
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    fl = {'tr','p','qu_tr','qu_p','qv_tr','qv_p'};
    for j = 1:numel(fl)
        P(ip).(fl{j}) = R(ip).(fl{j})(mlat_p, mlon_p); %#ok<SAGROW>
    end
end

%% ===== 色标：Precip4_diff_19lev，色阶分界与刻度对齐 =====
load Precip4_diff_19lev.mat;                        % 文件内变量名为小写 p：precip4_diff_19lev
cmap = double(precip4_diff_19lev);
if max(cmap(:)) > 1; cmap = cmap / 255; end         % 若为 0–255 的 RGB，转换为 0–1
% 去掉最中间的两个白色色阶
n0      = size(cmap, 1);
i_white = find(all(cmap >= white_thr, 2));
[~, ord] = sort(abs(i_white - (n0 + 1) / 2));       % 按离中心的距离排序
i_rm    = i_white(ord(1:min(2, numel(i_white))));
cmap(i_rm, :) = [];
fprintf('原色阶数 %d，去掉中间白色色阶 %d 个（第 %s 行），剩余 %d 个\n', ...
    n0, numel(i_rm), mat2str(sort(i_rm)'), size(cmap, 1));
if numel(i_rm) ~= 2
    warning('中间白色色阶不是 2 个，请检查 white_thr 或 colormap');
end

n_lev    = size(cmap, 1);
lev_w    = 2 * clim_tr / n_lev;                     % 色阶宽度（18 个色阶时为 0.002）
cb_edges = linspace(-clim_tr, clim_tr, n_lev + 1);  % 色阶分界
fprintf('色标：%d 个色阶，范围 ±%g，色阶宽度 %g\n', n_lev, clim_tr, lev_w);

% 刻度标签：按 tick_scale 换算，并按 tick_lbl_step 对称地隔开
k_edge  = round(cb_edges / lev_w * 2);              % 分界的半色阶序号（偶数色阶时为 0, ±2, ±4, …）
k_min   = min(abs(k_edge));                         % 最靠近 0 的分界序号
show_lb = mod((abs(k_edge) - k_min) / 2, tick_lbl_step) == 0;
% 只标注指定数值（10^-9 单位）；其余分界保留刻度线但不标数字
lbl_vals = -18:6:18;
show_lb  = show_lb & ismember(round(cb_edges * tick_scale, 2), lbl_vals);
cb_lbls = arrayfun(@(x) sprintf('%g', round(x * tick_scale, 2)), cb_edges, 'UniformOutput', false);
cb_lbls(~show_lb) = {''};

%% ===== 画图：4 行 2 列 =====
silver_rgb = [0.85 0.85 0.85];

cm2px    = 96/2.54;
fig_w_cm = 18;
% 保持长宽比：纬向每厘米对应的度数与原 40°S–40°N 版本（高 4.23 cm）相同
ph_ref = 4.23;   gap_cm = 5.17 - ph_ref;               % 原面板高度与行间空白
pw_cm  = 0.32*18;
ph_cm  = ph_ref * diff(plot_lat) / diff(plot_lat_ref);  % 60°/80° × 4.23 ≈ 3.17 cm
row_dy = ph_cm + gap_cm;
nrow  = 4;
y_bot = 2.6;                                         % 底部留给横向色标
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

% 矢量：模式网格约 2.8°，每隔 2 个格点画一个（约 5.6°，接近 Fig2 f 的 6°）
vec_skip = 2;
ref_len  = 5;                                        % 参考箭头长度（经度°），同 Fig2 f
iq = 1:vec_skip:numel(lon_p);  jq = 1:vec_skip:numel(lat_p);
[LONq, LATq] = meshgrid(lon_p(iq), lat_p(jq));
[LONs, LATs] = meshgrid(lon_p, lat_p);
land_img = repmat(reshape(silver_rgb,1,1,3), [size(land_p) 1]);

% 所有面板共用一个参考箭头
mag_all = [];
for ip = 1:n_panels
    if ~R(ip).ok; continue; end
    m = hypot(P(ip).qu_tr(jq,iq), P(ip).qv_tr(jq,iq));
    mag_all = [mag_all; m(isfinite(m))]; %#ok<AGROW>
end
ref = nice_ceil(prctile(mag_all, 90));
sc  = ref_len / ref;
fprintf('参考箭头 = %g g kg^-1 m s^-1 decade^-1\n', ref);

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

        % MFC 趋势显著 → 绿点
        sg = P(ip).p < sig_level & ~isnan(dat);
        plot(ax, LONs(sg), LATs(sg), 'o', 'MarkerSize', dot_sz, ...
             'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerEdgeColor', 'none');

        % 矢量：qu 或 qv 显著 → 紫色
        Uq = P(ip).qu_tr(jq,iq);  Vq = P(ip).qv_tr(jq,iq);
        S  = (P(ip).qu_p(jq,iq) < sig_level) | (P(ip).qv_p(jq,iq) < sig_level);
        Un = Uq;  Vn = Vq;  Un(S) = NaN;  Vn(S) = NaN;
        Us = Uq;  Vs = Vq;  Us(~S) = NaN; Vs(~S) = NaN;
        quiver(ax, LONq, LATq, Un*sc, Vn*sc, 0, 'k', 'LineWidth', 0.6, 'MaxHeadSize', 0.5);
        quiver(ax, LONq, LATq, Us*sc, Vs*sc, 0, 'Color', [0.75 0 0.95], ...
               'LineWidth', 0.9, 'MaxHeadSize', 0.5);
        ttl = sprintf('(%c) %s', letters(ip), panel_names{ip});
        if show_dipole_tag && R(ip).dipole && ~ismember(ip, no_dipole_tag_panels)
            ttl = [ttl ' (dipole)'];
        end
    else
        ttl = sprintf('(%c) %s (no data)', letters(ip), panel_names{ip});
    end

    set(ax, 'YDir', 'normal', 'Color', silver_rgb);
    colormap(ax, cmap);
    caxis(ax, [-clim_tr clim_tr]);
    contour(ax, topo_lon_v, topo_lat_v, topo_elev_v, [0 0], 'k-', 'LineWidth', 0.5);

    % 区域框（Fig. 1c）：WIO 橙色，EIO 蓝色（show_boxes = false 时不画）
    if show_boxes
        for ib = 1:size(boxes, 1)
            bx = boxes{ib, 1};  by = boxes{ib, 2};
            plot(ax, bx([1 2 2 1 1]), by([1 1 2 2 1]), '-', ...
                 'Color', boxes{ib, 3}, 'LineWidth', box_lw);
        end
    end
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
    set(ax, 'FontSize', fs, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
    ax.Position = pos;

    draw_ref_arrow(ax, plot_lon, plot_lat, ref_len, sprintf('%g', ref), fs);

    annotation(fig, 'textbox', [pos(1), pos(2)+ph+0.003, pw, 0.015], ...
               'String', ttl, 'FontSize', fs, 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
               'EdgeColor', 'none', 'FitBoxToText', 'off', 'Interpreter', 'tex');
    ax.Position = pos;
    ax_last = ax;
end

%% ===== 底部横向色标（横跨两列）：刻度放在色阶分界上 =====
cb = colorbar(ax_last, 'southoutside');
cb.Position = [xs(1)+0.1*pw, (y_bot - 1.3)/fig_h_cm, xs(2)+pw-xs(1)-0.2*pw, 0.3/fig_h_cm];
cb.Limits     = [-clim_tr clim_tr];
cb.Ticks      = cb_edges;
cb.TickLabels = cb_lbls;
cb.Label.String = sprintf('SON 850-hPa MFC trend, %d–%d (%s kg kg^{-1} s^{-1} decade^{-1})', ...
    yr_start, yr_end, tick_lbl_unit);
cb.FontSize       = fs;
cb.Label.FontSize = fs;
cb.TickDirection  = 'out';
ax_last.Position  = [xs(2), y_cm(end)/fig_h_cm, pw, ph];

set(findall(fig, 'Type', 'axes'), 'TickLength', [0.025 0.035]);

%% ===== 保存 =====
out_dir = 'C:\Users\yingy\Downloads\Rice_mac';
base = fullfile(out_dir, sprintf('FigS_MFC850_SON_trend_%d_%d_7models_MME', yr_start, yr_end));
exportgraphics(fig, [base '.png'], 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig, [base '.tif'], 'Resolution', 600, 'BackgroundColor', 'white');
fprintf('Saved: %s.png / .tif\n', base);

%% ===== 子函数 =====
function F = read_level(file, var, il)
% 读单个气压层，返回 lon x lat x time；-999 → NaN
F = ncread(file, var, [1 1 il 1], [Inf Inf 1 Inf]);
F = double(reshape(F, size(F,1), size(F,2), size(F,4)));
F(F <= -998) = NaN;
end

function [yr, mo] = model_time(file)
% months since YYYY-M-D：按逐月连续重建（防止时间轴末尾异常值）
info  = ncinfo(file);
names = {info.Variables.Name};
if any(strcmp(names, 'time')); tn = 'time'; else; tn = 'time_counter'; end
t_raw = double(ncread(file, tn));
units = ncreadatt(file, tn, 'units');
tok = regexp(units, '(\w+)\s+since\s+(\d{4})-(\d{1,2})', 'tokens', 'once');
if isempty(tok) || ~contains(lower(tok{1}), 'month')
    error('时间单位不是 months since：%s', units);
end
y0 = str2double(tok{2});  m0 = str2double(tok{3});
nt = numel(t_raw);
tm = round(t_raw(1)) + (0:nt-1)';
bad = find(abs(t_raw(:) - tm) > 0.5);
if ~isempty(bad)
    fprintf('  时间轴有 %d 个异常值（第 %s 个），已按逐月连续重建\n', numel(bad), mat2str(bad'));
end
yr = y0 + floor((m0 - 1 + tm) / 12);
mo = mod(m0 - 1 + tm, 12) + 1;
end

function MFC = calc_mfc(QU, QV, lon, lat)
% MFC = -∇·(qV)，球坐标中央差分（全球周期经度）；输入 lon x lat x t
Re   = 6.371e6;
phiM = reshape(deg2rad(lat), 1, []);
dlon = abs(median(diff(lon)));
dQU  = (circshift(QU,-1,1) - circshift(QU,1,1)) / deg2rad(2*dlon);
QVc  = QV .* cos(phiM);
dQV  = NaN(size(QV));
dQV(:,2:end-1,:) = (QVc(:,3:end,:) - QVc(:,1:end-2,:)) ./ (phiM(3:end) - phiM(1:end-2));
MFC  = -(dQU + dQV) ./ (Re * cos(phiM)) * 1e6;
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
x = years(:);  y = y(:);
ok = isfinite(y);  x = x(ok);  y = y(ok);  n = numel(y);
if n < 4; b = NaN; p = NaN; return; end
x = x - mean(x);
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

function draw_ref_arrow(ax, lon_range, lat_range, ref_len, label, fs)
% 面板右上角白底框内画参考箭头和数值
y0 = lat_range(2) - 6;
ht = text(ax, lon_range(2) - 1.5, y0, label, 'FontSize', fs, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
ext = ht.Extent;  delete(ht);
x0 = ext(1) - ref_len - 1.5;
yb = y0 - ext(4)/2 - 1.5;
rectangle(ax, 'Position', [x0-1.5, yb, lon_range(2)-(x0-1.5), lat_range(2)-yb], ...
          'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 0.5);
quiver(ax, x0, y0, ref_len, 0, 0, 'k', 'LineWidth', 0.6, 'MaxHeadSize', 0.8);
text(ax, lon_range(2) - 1.5, y0, label, 'FontSize', fs, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
