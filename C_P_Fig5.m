%% ========================================================================
%  LMDZ：提取 ts，以及最接近 850 hPa 的 q、u、v，各自单独保存，
%        并重插到 5°×5° 网格（lon 0:5:355，lat 90:-5:-90）
%  输入：C:\Users\yingy\Downloads\Rice_mac\Total.LMDZ_Monthly_Original.Grid.nc
%        原网格 lon 0–357.5 (2.5°, 144)，lat 90 → -90 (~1.27°, 143)
%  重插方法：原网格中心落在每个 5° 格子（中心 ± 2.5°）内的点做 cos(lat) 加权平均，
%            缺测（-999）自动跳过
%  输出（每个变量一个文件）：
%     原网格：LMDZ_<var>_OriginalGrid.nc   （save_original = true 时）
%     5° 网格：LMDZ_<var>_5deg.nc
%% ========================================================================
clear; clc;

%% ===== 设置 =====
fn_in   = 'C:\Users\yingy\Downloads\Rice_mac\Total.LMDZ_Monthly_Original.Grid.nc';
out_dir = 'C:\Users\yingy\Downloads\Rice_mac\LMDZ_5deg';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save_original = true;          % 是否同时保存原网格上的单变量文件

% 目标网格
lon_t = (0:5:355)';            % 72
lat_t = (90:-5:-90)';          % 37，北 → 南

%% ===== 坐标与层 =====
lon_s = double(ncread(fn_in, 'lon'));
lat_s = double(ncread(fn_in, 'lat'));
p     = double(ncread(fn_in, 'presnivs'));
[~, k] = min(abs(p - 850));
fprintf('最接近 850 hPa 的层：k = %d, presnivs = %.6f hPa\n', k, p(k));

%% ===== 时间轴（检查并修复异常值）=====
t_raw   = double(ncread(fn_in, 'time_counter'));
t_units = ncreadatt(fn_in, 'time_counter', 'units');     % months since 1979-1-15
t_cal   = ncreadatt(fn_in, 'time_counter', 'calendar');
nt  = numel(t_raw);
bad = find(t_raw(:) ~= (0:nt-1)');
if ~isempty(bad)
    fprintf('时间轴有 %d 个异常值（第 %s 个），按逐月连续重建为 0:%d\n', ...
        numel(bad), mat2str(bad'), nt-1);
end
time = (0:nt-1)';
fprintf('时间：1979-01 到 %d-%02d，共 %d 个月\n', 1979 + floor((nt-1)/12), mod(nt-1,12)+1, nt);

%% ===== 重插权重矩阵（只算一次）=====
[Wlon, Wlat] = box_weights(lon_s, lat_s, lon_t, lat_t);

%% ===== 变量设置 =====
% {输出名, 原变量名, 是否三维, long_name}
V = { ...
    'ts',   'ts', false, 'Temperature (surface)'; ...
    'q850', 'q',  true,  sprintf('Specific humidity at %.2f hPa (model level %d)', p(k), k); ...
    'u850', 'u',  true,  sprintf('u-wind at %.2f hPa (model level %d)', p(k), k); ...
    'v850', 'v',  true,  sprintf('v-wind at %.2f hPa (model level %d)', p(k), k)};

%% ===== 逐变量：读取 → 保存原网格 → 重插 → 保存 5° =====
for iv = 1:size(V, 1)
    oname = V{iv,1};  vname = V{iv,2};  is3d = V{iv,3};  lname = V{iv,4};
    units = ncreadatt(fn_in, vname, 'units');
    fprintf('\n===== %s =====\n', oname);

    if is3d
        X = ncread(fn_in, vname, [1 1 k 1], [Inf Inf 1 Inf]);   % lon x lat x 1 x t
    else
        X = ncread(fn_in, vname);                               % lon x lat x t
    end
    X = double(squeeze(X));
    X(X <= -998) = NaN;
    fprintf('原网格：size = %s, NaN 比例 = %.4f, 范围 = [%.4g, %.4g] %s\n', ...
        mat2str(size(X)), mean(isnan(X(:))), min(X(:)), max(X(:)), units);

    % --- 原网格单变量文件 ---
    if save_original
        fo = fullfile(out_dir, sprintf('LMDZ_%s_OriginalGrid.nc', oname));
        write_nc(fo, oname, X, lon_s, lat_s, time, t_units, t_cal, lname, units, ...
            sprintf('Extracted from %s', fn_in), p, k, is3d);
        fprintf('已保存：%s\n', fo);
    end

    % --- 重插到 5° ---
    Y = regrid_box(X, Wlon, Wlat);
    fprintf('5° 网格：size = %s, NaN 比例 = %.4f, 范围 = [%.4g, %.4g] %s\n', ...
        mat2str(size(Y)), mean(isnan(Y(:))), min(Y(:)), max(Y(:)), units);

    f5 = fullfile(out_dir, sprintf('LMDZ_%s_5deg.nc', oname));
    write_nc(f5, oname, Y, lon_t, lat_t, time, t_units, t_cal, lname, units, ...
        ['5x5 deg cos(lat)-weighted box average (lon 0:5:355, lat 90:-5:-90) of ' fn_in], ...
        p, k, is3d);
    fprintf('已保存：%s\n', f5);
end
fprintf('\n全部完成。输出目录：%s\n', out_dir);

%% ===== 子函数 =====
function [Wlon, Wlat] = box_weights(lon_s, lat_s, lon_t, lat_t)
% 原网格中心落在目标格子 [中心 ± 半格距] 内的点参与平均；纬向权重 cos(lat)
% 经度按 360° 周期处理；lat = ±90 的格子只包含极点一侧的半个格子
dlon = abs(median(diff(lon_t)));
dlat = abs(median(diff(lat_t)));
Wlon = zeros(numel(lon_t), numel(lon_s));
for i = 1:numel(lon_t)
    d = mod(lon_s(:)' - lon_t(i) + 180, 360) - 180;
    Wlon(i, d >= -dlon/2 & d < dlon/2) = 1;
end
w = cosd(lat_s(:))';  w(w < 0) = 0;
w(abs(lat_s(:)') >= 90) = 0;                 % 极点处 cos = 0
Wlat = zeros(numel(lat_t), numel(lat_s));
for j = 1:numel(lat_t)
    lo = lat_t(j) - dlat/2;  hi = lat_t(j) + dlat/2;
    sel = lat_s(:)' >= lo & lat_s(:)' < hi;
    if hi > 90;  sel = lat_s(:)' >= lo; end   % 北极格子含 90°N
    Wlat(j, sel) = w(sel);
    if sum(Wlat(j,:)) == 0 && any(sel)       % 只有极点时退化为等权
        Wlat(j, sel) = 1;
    end
end
end

function Y = regrid_box(X, Wlon, Wlat)
% X: lon_s x lat_s x t → Y: lon_t x lat_t x t（缺测自动跳过）
nt = size(X, 3);
Y = NaN(size(Wlon, 1), size(Wlat, 1), nt);
for it = 1:nt
    A = X(:,:,it);
    M = ~isnan(A);  A(~M) = 0;
    num = Wlon * A * Wlat';
    den = Wlon * double(M) * Wlat';
    B = num ./ den;  B(den == 0) = NaN;
    Y(:,:,it) = B;
end
end

function write_nc(fn, vname, X, lon, lat, time, t_units, t_cal, lname, units, note, p, k, is3d)
% 写单变量 NetCDF：dims lon x lat x time，缺测 -999
if exist(fn, 'file'); delete(fn); end
nccreate(fn, 'lon',  'Dimensions', {'lon',  numel(lon)}, 'Datatype', 'double');
nccreate(fn, 'lat',  'Dimensions', {'lat',  numel(lat)}, 'Datatype', 'double');
nccreate(fn, 'time', 'Dimensions', {'time', Inf},        'Datatype', 'double');
ncwrite(fn, 'lon',  lon);
ncwrite(fn, 'lat',  lat);
ncwrite(fn, 'time', time);
ncwriteatt(fn, 'lon',  'units', 'degrees_east');
ncwriteatt(fn, 'lat',  'units', 'degrees_north');
ncwriteatt(fn, 'time', 'units', t_units);
ncwriteatt(fn, 'time', 'calendar', t_cal);

Xs = single(X);  Xs(isnan(Xs)) = -999;
nccreate(fn, vname, 'Dimensions', {'lon', numel(lon), 'lat', numel(lat), 'time', Inf}, ...
    'Datatype', 'single', 'FillValue', single(-999), 'DeflateLevel', 4);
ncwrite(fn, vname, Xs);
ncwriteatt(fn, vname, 'long_name', lname);
ncwriteatt(fn, vname, 'units', units);
if is3d
    ncwriteatt(fn, vname, 'pressure_level_hPa', p(k));
    ncwriteatt(fn, vname, 'model_level_index', k);
end
ncwriteatt(fn, '/', 'source', 'LMDZ monthly output');
ncwriteatt(fn, '/', 'note', note);
ncwriteatt(fn, '/', 'time_note', 'time axis rebuilt as consecutive months (0..nt-1) since 1979-1-15');
end