lat = ncread('/Users/yingyingfeng/Downloads/Isotope-IOD/AIRS_RS_monthly_1deg.nc','lat');
lon =  ncread('/Users/yingyingfeng/Downloads/Isotope-IOD/AIRS_RS_monthly_1deg.nc','lon');
dD = ncread('/Users/yingyingfeng/Downloads/Isotope-IOD/AIRS_RS_monthly_1deg.nc','dD');
pressure = ncread('/Users/yingyingfeng/Downloads/Isotope-IOD/AIRS_RS_monthly_1deg.nc','Pressure');
pressure_per = permute(pressure,[2 1 3 4]);
pressure_fn = flipud(pressure_per);

dD_per = permute(dD,[2 1 3 4]);
dD_fn = flipud(dD_per);
lat_fn = flipud(lat);
lon_fn = lon;

clear lat; clear lon; clear dD;

lat = lat_fn;
lon = lon_fn;
dD  = dD_fn;

pressure_levels = [1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10];

dates  = datetime(2002, 9, 1) + calmonths(0:210);
years  = year(dates);
months = month(dates);

%% ============================================================
fprintf('计算气候态月平均...\n');
climatology = NaN(180, 360, 17, 12);
for mo = 1:12
    idx = (months == mo);
    climatology(:,:,:,mo) = mean(dD(:,:,:,idx), 4, 'omitnan');
end

dD_anom = NaN(size(dD));
for t = 1:211
    mo = months(t);
    dD_anom(:,:,:,t) = dD(:,:,:,t) - climatology(:,:,:,mo);
end
fprintf('季节性循环去除完成！\n');


%% Step 2: 参数设置
%% ============================================================
elnino_posIOD_yr = [2003 2009 2014];

%% Step 1: 定义经纬度网格
lat_vec = (89.5:-1:-89.5)';    % [180x1]，从北到南
lon_vec = (-179.5:1:179.5)';   % [360x1]，从西到东
[lon2d, lat2d] = meshgrid(lon_vec, lat_vec);  % 各[180x360]

%% Step 2: 生成陆地mask
fprintf('生成陆地mask...\n');
load coastlines  % 加载内置海岸线，变量名：coastlat, coastlon
land_mask = inpolygon(lon2d, lat2d, coastlon, coastlat);
land_mask = logical(land_mask);
fprintf('陆地格点比例: %.1f%%\n', sum(land_mask(:))/numel(land_mask)*100);

%% Step 3: 定义印度洋区域mask
% 印度洋范围：40°E–110°E, 30°S–30°N
lat_range = [-30 30];
lon_range = [40  110];

IO_mask = (lat2d >= lat_range(1) & lat2d <= lat_range(2)) & ...
          (lon2d >= lon_range(1) & lon2d <= lon_range(2));

% 最终mask：印度洋范围内去除陆地
IO_ocean_mask = IO_mask & ~land_mask;  % [180x360]
fprintf('印度洋海洋格点数: %d\n', sum(IO_ocean_mask(:)));

%% Step 4: 参数设置
elnino_posIOD_yr =  [2003 2009 2014];

month_order  = [3 4 5 6 7 8 9 10 11 12 1 2];
month_labels = {'Mar','Apr','May','Jun','Jul','Aug', ...
                'Sep','Oct','Nov','Dec','Jan','Feb'};

pressure_levels = [1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10];

% 只保留对流层（200 hPa及以下）
trop_idx      = pressure_levels >= 200;
pressure_trop = pressure_levels(trop_idx);  % [1000 925 850 700 600 500 400 300 250 200]
n_trop        = length(pressure_trop);
n_months      = length(month_order);
n_years       = length(elnino_posIOD_yr);

% 纬度余弦权重
lat_weights = cosd(lat2d);           % [180x360]
lat_weights(~IO_ocean_mask) = NaN;   % 只保留印度洋海洋格点

%% Step 5: 计算每个年份、月份、气压层的印度洋区域加权平均
dD_IO_avg = NaN(n_months, n_trop, n_years);

fprintf('计算印度洋区域加权平均...\n');
for yi = 1:n_years
    yr = elnino_posIOD_yr(yi);
    fprintf('  处理年份: %d\n', yr);

    for mi = 1:n_months
        mo = month_order(mi);

        % 跨年处理：1月和2月属于次年
        if mo <= 2
            target_yr = yr + 1;
        else
            target_yr = yr;
        end

        % 找时间索引
        t_idx = find(years == target_yr & months == mo);
        if isempty(t_idx)
            fprintf('  Warning: yr=%d mo=%d not found\n', target_yr, mo);
            continue;
        end

        % 取该时间步的dD anomaly [180x360x17]
        dD_t = dD_anom(:,:,:,t_idx(1));

        % 对每个对流层气压层计算加权平均
        for pi = 1:n_trop
            orig_pi = find(pressure_levels == pressure_trop(pi));

            dD_level = dD_t(:,:,orig_pi);      % [180x360]
            dD_masked = dD_level;
            dD_masked(~IO_ocean_mask) = NaN;   % 应用mask

            % 纬度加权
            w = lat_weights;
            w(isnan(dD_masked)) = NaN;

            dD_IO_avg(mi, pi, yi) = ...
                nansum(dD_masked(:) .* w(:)) / nansum(w(:));
        end
    end
end

%% Step 6: Composite mean（对所有年份取平均）
dD_hov = mean(dD_IO_avg, 3, 'omitnan');  % [12 x n_trop]
fprintf('计算完成！\n');

%% Step 7: 画Hovmoeller图
% 排列：X轴=月份，Y轴=气压层（1000hPa在下，200hPa在上）
% data_plot: [n_trop x 12]，flipud使1000hPa在最下行
data_plot     = flipud(dD_hov');                    % [n_trop x 12]
ytick_labels  = fliplr(arrayfun(@(p) sprintf('%d hPa', p), ...
                pressure_trop, 'UniformOutput', false));  % 200→1000 从上到下

figure('Position', [100 100 800 500], 'Color', 'w');
ax = axes;

imagesc(ax, 1:n_months, 1:n_trop, data_plot);
set(ax, 'YDir', 'normal');  % y=1在最下（对应1000hPa）

load BlueDarkRed18.mat;

% 蓝白红colormap
colormap(ax, BlueDarkRed18);
clim([-9 9]);

% Colorbar
cb = colorbar(ax, 'Location', 'eastoutside');
cb.Label.String = '\deltaD anomaly (‰)';
cb.Label.FontSize = 14;
cb.FontSize      = 14;

% X轴：月份
xticks(ax, 1:n_months);
xticklabels(ax, month_labels);
xlabel(ax, 'Month', 'FontSize', 13, 'FontWeight', 'bold');
ax.XAxis.TickLength = [0 0];

% Y轴：气压层
yticks(ax, 1:n_trop);
yticklabels(ax, ytick_labels);
ylabel(ax, 'Pressure (hPa)', 'FontSize', 13, 'FontWeight', 'bold');
ax.YAxis.TickLength = [0 0];

% 标题
title(ax, sprintf(['Independent El nino Events dD Anomaly\n' ...
    'Years: %s'], ...
    strjoin(arrayfun(@num2str, elnino_posIOD_yr, 'UniformOutput', false), ', ')), ...
    'FontSize', 14, 'FontWeight', 'bold');

% 细网格线
hold(ax, 'on');
for xi = 1.5 : 1 : n_months - 0.5
    plot(ax, [xi xi], [0.5 n_trop+0.5], '-', ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 0.4);
end
for yi = 1.5 : 1 : n_trop - 0.5
    plot(ax, [0.5 n_months+0.5], [yi yi], '-', ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 0.4);
end
hold(ax, 'off');

set(ax, 'FontSize', 14, 'Box', 'on', 'LineWidth', 1);

% 保存
exportgraphics(gcf, 'dD_Hovmoeller_IO_Independent_ElNino.png', 'Resolution', 600);


