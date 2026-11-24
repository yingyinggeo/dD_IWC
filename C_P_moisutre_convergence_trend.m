%% ========================================================================
% ERA5 850-hPa Moisture Convergence + Moisture Flux TREND
%
% Shading:
%   Trend of 850-hPa moisture convergence
%
% Vectors:
%   Trend of 850-hPa moisture flux (qu, qv)
%
% Significance:
%   MFC: p < 0.05 -> bright green dots
%
%   Moisture-flux vector:
%   p(QU) < 0.05 OR p(QV) < 0.05
%   -> bright purple vector
%
%   Non-significant vector -> black
%
% Trend period:
%   SON 2003-2020
%
% Trend unit:
%   per decade
%
% Domain:
%   0-180E, 40S-40N
%
% Land:
%   masked
%
% ========================================================================

clear;
clc;
close all;


%% ========================================================================
% 0. FILE PATHS
% ========================================================================

data_dir = 'C:\Users\yingy\Downloads\ERA5';

qfile = fullfile(data_dir, ...
    'era5_qs_90N90S_1deg_2000_2024.nc');

ufile = fullfile(data_dir, ...
    'era5_uwind_90N90S_1deg_2000_2024.nc');

vfile = fullfile(data_dir, ...
    'era5_vwind_90N90S_2000_2024.nc');


%% ========================================================================
% 1. VARIABLE NAMES
% ========================================================================

q_var = 'q';
u_var = 'u';
v_var = 'v';


%% ========================================================================
% 2. READ COORDINATES
% ========================================================================

lon  = double(ncread(ufile,'longitude'));
lat  = double(ncread(ufile,'latitude'));
lev  = double(ncread(ufile,'pressure_level'));
time = double(ncread(ufile,'valid_time'));

fprintf('\n============================================\n');
fprintf('ERA5 grid information\n');
fprintf('============================================\n');

fprintf('Longitude: %.1f to %.1f, N = %d\n', ...
    lon(1),lon(end),length(lon));

fprintf('Latitude : %.1f to %.1f, N = %d\n', ...
    lat(1),lat(end),length(lat));


%% ========================================================================
% 3. FIND 850 hPa
% ========================================================================

[~,ilev] = min(abs(lev-850));

fprintf('Selected pressure level = %.0f hPa\n',lev(ilev));

if abs(lev(ilev)-850) > 0.1
    error('850 hPa level was not found.');
end


%% ========================================================================
% 4. TIME
% ========================================================================

dates = datetime(1970,1,1) + seconds(time);

yr = year(dates);
mo = month(dates);

fprintf('First time = %s\n',datestr(dates(1),'yyyy-mm-dd'));
fprintf('Last time  = %s\n',datestr(dates(end),'yyyy-mm-dd'));
fprintf('Number of months = %d\n',length(dates));


%% ========================================================================
% 5. READ 850-hPa U
% ========================================================================

fprintf('\nReading U...\n');

u850 = squeeze(ncread(ufile,u_var, ...
    [1 1 ilev 1], ...
    [Inf Inf 1 Inf]));

u850 = double(u850);


%% ========================================================================
% 6. READ 850-hPa V
% ========================================================================

fprintf('Reading V...\n');

v850 = squeeze(ncread(vfile,v_var, ...
    [1 1 ilev 1], ...
    [Inf Inf 1 Inf]));

v850 = double(v850);


%% ========================================================================
% 7. READ 850-hPa Q
% ========================================================================

fprintf('Reading Q...\n');

q850 = squeeze(ncread(qfile,q_var, ...
    [1 1 ilev 1], ...
    [Inf Inf 1 Inf]));

q850 = double(q850);


%% ========================================================================
% 8. DIMENSION CHECK
% ========================================================================

if ~isequal(size(u850),size(v850))
    error('U and V dimensions differ.');
end

if ~isequal(size(u850),size(q850))
    error('U and Q dimensions differ.');
end

nlon = length(lon);
nlat = length(lat);
nt   = length(time);


%% ========================================================================
% 9. MOISTURE FLUX
%
% QU = q*u
% QV = q*v
% ========================================================================

fprintf('Calculating moisture flux...\n');

QU = q850 .* u850;

QV = q850 .* v850;


%% ========================================================================
% 10. EARTH GEOMETRY
% ========================================================================

Re = 6.371e6;

phi = deg2rad(lat);

cosphi = cos(phi);


% Actual longitude spacing
dlon = abs(median(diff(lon)));

% centered difference spans two grid cells
dlambda_center = deg2rad(2*dlon);


%% ========================================================================
% 11. MONTHLY MOISTURE CONVERGENCE
%
% MFC = -div(qV)
% ========================================================================

fprintf('Calculating monthly moisture convergence...\n');

MFC = NaN(nlon,nlat,nt);


for it = 1:nt

    if mod(it,12)==0
        fprintf('Processing month %d / %d\n',it,nt);
    end

    qu = QU(:,:,it);
    qv = QV(:,:,it);


    %% ZONAL DERIVATIVE

    dQUdlambda = NaN(nlon,nlat);


    dQUdlambda(2:nlon-1,:) = ...
        (qu(3:nlon,:) - qu(1:nlon-2,:)) ...
        ./ dlambda_center;


    % periodic boundaries
    dQUdlambda(1,:) = ...
        (qu(2,:) - qu(nlon,:)) ...
        ./ dlambda_center;


    dQUdlambda(nlon,:) = ...
        (qu(1,:) - qu(nlon-1,:)) ...
        ./ dlambda_center;


    %% MERIDIONAL DERIVATIVE

    QVcos = qv .* reshape(cosphi,1,nlat);

    dQVcosdphi = NaN(nlon,nlat);


    for j = 2:nlat-1

        dphi_center = phi(j+1) - phi(j-1);

        dQVcosdphi(:,j) = ...
            (QVcos(:,j+1) - QVcos(:,j-1)) ...
            ./ dphi_center;

    end


    %% DIVERGENCE

    divQV = NaN(nlon,nlat);


    for j = 2:nlat-1

        divQV(:,j) = ...
            dQUdlambda(:,j) ./ ...
            (Re*cosphi(j)) ...
            + ...
            dQVcosdphi(:,j) ./ ...
            (Re*cosphi(j));

    end


    %% MOISTURE CONVERGENCE

    MFC(:,:,it) = -divQV;

end


fprintf('Monthly MFC finished.\n');


%% ========================================================================
% 12. CALCULATE ANNUAL SON MEANS
%
% One SON value per year:
%
% 2003
% 2004
% ...
% 2020
%
% ========================================================================

yr_start = 2003;
yr_end   = 2020;

years = (yr_start:yr_end)';

nyear = length(years);


MFC_year = NaN(nlon,nlat,nyear);

QU_year  = NaN(nlon,nlat,nyear);

QV_year  = NaN(nlon,nlat,nyear);


fprintf('\nCalculating annual SON means...\n');


for iy = 1:nyear

    yy = years(iy);

    idx = ...
        yr == yy & ...
        ismember(mo,[9 10 11]);


    fprintf('%d: %d months\n',yy,sum(idx));


    if sum(idx) ~= 3
        warning('Year %d does not contain exactly 3 SON months.',yy);
    end


    MFC_year(:,:,iy) = ...
        mean(MFC(:,:,idx),3,'omitnan');


    QU_year(:,:,iy) = ...
        mean(QU(:,:,idx),3,'omitnan');


    QV_year(:,:,iy) = ...
        mean(QV(:,:,idx),3,'omitnan');

end


%% ========================================================================
% 13. LINEAR TREND + SIGNIFICANCE
%
% Model:
%
% Y = beta0 + beta1 * year
%
% beta1 = trend per year
%
% p value:
% two-sided t test for regression slope
%
% ========================================================================

fprintf('\nCalculating trends and significance...\n');


MFC_trend = NaN(nlon,nlat);
MFC_p     = NaN(nlon,nlat);

QU_trend  = NaN(nlon,nlat);
QU_p      = NaN(nlon,nlat);

QV_trend  = NaN(nlon,nlat);
QV_p      = NaN(nlon,nlat);


% Center time for numerical stability
x = double(years);

x = x - mean(x);


for ii = 1:nlon

    if mod(ii,30)==0
        fprintf('Longitude %d / %d\n',ii,nlon);
    end


    for jj = 1:nlat


        %% ================================================================
        % MFC TREND
        % ================================================================

        y = squeeze(MFC_year(ii,jj,:));

        good = isfinite(y);


        if sum(good) >= 10

            xx = x(good);
            yy1 = y(good);


            X = [ones(length(xx),1) xx];

            b = X \ yy1;


            % slope
            MFC_trend(ii,jj) = b(2);


            % residuals
            residual = yy1 - X*b;

            n = length(yy1);

            df = n - 2;


            SSE = sum(residual.^2);

            s2 = SSE/df;


            C = inv(X'*X);

            SE = sqrt(s2*C(2,2));


            tstat = b(2)/SE;


            MFC_p(ii,jj) = ...
                2*(1-tcdf(abs(tstat),df));

        end


        %% ================================================================
        % QU TREND
        % ================================================================

        y = squeeze(QU_year(ii,jj,:));

        good = isfinite(y);


        if sum(good) >= 10

            xx = x(good);
            yy1 = y(good);


            X = [ones(length(xx),1) xx];

            b = X \ yy1;


            QU_trend(ii,jj) = b(2);


            residual = yy1-X*b;

            n = length(yy1);

            df = n-2;


            SSE = sum(residual.^2);

            s2 = SSE/df;


            C = inv(X'*X);

            SE = sqrt(s2*C(2,2));


            tstat = b(2)/SE;


            QU_p(ii,jj) = ...
                2*(1-tcdf(abs(tstat),df));

        end


        %% ================================================================
        % QV TREND
        % ================================================================

        y = squeeze(QV_year(ii,jj,:));

        good = isfinite(y);


        if sum(good) >= 10

            xx = x(good);
            yy1 = y(good);


            X = [ones(length(xx),1) xx];

            b = X \ yy1;


            QV_trend(ii,jj) = b(2);


            residual = yy1-X*b;

            n = length(yy1);

            df = n-2;


            SSE = sum(residual.^2);

            s2 = SSE/df;


            C = inv(X'*X);

            SE = sqrt(s2*C(2,2));


            tstat = b(2)/SE;


            QV_p(ii,jj) = ...
                2*(1-tcdf(abs(tstat),df));

        end

    end

end


%% ========================================================================
% 14. CONVERT TREND TO PER DECADE
%
% MFC:
% kg kg^-1 s^-1 decade^-1
%
% QU/QV:
% kg kg^-1 m s^-1 decade^-1
% ========================================================================

MFC_trend = MFC_trend * 10;

QU_trend = QU_trend * 10;

QV_trend = QV_trend * 10;


%% ========================================================================
% 15. MFC PLOT UNIT
%
% -> 10^-6 kg kg^-1 s^-1 decade^-1
% ========================================================================

MFC_trend = MFC_trend * 1e6;


%% ========================================================================
% 16. VECTOR SIGNIFICANCE
%
% User definition:
%
% QU significant OR QV significant
%
% p < 0.05
% ========================================================================

MFC_sig = MFC_p < 0.05;

VECTOR_sig = ...
    QU_p < 0.05 | ...
    QV_p < 0.05;


%% ========================================================================
% 17. LONGITUDE -> 0-359
% ========================================================================

lon360 = mod(lon,360);

[lon360,ix] = sort(lon360);


MFC_trend = MFC_trend(ix,:);
MFC_sig   = MFC_sig(ix,:);

QU_trend = QU_trend(ix,:);
QV_trend = QV_trend(ix,:);

QU_p = QU_p(ix,:);
QV_p = QV_p(ix,:);

VECTOR_sig = VECTOR_sig(ix,:);


%% ========================================================================
% 18. TRANSPOSE TO LAT x LON
% ========================================================================

MFC_plot = MFC_trend';

SIG_plot = MFC_sig';


QU_plot = QU_trend';

QV_plot = QV_trend';


VSIG_plot = VECTOR_sig';


%% ========================================================================
% 19. DOMAIN
% ========================================================================

lon_range = [0 180];

lat_range = [-40 40];


ilon = ...
    lon360 >= lon_range(1) & ...
    lon360 <= lon_range(2);


ilat = ...
    lat >= lat_range(1) & ...
    lat <= lat_range(2);


lon_sub = lon360(ilon);

lat_sub = lat(ilat);


MFC_sub = MFC_plot(ilat,ilon);

SIG_sub = SIG_plot(ilat,ilon);


QU_sub = QU_plot(ilat,ilon);

QV_sub = QV_plot(ilat,ilon);


VSIG_sub = VSIG_plot(ilat,ilon);


%% ========================================================================
% 20. LATITUDE SOUTH -> NORTH
% ========================================================================

if lat_sub(1) > lat_sub(end)

    lat_sub = flipud(lat_sub);

    MFC_sub = flipud(MFC_sub);

    SIG_sub = flipud(SIG_sub);

    QU_sub = flipud(QU_sub);

    QV_sub = flipud(QV_sub);

    VSIG_sub = flipud(VSIG_sub);

end


%% ========================================================================
% 21. TOPOGRAPHY
% ========================================================================

load topo;


topo_lon = 0.5:359.5;

topo_lat = -89.5:89.5;


topo_ext = [topo topo(:,1)];

topo_lon_ext = [topo_lon 360.5];


buf = 5;


tlon_idx = ...
    topo_lon_ext >= lon_range(1)-buf & ...
    topo_lon_ext <= lon_range(2)+buf;


tlat_idx = ...
    topo_lat >= lat_range(1)-buf & ...
    topo_lat <= lat_range(2)+buf;


topo_lon_v = topo_lon_ext(tlon_idx);

topo_lat_v = topo_lat(tlat_idx);

topo_elev_v = ...
    topo_ext(tlat_idx,tlon_idx);


%% ========================================================================
% 22. LAND MASK
% ========================================================================

[LON_sub,LAT_sub] = ...
    meshgrid(lon_sub,lat_sub);


land_elev = interp2( ...
    topo_lon_ext, ...
    topo_lat, ...
    double(topo_ext), ...
    LON_sub, ...
    LAT_sub, ...
    'nearest');


land_mask = land_elev > 0;


% mask all plotted quantities over land

MFC_sub(land_mask) = NaN;

SIG_sub(land_mask) = false;


QU_sub(land_mask) = NaN;

QV_sub(land_mask) = NaN;

VSIG_sub(land_mask) = false;


%% ========================================================================
% 23. COLORMAP
% ========================================================================

if exist('BlueWhiteOrangeRed.mat','file')

    C = load('BlueWhiteOrangeRed.mat');

    fn = fieldnames(C);

    cmap = C.(fn{1});

else

    warning('BlueWhiteOrangeRed.mat not found. Using turbo.');

    cmap = turbo(256);

end


%% ========================================================================
% 24. FIGURE SIZE
%
% Actual geographic ratio:
%
% 180 / 80 = 2.25
% ========================================================================

cm2px = 96/2.54;


fig_w_cm = 18;

fig_w = round(fig_w_cm*cm2px);


lon_span = diff(lon_range);

lat_span = diff(lat_range);

map_aspect = lon_span/lat_span;


pw_norm = 0.84;


panel_w_px = round(pw_norm*fig_w);

panel_h_px = round(panel_w_px/map_aspect);


mt_cm = 1.3;

mb_cm = 2.8;


mt_px = round(mt_cm*cm2px);

mb_px = round(mb_cm*cm2px);


fig_h = ...
    mt_px + ...
    panel_h_px + ...
    mb_px;


fig_h_cm = fig_h/cm2px;


mt_norm = mt_px/fig_h;

ph_norm = panel_h_px/fig_h;


pos_y = 1-mt_norm-ph_norm;


%% ========================================================================
% 25. CREATE FIGURE
% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Position',[50 50 fig_w fig_h], ...
    'PaperUnits','centimeters', ...
    'PaperSize',[fig_w_cm fig_h_cm]);


%% ========================================================================
% 26. AXES
% ========================================================================

ax = axes( ...
    'Position',[ ...
    (1-pw_norm)/2, ...
    pos_y, ...
    pw_norm, ...
    ph_norm]);


hold(ax,'on');


%% ========================================================================
% 27. MFC TREND SHADING
% ========================================================================

imagesc( ...
    ax, ...
    lon_sub, ...
    lat_sub, ...
    MFC_sub, ...
    'AlphaData',~isnan(MFC_sub));


set(ax, ...
    'YDir','normal', ...
    'Color',[0.90 0.90 0.90]);


colormap(ax,cmap);


%% ========================================================================
% 28. COLOR RANGE
%
% Adjust this after seeing the trend magnitude.
%
% Current starting choice:
% -0.02 to +0.02 per decade
% ========================================================================

clim(ax,[-0.02 0.02]);


%% ========================================================================
% 29. SIGNIFICANT MFC TREND
%
% Bright green dots
%
% Plot every 3 degrees to avoid excessive dots.
% ========================================================================

sig_skip = 3;


ilon_sig = 1:sig_skip:length(lon_sub);

ilat_sig = 1:sig_skip:length(lat_sub);


SIG_sparse = ...
    SIG_sub(ilat_sig,ilon_sig);


LON_sig = ...
    LON_sub(ilat_sig,ilon_sig);

LAT_sig = ...
    LAT_sub(ilat_sig,ilon_sig);


sig_idx = ...
    SIG_sparse & ...
    isfinite(MFC_sub(ilat_sig,ilon_sig));


% Bright green
green_sig = [0.1 0.95 0.1];


scatter( ...
    ax, ...
    LON_sig(sig_idx), ...
    LAT_sig(sig_idx), ...
    7, ...
    green_sig, ...
    'filled', ...
    'MarkerEdgeColor','none');


%% ========================================================================
% 30. MOISTURE-FLUX TREND VECTORS
%
% Every 6 degrees
%
% Black:
% non-significant
%
% Bright purple:
% QU OR QV significant at p < 0.05
% ========================================================================

skip = 6;


ilon_q = 1:skip:length(lon_sub);

ilat_q = 1:skip:length(lat_sub);


lon_q = lon_sub(ilon_q);

lat_q = lat_sub(ilat_q);


QU_q = QU_sub(ilat_q,ilon_q);

QV_q = QV_sub(ilat_q,ilon_q);


VSIG_q = VSIG_sub(ilat_q,ilon_q);


[LON_q,LAT_q] = ...
    meshgrid(lon_q,lat_q);


%% ------------------------------------------------------------------------
% Non-significant vectors
% -------------------------------------------------------------------------

nonsig_vector = ...
    ~VSIG_q & ...
    isfinite(QU_q) & ...
    isfinite(QV_q);


QU_nonsig = QU_q;

QV_nonsig = QV_q;


QU_nonsig(~nonsig_vector) = NaN;

QV_nonsig(~nonsig_vector) = NaN;


quiver( ...
    ax, ...
    LON_q, ...
    LAT_q, ...
    QU_nonsig, ...
    QV_nonsig, ...
    'k', ...
    'LineWidth',0.7, ...
    'AutoScale','on', ...
    'AutoScaleFactor',1.1);


%% ------------------------------------------------------------------------
% Significant vectors
%
% Significant if:
%
% p(QU)<0.05 OR p(QV)<0.05
% -------------------------------------------------------------------------

sig_vector = ...
    VSIG_q & ...
    isfinite(QU_q) & ...
    isfinite(QV_q);


QU_sig = QU_q;

QV_sig = QV_q;


QU_sig(~sig_vector) = NaN;

QV_sig(~sig_vector) = NaN;


% Bright purple
purple_sig = [0.75 0.00 0.95];


quiver( ...
    ax, ...
    LON_q, ...
    LAT_q, ...
    QU_sig, ...
    QV_sig, ...
    'Color',purple_sig, ...
    'LineWidth',1.0, ...
    'AutoScale','on', ...
    'AutoScaleFactor',1.1);


%% ========================================================================
% 31. COASTLINE
% ========================================================================

contour( ...
    ax, ...
    topo_lon_v, ...
    topo_lat_v, ...
    topo_elev_v, ...
    [0 0], ...
    'k-', ...
    'LineWidth',0.5);


%% ========================================================================
% 32. AXIS
% ========================================================================

xlim(ax,lon_range);

ylim(ax,lat_range);


daspect(ax,[1 1 1]);


set(ax, ...
    'XTick',0:30:180, ...
    'XTickLabel', ...
    {'0°','30°E','60°E','90°E', ...
     '120°E','150°E','180°'}, ...
    'YTick',-40:20:40, ...
    'YTickLabel', ...
    {'40°S','20°S','0°','20°N','40°N'}, ...
    'FontSize',10, ...
    'TickDir','out', ...
    'Box','off');


%% ========================================================================
% 33. TOP / RIGHT BORDER
% ========================================================================

xl = ax.XLim;

yl = ax.YLim;


plot(ax, ...
    [xl(2) xl(2)], ...
    yl, ...
    'k-', ...
    'LineWidth',0.5);


plot(ax, ...
    xl, ...
    [yl(2) yl(2)], ...
    'k-', ...
    'LineWidth',0.5);


%% ========================================================================
% 34. TITLE
% ========================================================================

title(ax, ...
    sprintf( ...
    'Trend: 850-hPa Moisture Convergence and Moisture Flux (%d–%d)', ...
    yr_start,yr_end), ...
    'FontSize',10, ...
    'FontWeight','bold');


%% ========================================================================
% 35. COLORBAR
% ========================================================================

cb = colorbar(ax,'southoutside');


cb.FontSize = 10;

cb.TickDirection = 'out';


cb.Ticks = -0.02:0.01:0.02;


cb.TickLabels = ...
    {'-0.02','-0.01','0','0.01','0.02'};


cb.Label.String = ...
    ['Moisture convergence trend ' ...
     '(10^{-6} kg kg^{-1} s^{-1} decade^{-1})'];


cb.Label.FontSize = 10;


%% ========================================================================
% 36. SAVE DATA
% ========================================================================

savefile = fullfile( ...
    data_dir, ...
    'ERA5_850hPa_MFC_MF_SON_trend_2003_2020.mat');


save( ...
    savefile, ...
    'MFC_trend', ...
    'MFC_p', ...
    'QU_trend', ...
    'QV_trend', ...
    'QU_p', ...
    'QV_p', ...
    'VECTOR_sig', ...
    'lon360', ...
    'lat', ...
    'years');


fprintf('\nMAT file saved:\n%s\n',savefile);


%% ========================================================================
% 37. EXPORT FIGURE
%
% exportgraphics
% 600 dpi
% ========================================================================

figfile = fullfile( ...
    data_dir, ...
    'ERA5_850hPa_MFC_MF_SON_trend_2003_2020_0E180E.tif');


exportgraphics( ...
    fig, ...
    figfile, ...
    'Resolution',600, ...
    'BackgroundColor','white');


fprintf('\nFigure saved:\n%s\n',figfile);


%% ========================================================================
% DONE
% ========================================================================

fprintf('\n============================================\n');
fprintf('DONE\n');
fprintf('Trend period = 2003-2020 SON\n');
fprintf('Trend unit = per decade\n');
fprintf('Green dots = significant MFC trend (p < 0.05)\n');
fprintf('Purple vector = QU or QV significant (p < 0.05)\n');
fprintf('Black vector = neither QU nor QV significant\n');
fprintf('Land = masked\n');
fprintf('============================================\n');