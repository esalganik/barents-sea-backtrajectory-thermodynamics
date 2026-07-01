%% plot_draft_pdf_with_model_from_exported_netcdf.m
%
% Recreate draft PDF figure using:
%   1. Original ULS draft observations
%   2. Exported back-trajectory NetCDF files containing model draft estimates
%
% The script no longer loads the intermediate FDD or thermal-resistance
% MATLAB .mat files. Model drafts are read directly from the exported
% NetCDF products:
%
%   - back_trajectories_M1_with_forcing_and_draft.nc
%   - back_trajectories_M2_with_forcing_and_draft.nc
%
% Required variables in exported NetCDF files:
%   - arrival_time
%   - draft_fdd_arrival_m
%   - draft_therm_arrival_m

clear; clc; close all

% Paths

projectDir = 'C:\Users\evsalg001\Documents\MATLAB\Sonar Barents Sea';

folder_moor = fullfile(projectDir,'data','raw','Nansen_Legacy_ULS_data');
folder_out  = fullfile(projectDir,'figures');
processedDir = fullfile(projectDir,'data','processed');

modelFile_M1 = fullfile(processedDir, ...
    'back_trajectories_M1_with_forcing_and_draft.nc');

modelFile_M2 = fullfile(processedDir, ...
    'back_trajectories_M2_with_forcing_and_draft.nc');

files_M1 = {
    fullfile(folder_moor,'M1_2018_2019_sig500_812_v1.nc')
    fullfile(folder_moor,'M1_2019_2020_sig500_809_v1.nc')
    };

files_M2 = {
    fullfile(folder_moor,'M2_2018_2019_sig500_802_v1.nc')
    fullfile(folder_moor,'M2_2019_2020_sig500_812_v1.nc')
    fullfile(folder_moor,'M2_2020_2021_sig500_809_v1.nc')
    };

% Plot settings

draft_min = 0;
draft_max = 5;
plot_max_draft = 3;
bin_size = 0.02;

smooth_days_pdf  = 3;
smooth_days_line = 10;

pdf_max = 2.0;
remove_no_ice = true;

cEarly = [0.55 0.82 0.45];
cLate  = [0.95 0.75 0.10];
cMean  = [213 94 0]/255;
cFDD   = [58 174 140]/255;
cTherm = [117 112 179]/255;

% Early/late regimes

M1_early = [
    datetime(2019,1,1)  datetime(2019,4,30)
    datetime(2020,1,1)  datetime(2020,3,31)
    ];

M1_late = [
    datetime(2019,5,1)  datetime(2019,9,1)
    datetime(2020,4,1)  datetime(2020,7,1)
    ];

M2_early = [
    datetime(2019,1,1)  datetime(2019,4,30)
    datetime(2020,1,1)  datetime(2020,3,31)
    datetime(2021,1,1)  datetime(2021,5,31)
    ];

M2_late = [
    datetime(2019,5,1)  datetime(2019,9,1)
    datetime(2020,4,1)  datetime(2020,7,1)
    datetime(2021,6,1)  datetime(2021,8,1)
    ];

% Read model draft estimates from exported NetCDF files

arrival_time_M1 = read_arrival_time_from_t(modelFile_M1);
arrival_time_M2 = read_arrival_time_from_t(modelFile_M2);

draft_fdd_M1 = double(ncread(modelFile_M1,'draft_fdd_arrival_m'));
draft_fdd_M2 = double(ncread(modelFile_M2,'draft_fdd_arrival_m'));

draft_therm_M1 = double(ncread(modelFile_M1,'draft_therm_arrival_m'));
draft_therm_M2 = double(ncread(modelFile_M2,'draft_therm_arrival_m'));

% Read ULS observations

[t1, d1] = read_all(files_M1, draft_min, draft_max, remove_no_ice);
[t2, d2] = read_all(files_M2, draft_min, draft_max, remove_no_ice);

% Make daily draft PDFs

[u1, Z1] = make_pdf(t1, d1, bin_size, plot_max_draft);
[u2, Z2] = make_pdf(t2, d2, bin_size, plot_max_draft);

tmin = min([u1(1), u2(1)]);
tmax = max([u1(end), u2(end)]);

if smooth_days_pdf > 1
    Z1 = movmean(Z1, smooth_days_pdf, 2, 'omitnan');
    Z2 = movmean(Z2, smooth_days_pdf, 2, 'omitnan');
end

mean1 = daily_mean(t1, d1, u1);
mean2 = daily_mean(t2, d2, u2);

if smooth_days_line > 1
    mean1 = movmean(mean1, smooth_days_line, 'omitnan');
    mean2 = movmean(mean2, smooth_days_line, 'omitnan');
end

z = draft_min:bin_size:plot_max_draft;
zc = z(1:end-1) + bin_size/2;

%% Plot

out_png = fullfile(folder_out, ...
    'draft_pdf_M1_M2_with_FDD_and_thermal_draft_from_exported_netcdf.png');

figure('Color','w','Position',[100 100 1200 800])

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% M1 panel

ax1 = nexttile;

imagesc(u1, zc, Z1)
set(gca,'YDir','reverse')
ylabel('Sea-ice draft (m)')
colormap(flipud(gray))

cb = colorbar;
ylabel(cb,'Probability density (m^{-1})')

clim([0 pdf_max])
ylim([draft_min plot_max_draft])
xlim([tmin tmax])
hold on

plot(u1, mean1, ...
    'Color', cMean, ...
    'LineWidth', 1.5)

plot(arrival_time_M1, draft_fdd_M1, '.', ...
    'Color', cFDD, ...
    'MarkerSize', 8)

plot(arrival_time_M1, draft_therm_M1, '.', ...
    'Color', cTherm, ...
    'MarkerSize', 8)

text(datetime(2019,1,10), 2.8, 'M1', ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'Color','k')

add_regime_bars(gca, M1_early, M1_late, cEarly, cLate)

legend('Obs.: Mean','Model: FDD + OHF','Model: thermal resistance + SM-LG snow', ...
    'Location','southeast')

% M2 panel

ax2 = nexttile;

imagesc(u2, zc, Z2)
set(gca,'YDir','reverse')
ylabel('Sea-ice draft (m)')
colormap(flipud(gray))

cb = colorbar;
ylabel(cb,'Probability density (m^{-1})')

clim([0 pdf_max])
ylim([draft_min plot_max_draft])
xlim([tmin tmax])
hold on

plot(u2, mean2, ...
    'Color', cMean, ...
    'LineWidth', 1.5)

plot(arrival_time_M2, draft_fdd_M2, '.', ...
    'Color', cFDD, ...
    'MarkerSize', 8)

plot(arrival_time_M2, draft_therm_M2, '.', ...
    'Color', cTherm, ...
    'MarkerSize', 8)

text(datetime(2019,1,10), 2.8, 'M2', ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'Color','k')

add_regime_bars(gca, M2_early, M2_late, cEarly, cLate)

legend('Obs.: Mean','Model: FDD + OHF','Model: thermal resistance + SM-LG snow', ...
    'Location','southeast')

% Formatting

linkaxes([ax1 ax2],'x')

for ax = [ax1 ax2]
    xticks(ax, datetime(2019,1,1):calmonths(4):datetime(2021,9,1))
    xtickformat(ax,'MMM yyyy')
end

% exportgraphics(gcf, out_png, 'Resolution', 300)
% disp(['Saved figure: ' out_png])

%% Helpers

function t_dt = read_netcdf_time(ncfile,varname)

    t_raw = double(ncread(ncfile,varname));

    try
        units = ncreadatt(ncfile,varname,'units');
    catch
        units = 'days since 1970-01-01 00:00:00';
    end

    units = char(units);

    if contains(units,'days since 1970-01-01')
        t_dt = datetime(1970,1,1) + days(t_raw);
    elseif contains(units,'seconds since 1970-01-01')
        t_dt = datetime(1970,1,1) + seconds(t_raw);
    elseif contains(units,'hours since 1970-01-01')
        t_dt = datetime(1970,1,1) + hours(t_raw);
    else
        error('Unsupported time units in %s: %s', varname, units)
    end

end

function add_regime_bars(ax, early_periods, late_periods, cEarly, cLate)

axes(ax)

yl = ylim;
yr = yl(2) - yl(1);

y1 = yl(1) + 0.005*yr;
y2 = yl(1) + 0.025*yr;
ytxt = yl(1) - 0.010*yr;

for k = 1:size(early_periods,1)
    x1 = early_periods(k,1);
    x2 = early_periods(k,2);

    patch([x1 x2 x2 x1], [y1 y1 y2 y2], cEarly, ...
        'EdgeColor','none', ...
        'FaceAlpha',1, ...
        'HandleVisibility','off')

    text(x1 + (x2-x1)/2, ytxt, sprintf('EARLY %02d', mod(year(x1),100)), ...
        'Color', cEarly, ...
        'FontWeight','bold', ...
        'FontSize',11, ...
        'FontAngle','italic', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'HandleVisibility','off')
end

for k = 1:size(late_periods,1)
    x1 = late_periods(k,1);
    x2 = late_periods(k,2);

    patch([x1 x2 x2 x1], [y1 y1 y2 y2], cLate, ...
        'EdgeColor','none', ...
        'FaceAlpha',1, ...
        'FaceAlpha',1, ...
        'HandleVisibility','off')

    text(x1 + (x2-x1)/2, ytxt, sprintf('LATE %02d', mod(year(x1),100)), ...
        'Color', cLate, ...
        'FontWeight','bold', ...
        'FontSize',12, ...
        'FontAngle','italic', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'HandleVisibility','off')
end

end

function [t_all, d_all] = read_all(files, dmin, dmax, remove_no_ice)

t_all = [];
d_all = [];

for k = 1:numel(files)
    f = files{k};

    time = ncread(f,'TIME');
    d = ncread(f,'SEA_ICE_DRAFT_PING');
    icef = ncread(f,'SEA_ICE_FRACTION');

    t = datetime(1970,1,1) + days(time);

    d = double(d);
    icef = double(icef(:));

    d(d < dmin | d > dmax) = NaN;

    if remove_no_ice
        for i = 1:length(t)
            if icef(i) <= 0 || isnan(icef(i))
                d(i,:) = NaN;
            end
        end
    end

    t_all = [t_all; t(:)];
    d_all = [d_all; d];
end

[t_all, idx] = sort(t_all);
d_all = d_all(idx,:);

[t_all, ia] = unique(t_all, 'stable');
d_all = d_all(ia,:);

end

function [u_day, Z] = make_pdf(t, d, binsize, maxdraft)

t_day = dateshift(t,'start','day');
u_day = unique(t_day);

z = 0:binsize:maxdraft;
Z = nan(length(z)-1, length(u_day));

for i = 1:length(u_day)
    ii = (t_day == u_day(i));
    x = d(ii,:);
    x = x(:);
    x = x(~isnan(x));

    if ~isempty(x)
        Z(:,i) = histcounts(x, z, 'Normalization', 'pdf')';
    end
end

end

function m = daily_mean(t, d, u_day)

t_day = dateshift(t,'start','day');
m = nan(size(u_day));

for i = 1:length(u_day)
    ii = (t_day == u_day(i));
    x = d(ii,:);
    x = x(:);
    x = x(~isnan(x));

    if ~isempty(x)
        m(i) = mean(x,'omitnan');
    end
end

end

function arrival_time = read_arrival_time_from_t(ncfile)

    t = double(ncread(ncfile,'t'));

    % In the trajectory products, row 1 is the arrival point.
    t_arrival = t(1,:);

    try
        units = ncreadatt(ncfile,'t','units');
    catch
        units = 'days since 1970-01-01 00:00:00';
    end

    units = char(units);

    if contains(units,'days since 1970-01-01')
        arrival_time = datetime(1970,1,1) + days(t_arrival);
    elseif contains(units,'seconds since 1970-01-01')
        arrival_time = datetime(1970,1,1) + seconds(t_arrival);
    elseif contains(units,'hours since 1970-01-01')
        arrival_time = datetime(1970,1,1) + hours(t_arrival);
    else
        error('Unsupported time units in t: %s', units)
    end

end