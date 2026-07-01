% plot_combined_obs_model_draft_correlations.m
%
% Compare combined M1 + M2 observed ULS draft with thermodynamic model
% draft estimates.
%
% Observations:
%   - daily median draft
%   - daily modal draft
%
% Models:
%   - FDD model
%   - resistive model
%
% Markers are colour-coded by month of year.
% No-ice periods are ignored using SEA_ICE_FRACTION.
% Near-zero model draft values are removed.
% Thick/ridged outliers are excluded from the correlation plots.
% Nothing is saved.

clear; clc; close all

% Paths

projectDir = 'C:\Users\evsalg001\Documents\MATLAB\Sonar Barents Sea';

folder_moor = fullfile(projectDir,'data','raw','Nansen_Legacy_ULS_data');
exportDir   = fullfile(projectDir,'export');

modelFile_M1 = fullfile(exportDir, ...
    'back_trajectories_M1_with_forcing_and_draft.nc');

modelFile_M2 = fullfile(exportDir, ...
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

% Settings

draft_min = 0.20;         % ignore no-ice / very thin values [m]
draft_max = 5.00;         % maximum accepted sonar draft [m]
draft_corr_max = 1.80;    % maximum draft used in correlations [m]
mode_bin = 0.05;          % bin size for modal draft [m]
marker_size = 25;         % scatter marker size

% Read observations

[t_obs_M1, d_obs_M1] = read_all_sonar(files_M1, draft_min, draft_max);
[t_obs_M2, d_obs_M2] = read_all_sonar(files_M2, draft_min, draft_max);

[daily_M1, obs_median_M1, obs_mode_M1] = daily_median_mode(t_obs_M1, d_obs_M1, mode_bin);
[daily_M2, obs_median_M2, obs_mode_M2] = daily_median_mode(t_obs_M2, d_obs_M2, mode_bin);

% Read model arrival drafts

t_mod_M1 = read_arrival_time_from_t(modelFile_M1);
t_mod_M2 = read_arrival_time_from_t(modelFile_M2);

draft_fdd_M1 = double(ncread(modelFile_M1,'draft_fdd_arrival_m'));
draft_fdd_M2 = double(ncread(modelFile_M2,'draft_fdd_arrival_m'));

draft_resistive_M1 = double(ncread(modelFile_M1,'draft_therm_arrival_m'));
draft_resistive_M2 = double(ncread(modelFile_M2,'draft_therm_arrival_m'));

% Convert model estimates to daily medians

[mod_day_M1, fdd_day_M1] = daily_model_median(t_mod_M1, draft_fdd_M1, draft_min);
[~, resistive_day_M1] = daily_model_median(t_mod_M1, draft_resistive_M1, draft_min);

[mod_day_M2, fdd_day_M2] = daily_model_median(t_mod_M2, draft_fdd_M2, draft_min);
[~, resistive_day_M2] = daily_model_median(t_mod_M2, draft_resistive_M2, draft_min);

% Match observations and model estimates separately for M1 and M2

[fdd_obs_median_M1, fdd_mod_median_M1, fdd_time_median_M1] = match_days( ...
    daily_M1, obs_median_M1, mod_day_M1, fdd_day_M1, draft_min);

[fdd_obs_mode_M1, fdd_mod_mode_M1, fdd_time_mode_M1] = match_days( ...
    daily_M1, obs_mode_M1, mod_day_M1, fdd_day_M1, draft_min);

[resistive_obs_median_M1, resistive_mod_median_M1, resistive_time_median_M1] = match_days( ...
    daily_M1, obs_median_M1, mod_day_M1, resistive_day_M1, draft_min);

[resistive_obs_mode_M1, resistive_mod_mode_M1, resistive_time_mode_M1] = match_days( ...
    daily_M1, obs_mode_M1, mod_day_M1, resistive_day_M1, draft_min);

[fdd_obs_median_M2, fdd_mod_median_M2, fdd_time_median_M2] = match_days( ...
    daily_M2, obs_median_M2, mod_day_M2, fdd_day_M2, draft_min);

[fdd_obs_mode_M2, fdd_mod_mode_M2, fdd_time_mode_M2] = match_days( ...
    daily_M2, obs_mode_M2, mod_day_M2, fdd_day_M2, draft_min);

[resistive_obs_median_M2, resistive_mod_median_M2, resistive_time_median_M2] = match_days( ...
    daily_M2, obs_median_M2, mod_day_M2, resistive_day_M2, draft_min);

[resistive_obs_mode_M2, resistive_mod_mode_M2, resistive_time_mode_M2] = match_days( ...
    daily_M2, obs_mode_M2, mod_day_M2, resistive_day_M2, draft_min);

% Combine M1 and M2

fdd_obs_median_all = [fdd_obs_median_M1; fdd_obs_median_M2];
fdd_mod_median_all = [fdd_mod_median_M1; fdd_mod_median_M2];
fdd_time_median_all = [fdd_time_median_M1; fdd_time_median_M2];

fdd_obs_mode_all = [fdd_obs_mode_M1; fdd_obs_mode_M2];
fdd_mod_mode_all = [fdd_mod_mode_M1; fdd_mod_mode_M2];
fdd_time_mode_all = [fdd_time_mode_M1; fdd_time_mode_M2];

resistive_obs_median_all = [resistive_obs_median_M1; resistive_obs_median_M2];
resistive_mod_median_all = [resistive_mod_median_M1; resistive_mod_median_M2];
resistive_time_median_all = [resistive_time_median_M1; resistive_time_median_M2];

resistive_obs_mode_all = [resistive_obs_mode_M1; resistive_obs_mode_M2];
resistive_mod_mode_all = [resistive_mod_mode_M1; resistive_mod_mode_M2];
resistive_time_mode_all = [resistive_time_mode_M1; resistive_time_mode_M2];

% Plot correlations

figure('Color','w','Position',[100 100 1050 850])

tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

nexttile
plot_correlation(fdd_obs_median_all, fdd_mod_median_all, fdd_time_median_all, ...
    'Daily median draft vs FDD model', ...
    draft_min, draft_corr_max, marker_size)

nexttile
plot_correlation(resistive_obs_median_all, resistive_mod_median_all, resistive_time_median_all, ...
    'Daily median draft vs resistive model', ...
    draft_min, draft_corr_max, marker_size)

nexttile
plot_correlation(fdd_obs_mode_all, fdd_mod_mode_all, fdd_time_mode_all, ...
    'Daily modal draft vs FDD model', ...
    draft_min, draft_corr_max, marker_size)

nexttile
plot_correlation(resistive_obs_mode_all, resistive_mod_mode_all, resistive_time_mode_all, ...
    'Daily modal draft vs resistive model', ...
    draft_min, draft_corr_max, marker_size)

sgtitle('Combined M1 and M2: observed draft vs model draft')

outFile = fullfile(projectDir, 'figures', ...
    'combined_M1_M2_obs_model_draft_correlations.png');
exportgraphics(gcf, outFile, 'Resolution', 300)
fprintf('Saved figure: %s\n', outFile)

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Ticks = 1:12;
cb.TickLabels = {'Jan','Feb','Mar','Apr','May','Jun', ...
                 'Jul','Aug','Sep','Oct','Nov','Dec'};
cb.Label.String = 'Month of year';

%% Helpers

function [t_all, d_all] = read_all_sonar(files, dmin, dmax)

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

    for i = 1:length(t)
        if icef(i) <= 0 || isnan(icef(i))
            d(i,:) = NaN;
        end
    end

    t_all = [t_all; t(:)];
    d_all = [d_all; d];

end

[t_all, idx] = sort(t_all);
d_all = d_all(idx,:);

[t_all, ia] = unique(t_all,'stable');
d_all = d_all(ia,:);

end

function [u_day, daily_median, daily_mode] = daily_median_mode(t, d, mode_bin)

t_day = dateshift(t,'start','day');
u_day = unique(t_day);

daily_median = nan(size(u_day));
daily_mode = nan(size(u_day));

for i = 1:numel(u_day)

    use = t_day == u_day(i);

    x = d(use,:);
    x = x(:);
    x = x(~isnan(x));

    if isempty(x)
        continue
    end

    daily_median(i) = median(x,'omitnan');
    daily_mode(i) = binned_mode(x, mode_bin);

end

end

function [u_day, model_daily] = daily_model_median(t, x, draft_min)

t = t(:);
x = x(:);

valid = ~isnat(t) & ~isnan(x) & x >= draft_min;

t = t(valid);
x = x(valid);

t_day = dateshift(t,'start','day');
u_day = unique(t_day);

model_daily = nan(size(u_day));

for i = 1:numel(u_day)

    use = t_day == u_day(i);
    model_daily(i) = median(x(use),'omitnan');

end

end

function [obs_match, mod_match, time_match] = match_days(obs_day, obs_value, mod_day, mod_value, draft_min)

[common_day, ia, ib] = intersect(obs_day, mod_day);

obs_match = obs_value(ia);
mod_match = mod_value(ib);
time_match = common_day;

valid = ~isnan(obs_match) & ~isnan(mod_match) & ...
        obs_match >= draft_min & mod_match >= draft_min;

obs_match = obs_match(valid);
mod_match = mod_match(valid);
time_match = time_match(valid);

end

function m = binned_mode(x, bin_size)

x = x(:);
x = x(~isnan(x));

if isempty(x)
    m = NaN;
    return
end

edges = 0:bin_size:(ceil(max(x)/bin_size)*bin_size + bin_size);
counts = histcounts(x, edges);

if all(counts == 0)
    m = NaN;
    return
end

[~, imax] = max(counts);
m = edges(imax) + bin_size/2;

end

function plot_correlation(obs, model, t, ttl, draft_min, draft_corr_max, marker_size)

obs = obs(:);
model = model(:);
t = t(:);

valid = ~isnan(obs) & ~isnan(model) & ~isnat(t) & ...
        obs >= draft_min & model >= draft_min & ...
        obs <= draft_corr_max & model <= draft_corr_max;

obs = obs(valid);
model = model(valid);
t = t(valid);

month_num = month(t);

scatter(obs, model, marker_size, month_num, ...
    'filled', ...
    'MarkerFaceAlpha', 0.55, ...
    'MarkerEdgeColor', 'none')

hold on
box on

load('roma.mat');
idx = round(linspace(1, size(roma,1), 12));
colormap(roma(idx,:))

clim([0.5 12.5])

xlabel('Observed draft (m)')
ylabel('Model draft (m)')
title(ttl, 'Interpreter','none')

lims = [draft_min draft_corr_max];

plot(lims, lims, 'k--', 'LineWidth', 1.0)

xlim(lims)
ylim(lims)
axis square

if numel(obs) < 3

    text(0.05,0.95,'Not enough points', ...
        'Units','normalized', ...
        'VerticalAlignment','top')

    return

end

p = polyfit(obs, model, 1);
xfit = linspace(lims(1), lims(2), 100);
yfit = polyval(p, xfit);

plot(xfit, yfit, 'r-', 'LineWidth', 1.3)

R_pearson = corr(obs, model, ...
    'Rows','complete', ...
    'Type','Pearson');

RMSE = sqrt(mean((model - obs).^2, 'omitnan'));
bias = mean(model - obs, 'omitnan');

txt = sprintf(['R = %.2f\n' ...
               'RMSE = %.2f m\n' ...
               'Bias = %.2f m\n' ...
               'N = %d'], ...
               R_pearson, RMSE, bias, numel(obs));

text(0.05, 0.95, txt, ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontSize',9, ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.7 0.7 0.7])

end

function arrival_time = read_arrival_time_from_t(ncfile)

t = double(ncread(ncfile,'t'));

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

arrival_time = arrival_time(:);

end