%% a4_hi_from_hs_Ta.m
%
% Estimate thermodynamic sea-ice growth along back-trajectories using a
% simple thermal-resistance model.
%
% The script combines ERA5 2-m air temperature along the trajectories with
% SM-LG snow depth and a prescribed monthly ocean heat flux. Ice growth is
% calculated when air temperature is below the seawater freezing point,
% while basal melt is estimated from the ocean heat flux.
%
% Outputs:
%   - trajectory_thermal_resistance_SMLG_M1_M2.mat
%
% The output file contains trajectory variables, model parameters,
% stepwise growth/melt diagnostics, final thermodynamic ice thickness,
% estimated draft, snow and temperature summaries, and a summary table.
%
% Required input files:
%   - ERA5_T2m_on_backtrajectories.mat
%   - SnowDepth_LG_on_backtrajectories.mat
%
% Developed for Barents Sea mooring back-trajectory analysis.

clear; clc; close all

% Paths

processedDir = 'C:\Users\evsalg001\Documents\MATLAB\Sonar Barents Sea\data\processed';

era5File = fullfile(processedDir, ...
    'ERA5_T2m_on_backtrajectories.mat');

snowFile = fullfile(processedDir, ...
    'SnowDepth_LG_on_backtrajectories.mat');

outputFile = fullfile(processedDir, ...
    'trajectory_thermal_resistance_SMLG_M1_M2.mat');

% Load input data

load(era5File, ...
    'Ttraj_C','Ttraj_K','traj_time','lat_all','lon_all','t_all', ...
    'mooring_id','mooring_num','idx_M1','idx_M2','nM1','nM2');

load(snowFile, 'Hs_traj');

% Model parameters

Tf = -1.8;          % seawater freezing temperature [deg C]

rho_i = 900;        % sea-ice density [kg m^-3]
Lf = 334000;        % latent heat of fusion [J kg^-1]

ki = 2.03;          % thermal conductivity of sea ice [W m^-1 K^-1]
ks = 0.31;          % thermal conductivity of snow [W m^-1 K^-1]

h0 = 0.05;          % initial ice thickness at trajectory origin [m]
h_min = 0.02;       % minimum thickness used for conductive resistance [m]

% Monthly ocean heat flux [W m^-2]
ohf_monthly = [ ...
    3.6
    4.1
    3.6
    3.3
    3.4
    6.9
    17.1
    18.7
    7.0
    4.7
    3.7
    3.6];

% Prepare time stepping

arrival_time = traj_time(1,:);
nStep = size(Ttraj_C,1);
nTraj = size(Ttraj_C,2);

dt_days = nan(nStep,nTraj);
dt_days(2:end,:) = abs(days(traj_time(1:end-1,:) - traj_time(2:end,:)));
dt_days(dt_days <= 0) = nan;
dt_days(isnan(Ttraj_C)) = nan;

dt_sec = dt_days * 86400;

% Ocean heat flux and basal melt

month_index = month(traj_time);
OHF = nan(size(Ttraj_C));

for m = 1:12
    OHF(month_index == m) = ohf_monthly(m);
end

OHF(isnan(Ttraj_C)) = nan;

dmelt_ocean_m = OHF .* dt_sec ./ (rho_i * Lf);
dmelt_ocean_m(isnan(Ttraj_C)) = nan;

% Run thermal-resistance model

h_therm_step_m = nan(nStep,nTraj);
growth_atm_step_m = nan(nStep,nTraj);
melt_ocean_step_m = nan(nStep,nTraj);
Q_atm_step = nan(nStep,nTraj);

season_start = NaT(1,nTraj);

for j = 1:nTraj

    if isnat(arrival_time(j))
        continue
    end

    ay = year(arrival_time(j));

    if month(arrival_time(j)) < 8
        season_start(j) = datetime(ay-1,8,1);
    else
        season_start(j) = datetime(ay,8,1);
    end

    h_j = h0;

    idx_valid = find(~isnan(Ttraj_C(:,j)) & ~isnat(traj_time(:,j)),1,'last');

    if isempty(idx_valid)
        continue
    end

    h_therm_step_m(idx_valid,j) = h_j;

    for i = idx_valid:-1:2

        if isnan(Ttraj_C(i,j)) || isnan(dt_sec(i,j))
            continue
        end

        Ta = Ttraj_C(i,j);

        hs = Hs_traj(i,j);
        if isnan(hs) || hs < 0
            hs = 0;
        end

        h_cond = max(h_j,h_min);

        if Ta < Tf
            R = h_cond / ki + hs / ks;
            Q_atm = (Tf - Ta) / R;
            dh_growth = Q_atm * dt_sec(i,j) / (rho_i * Lf);
        else
            Q_atm = 0;
            dh_growth = 0;
        end

        dh_melt = dmelt_ocean_m(i,j);
        if isnan(dh_melt)
            dh_melt = 0;
        end

        h_j = max(h_j + dh_growth - dh_melt,0);

        h_therm_step_m(i-1,j) = h_j;
        growth_atm_step_m(i,j) = dh_growth;
        melt_ocean_step_m(i,j) = dh_melt;
        Q_atm_step(i,j) = Q_atm;
    end
end

% Final thickness and draft estimates

h_therm_m = h_therm_step_m(1,:);
draft_therm_m = 0.85 .* h_therm_m;

% Snow-depth summaries

Hs_arrival = Hs_traj(1,:);
Hs_mean_all = mean(Hs_traj,1,'omitnan');

Hs_mean_lastfreeze = nan(1,nTraj);
Hs_int_lastfreeze = nan(1,nTraj);

for j = 1:nTraj

    tt = traj_time(:,j);
    hs = Hs_traj(:,j);

    use = tt >= season_start(j) & tt <= arrival_time(j) & ...
        ~isnat(tt) & ~isnan(hs);

    if sum(use) > 0
        Hs_mean_lastfreeze(j) = mean(hs(use),'omitnan');
    end

    if sum(use) > 1
        tt_use = tt(use);
        hs_use = hs(use);

        [tt_use,idx] = sort(tt_use);
        hs_use = hs_use(idx);

        Hs_int_lastfreeze(j) = trapz(days(tt_use - tt_use(1)),hs_use);
    end
end

% Temperature and trajectory summaries

Tmean_C = mean(Ttraj_C,1,'omitnan');
Tmedian_C = median(Ttraj_C,1,'omitnan');
Tmin_C = min(Ttraj_C,[],1,'omitnan');
Tmax_C = max(Ttraj_C,[],1,'omitnan');

cold_fraction = mean(Ttraj_C < Tf,1,'omitnan');
traj_duration_days = max(abs(days(traj_time - arrival_time)),[],1,'omitnan');

valid_points = ~isnan(Ttraj_C) & ~isnan(lat_all) & ~isnan(lon_all);
traj_age_points = sum(valid_points,1);

lat_arrival = lat_all(1,:);
lon_arrival = lon_all(1,:);

lat_origin = nan(1,nTraj);
lon_origin = nan(1,nTraj);

for j = 1:nTraj
    idx = find(valid_points(:,j),1,'last');

    if ~isempty(idx)
        lat_origin(j) = lat_all(idx,j);
        lon_origin(j) = lon_all(idx,j);
    end
end

% Create summary table

summary = table();

summary.mooring_id = mooring_id(:);
summary.mooring_num = mooring_num(:);
summary.arrival_time = arrival_time(:);
summary.season_start = season_start(:);

summary.h_therm_m = h_therm_m(:);
summary.draft_therm_m = draft_therm_m(:);

summary.Hs_arrival = Hs_arrival(:);
summary.Hs_mean_all = Hs_mean_all(:);
summary.Hs_mean_lastfreeze = Hs_mean_lastfreeze(:);
summary.Hs_int_lastfreeze = Hs_int_lastfreeze(:);

summary.Tmean_C = Tmean_C(:);
summary.Tmedian_C = Tmedian_C(:);
summary.Tmin_C = Tmin_C(:);
summary.Tmax_C = Tmax_C(:);
summary.cold_fraction = cold_fraction(:);
summary.traj_duration_days = traj_duration_days(:);
summary.traj_age_points = traj_age_points(:);

summary.lat_arrival = lat_arrival(:);
summary.lon_arrival = lon_arrival(:);
summary.lat_origin = lat_origin(:);
summary.lon_origin = lon_origin(:);

% Save output

save(outputFile, ...
    'Ttraj_C','Ttraj_K','Hs_traj','traj_time','arrival_time', ...
    'lat_all','lon_all','t_all', ...
    'mooring_id','mooring_num','idx_M1','idx_M2','nM1','nM2', ...
    'OHF','dmelt_ocean_m', ...
    'h_therm_step_m','h_therm_m','draft_therm_m', ...
    'growth_atm_step_m','melt_ocean_step_m','Q_atm_step', ...
    'Hs_arrival','Hs_mean_all','Hs_mean_lastfreeze','Hs_int_lastfreeze', ...
    'summary','Tf','rho_i','Lf','ki','ks','h0','h_min','ohf_monthly','-v7.3');

fprintf('Saved %s\n', outputFile);

%% Plot thickness time series

figure
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
plot(arrival_time(idx_M1),h_therm_m(idx_M1),'.','MarkerSize',8)
ylabel('Ice thickness (m)')
title('M1: thermal-resistance model with SM-LG snow')
grid on

nexttile
plot(arrival_time(idx_M2),h_therm_m(idx_M2),'.','MarkerSize',8)
xlabel('Observation / arrival time')
ylabel('Ice thickness (m)')
title('M2: thermal-resistance model with SM-LG snow')
grid on

%% Plot draft time series

figure
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
plot(arrival_time(idx_M1),draft_therm_m(idx_M1),'.','MarkerSize',8)
ylabel('Estimated draft (m)')
title('M1')
set(gca,'YDir','reverse')
ylim([0 3])
grid on

nexttile
plot(arrival_time(idx_M2),draft_therm_m(idx_M2),'.','MarkerSize',8)
xlabel('Observation / arrival time')
ylabel('Estimated draft (m)')
title('M2')
set(gca,'YDir','reverse')
ylim([0 3])
grid on

%% Plot snow-depth relation

figure
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
c1 = datenum(arrival_time(idx_M1));
scatter(Hs_mean_lastfreeze(idx_M1),h_therm_m(idx_M1),12,c1,'filled')
xlabel('Mean snow depth, last freezing season (m)')
ylabel('Thermal model ice thickness (m)')
title('M1')
cb = colorbar;
cb.Label.String = 'Arrival time';
datetick(cb,'y','mmm yyyy')
grid on

nexttile
c2 = datenum(arrival_time(idx_M2));
scatter(Hs_mean_lastfreeze(idx_M2),h_therm_m(idx_M2),12,c2,'filled')
xlabel('Mean snow depth, last freezing season (m)')
ylabel('Thermal model ice thickness (m)')
title('M2')
cb = colorbar;
cb.Label.String = 'Arrival time';
datetick(cb,'y','mmm yyyy')
grid on