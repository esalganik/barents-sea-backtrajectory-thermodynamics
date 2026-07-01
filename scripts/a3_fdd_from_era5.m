%% a3_fdd_from_era5.m
%
% Estimate sea-ice thickness along back-trajectories using freezing degree
% days (FDD) and a simple ocean-heat-flux correction.
%
% The script applies three approaches:
%   1. Full-trajectory FDD thickness estimate.
%   2. Seasonal FDD thickness estimate, starting from 1 August of the
%      relevant freezing season.
%   3. FDD thickness estimate with monthly ocean heat flux (OHF) used to
%      reduce ice thickness by basal melt.
%
% Outputs:
%   - trajectory_FDD_three_approaches_M1_M2.mat
%
% The output file contains the original trajectory variables, FDD fields,
% OHF-derived melt, estimated ice thickness and draft, trajectory summary
% statistics, and a summary table.
%
% Required input file:
%   - ERA5_T2m_on_backtrajectories.mat
%
% Developed for Barents Sea mooring back-trajectory analysis.

clear; clc; close all

% Paths

processedDir = 'C:\Users\evsalg001\Documents\MATLAB\Sonar Barents Sea\data\processed';

era5File = fullfile(processedDir, ...
    'ERA5_T2m_on_backtrajectories.mat');

outputFile = fullfile(processedDir, ...
    'trajectory_FDD_three_approaches_M1_M2.mat');

% Load input data

load(era5File, ...
    'Ttraj_C','Ttraj_K','traj_time','lat_all','lon_all','t_all', ...
    'mooring_id','mooring_num','idx_M1','idx_M2','nM1','nM2');

% Model parameters

Tf   = -1.8;       % seawater freezing temperature [deg C]
a_cm = 1.33;       % FDD coefficient [cm degC^-p days^-p]
p    = 0.58;       % FDD exponent [-]

rho_i = 900;       % sea-ice density [kg m^-3]
Lf    = 334000;    % latent heat of fusion [J kg^-1]

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

dt_days = abs(days(diff(traj_time,1,1)));
dt_first = median(dt_days,1,'omitnan');
dt_days = [dt_first; dt_days];

dt_days(dt_days <= 0) = nan;
dt_days(isnan(Ttraj_C)) = nan;

dt_sec = dt_days * 86400;

% Freezing degree days

dFDD = max(Tf - Ttraj_C,0) .* dt_days;
dFDD(isnan(Ttraj_C)) = nan;

% Ocean heat flux and basal melt

month_index = month(traj_time);
OHF = nan(size(Ttraj_C));

for m = 1:12
    OHF(month_index == m) = ohf_monthly(m);
end

OHF(isnan(Ttraj_C)) = nan;

dmelt_ocean_m = OHF .* dt_sec ./ (rho_i * Lf);
dmelt_ocean_m(isnan(Ttraj_C)) = nan;

% Approach 1: full-trajectory FDD

FDD_full = sum(dFDD,1,'omitnan');
h_full_m = a_cm .* (FDD_full.^p) ./ 100;

% Approach 2: seasonal FDD since 1 August

FDD_season = nan(1,nTraj);
h_season_m = nan(1,nTraj);
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

    use = traj_time(:,j) >= season_start(j) & ...
          traj_time(:,j) <= arrival_time(j) & ...
          ~isnan(dFDD(:,j));

    FDD_season(j) = sum(dFDD(use,j),'omitnan');
    h_season_m(j) = a_cm .* (FDD_season(j).^p) ./ 100;
end

% Approach 3: FDD with ocean-heat-flux melt correction

h_ohf_step_m = nan(nStep,nTraj);
FDD_ohf_equiv = nan(nStep,nTraj);

for j = 1:nTraj

    FDD_j = 0;
    h_j = 0;

    for i = nStep:-1:1

        if isnan(Ttraj_C(i,j)) || isnan(dt_days(i,j))
            h_ohf_step_m(i,j) = nan;
            FDD_ohf_equiv(i,j) = nan;
            continue
        end

        if Ttraj_C(i,j) < Tf
            FDD_j = FDD_j + dFDD(i,j);
        end

        h_growth_m = a_cm .* (FDD_j.^p) ./ 100;

        if h_growth_m > h_j
            h_j = h_growth_m;
        end

        melt_i = dmelt_ocean_m(i,j);

        if isnan(melt_i)
            melt_i = 0;
        end

        h_j = max(h_j - melt_i,0);

        if h_j <= 0
            FDD_j = 0;
        else
            FDD_j = ((h_j * 100) / a_cm)^(1/p);
        end

        h_ohf_step_m(i,j) = h_j;
        FDD_ohf_equiv(i,j) = FDD_j;
    end
end

h_ohf_m = h_ohf_step_m(1,:);
FDD_ohf_equiv_final = FDD_ohf_equiv(1,:);

draft_ohf_m = 0.85 .* h_ohf_m;

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

summary.FDD_full = FDD_full(:);
summary.h_full_m = h_full_m(:);

summary.FDD_season = FDD_season(:);
summary.h_season_m = h_season_m(:);

summary.FDD_ohf_equiv = FDD_ohf_equiv_final(:);
summary.h_ohf_m = h_ohf_m(:);
summary.draft_ohf_m = draft_ohf_m(:);

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
    'Ttraj_C','Ttraj_K','traj_time','arrival_time', ...
    'lat_all','lon_all','t_all', ...
    'mooring_id','mooring_num','idx_M1','idx_M2','nM1','nM2', ...
    'dFDD','FDD_full','h_full_m', ...
    'FDD_season','h_season_m','season_start', ...
    'OHF','dmelt_ocean_m', ...
    'h_ohf_step_m','h_ohf_m','draft_ohf_m', ...
    'FDD_ohf_equiv','FDD_ohf_equiv_final', ...
    'summary','Tf','a_cm','p','rho_i','Lf','ohf_monthly','-v7.3');

fprintf('Saved %s\n', outputFile);

%% Plot thickness estimates

figure
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
plot(arrival_time(idx_M1),h_full_m(idx_M1),'.','MarkerSize',8)
hold on
plot(arrival_time(idx_M1),h_season_m(idx_M1),'.','MarkerSize',8)
plot(arrival_time(idx_M1),h_ohf_m(idx_M1),'.','MarkerSize',8)
ylabel('Ice thickness (m)')
title('M1')
legend('FDD','Seasonal FDD','FDD + OHF','Location','best')
grid on

nexttile
plot(arrival_time(idx_M2),h_full_m(idx_M2),'.','MarkerSize',8)
hold on
plot(arrival_time(idx_M2),h_season_m(idx_M2),'.','MarkerSize',8)
plot(arrival_time(idx_M2),h_ohf_m(idx_M2),'.','MarkerSize',8)
xlabel('Observation / arrival time')
ylabel('Ice thickness (m)')
title('M2')
legend('FDD','Seasonal FDD','FDD + OHF','Location','best')
grid on

%% Plot FDD estimates

figure
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
plot(arrival_time(idx_M1),FDD_full(idx_M1),'.','MarkerSize',8)
hold on
plot(arrival_time(idx_M1),FDD_season(idx_M1),'.','MarkerSize',8)
plot(arrival_time(idx_M1),FDD_ohf_equiv_final(idx_M1),'.','MarkerSize',8)
ylabel('FDD (^oC days)')
title('M1')
legend('FDD','Seasonal FDD','FDD + OHF','Location','best')
grid on

nexttile
plot(arrival_time(idx_M2),FDD_full(idx_M2),'.','MarkerSize',8)
hold on
plot(arrival_time(idx_M2),FDD_season(idx_M2),'.','MarkerSize',8)
plot(arrival_time(idx_M2),FDD_ohf_equiv_final(idx_M2),'.','MarkerSize',8)
xlabel('Observation / arrival time')
ylabel('FDD (^oC days)')
title('M2')
legend('FDD','Seasonal FDD','FDD + OHF','Location','best')
grid on