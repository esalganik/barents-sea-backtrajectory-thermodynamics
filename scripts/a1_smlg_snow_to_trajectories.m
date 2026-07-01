%% a1_smlg_snow_to_trajectories.m
%
% Interpolate SM-LG snow depth onto Barents Sea back-trajectories.
%
% The script reads M1 and M2 back-trajectory NetCDF files, combines them
% into one trajectory array, loads gridded SM-LG snow-depth data, projects
% trajectory positions to the EASE-grid projection, and interpolates snow
% depth linearly in space and time onto each trajectory position.
%
% Outputs:
%   - SnowDepth_LG_on_backtrajectories.mat
%
% The output file contains snow depth along the trajectories, snow-depth
% summary metrics, trajectory time and position, projected coordinates,
% mooring identifiers, and indexing variables for M1 and M2.
%
% Required input files:
%   - back_trajectories_M1_v3_withCS2_v3.0.nc
%   - back_trajectories_M2_v3_withCS2_v3.0.nc
%   - SM_snod_MERRA2_ease_01Aug2018-31Jul2021.nc
%
% Developed for Barents Sea mooring back-trajectory analysis.

clear; clc; close all

% Paths

projectDir = 'C:\Users\evsalg001\Documents\MATLAB\Sonar Barents Sea';

rawDir = fullfile(projectDir,'data','raw');
processedDir = fullfile(projectDir,'data','processed');

trajDir = fullfile(rawDir,'Trajectories');
snowDir = 'C:\Users\evsalg001\Documents\MATLAB\datasets\SnowModel-LG';

traj_files = {
    fullfile(trajDir,'back_trajectories_M1_v3_withCS2_v3.0.nc')
    fullfile(trajDir,'back_trajectories_M2_v3_withCS2_v3.0.nc')
};

snow_file = fullfile(snowDir,'SM_snod_MERRA2_ease_01Aug2018-31Jul2021.nc');

out_file = fullfile(processedDir,'SnowDepth_LG_on_backtrajectories.mat');

% Load trajectory files

vars = {'t','lat','lon'};

data = struct();
Tmax = 0;

for i = 1:length(traj_files)

    file = traj_files{i};
    fprintf('Loading trajectory file: %s\n', file);

    for v = 1:length(vars)
        name = vars{v};
        data(i).(name) = double(ncread(file,name));
    end

    Tmax = max(Tmax,size(data(i).t,1));
end

% Pad trajectory arrays to common length

for i = 1:length(traj_files)

    for v = 1:length(vars)

        name = vars{v};
        A = data(i).(name);

        if isvector(A)
            A = A(:);
        end

        [Ti,Ni] = size(A);

        if Ti < Tmax
            A = [A; nan(Tmax-Ti,Ni)];
        end

        data(i).(name) = A;
    end
end

% Combine M1 and M2 trajectories

nM1 = size(data(1).t,2);
nM2 = size(data(2).t,2);

idx_M1 = 1:nM1;
idx_M2 = nM1+1:nM1+nM2;

mooring_id = [repmat("M1",1,nM1), repmat("M2",1,nM2)];
mooring_num = [ones(1,nM1), 2*ones(1,nM2)];

t_all   = [data(1).t,   data(2).t];
lat_all = [data(1).lat, data(2).lat];
lon_all = [data(1).lon, data(2).lon];

traj_time = datetime(1970,1,1) + days(t_all);

% Load SM-LG snow grid and time

fprintf('Loading snow-model grid and time...\n')

x = double(ncread(snow_file,'x'));
y = double(ncread(snow_file,'y'));

time_raw = double(ncread(snow_file,'time'));

units = '';
try
    units = char(ncreadatt(snow_file,'time','units'));
catch
end

time_snow_all = convert_model_time(time_raw,units);

[time_snow_all,idx_sort] = sort(time_snow_all);

% Project trajectory coordinates to EASE grid

lon_for_proj = lon_all;
lon_for_proj(lon_for_proj > 180) = lon_for_proj(lon_for_proj > 180) - 360;

[x_traj,y_traj] = projfwd(projcrs(3408),lat_all,lon_for_proj);

% Interpolate snow depth onto trajectory positions

Hs_traj = nan(size(traj_time));

fprintf('Interpolating snow depth onto back trajectories...\n')

for it = 2:length(time_snow_all)

    t1 = time_snow_all(it-1);
    t2 = time_snow_all(it);

    mask = traj_time >= t1 & traj_time <= t2 & ...
        ~isnan(x_traj) & ~isnan(y_traj) & ~isnan(t_all);

    if ~any(mask(:))
        continue
    end

    i1 = idx_sort(it-1);
    i2 = idx_sort(it);

    S1 = double(ncread(snow_file,'snod',[1 1 i1],[Inf Inf 1]));
    S2 = double(ncread(snow_file,'snod',[1 1 i2],[Inf Inf 1]));

    S1 = squeeze(S1);
    S2 = squeeze(S2);

    F1 = griddedInterpolant({x,y},S1,'linear','none');
    F2 = griddedInterpolant({x,y},S2,'linear','none');

    h1 = F1(x_traj(mask),y_traj(mask));
    h2 = F2(x_traj(mask),y_traj(mask));

    w = seconds(traj_time(mask) - t1) ./ seconds(t2 - t1);

    Hs_traj(mask) = (1-w).*h1 + w.*h2;
end

% Calculate snow-depth metrics

fprintf('Calculating snow-depth metrics...\n')

Hs_arrival = Hs_traj(1,:);
Hs_mean_all = mean(Hs_traj,1,'omitnan');

Hs_mean_lastfreeze = nan(1,size(Hs_traj,2));
Hs_int_lastfreeze = nan(1,size(Hs_traj,2));

for j = 1:size(Hs_traj,2)

    t_arr = traj_time(1,j);

    if isnat(t_arr) || isnan(t_all(1,j))
        continue
    end

    if month(t_arr) >= 9
        freeze_start = datetime(year(t_arr),9,1);
    else
        freeze_start = datetime(year(t_arr)-1,9,1);
    end

    freeze_end = t_arr;

    tt = traj_time(:,j);
    hs = Hs_traj(:,j);

    mask = tt >= freeze_start & tt <= freeze_end & ...
        ~isnat(tt) & ~isnan(hs);

    if sum(mask) > 0
        Hs_mean_lastfreeze(j) = mean(hs(mask),'omitnan');
    end

    if sum(mask) > 1
        tt_use = tt(mask);
        hs_use = hs(mask);

        [tt_use,idx] = sort(tt_use);
        hs_use = hs_use(idx);

        Hs_int_lastfreeze(j) = trapz(days(tt_use - tt_use(1)),hs_use);
    end
end

% Save output

save(out_file, ...
    'Hs_traj','Hs_arrival','Hs_mean_all','Hs_mean_lastfreeze','Hs_int_lastfreeze', ...
    'traj_time','lat_all','lon_all','x_traj','y_traj','t_all', ...
    'mooring_id','mooring_num','idx_M1','idx_M2','nM1','nM2', ...
    'traj_files','snow_file','-v7.3');

fprintf('Done. Output saved to %s\n',out_file);
fprintf('M1 trajectories: %d columns\n',nM1);
fprintf('M2 trajectories: %d columns\n',nM2);

%% Plot snow depth at arrival points

figure('Units','inches','Position',[1 5 10 3],'Color','w')
plot(traj_time(1,idx_M1),Hs_arrival(idx_M1),'.')
hold on
plot(traj_time(1,idx_M2),Hs_arrival(idx_M2),'.')
xlabel('Arrival time')
ylabel('Snow depth at arrival (m)')
legend('M1','M2','Location','best')
title('MERRA2 SM-LG snow depth at trajectory arrival points')
xlim([datetime(2018,9,1) datetime(2021,10,1)])
xticks(datetime(2019,1,1):calmonths(4):datetime(2021,9,1))
xtickformat('MMM yyyy')
grid on

%% Plot mean snow depth along full trajectories

figure('Units','inches','Position',[1 5 10 3],'Color','w')
plot(traj_time(1,idx_M1),Hs_mean_all(idx_M1),'.')
hold on
plot(traj_time(1,idx_M2),Hs_mean_all(idx_M2),'.')
xlabel('Arrival time')
ylabel('Mean snow depth along trajectory (m)')
legend('M1','M2','Location','best')
title('Mean MERRA2 SM-LG snow depth along full back trajectories')
xlim([datetime(2018,9,1) datetime(2021,10,1)])
xticks(datetime(2019,1,1):calmonths(4):datetime(2021,9,1))
xtickformat('MMM yyyy')
grid on

%% Plot mean snow depth during last freezing season

figure('Units','inches','Position',[1 5 10 3],'Color','w')
plot(traj_time(1,idx_M1),Hs_mean_lastfreeze(idx_M1),'.')
hold on
plot(traj_time(1,idx_M2),Hs_mean_lastfreeze(idx_M2),'.')
xlabel('Arrival time')
ylabel('Mean snow depth, last freezing season (m)')
legend('M1','M2','Location','best')
title('Mean MERRA2 SM-LG snow depth during last freezing season')
xlim([datetime(2018,9,1) datetime(2021,10,1)])
xticks(datetime(2019,1,1):calmonths(4):datetime(2021,9,1))
xtickformat('MMM yyyy')
grid on

%% Helper functions

function model_time = convert_model_time(time_raw,units)

    units = char(units);

    tok = regexp(units, ...
        '(\w+)\s+since\s+(\d{1,4})-(\d{1,2})-(\d{1,2})(?:[ T](\d{1,2}:\d{2}:\d{2}))?', ...
        'tokens','once');

    if isempty(tok)
        error('Could not parse time units: %s',units)
    end

    base_unit = lower(tok{1});
    yyyy = str2double(tok{2});
    mm   = str2double(tok{3});
    dd   = str2double(tok{4});

    if numel(tok) >= 5 && ~isempty(tok{5})
        base_time = tok{5};
    else
        base_time = '00:00:00';
    end

    hhmmss = sscanf(base_time,'%d:%d:%d');

    t0 = datetime(yyyy,mm,dd,hhmmss(1),hhmmss(2),hhmmss(3));

    switch base_unit
        case {'second','seconds'}
            model_time = t0 + seconds(time_raw);
        case {'hour','hours'}
            model_time = t0 + hours(time_raw);
        case {'day','days'}
            model_time = t0 + days(time_raw);
        otherwise
            error('Unsupported time unit: %s',base_unit)
    end
end