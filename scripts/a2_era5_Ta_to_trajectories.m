%% a2_era5_Ta_to_trajectories.m
%
% Interpolate ERA5 2-m air temperature onto Barents Sea back-trajectories.
%
% The script reads M1 and M2 back-trajectory NetCDF files, combines them
% into one trajectory array, loads ERA5 2-m air temperature fields, and
% interpolates ERA5 T2m linearly in space and time onto each trajectory
% position.
%
% Outputs:
%   - ERA5_T2m_on_backtrajectories.mat
%
% The output file contains interpolated air temperature in Kelvin and
% degrees Celsius, trajectory time and position, mooring identifiers, and
% indexing variables for M1 and M2.
%
% Required input files:
%   - back_trajectories_M1_v3_withCS2_v3.0.nc
%   - back_trajectories_M2_v3_withCS2_v3.0.nc
%   - ta2018.nc
%   - ta2019.nc
%   - ta2020.nc
%   - ta2021.nc
%
% Developed for Barents Sea mooring back-trajectory analysis.

clear; clc; close all

% Paths

projectDir = 'C:\Users\evsalg001\Documents\MATLAB\Sonar Barents Sea';

rawDir = fullfile(projectDir,'data','raw');
processedDir = fullfile(projectDir,'data','processed');

trajDir = fullfile(rawDir,'Trajectories');
era5Dir = fullfile(rawDir,'ERA5');

traj_files = {
    fullfile(trajDir,'back_trajectories_M1_v3_withCS2_v3.0.nc')
    fullfile(trajDir,'back_trajectories_M2_v3_withCS2_v3.0.nc')
};

era_files = {
    fullfile(era5Dir,'ta2018.nc')
    fullfile(era5Dir,'ta2019.nc')
    fullfile(era5Dir,'ta2020.nc')
    fullfile(era5Dir,'ta2021.nc')
};

out_file = fullfile(processedDir,'ERA5_T2m_on_backtrajectories.mat');

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

% Convert trajectory longitude to 0-360 degrees for ERA5 interpolation.
lon_interp = lon_all;
lon_interp(lon_interp < 0) = lon_interp(lon_interp < 0) + 360;

% Load ERA5 files

lat_era = [];
lon_era = [];
time_era = [];
t2m_era = [];

for i = 1:length(era_files)

    file = era_files{i};
    fprintf('Loading ERA5 file: %s\n', file);

    lat = double(ncread(file,'latitude'));
    lon = double(ncread(file,'longitude'));
    time = double(ncread(file,'valid_time'));
    t2m = double(ncread(file,'t2m'));

    info = ncinfo(file);
    k = strcmp({info.Variables.Name},'valid_time');
    attrs = info.Variables(k).Attributes;
    units = attrs(strcmp({attrs.Name},'units')).Value;

    this_time = convert_era_time(time,units);

    if isempty(lat_era)
        lat_era = lat;
        lon_era = lon;
    end

    time_era = [time_era; this_time(:)];
    t2m_era = cat(3,t2m_era,t2m);
end

% Sort ERA5 data by time and latitude

[time_era,idx_sort] = sort(time_era);
t2m_era = t2m_era(:,:,idx_sort);

if lat_era(1) > lat_era(end)
    lat_era = flipud(lat_era);
    t2m_era = t2m_era(:,end:-1:1,:);
end

% Interpolate ERA5 T2m onto trajectory positions

Ttraj_K = nan(size(traj_time));

fprintf('Interpolating ERA5 T2m onto back trajectories...\n');

for n = 1:numel(traj_time)

    if isnan(t_all(n)) || isnan(lat_all(n)) || isnan(lon_interp(n))
        continue
    end

    tt = traj_time(n);

    if tt < time_era(1) || tt > time_era(end)
        continue
    end

    i2 = find(time_era >= tt,1,'first');

    if isempty(i2)
        continue

    elseif i2 == 1

        Tmap = t2m_era(:,:,1);
        Ttraj_K(n) = interp2(lon_era,lat_era,Tmap', ...
            lon_interp(n),lat_all(n),'linear');

    else

        i1 = i2 - 1;

        t1 = time_era(i1);
        t2 = time_era(i2);

        if t1 == t2
            w = 0;
        else
            w = seconds(tt - t1) / seconds(t2 - t1);
        end

        Tmap1 = t2m_era(:,:,i1);
        Tmap2 = t2m_era(:,:,i2);

        T1 = interp2(lon_era,lat_era,Tmap1', ...
            lon_interp(n),lat_all(n),'linear');

        T2 = interp2(lon_era,lat_era,Tmap2', ...
            lon_interp(n),lat_all(n),'linear');

        Ttraj_K(n) = (1-w)*T1 + w*T2;
    end
end

Ttraj_C = Ttraj_K - 273.15;

% Save output

save(out_file, ...
    'Ttraj_C','Ttraj_K','traj_time','lat_all','lon_all','lon_interp','t_all', ...
    'mooring_id','mooring_num','idx_M1','idx_M2','nM1','nM2', ...
    'traj_files','era_files','-v7.3');

fprintf('Done. Output saved to %s\n',out_file);
fprintf('M1 trajectories: %d columns\n',nM1);
fprintf('M2 trajectories: %d columns\n',nM2);

%% Plot ERA5 T2m at arrival points

figure
plot(traj_time(1,idx_M1),Ttraj_C(1,idx_M1),'.')
hold on
plot(traj_time(1,idx_M2),Ttraj_C(1,idx_M2),'.')
xlabel('Arrival time')
ylabel('ERA5 T2m at arrival (^oC)')
legend('M1','M2','Location','best')
title('ERA5 T2m interpolated to trajectory arrival points')
grid on

%% Helper functions

function era_time = convert_era_time(time,units)

    units = char(units);

    if contains(units,'seconds since')
        ref_txt = strtrim(extractAfter(units,'seconds since'));
        ref_date = parse_ref_date(ref_txt);
        era_time = ref_date + seconds(time);

    elseif contains(units,'hours since')
        ref_txt = strtrim(extractAfter(units,'hours since'));
        ref_date = parse_ref_date(ref_txt);
        era_time = ref_date + hours(time);

    elseif contains(units,'days since')
        ref_txt = strtrim(extractAfter(units,'days since'));
        ref_date = parse_ref_date(ref_txt);
        era_time = ref_date + days(time);

    else
        error('Unknown ERA5 time units: %s',units);
    end
end

function ref_date = parse_ref_date(ref_txt)

    ref_txt = char(ref_txt);

    fmts = {'yyyy-MM-dd HH:mm:ss','yyyy-MM-dd HH:mm','yyyy-MM-dd'};

    for i = 1:length(fmts)
        try
            ref_date = datetime(ref_txt,'InputFormat',fmts{i});
            return
        catch
        end
    end

    error('Could not parse reference date: %s',ref_txt);
end