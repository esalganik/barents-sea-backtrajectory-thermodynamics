# Barents Sea Back-Trajectory Draft Estimates

This repository contains MATLAB scripts for estimating sea-ice draft along Barents Sea back-trajectories using ERA5 air temperature, SM-LG snow depth, and simple thermodynamic growth models.

The workflow interpolates atmospheric and snow forcing onto M1 and M2 back-trajectories, estimates ice thickness and draft using two approaches, compares the estimates with ULS draft observations, and exports a compact NetCDF product for sharing.

## Repository structure

```text
.
├── scripts/
│   ├── a1_smlg_snow_to_trajectories.m
│   ├── a2_era5_Ta_to_trajectories.m
│   ├── a3_fdd_from_era5.m
│   ├── a4_hi_from_hs_Ta.m
│   ├── a5_plot_FDD_and_resistive_models.m
│   ├── a6_netcdf_export.m
│   └── a7_netcdf_import.m
├── data/
│   ├── raw/
│   │   ├── Trajectories/
│   │   ├── ERA5/
│   │   └── Nansen_Legacy_ULS_data/
│   └── processed/
├── figures/
└── README.md
```

The SM-LG snow-model file is expected separately, for example:

```text
C:\Users\evsalg001\Documents\MATLAB\datasets\SnowModel-LG
```

## Scripts

### `scripts/a1_smlg_snow_to_trajectories.m`

Interpolates SM-LG snow depth onto the M1 and M2 back-trajectories.

Output:

```text
data/processed/SnowDepth_LG_on_backtrajectories.mat
```

### `scripts/a2_era5_Ta_to_trajectories.m`

Interpolates ERA5 2-m air temperature onto the M1 and M2 back-trajectories.

Output:

```text
data/processed/ERA5_T2m_on_backtrajectories.mat
```

### `scripts/a3_fdd_from_era5.m`

Estimates sea-ice thickness from freezing degree days using three approaches: full-trajectory FDD, seasonal FDD, and FDD with monthly ocean heat flux correction.

Output:

```text
data/processed/trajectory_FDD_three_approaches_M1_M2.mat
```

### `scripts/a4_hi_from_hs_Ta.m`

Estimates sea-ice thickness using a thermal-resistance model forced by ERA5 air temperature and SM-LG snow depth.

Output:

```text
data/processed/trajectory_thermal_resistance_SMLG_M1_M2.mat
```

### `scripts/a5_plot_FDD_and_resistive_models.m`

Plots observed ULS draft distributions together with modelled draft estimates from the FDD and thermal-resistance models.

Output:

```text
figures/draft_pdf_M1_M2_with_FDD_and_thermal_draft.png
```

### `scripts/a6_netcdf_export.m`

Creates compact sharable NetCDF files containing the back-trajectories, interpolated forcing, and estimated draft.

Outputs:

```text
data/processed/back_trajectories_M1_with_forcing_and_draft.nc
data/processed/back_trajectories_M2_with_forcing_and_draft.nc
```

### `scripts/a7_netcdf_import.m`

Recreates the final draft-comparison figure using the exported NetCDF files instead of the intermediate MATLAB `.mat` files.

## Required input data

The workflow requires back-trajectory NetCDF files, ERA5 2-m air-temperature files, SM-LG snow-depth data, and Nansen Legacy ULS draft observations.

Example expected input structure:

```text
data/raw/Trajectories/
data/raw/ERA5/
data/raw/Nansen_Legacy_ULS_data/
```

The SM-LG snow-depth file is expected at:

```text
C:\Users\evsalg001\Documents\MATLAB\datasets\SnowModel-LG\SM_snod_MERRA2_ease_01Aug2018-31Jul2021.nc
```

## Running the workflow

Run the scripts in order:

```matlab
scripts/a1_smlg_snow_to_trajectories
scripts/a2_era5_Ta_to_trajectories
scripts/a3_fdd_from_era5
scripts/a4_hi_from_hs_Ta
scripts/a5_plot_FDD_and_resistive_models
scripts/a6_netcdf_export
scripts/a7_netcdf_import
```

## Exported NetCDF variables

The exported NetCDF files include:

```text
t
lat
lon
air_temperature_2m
snow_depth_smlg
draft_fdd_m
draft_therm_m
draft_fdd_arrival_m
draft_therm_arrival_m
```

Row 1 corresponds to the arrival location. Increasing row index follows the back-trajectory backward in time.

## MATLAB requirements

The scripts require MATLAB with NetCDF support. Snow-depth interpolation also requires coordinate projection support through `projcrs` and `projfwd`.

## Notes

The original trajectory NetCDF files are not modified. Derived variables are written to separate processed `.mat` files and compact sharable NetCDF files.

The draft estimates are model-derived and should be interpreted as simple thermodynamic estimates rather than direct observations.
