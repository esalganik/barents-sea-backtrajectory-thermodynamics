# Barents Sea Back-Trajectory Draft Estimates

This repository contains MATLAB scripts for estimating sea-ice draft along Barents Sea back-trajectories using ERA5 air temperature, SM-LG snow depth, and simple thermodynamic growth models.

The workflow interpolates atmospheric and snow forcing onto M1 and M2 back-trajectories, estimates ice thickness and draft using two approaches, compares the estimates with ULS draft observations, and exports a compact NetCDF product for sharing.

## Repository structure

```text
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
├── export/
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

## Thermodynamic models

Two simple thermodynamic approaches are used to estimate sea-ice thickness and draft along the back-trajectories.

### 1. Freezing-degree-day model with ocean heat flux correction

The freezing-degree-day (FDD) model estimates ice growth from accumulated freezing degree days:

$$
\mathrm{FDD} = \sum \max(T_f - T_a, 0)\Delta t
$$

where (T_f) is the seawater freezing temperature, (T_a) is ERA5 2-m air temperature, and (\Delta t) is the trajectory time step in days.

Ice thickness is estimated using an empirical power-law relation:

$$
h_i = a , \mathrm{FDD}^{p}
$$

where (a = 1.33) cm ((^\circ\mathrm{C},\mathrm{days})^{-p}) and (p = 0.58). The result is converted from cm to m.

Basal melt from ocean heat flux is calculated as:

$$
\Delta h_\mathrm{melt} = \frac{Q_o \Delta t}{\rho_i L_f}
$$

where (Q_o) is the prescribed monthly ocean heat flux, (\rho_i = 900) kg m(^{-3}) is sea-ice density, and (L_f = 334000) J kg(^{-1}) is the latent heat of fusion.

In the FDD + OHF model, freezing increases the FDD-equivalent ice thickness, while ocean heat flux reduces the thickness by basal melt.

### 2. Thermal-resistance model with SM-LG snow depth

The thermal-resistance model estimates conductive atmospheric heat flux through snow and sea ice:

$$
R = \frac{h_i}{k_i} + \frac{h_s}{k_s}
$$

$$
Q_\mathrm{atm} = \frac{T_f - T_a}{R}
$$

where (h_i) is ice thickness, (h_s) is SM-LG snow depth, (k_i = 2.03) W m(^{-1}) K(^{-1}) is the thermal conductivity of sea ice, and (k_s = 0.31) W m(^{-1}) K(^{-1}) is the thermal conductivity of snow.

When air temperature is below the freezing point, ice growth is calculated as:

$$
\Delta h_\mathrm{growth} = \frac{Q_\mathrm{atm} \Delta t}{\rho_i L_f}
$$

Ocean heat flux melt is calculated as:

$$
\Delta h_\mathrm{melt} = \frac{Q_o \Delta t}{\rho_i L_f}
$$

Ice thickness is updated as:

$$
h_i(t+\Delta t) = \max(h_i(t) + \Delta h_\mathrm{growth} - \Delta h_\mathrm{melt}, 0)
$$

The initial ice thickness is set to (h_0 = 0.05) m, and a minimum conductive thickness of (h_\mathrm{min} = 0.02) m is used to avoid unrealistically small thermal resistance.

### Draft conversion

For both models, ice draft is estimated from ice thickness as:

$$
d = 0.85 h_i
$$

The exported NetCDF files contain full-trajectory draft estimates and arrival-only draft estimates.

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
