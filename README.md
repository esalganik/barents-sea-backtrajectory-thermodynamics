# Barents Sea Back-Trajectory Draft Estimates

MATLAB workflow for estimating sea-ice draft along Barents Sea back-trajectories using ERA5 air temperature, SM-LG snow depth, and simple thermodynamic growth models.

The workflow interpolates forcing onto M1 and M2 back-trajectories, estimates ice thickness and draft, compares the estimates with ULS observations, and exports compact NetCDF files.

## Associated manuscript

This repository was developed to support the manuscript:

**Springtime shifts from thin, locally formed sea ice to thick, imported ice in the northwestern Barents Sea**  
Øyvind Foss, Samuel Brenner, Evgenii Salganik, Jack C. Landy, Arild Sundfjord, Sebastian Gerland, Hiroshi Sumata, and Mats A. Granskog.

The scripts provide supporting analysis for estimating thermodynamic ice growth along sea-ice back-trajectories and comparing these estimates with moored ULS sea-ice draft observations from M1 and M2 in the northwestern Barents Sea.

## Repository structure

```
├── scripts/
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

## Workflow

Run the scripts in order:

| Step | Script                                       | Purpose                                                                |
| ---- | -------------------------------------------- | ---------------------------------------------------------------------- |
| 1    | `scripts/a1_smlg_snow_to_trajectories.m`     | Interpolate SM-LG snow depth onto back-trajectories                    |
| 2    | `scripts/a2_era5_Ta_to_trajectories.m`       | Interpolate ERA5 2-m air temperature onto back-trajectories            |
| 3    | `scripts/a3_fdd_from_era5.m`                 | Estimate ice thickness and draft using the FDD + OHF model             |
| 4    | `scripts/a4_hi_from_hs_Ta.m`                 | Estimate ice thickness and draft using the thermal-resistance model    |
| 5    | `scripts/a5_plot_FDD_and_resistive_models.m` | Compare modelled draft with ULS draft observations                     |
| 6    | `scripts/a6_netcdf_export.m`                 | Export compact sharable NetCDF files                                   |
| 7    | `scripts/a7_netcdf_import.m`                 | Recreate the final plot from exported NetCDF files                     |
| 8    | `scripts/a8_model_vs_uls.m`                  | Compare observed ULS draft with modelled draft using correlation plots |

Main outputs:

```
data/processed/SnowDepth_LG_on_backtrajectories.mat
data/processed/ERA5_T2m_on_backtrajectories.mat
data/processed/trajectory_FDD_three_approaches_M1_M2.mat
data/processed/trajectory_thermal_resistance_SMLG_M1_M2.mat

export/back_trajectories_M1_with_forcing_and_draft.nc
export/back_trajectories_M2_with_forcing_and_draft.nc

figures/draft_pdf_M1_M2_with_FDD_and_thermal_draft.png
figures/combined_M1_M2_obs_model_draft_correlations.png
```

## Thermodynamic models

Two simple thermodynamic approaches are used to estimate sea-ice thickness and draft along the back-trajectories.

### 1. Freezing-degree-day model with ocean heat flux correction

The freezing-degree-day model estimates ice growth from accumulated freezing degree days:

$$
\mathrm{FDD} = \sum \max(T_f - T_a, 0)\Delta t
$$

where `T_f` is the seawater freezing temperature, `T_a` is ERA5 2-m air temperature, and `Δt` is the trajectory time step in days.

Ice thickness is estimated using the empirical relation of Lebedev (1938):

$$
h_i = 1.33 \times \mathrm{FDD}^{0.58} \quad [\mathrm{cm}]
$$

where `FDD` is given in degree-Celsius days. The result is converted from cm to m before further calculations.

Basal melt from ocean heat flux is calculated as:

$$
\Delta h_\mathrm{melt} = \frac{Q_o \Delta t}{\rho_i L_f}
$$

where `Q_o` is the prescribed monthly ocean heat flux, taken from Krishfield et al. (2005; https://doi.org/10.1029/2004JC002293), `ρ_i = 900 kg m^-3` is sea-ice density, and `L_f = 334000 J kg^-1` is the latent heat of fusion.

### 2. Thermal-resistance model with SM-LG snow depth

The thermal-resistance model estimates conductive heat flux through snow and sea ice:

$$
R = \frac{h_i}{k_i} + \frac{h_s}{k_s}
$$

$$
Q_\mathrm{atm} = \frac{T_f - T_a}{R}
$$

where `h_i` is ice thickness, `h_s` is SM-LG snow depth, `k_i = 2.03 W m^-1 K^-1`, and `k_s = 0.31 W m^-1 K^-1`.

Ice growth is calculated as:

$$
\Delta h_\mathrm{growth} = \frac{Q_\mathrm{atm} \Delta t}{\rho_i L_f}
$$

Ice thickness is updated as:

$$
h_i(t+\Delta t) = \max(h_i(t) + \Delta h_\mathrm{growth} - \Delta h_\mathrm{melt}, 0)
$$

The initial ice thickness is `h0 = 0.05 m`, and the minimum conductive thickness is `h_min = 0.02 m`.

### Draft conversion

For both models, ice draft is estimated from ice thickness as:

$$
d = 0.85 h_i
$$

The exported NetCDF files contain full-trajectory draft estimates and arrival-only draft estimates.

## Required input data

Expected raw-data folders:

```
data/raw/Trajectories/
data/raw/ERA5/
data/raw/Nansen_Legacy_ULS_data/
```

The SM-LG snow-depth file is stored separately, for example:

```
C:\Users\evsalg001\Documents\MATLAB\datasets\SnowModel-LG\SM_snod_MERRA2_ease_01Aug2018-31Jul2021.nc
```

## Exported NetCDF variables

The exported files in `export/` include:

```
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

## Requirements

MATLAB with NetCDF support is required. Snow-depth interpolation also requires coordinate projection support through `projcrs` and `projfwd`.

## Notes

The original trajectory NetCDF files are not modified. Intermediate products are written to `data/processed/`, while final sharable NetCDF files are written to `export/`.

The draft estimates are model-derived and should be interpreted as simple thermodynamic estimates rather than direct observations.
