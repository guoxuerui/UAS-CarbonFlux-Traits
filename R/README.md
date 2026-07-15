# Modular PROSAIL inversion scripts

This folder contains a cleaned and modularized version of the original
`LAI and Cab inversion` R Markdown workflow.

## Structure

```text
R/
  00_setup.R
  01_generate_PROSAIL_LUT.R
  02_resample_sensor_response.R
  03_invert_UAS_traits.R
  04_run_all_dates.R

config/
  config.yml

notebooks/
  PROSAIL_inversion_workflow.Rmd
```

## Main improvements

- Removed hard-coded Windows working directories.
- Replaced five repeated campaign blocks with one configuration-driven loop.
- Separated LUT generation, sensor resampling, and raster inversion.
- Added argument checks and clearer function names.
- Preserved the original retrieval targets: LAI, Cab, and LIDFa.
- Added RMSE as an output raster and an inversion table as an RDS file.

## Required packages

```r
install.packages(c(
  "hsdar",
  "raster",
  "dplyr",
  "snowfall",
  "truncnorm",
  "yaml"
))
```

Install `prosail` according to the package author's instructions.

## Before running

1. Place the required spectral helper files in `R/`:
   - `SpecPROSPECT_FullRange-data.R`
   - `SpecSOIL-data.R`
   - `SpecATM-data.R`

2. Edit `config/config.yml`.

3. Ensure every input raster has five reflectance bands in the configured
   order.

4. Run from the repository root:

```r
source("R/04_run_all_dates.R")
```

## Important scientific note

The sensor center wavelengths, FWHM values, PROSAIL parameter distributions,
and viewing/illumination geometry are study-specific assumptions. Verify them
before applying the workflow to another sensor, site, or campaign.
