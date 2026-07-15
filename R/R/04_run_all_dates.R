# Run trait inversion for all campaign dates -----------------------------

source(file.path("R", "00_setup.R"))
source(file.path("R", "01_generate_PROSAIL_LUT.R"))
source(file.path("R", "02_resample_sensor_response.R"))
source(file.path("R", "03_invert_UAS_traits.R"))

config <- load_config()

# These objects must be supplied by the repository or loaded here.
# Expected objects:
#   SpecPROSPECT_FullRange
#   SpecSOIL
#   SpecATM
source(config$inputs$spectral_prospect_r)
source(config$inputs$spectral_soil_r)
source(config$inputs$spectral_atmosphere_r)

lut <- generate_prosail_lut(
  n_samples = config$prosail$n_samples,
  spectral_soil = SpecSOIL,
  spectral_atmosphere = SpecATM,
  spectral_prospect = SpecPROSPECT_FullRange,
  seed = config$prosail$seed
)

lut_bands <- resample_lut_to_sensor(
  lut_spectra = lut$spectra,
  wavelengths = seq(400, 2500),
  centers = unlist(config$sensor$centers_nm),
  fwhm = unlist(config$sensor$fwhm_nm),
  band_names = unlist(config$sensor$band_names)
)

for (campaign in config$campaigns) {
  message("Processing campaign: ", campaign$date)

  result <- invert_multispectral_raster(
    input_raster = campaign$input_raster,
    lut_bands = lut_bands,
    lut_parameters = lut$parameters,
    output_directory = campaign$output_directory,
    output_prefix = paste0("Mica_", campaign$date, "_1m"),
    band_order = unlist(config$sensor$band_names),
    n_cores = config$processing$n_cores,
    chunk_size = config$processing$chunk_size
  )

  saveRDS(
    result$inversion_table,
    file.path(
      campaign$output_directory,
      paste0("inversion_table_", campaign$date, ".rds")
    )
  )
}
