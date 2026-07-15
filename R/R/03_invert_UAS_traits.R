# Invert LAI, Cab, and LIDFa from a multispectral raster -----------------

rmse_vector <- function(simulated, observed) {
  sqrt(mean((simulated - observed)^2, na.rm = TRUE))
}

invert_single_pixel <- function(observed, lut_bands, lut_parameters) {
  if (any(!is.finite(observed))) {
    return(NULL)
  }

  errors <- apply(
    lut_bands,
    MARGIN = 1,
    FUN = rmse_vector,
    observed = observed
  )

  best_index <- which.min(errors)

  result <- lut_parameters[best_index, , drop = FALSE]
  result$RMSE <- errors[best_index]
  result
}

invert_multispectral_raster <- function(
    input_raster,
    lut_bands,
    lut_parameters,
    output_directory,
    output_prefix,
    band_order = c("B1", "B2", "B3", "B5", "B4"),
    n_cores = max(1, parallel::detectCores() - 2),
    chunk_size = 1000
) {
  if (!file.exists(input_raster)) {
    stop("Input raster not found: ", input_raster)
  }

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  raster_stack <- raster::brick(input_raster)

  if (raster::nlayers(raster_stack) != length(band_order)) {
    stop(
      "Input raster has ", raster::nlayers(raster_stack),
      " layers, but ", length(band_order), " bands were expected."
    )
  }

  names(raster_stack) <- band_order

  observed_values <- as.data.frame(raster::values(raster_stack))
  names(observed_values) <- band_order

  valid_rows <- stats::complete.cases(observed_values)
  valid_indices <- which(valid_rows)

  if (length(valid_indices) == 0) {
    stop("No valid pixels were found in the input raster.")
  }

  snowfall::sfInit(parallel = TRUE, cpus = n_cores)
  on.exit(snowfall::sfStop(), add = TRUE)

  snowfall::sfExport(
    list = c(
      "observed_values", "valid_indices",
      "lut_bands", "lut_parameters",
      "invert_single_pixel", "rmse_vector"
    )
  )

  snowfall::sfLibrary(dplyr)

  process_index <- function(index) {
    observed <- as.numeric(observed_values[index, , drop = TRUE])
    result <- invert_single_pixel(observed, lut_bands, lut_parameters)
    if (is.null(result)) {
      return(NULL)
    }
    result$cell_index <- index
    result
  }

  chunks <- split(
    valid_indices,
    ceiling(seq_along(valid_indices) / chunk_size)
  )

  chunk_results <- vector("list", length(chunks))

  for (i in seq_along(chunks)) {
    chunk_results[[i]] <- snowfall::sfLapply(
      chunks[[i]],
      process_index
    )
  }

  inversion <- dplyr::bind_rows(unlist(chunk_results, recursive = FALSE))

  template <- raster::raster(raster_stack)
  output_layers <- list()

  traits <- c(LAI = "lai", Cab = "CHL", LIDFa = "LIDFa", RMSE = "RMSE")

  for (output_name in names(traits)) {
    source_column <- traits[[output_name]]
    output_raster <- template
    values <- rep(NA_real_, raster::ncell(template))
    values[inversion$cell_index] <- inversion[[source_column]]
    raster::values(output_raster) <- values

    output_path <- file.path(
      output_directory,
      paste0(output_prefix, "_", output_name, ".tif")
    )

    raster::writeRaster(
      output_raster,
      filename = output_path,
      options = c("COMPRESS=DEFLATE"),
      overwrite = TRUE
    )

    output_layers[[output_name]] <- output_path
  }

  list(
    outputs = output_layers,
    inversion_table = inversion
  )
}
