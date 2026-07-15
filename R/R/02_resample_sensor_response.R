# Resample PROSAIL spectra to multispectral sensor bands -----------------

build_gaussian_response <- function(
    wavelengths = 400:2500,
    centers = c(475, 560, 668, 840, 717),
    fwhm = c(20, 20, 10, 40, 10)
) {
  if (length(centers) != length(fwhm)) {
    stop("centers and fwhm must have the same length.")
  }

  response_matrix <- t(vapply(
    seq_along(centers),
    function(i) {
      response <- stats::dnorm(
        wavelengths,
        mean = centers[i],
        sd = fwhm[i] / 2
      )

      range_response <- range(response, na.rm = TRUE)
      if (diff(range_response) == 0) {
        return(rep(0, length(response)))
      }

      (response - range_response[1]) / diff(range_response)
    },
    numeric(length(wavelengths))
  ))

  hsdar::speclib(response_matrix, w = wavelengths)
}

resample_lut_to_sensor <- function(
    lut_spectra,
    wavelengths = 400:2500,
    centers = c(475, 560, 668, 840, 717),
    fwhm = c(20, 20, 10, 40, 10),
    band_names = c("B1", "B2", "B3", "B5", "B4")
) {
  if (ncol(lut_spectra) != length(wavelengths)) {
    stop("Number of LUT spectral columns must match wavelengths.")
  }

  sensor_definition <- data.frame(center = centers)
  sensor_definition$fwhm <- fwhm

  spectral_library <- hsdar::speclib(
    as.matrix(lut_spectra),
    w = as.numeric(wavelengths)
  )

  response <- build_gaussian_response(
    wavelengths = wavelengths,
    centers = centers,
    fwhm = fwhm
  )

  resampled <- hsdar::spectralResampling(
    spectral_library,
    sensor = sensor_definition,
    response_function = response,
    rm.NA = TRUE,
    continuousdata = "auto"
  )

  output <- as.data.frame(resampled@spectra@spectra_ma)
  names(output) <- band_names
  output
}
