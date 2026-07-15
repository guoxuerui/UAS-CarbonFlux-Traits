# Generate a PROSAIL look-up table ---------------------------------------

generate_prosail_lut <- function(
    n_samples = 5000,
    spectral_soil,
    spectral_atmosphere,
    spectral_prospect,
    seed = 42
) {
  set.seed(seed)

  parameter_bounds <- list(
    min = data.frame(
      N = 1.0, CHL = 0, CAR = 0, ANT = 0.56,
      EWT = 0.004, LMA = 0.00166, lai = 0,
      LIDFa = 25, fraction_brown = 0,
      q = 0.1, psoil = 0.5, TypeLidf = 2,
      tts = 30, tto = 0, psi = 0
    ),
    max = data.frame(
      N = 2.5, CHL = 80, CAR = 25.28, ANT = 2.81,
      EWT = 0.034, LMA = 0.0331, lai = 8,
      LIDFa = 78, fraction_brown = 1,
      q = 0.1, psoil = 0.5, TypeLidf = 2,
      tts = 30, tto = 0, psi = 0
    ),
    mean = data.frame(
      N = 1.6, CHL = 44.8, CAR = 8.84, ANT = 1.23,
      EWT = 0.0203, LMA = 0.005575, lai = 3.96,
      LIDFa = 61, fraction_brown = 0.17,
      q = 0.1, psoil = 0.5, TypeLidf = 2,
      tts = 30, tto = 0, psi = 0
    ),
    sd = data.frame(
      N = 0.3, CHL = 12.2, CAR = 5.14, ANT = 0.36,
      EWT = 0.0062, LMA = 0.0008, lai = 2.0,
      LIDFa = 13, fraction_brown = 0.31,
      q = 0, psoil = 0, TypeLidf = 0,
      tts = 0, tto = 0, psi = 0
    )
  )

  input_prosail <- list()

  for (parameter_name in names(parameter_bounds$min)) {
    lower <- parameter_bounds$min[[parameter_name]]
    upper <- parameter_bounds$max[[parameter_name]]
    mean_value <- parameter_bounds$mean[[parameter_name]]
    sd_value <- parameter_bounds$sd[[parameter_name]]

    if (sd_value == 0 || lower == upper) {
      input_prosail[[parameter_name]] <- rep(lower, n_samples)
    } else {
      input_prosail[[parameter_name]] <- truncnorm::rtruncnorm(
        n = n_samples,
        a = lower,
        b = upper,
        mean = mean_value,
        sd = sd_value
      )
    }
  }

  input_prosail <- as.data.frame(input_prosail)
  input_prosail$PROT <- 0
  input_prosail$CBC <- 0
  input_prosail$TypeLidf <- 2
  input_prosail$q <- 0.1
  input_prosail$tts <- 30
  input_prosail$tto <- 0
  input_prosail$psi <- 0
  input_prosail$psoil <- 0.5

  lut_result <- prosail::Generate_LUT_PROSAIL(
    InputPROSAIL = input_prosail,
    SpecPROSPECT = spectral_prospect,
    BandNames = NULL,
    SpecSOIL = spectral_soil,
    SpecATM = spectral_atmosphere,
    SAILversion = "4SAIL",
    BrownLOP = NULL
  )

  list(
    parameters = input_prosail,
    spectra = lut_result$BRF,
    raw = lut_result
  )
}
