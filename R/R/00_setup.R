# Project setup -----------------------------------------------------------

required_packages <- c(
  "prosail", "hsdar", "raster", "dplyr", "snowfall",
  "truncnorm", "yaml"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the workflow."
  )
}

options(stringsAsFactors = FALSE)

# Use a reproducible random seed for LUT generation.
set.seed(42)

# Resolve project-relative paths.
project_path <- function(...) {
  file.path(normalizePath(".", winslash = "/", mustWork = FALSE), ...)
}

# Load project configuration.
load_config <- function(path = project_path("config", "config.yml")) {
  if (!file.exists(path)) {
    stop("Configuration file not found: ", path)
  }
  yaml::read_yaml(path)
}
