# 00_install_packages.R
# Install all required packages for this project.
#
# Run this script ONCE before running any other scripts.
# Only installs what's missing — safe to re-run.

required_packages <- c(
  "dplyr", "tidyr", "stringr", "lubridate", "purrr",
  "readr", "arrow", "fst", "data.table", "glue",
  "janitor", "jsonlite"
)

missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing)
} else {
  message("All required packages are already installed.")
}

# Verify
check <- vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)
if (all(check)) {
  message("\n=== Setup complete. You can now run 01_cohort_identification.R ===")
} else {
  warning(
    "Failed to install: ",
    paste(required_packages[!check], collapse = ", "),
    "\nPlease install them manually."
  )
}
