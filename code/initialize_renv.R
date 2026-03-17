# Setup R environment using renv
# Install renv if not already installed:
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Initialize renv for the project:
renv::init()

# All required packages (deduped)
pkgs <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "lubridate",
  "purrr",
  "readr",
  "arrow",
  "fst",
  "data.table",
  "glue",
  "janitor"
)
pkgs <- unique(pkgs)

# Install and lock
renv::install(pkgs)
renv::settings$snapshot.type("all")
renv::snapshot(prompt = FALSE) # Then run snapshot again
