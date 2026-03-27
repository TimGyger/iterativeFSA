#######################################################
### Packages
#######################################################

# Package names
packages <- c("fields", "ggplot2", "dplyr", "ggpubr", "gpboost", "here")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Ensure remotes is installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)

# Relative paths inside the project
here::i_am("Packages.R")
rfu_path <- here::here("external", "RandomFieldsUtils")
rf_path  <- here::here("external", "RandomFields")

patch_rfu <- here::here("Patches", "RandomFieldsUtils.patch")
patch_rf  <- here::here("Patches", "RandomFields.patch")

# Checks
stopifnot(dir.exists(rfu_path))
stopifnot(dir.exists(rf_path))
stopifnot(file.exists(patch_rfu))
stopifnot(file.exists(patch_rf))

# Apply patch only if needed
apply_patch_if_needed <- function(repo_path, patch_file) {
  check <- system2("git", c("-C", repo_path, "apply", "--check", patch_file))
  if (check == 0) {
    status <- system2("git", c("-C", repo_path, "apply", patch_file))
    if (status != 0) {
      stop("Applying patch failed for: ", repo_path)
    }
  } else {
    message("Patch not applied for ", repo_path,
            " (possibly already applied or repo version does not match patch).")
  }
}

# Install package from local source directory
install_local_if_needed <- function(pkg, path, force = FALSE) {
  needs_install <- !requireNamespace(pkg, quietly = TRUE) || force
  if (needs_install) {
    remotes::install_local(
      path = path,
      upgrade = "never",
      dependencies = FALSE,
      force = TRUE
    )
  }
}

apply_patch_if_needed(rfu_path, patch_rfu)
apply_patch_if_needed(rf_path, patch_rf)

# Install in the correct order
install_local_if_needed("RandomFieldsUtils", rfu_path)
install_local_if_needed("RandomFields", rf_path)

# Load
library(RandomFieldsUtils)
library(RandomFields)
