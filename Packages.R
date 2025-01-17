#######################################################
### Packages
#######################################################

### Install/Load Packages
# Package names
packages <- c("fields","ggplot2", "dplyr", "ggpubr", "gpboost")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Random Fields
# Function to check and install a package if not already installed
install_if_missing <- function(package, version = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (!is.null(version)) {
      remotes::install_version(package, version)
    } else {
      install.packages(package)
    }
  }
}

# Ensure remotes is installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Load remotes library
library(remotes)

# Install specific versions if not yet installed
install_if_missing("RandomFieldsUtils", "1.2.5")
install_if_missing("RandomFields", "3.3.14")

# Load the required libraries
library(RandomFields)
library(RandomFieldsUtils)