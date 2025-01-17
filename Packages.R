#######################################################
### Packages
#######################################################

### Install/Load Packages
# Package names
packages <- c("fields","ggplot2", "dplyr", "ggpubr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

install.packages("remotes")
library(remotes)
install_version("RandomFieldsUtils", "1.2.5")
install_version("RandomFields", "3.3.14")
library(RandomFields)