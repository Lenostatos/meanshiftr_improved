# Set up renv for package management
install.packages("remotes")
remotes::install_github("rstudio/renv")
renv::init()

# Install packages needed for development
install.packages(c("roxygen2", "testthat", "knitr"))
remotes::install_github("r-lib/devtools")
remotes::install_github("r-lib/usethis")
install.packages("styler")

# Make usethis functions available on the command line without package prefix.
# usethis::use_usethis()

# Create bare-bones package structure
usethis::create_package(".")

# Make usethis functions available on the command line without package prefix.
# usethis::use_usethis()

# Ignore this file
usethis::use_build_ignore("package_development_code.R")

# Add a license
usethis::use_gpl3_license(name = "Nikolai Knapp")

# Install packages that this package uses
usethis::use_package("assertthat")

usethis::use_rcpp()
usethis::use_r(name = "crownsegmentr-package.R")
# The following stackoverflow post provides valuable info on the Makevars file:
# https://stackoverflow.com/questions/43597632/understanding-the-contents-of-the-makevars-file-in-r-macros-variables-r-ma
