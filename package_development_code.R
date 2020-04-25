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
usethis::use_gpl3_license(name = "Leon Steinmeier")

# Install packages that this package uses
usethis::use_package("assertthat")

install.packages("dbscan")
usethis::use_package("dbscan")

install.packages("data.table")
usethis::use_package("data.table")

usethis::use_package("BH", type = "LinkingTo")

usethis::use_rcpp()
usethis::use_r(name = "crownsegmentr-package.R")
# The following stackoverflow post provides valuable info on the Makevars file:
# https://stackoverflow.com/questions/43597632/understanding-the-contents-of-the-makevars-file-in-r-macros-variables-r-ma

devtools::document()
# Or use keyboard shortcut Ctrl + Shift + D

usethis::use_r(name = "try_to_extract_coordinate_data")
# Execute usethis::use_test() interactively when the newly created R file is
# active in order to create an R file with tests.

usethis::use_r(name = "segment_tree_crowns")
