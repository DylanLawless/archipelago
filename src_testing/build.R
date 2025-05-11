library(devtools)
setwd("~/web/archipelago")

# Update documentation
devtools::document()

# Run tests
devtools::test()

# Check the package with CRAN-level strictness
devtools::check()

# Optionally build and check the site (if using pkgdown)
if (requireNamespace("pkgdown", quietly = TRUE)) {
  pkgdown::build_site()
}

# Build CRAN tarball (goes to parent dir)
devtools::build(path = "..")

# Then check the tarball manually:

ls -lh ../archipelago_0.1.0.tar.gz
R CMD check ../archipelago_0.1.0.tar.gz --as-cran

# If that’s all clean and you’re ready to submit:

# From inside the package directory
setwd("~/web/archipelago")
devtools::release()


