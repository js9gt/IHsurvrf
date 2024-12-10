

.onLoad <- function(libname, pkgname) {
  # Automatically load required packages
  library(survival, character.only = TRUE)
  library(tidyr, character.only = TRUE)
  library(dplyr, character.only = TRUE)
}
