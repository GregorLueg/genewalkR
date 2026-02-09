.onLoad <- function(...) {
  # S7
  S7::methods_register()

  # Create torch cache directories
  cache_dir <- file.path(tempdir(), ".torch")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(TORCH_HOME = cache_dir)
}
