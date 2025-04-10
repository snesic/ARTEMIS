#Source python and cython functions on load
#' @export
.onLoad <- function(libname, pkgname) {
  
  # Python environment details
  py_config <- reticulate::py_config()
  python_env_path <- py_config$python
  python_site_packages <- file.path(dirname(python_env_path), "lib", "python3.12", "site-packages")

  # Check if compiled .so or .pyd exists
  compiled_so <- file.path(python_site_packages, "TSW_scoreMat.so")
  compiled_pyd <- file.path(python_site_packages, "TSW_scoreMat.pyd")

  # If neither .so nor .pyd exists, compile the Cython code
  if (!file.exists(compiled_so) && !file.exists(compiled_pyd)) {
    # Compile the Cython code
    cython_dir <- system.file("python/cython", package = "ARTEMIS")
    system(paste("python", file.path(cython_dir, "setup.py"), "build_ext --inplace"))
  }

  reticulate::source_python(system.file("./Python/init.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./Python/score.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./Python/align.py",package="ARTEMIS"),envir=globalenv())
  reticulate::source_python(system.file("./Python/main.py",package="ARTEMIS"),envir=globalenv())
}