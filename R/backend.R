# Management of the two Stan backends.
#
# The Stan programs live in inst/stan.
# - "rstan" (default): the programs are compiled when the package is installed,
#   by rstantools, and are available in the `stanmodels` list. No compilation
#   is needed at run time.
# - "cmdstanr": the programs are compiled the first time they are used and the
#   executables are saved in the user cache directory. Each executable is
#   identified by the CmdStan version and by a hash of the Stan code, so it is
#   reused across sessions and package reinstallations, and recompiled only
#   when the Stan code or the CmdStan version changes.

#' Names of the Stan programs shipped with the package
#' @noRd
stan_programs <- function() {
  files <- list.files(system.file("stan", package = "baysesp"), pattern = "\\.stan$")
  sub("\\.stan$", "", files)
}

#' Path of a Stan program shipped with the package
#' @noRd
stan_file <- function(program) {
  path <- system.file("stan", paste0(program, ".stan"), package = "baysesp")
  if (!nzchar(path)) {
    stop("Stan program '", program, "' not found in the package.", call. = FALSE)
  }
  path
}

#' Is the cmdstanr backend usable?
#' @noRd
cmdstan_available <- function() {
  requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                      error = function(e) NULL))
}

#' Choose the Stan backend, falling back to rstan when cmdstanr is unavailable
#' @noRd
resolve_backend <- function(backend = c("rstan", "cmdstanr")) {
  backend <- match.arg(backend)
  if (backend == "cmdstanr" && !cmdstan_available()) {
    message("'cmdstanr' or CmdStan is not available: using the 'rstan' backend.")
    backend <- "rstan"
  }
  backend
}

#' Cache directory for cmdstanr executables, one per CmdStan version
#' @noRd
cache_dir <- function() {
  dir <- file.path(tools::R_user_dir("baysesp", which = "cache"),
                   paste0("cmdstan-", cmdstanr::cmdstan_version()))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

#' Hash of a Stan program, including the files in inst/stan/include
#' @noRd
stan_hash <- function(program) {
  include_files <- list.files(system.file("stan", "include", package = "baysesp"),
                              full.names = TRUE)
  files <- c(stan_file(program), sort(include_files))
  code <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  hash_file <- tempfile()
  on.exit(unlink(hash_file))
  writeLines(code, hash_file)
  substr(unname(tools::md5sum(hash_file)), 1, 12)
}

#' Compiled cmdstanr model, reused across sessions and reinstallations
#' @noRd
cmdstan_program <- function(program) {
  # The Stan file is copied into the cache under a name that contains its hash:
  # a reinstalled but unchanged program maps to the same, already compiled, file.
  cached_file <- file.path(cache_dir(), paste0(program, "_", stan_hash(program), ".stan"))
  if (!file.exists(cached_file)) {
    file.copy(stan_file(program), cached_file)
  }
  cmdstanr::cmdstan_model(
    cached_file,
    include_paths = system.file("stan", package = "baysesp"),
    quiet = TRUE
  )
}

#' Compile the Stan programs for the cmdstanr backend
#'
#' The default backend, `rstan` (Stan Development Team, 2025), uses programs
#' written in Stan (Carpenter et al., 2017) and compiled when the package is
#' installed. With the optional `cmdstanr` backend each Stan program is
#' compiled the first time it is used and the executable is saved in the user
#' cache directory. Call this function once after installing or updating the
#' package to compile all programs in advance.
#'
#' @param programs Names of the Stan programs to compile. By default, all.
#' @return The paths of the compiled executables, invisibly.
#' @references
#' Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1), 1–32. \doi{10.18637/jss.v076.i01}
#'
#' Stan Development Team (2025). RStan: The R interface to Stan. R package. <https://CRAN.R-project.org/package=rstan>
#' @export
compile_models <- function(programs = stan_programs()) {
  if (!cmdstan_available()) {
    stop("'cmdstanr' and CmdStan are required to compile the models.", call. = FALSE)
  }
  paths <- vapply(programs, function(p) cmdstan_program(p)$exe_file(), character(1))
  invisible(paths)
}
