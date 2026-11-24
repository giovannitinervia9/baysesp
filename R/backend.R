# Gestione dei due backend di Stan.
#
# I programmi Stan stanno in inst/stan. Con il backend "rstan" sono compilati
# all'installazione da rstantools e disponibili nella lista `stanmodels`.
# Con il backend "cmdstanr" sono compilati al primo uso e salvati nella cache
# dell'utente, in una sottocartella per ogni versione di CmdStan.

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
resolve_backend <- function(backend = c("cmdstanr", "rstan")) {
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

#' Compiled cmdstanr model, reused across sessions
#' @noRd
cmdstan_program <- function(program) {
  cmdstanr::cmdstan_model(
    stan_file(program),
    dir = cache_dir(),
    include_paths = system.file("stan", "include", package = "baysesp"),
    quiet = TRUE
  )
}

#' Compile the Stan programs for the cmdstanr backend
#'
#' With the `cmdstanr` backend each Stan program is compiled the first time it
#' is used and the executable is saved in the user cache directory. Call this
#' function once after installing or updating the package to compile all
#' programs in advance.
#'
#' @param programs Names of the Stan programs to compile. By default, all.
#' @return The paths of the compiled executables, invisibly.
#' @export
compile_models <- function(programs = stan_programs()) {
  if (!cmdstan_available()) {
    stop("'cmdstanr' and CmdStan are required to compile the models.", call. = FALSE)
  }
  paths <- vapply(programs, function(p) cmdstan_program(p)$exe_file(), character(1))
  invisible(paths)
}
