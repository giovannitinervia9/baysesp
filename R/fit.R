#' @include models.R data.R
NULL

#' Bayesian estimation of the accuracy of diagnostic tests
#'
#' Fits a model for the accuracy of one or more binary tests under partial
#' verification. The model object holds the priors and the verification
#' mechanism; see [saturated()].
#'
#' @param data Study data, created by [baysesp_data()].
#' @param model Model, by default [saturated()] with MAR verification.
#' @param backend Stan interface: `"rstan"` (default; the programs are compiled
#'   when the package is installed) or `"cmdstanr"` (each program is compiled
#'   the first time it is used, see [compile_models()]). If `cmdstanr` or
#'   CmdStan are not available, `"rstan"` is used with a message.
#' @param chains Number of Markov chains.
#' @param iter_warmup,iter_sampling Number of warmup and sampling iterations per
#'   chain.
#' @param seed Random seed, or `NULL`.
#' @param ... Further arguments passed to [rstan::sampling()] or to the
#'   `$sample()` method of `cmdstanr`, such as `cores`/`parallel_chains`.
#'
#' @return An object of class `baysesp_fit`, with properties `@model`, `@data`,
#'   `@draws` (a `posterior::draws_array`), `@backend` and `@stan_fit` (the
#'   object returned by the backend). The draws contain `prevalence`,
#'   `theta_healthy` and `theta_diseased` (the distributions of the profiles in
#'   the healthy and in the diseased) and, for the saturated model, the
#'   parameters of the transparent parametrization: `omega`, `v` and `phi1`
#'   (0 for never-verified profiles), `q` if some profile is never verified and
#'   `beta` with [mnar()] verification.
#' @seealso [test_accuracy()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_study(n = 1000, prevalence = 0.2,
#'                     sensitivity = c(A = 0.9, B = 0.8),
#'                     specificity = c(A = 0.85, B = 0.95),
#'                     verification_prob = c(0.2, 0.9, 0.9, 1))
#' d <- baysesp_data(d, tests = c("A", "B"), verified = "verified", disease = "disease")
#' fit <- baysesp(d, seed = 1)
#' test_accuracy(fit)
#' }
#' @export
baysesp <- function(data, model = saturated(), backend = c("rstan", "cmdstanr"),
                    chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                    seed = NULL, ...) {
  if (!S7::S7_inherits(data, baysesp_data)) {
    stop("`data` must be created with baysesp_data().", call. = FALSE)
  }
  if (!S7::S7_inherits(model, accuracy_model)) {
    stop("`model` must be a model, such as saturated().", call. = FALSE)
  }
  backend <- resolve_backend(backend)
  sdata <- stan_data(model, model@verification, data)
  program <- stan_program(model)
  variables <- c("prevalence", "theta_healthy", "theta_diseased",
                 "omega", "v", "phi1", if (sdata$J0 > 0L) "q", if (sdata$n_beta > 0L) "beta")

  if (backend == "rstan") {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
    stan_fit <- rstan::sampling(
      stanmodels[[program]], data = sdata, pars = variables,
      chains = chains, iter = iter_warmup + iter_sampling, warmup = iter_warmup,
      seed = seed, refresh = 0, ...
    )
    draws <- posterior::as_draws_array(rstan::extract(stan_fit, pars = variables,
                                                      permuted = FALSE))
  } else {
    stan_fit <- cmdstan_program(program)$sample(
      data = sdata, chains = chains, iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, seed = seed, refresh = 0,
      show_messages = FALSE, ...
    )
    draws <- stan_fit$draws(variables)
  }

  baysesp_fit(model = model, data = data, draws = draws, backend = backend,
              stan_fit = stan_fit)
}

#' @rdname baysesp
#' @usage NULL
#' @export
baysesp_fit <- S7::new_class(
  "baysesp_fit",
  properties = list(
    model = accuracy_model,
    data = baysesp_data,
    draws = S7::class_any,
    backend = S7::class_character,
    stan_fit = S7::class_any
  )
)

S7::method(print, baysesp_fit) <- function(x, ...) {
  cat("<baysesp_fit> ", sub("^baysesp::", "", class(x@model)[1]), " model, ",
      toupper(sub("^baysesp::", "", class(x@model@verification)[1])),
      " verification, ", posterior::nchains(x@draws), " chains x ",
      posterior::niterations(x@draws), " draws (", x@backend, ")\n\n", sep = "")
  print(test_accuracy(x), row.names = FALSE, digits = 3)
  invisible(x)
}

# Prevalence and distributions of the profiles, the parameters shared by all
# models, as a vector and two draws x profiles matrices.
core_draws <- function(fit) {
  m <- unclass(posterior::as_draws_matrix(fit@draws))
  J <- nrow(fit@data@profiles)
  list(
    prevalence = unname(m[, "prevalence"]),
    theta_healthy = unname(m[, paste0("theta_healthy[", seq_len(J), "]"), drop = FALSE]),
    theta_diseased = unname(m[, paste0("theta_diseased[", seq_len(J), "]"), drop = FALSE])
  )
}
